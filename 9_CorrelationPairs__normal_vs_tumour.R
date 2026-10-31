# 1. Load required libraries
library(tidyverse) # Includes dplyr, tidyr, readr, ggplot2, purrr

##############################################################
## INPUT & CONFIGURATION
##############################################################

# 2. Load your RNA-seq TPM matrix
tpm_data <- read_csv("SRA_TPM_CHOL_lncRNA_mRNA.csv")

# 3. Define target genes and sample cohorts
lncRNAs <- c("lnc-ADH4-3:2", "lnc-ADH4-3:3", "lnc-CYP2A6-2:1", "LINC01093:1")
mRNAs   <- c("LINC01093", "ADH4", "CYP2A6")
target_genes <- c(lncRNAs, mRNAs)

# Define control and tumour sample identifiers
control_samples <- paste0("C", 1:15)
tumor_samples   <- paste0("T", 1:15)
all_samples     <- c(control_samples, tumor_samples)

# Create output directory
if(!dir.exists("CorrelationPlots")) dir.create("CorrelationPlots")

##############################################################
## DATA CLEANING & SHAPING
##############################################################

# 4. Filter, clean, and shape the data

expression_clean <- tpm_data %>%
  filter(TranscriptSymbol %in% target_genes) %>%
  pivot_longer(
    cols = all_of(all_samples),
    names_to = "Patient",
    values_to = "TPM"
  ) %>%
  pivot_wider(
    names_from = TranscriptSymbol,
    values_from = TPM
  ) %>%
  mutate(
    Cohort = if_else(
      Patient %in% control_samples,
      "Control",
      "Tumor"
    )
  )

# Generate all possible combination pairs dynamically

analysis_pairs <- expand.grid(
  lncRNA = lncRNAs,
  mRNA = mRNAs,
  stringsAsFactors = FALSE
) %>%
  split(seq(nrow(.))) %>%
  map(~ c(.$lncRNA, .$mRNA))

##############################################################
## CORRELATION ANALYSIS (Control Only vs. Tumor Only)
##############################################################

# 5. Function to run correlation tests safely

run_correlation <- function(df, pairs_list, analysis_label) {
  
  map_df(pairs_list, function(pair) {
    
    lnc <- pair[1]
    pcg <- pair[2]
    
    x <- df[[lnc]]
    y <- df[[pcg]]
    
    spearman_test <- cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
    
    tibble(
      Analysis        = analysis_label,
      lncRNA          = lnc,
      Target_Gene     = pcg,
      Spearman_rho    = spearman_test$estimate,
      Spearman_pvalue = spearman_test$p.value
    )
  })
}

# 6. Execute subgroup analysis (Control vs Tumor)

results_control <- run_correlation(
  filter(expression_clean, Cohort == "Control"),
  analysis_pairs,
  "Control"
)

results_tumor <- run_correlation(
  filter(expression_clean, Cohort == "Tumor"),
  analysis_pairs,
  "Tumor"
)

# Combine results and run FDR correction across all tests

final_results <- bind_rows(
  results_control,
  results_tumor
) %>%
  mutate(
    Spearman_FDR = p.adjust(
      Spearman_pvalue,
      method = "BH"
    )
  )

# Save the master results table

write_csv(
  final_results,
  "Correlation_Results.csv"
)

##############################################################
## FILTERING TARGET PAIRS FOR PLOTTING (EITHER MUST SATISFY)
##############################################################

# 7. Identify pairs meeting thresholds in either subgroup

significant_pairs_df <- final_results %>%
  mutate(
    meets_cutoff = (
      Spearman_rho > 0.4 |
        Spearman_rho < -0.29
    ) & Spearman_FDR < 0.05
  ) %>%
  group_by(lncRNA, Target_Gene) %>%
  filter(any(meets_cutoff)) %>%
  distinct(lncRNA, Target_Gene)

##############################################################
## PLOTTING & PRINTING
##############################################################

if(nrow(significant_pairs_df) == 0) {
  
  cat("\n=======================================================\n")
  cat("No pairs matched the criteria in EITHER 'Control' or 'Tumor' subgroups.\n")
  cat("=======================================================\n")
  
} else {
  
  # 8. Create text labels with dynamically computed layout positions
  
  annotation_labels <- final_results %>%
    semi_join(
      significant_pairs_df,
      by = c("lncRNA", "Target_Gene")
    ) %>%
    mutate(
      FDR_fmt = if_else(
        Spearman_FDR < 0.001,
        formatC(
          Spearman_FDR,
          format = "e",
          digits = 2
        ),
        sprintf("%.3f", Spearman_FDR)
      ),
      Label_Line = sprintf(
        "%s ρ = %.2f\nFDR = %s",
        Analysis,
        Spearman_rho,
        FDR_fmt
      )
    ) %>%
    group_by(lncRNA, Target_Gene) %>%
    summarize(
      Stats_Label = paste(
        Label_Line,
        collapse = "\n"
      ),
      .groups = "drop"
    )
  
  # 9. Isolate only the plot data needed for target combinations
  
  significant_pairs_list <- split(
    significant_pairs_df,
    seq(nrow(significant_pairs_df))
  ) %>%
    map(~ c(.$lncRNA, .$Target_Gene))
  
  plot_data <- map_df(
    significant_pairs_list,
    function(pair) {
      
      expression_clean %>%
        select(
          Patient,
          Cohort,
          lncRNA_val = all_of(pair[1]),
          Target_val = all_of(pair[2])
        ) %>%
        mutate(
          lncRNA = pair[1],
          Target_Gene = pair[2]
        )
    }
  )
  
  ############################################################
  ## 10. GENERATE FILTERED MULTI-PANEL FACET GRID
  ############################################################
  
  # Exact log10(TPM + 1) transformation for both axes
  
  log10_tpm_trans <- scales::trans_new(
    name = "log10_TPM_plus_1",
    transform = function(x) log10(x + 1),
    inverse = function(x) (10^x) - 1
  )
  
  scatter_plots <- ggplot(
    plot_data,
    aes(
      x = lncRNA_val,
      y = Target_val
    )
  ) +
    
    # Control linear regression trend
    geom_smooth(
      data = filter(plot_data, Cohort == "Control"),
      aes(color = "Control Only (C1-C15)"),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      alpha = 0.1
    ) +
    
    # Tumour linear regression trend
    geom_smooth(
      data = filter(plot_data, Cohort == "Tumor"),
      aes(color = "Tumor Only (T1-T15)"),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      alpha = 0.15
    ) +
    
    geom_point(
      aes(
        shape = Cohort,
        fill = Cohort
      ),
      size = 2.5,
      alpha = 0.7
    ) +
    
    geom_text(
      data = annotation_labels,
      aes(label = Stats_Label),
      x = -Inf,
      y = Inf,
      hjust = -0.05,
      vjust = 1.05,
      size = 3.0,
      fontface = "plain",
      inherit.aes = FALSE
    ) +
    
    facet_grid(
      lncRNA ~ Target_Gene,
      scales = "free"
    ) +
    
    # CORRECTED: Exact log10(TPM + 1) transformation
    
    scale_x_continuous(
      trans = log10_tpm_trans,
      labels = function(x) round(log10(x + 1), 2)
    ) +
    
    scale_y_continuous(
      trans = log10_tpm_trans,
      labels = function(x) round(log10(x + 1), 2)
    ) +
    
    scale_shape_manual(
      values = c(
        "Control" = 24,
        "Tumor" = 21
      )
    ) +
    
    scale_fill_manual(
      values = c(
        "Control" = "#2c3e50",
        "Tumor" = "#e74c3c"
      )
    ) +
    
    scale_color_manual(
      values = c(
        "Control Only (C1-C15)" = "#2c3e50",
        "Tumor Only (T1-T15)" = "#e74c3c"
      )
    ) +
    
    labs(
      title = "Filtered lncRNA Expression Correlations Across Patients",
      
      subtitle = "Panels filtered by: Either 'Control' or 'Tumor' satisfies (rho > 0.4 OR rho < -0.29) & FDR < 0.05",
      
      # CORRECTED: Labels now match the actual transformation
      
      x = "lncRNA Expression Level log10(TPM + 1)",
      
      y = "Target PCG Expression Level log10(TPM + 1)",
      
      fill = "Patient Cohort",
      shape = "Patient Cohort",
      
      # Clarified that the displayed trends are linear regression trends
      
      color = "Linear Regression Trend"
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      strip.text = element_text(
        face = "bold",
        size = 10
      ),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    )
  
  # Output image to session
  
  print(scatter_plots)
  
  # Save image file
  
  ggsave(
    "CorrelationPlots/Filtered_lncRNA_correlation_grid.png",
    plot = scatter_plots,
    width = 15,
    height = 12,
    dpi = 1200
  )
}

##############################################################
## FINAL PRINTOUT
##############################################################

cat("\n=======================================================\n")
cat("FINAL FILTERED RESULTS SUMMARY (rho > 0.4 OR rho < -0.29):\n")
cat("=======================================================\n")

final_filtered_print <- final_results %>%
  filter(
    (Spearman_rho > 0.4 | Spearman_rho < -0.29) &
      Spearman_FDR < 0.05
  )

if(nrow(final_filtered_print) == 0) {
  
  print("Empty dataset (0 individual analysis lines met thresholds).")
  
} else {
  
  print(
    as.data.frame(final_filtered_print),
    row.names = FALSE
  )
}

cat("=======================================================\n")

write.csv(
  final_filtered_print,
  "Correlation_Results_significant_data.csv",
  row.names = FALSE
)

