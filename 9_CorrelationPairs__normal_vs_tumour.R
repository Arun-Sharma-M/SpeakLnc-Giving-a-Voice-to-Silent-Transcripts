
# 1. Load required libraries
library(tidyverse) # Includes dplyr, tidyr, readr, ggplot2, purrr

##############################################################
## INPUT & CONFIGURATION
##############################################################
# 2. Load your RNA-seq TPM matrix
tpm_data <- read_csv("SRA_TPM_CHOL_lncRNA_mRNA.csv")

# 3. Define target genes and sample cohorts
lncRNAs <- c("lnc-ADH4-3:2", "lnc-ADH4-3:3", "lnc-CYP2A6-2:1", "LINC01093:1")
mRNAs   <- c("LINC01093","ADH4", "CYP2A6")
target_genes <- c(lncRNAs, mRNAs)

# >>> ADDED COMMENT: Define the control and tumour sample identifiers expected in the TPM matrix and combine them into a single vector for subsequent data reshaping.
control_samples <- paste0("C", 1:15)
tumor_samples   <- paste0("T", 1:15)
all_samples     <- c(control_samples, tumor_samples)

# >>> ADDED COMMENT: Create the output directory for correlation plots if it does not already exist.
if(!dir.exists("CorrelationPlots")) dir.create("CorrelationPlots")

##############################################################
## DATA CLEANING & SHAPING
##############################################################
# 4. Filter, clean, and shape the data

# >>> ADDED COMMENT: Retain the specified lncRNAs and mRNAs, reshape the expression matrix into long format and then into one row per patient, with each target transcript represented by a separate column -. 
# - Assign each patient to the control or tumour cohort based on the sample identifier.
expression_clean <- tpm_data %>%
  filter(TranscriptSymbol %in% target_genes) %>%
  pivot_longer(cols = all_of(all_samples), names_to = "Patient", values_to = "TPM") %>%
  pivot_wider(names_from = TranscriptSymbol, values_from = TPM) %>%
  mutate(Cohort = if_else(Patient %in% control_samples, "Control", "Tumor"))

# Generate all possible combination pairs dynamically

# >>> ADDED COMMENT: Generate every possible lncRNA–mRNA combination so that each specified pair can be tested systematically without manually listing the combinations.
analysis_pairs <- expand.grid(lncRNA = lncRNAs, mRNA = mRNAs, stringsAsFactors = FALSE) %>%
  split(seq(nrow(.))) %>%
  map(~ c(.$lncRNA, .$mRNA))

##############################################################
## CORRELATION ANALYSIS (Control Only vs. Tumor Only)
##############################################################

# 5. Function to run correlation tests safely
# >>> ADDED COMMENT: Define a function that performs a separate Spearman rank-correlation test for each lncRNA–mRNA pair within a supplied cohort and returns the correlation coefficient, raw p-value and analysis-group label.
run_correlation <- function(df, pairs_list, analysis_label) {
  map_df(pairs_list, function(pair) {
    lnc <- pair[1]
    pcg <- pair[2]
    
    x <- df[[lnc]]
    y <- df[[pcg]]
    
    spearman_test <- cor.test(x, y, method = "spearman", exact = FALSE)
    
    tibble(
      Analysis        = analysis_label,
      lncRNA          = lnc,
      Target_Gene     = pcg,
      Spearman_rho    = spearman_test$estimate,
      Spearman_pvalue = spearman_test$p.value
    )
  })
}

# 6. Execute Subgroup Analysis (Control vs Tumor)
# >>> ADDED COMMENT: Apply the correlation function independently to control and tumour samples, allowing the direction and strength of each lncRNA–mRNA association to be assessed separately in the two cohorts.
results_control <- run_correlation(filter(expression_clean, Cohort == "Control"), analysis_pairs, "Control")
results_tumor   <- run_correlation(filter(expression_clean, Cohort == "Tumor"), analysis_pairs, "Tumor")

# >>> ADDED COMMENT: Combine the control and tumour correlation results into one table and apply Benjamini–Hochberg adjustment across all correlation tests to obtain FDR-adjusted p-values.
# Combine results and run FDR correction across all tests
final_results <- bind_rows(results_control, results_tumor) %>%
  mutate(Spearman_FDR = p.adjust(Spearman_pvalue, method = "BH"))

# Save the master text results table
write_csv(final_results, "Correlation_Results.csv")

##############################################################
## FILTERING TARGET PAIRS FOR PLOTTING (EITHER MUST SATISFY)
##############################################################
# >>> ADDED COMMENT: Flag correlation results satisfying the specified effect-size and FDR criteria, then retain each lncRNA–mRNA pair if at least one of its cohort-specific results meets both thresholds.
# 7. Identify pairs meeting thresholds in EITHER subgroup (Control OR Tumor)
significant_pairs_df <- final_results %>%
  mutate(meets_cutoff = (Spearman_rho > 0.4 | Spearman_rho < -0.29) & Spearman_FDR < 0.05) %>%
  group_by(lncRNA, Target_Gene) %>%
  # Check if AT LEAST ONE entry ('Control' or 'Tumor') for this pair passed the cutoff
  filter(any(meets_cutoff)) %>% 
  distinct(lncRNA, Target_Gene)

##############################################################
## PLOTTING & PRINTING
##############################################################
# >>> ADDED COMMENT: Check whether any lncRNA–mRNA pairs satisfy the selection criteria before attempting to construct the correlation plots.
if(nrow(significant_pairs_df) == 0) {
  
  cat("\n=======================================================\n")
  cat("No pairs matched the criteria in EITHER 'Control' or 'Tumor' subgroups.\n")
  cat("=======================================================\n")
  
} else {
  
  # 8. Create text labels with dynamically computed layout positions per grid cell
  # >>> ADDED COMMENT: Format the FDR values and construct annotation labels containing the cohort name, Spearman correlation coefficient and adjusted p-value for each selected pair.
  annotation_labels <- final_results %>%
    semi_join(significant_pairs_df, by = c("lncRNA", "Target_Gene")) %>%
    mutate(
      FDR_fmt = if_else(Spearman_FDR < 0.001, formatC(Spearman_FDR, format = "e", digits = 2), sprintf("%.3f", Spearman_FDR)),
      Label_Line = sprintf("%s ρ = %.2f\nFDR = %s", Analysis, Spearman_rho, FDR_fmt)
    ) %>%
    group_by(lncRNA, Target_Gene) %>%
    summarize(
      Stats_Label = paste(Label_Line, collapse = "\n"), # Double break between groups for clarity
      .groups = "drop"
    )
  
  # 9. Isolate only the plot data needed for target combinations
  # >>> ADDED COMMENT: Create a list of the selected lncRNA–mRNA combinations and extract the corresponding expression values and cohort labels for plotting.
  significant_pairs_list <- split(significant_pairs_df, seq(nrow(significant_pairs_df))) %>%
    map(~ c(.$lncRNA, .$Target_Gene))
  
  # >>> ADDED COMMENT: Assemble the expression data for each selected pair into a single plotting table, retaining the patient identifier, cohort, lncRNA expression and target-gene expression.
  plot_data <- map_df(significant_pairs_list, function(pair) {
    expression_clean %>%
      select(Patient, Cohort, lncRNA_val = all_of(pair[1]), Target_val = all_of(pair[2])) %>%
      mutate(lncRNA = pair[1], Target_Gene = pair[2])
  })
  
  # 10. Generate the filtered Multi-Panel Facet Grid
  # >>> ADDED COMMENT: Construct a faceted scatter-plot grid for the selected lncRNA–mRNA pairs. Display control and tumour observations separately, overlay cohort-specific linear regression trends with confidence intervals, annotate the correlation statistics, and use pseudo-logarithmic axes to accommodate expression values spanning a wide range.
  scatter_plots <- ggplot(plot_data, aes(x = lncRNA_val, y = Target_val)) +
    geom_smooth(data = filter(plot_data, Cohort == "Control"), aes(color = "Control Only (C1-C8)"), method = "lm", formula = y ~ x, se = TRUE, alpha = 0.1) +
    geom_smooth(data = filter(plot_data, Cohort == "Tumor"), aes(color = "Tumor Only (T1-T58)"), method = "lm", formula = y ~ x, se = TRUE, alpha = 0.15) +
    geom_point(aes(shape = Cohort, fill = Cohort), size = 2.5, alpha = 0.7) +
    geom_text(
      data = annotation_labels,
      aes(label = Stats_Label),
      x = -Inf, y = Inf,   # Snaps label to the left wall
      hjust = -0.05,       # Pushes label slightly right to avoid clipping the axis line
      vjust = 1.05,        # Keeps label just below the top ceiling line
      size = 3.0,          # Slightly smaller text size to fit both stats blocks cleanly
      fontface = "plain",
      inherit.aes = FALSE 
    ) +
    facet_grid(lncRNA ~ Target_Gene, scales = "free") +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    scale_shape_manual(values = c("Control" = 24, "Tumor" = 21)) + 
    scale_fill_manual(values = c("Control" = "#2c3e50", "Tumor" = "#e74c3c")) +
    scale_color_manual(values = c("Control Only (C1-C8)" = "#2c3e50", "Tumor Only (T1-T58)" = "#e74c3c")) +
    labs(
      title = "Filtered lncRNA Expression Correlations Across Patients",
      subtitle = "Panels filtered by: Either 'Control' or 'Tumor' satisfies (rho > 0.4 OR rho < -0.29) & FDR < 0.05",
      x = "lncRNA Expression Level log10(TPM + 1)",
      y = "Target PCG Expression Level log10(TPM + 1)",
      fill = "Patient Cohort", 
      shape = "Patient Cohort",
      color = "Regression Trend"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # Output image to session
  print(scatter_plots)
  
  # Save image file
  ggsave("CorrelationPlots/Filtered_lncRNA_correlation_grid.png", plot = scatter_plots, width = 11, height = 12, dpi = 300)
}

##############################################################
## FINAL PRINTOUT
##############################################################

# >>> ADDED COMMENT: Print a summary of the correlation results that satisfy the specified Spearman coefficient and FDR thresholds, then export those qualifying cohort-specific results to a separate CSV file.
cat("\n=======================================================\n")
cat("FINAL FILTERED RESULTS SUMMARY (rho > 0.4 OR rho < -0.29):\n")
cat("=======================================================\n")
final_filtered_print <- final_results %>% 
  filter((Spearman_rho > 0.4 | Spearman_rho < -0.29) & Spearman_FDR < 0.05)

if(nrow(final_filtered_print) == 0) {
  print("Empty dataset (0 individual analysis lines met thresholds).")
} else {
  print(as.data.frame(final_filtered_print), row.names = FALSE)
}
cat("=======================================================\n")
write.csv( final_filtered_print, "Correlation_Results_significant_data.csv", row.names=FALSE )
