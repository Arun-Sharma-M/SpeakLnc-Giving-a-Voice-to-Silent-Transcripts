############################################################
# 1. LOAD LIBRARIES
############################################################

library(survival)
library(survminer)
library(tidyverse)
library(stringr)


############################################################
# 2. LOAD CLINICAL DATA
############################################################

clinical_CHOL <- read.csv(
  "clinical_CHOL_TCGA.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Check survival status
print(table(clinical_CHOL$vital_status, useNA = "ifany"))

# Event: 1 = dead; 0 = alive/censored
clinical_CHOL$deceased <- ifelse(
  clinical_CHOL$vital_status == "Alive", 0,
  ifelse(clinical_CHOL$vital_status == "Dead", 1, NA)
)

# Convert days to months
clinical_CHOL$Months_to_last_follow_up <-
  as.numeric(clinical_CHOL$days_to_last_follow_up) / 30.4375

clinical_CHOL$Months_to_death <-
  as.numeric(clinical_CHOL$days_to_death) / 30.4375

# Calculate overall survival
clinical_CHOL$overall_survival <- ifelse(
  clinical_CHOL$deceased == 1,
  clinical_CHOL$Months_to_death,
  clinical_CHOL$Months_to_last_follow_up
)

# Retain valid clinical records
clinical_CHOL <- clinical_CHOL %>%
  filter(
    !is.na(samplename),
    trimws(samplename) != "",
    !is.na(deceased),
    !is.na(overall_survival),
    is.finite(overall_survival),
    overall_survival >= 0
  )

cat("Valid clinical records:", nrow(clinical_CHOL), "\n")
print(table(clinical_CHOL$deceased))


############################################################
# 3. LOAD EXPRESSION DATA
############################################################

TCGA_CHOL_TPM_Data_full <- read.csv(
  "TPM_CHOL_lncRNA_Tumour_Deseq.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Select the 35 tumour samples
TCGA_CHOL_TPM_Data <- TCGA_CHOL_TPM_Data_full[
  , 7:41, drop = FALSE
]


############################################################
# 4. LOAD CANDIDATE GENES
############################################################

transcript_data <- read.csv(
  "1.Merge_SRA_TCGA_lncRNA_Pvalue_L2FC_dwnr.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

genes_to_analyze <- unique(
  na.omit(transcript_data$Merge_SRA_TCGA_lncRNA_Pvalue_L2FC)
)

# Keep genes present in the expression matrix
genes_to_analyze <- intersect(
  genes_to_analyze,
  rownames(TCGA_CHOL_TPM_Data)
)

cat("Genes to analyse:", length(genes_to_analyze), "\n")


############################################################
# 5. RESHAPE EXPRESSION DATA
############################################################

exp_long <- TCGA_CHOL_TPM_Data %>%
  as.data.frame() %>%
  rownames_to_column(var = "GeneSymbol") %>%
  pivot_longer(
    cols = -GeneSymbol,
    names_to = "samplename",
    values_to = "TPM"
  )

exp_long$TPM <- as.numeric(exp_long$TPM)


############################################################
# 6. CHECK SAMPLE MATCHING
############################################################

cat(
  "Expression samples:",
  length(unique(exp_long$samplename)), "\n"
)

cat(
  "Clinical samples:",
  length(unique(clinical_CHOL$samplename)), "\n"
)

cat(
  "Common samples:",
  length(intersect(
    unique(exp_long$samplename),
    unique(clinical_CHOL$samplename)
  )), "\n"
)


############################################################
# 7. INITIALISE RESULTS AND PLOT LISTS
############################################################

survival_results <- list()
median_plots <- list()
quantile_plots <- list()


############################################################
# 8. SURVIVAL ANALYSIS FOR EACH GENE
############################################################

for (gene in genes_to_analyze) {
  
  cat("Analysing:", gene, "\n")
  
  # Match expression with clinical data
  gene_data <- exp_long %>%
    filter(GeneSymbol == gene) %>%
    inner_join(clinical_CHOL, by = "samplename") %>%
    filter(
      !is.na(TPM),
      is.finite(TPM),
      !is.na(overall_survival),
      is.finite(overall_survival),
      !is.na(deceased),
      overall_survival >= 0
    )
  
  # Avoid counting duplicate patient records
  if (anyDuplicated(gene_data$samplename)) {
    warning("Duplicate sample IDs for gene: ", gene)
    gene_data <- gene_data %>%
      distinct(samplename, .keep_all = TRUE)
  }
  
  if (nrow(gene_data) < 2 ||
      length(unique(gene_data$deceased)) < 2) {
    next
  }
  
  
  ##########################################################
  # 8A. MEDIAN CUTOFF: HIGH VS LOW
  ##########################################################
  
  median_cutoff <- median(
    gene_data$TPM,
    na.rm = TRUE
  )
  
  gene_data$median_group <- ifelse(
    gene_data$TPM >= median_cutoff,
    "High",
    "Low"
  )
  
  gene_data$median_group <- factor(
    gene_data$median_group,
    levels = c("Low", "High")
  )
  
  if (length(unique(gene_data$median_group)) == 2) {
    
    # Kaplan-Meier model
    km_med <- survfit(
      Surv(overall_survival, deceased) ~ median_group,
      data = gene_data
    )
    
    # Log-rank test
    logrank_med <- survdiff(
      Surv(overall_survival, deceased) ~ median_group,
      data = gene_data
    )
    
    p_logrank_med <- pchisq(
      logrank_med$chisq,
      df = length(logrank_med$n) - 1,
      lower.tail = FALSE
    )
    
    # Cox regression: High vs Low
    cox_med <- coxph(
      Surv(overall_survival, deceased) ~ median_group,
      data = gene_data
    )
    
    sm_med <- summary(cox_med)
    
    HR_med <- sm_med$coefficients[1, "exp(coef)"]
    p_cox_med <- sm_med$coefficients[1, "Pr(>|z|)"]
    
    CI_med <- sm_med$conf.int[
      1, c("lower .95", "upper .95")
    ]
    
    # Store results
    survival_results[[paste0(gene, "_Median")]] <- data.frame(
      Gene = gene,
      Cutoff = "Median",
      Cutoff_Value = median_cutoff,
      
      Low_N = sum(gene_data$median_group == "Low"),
      High_N = sum(gene_data$median_group == "High"),
      
      Deaths_Low = sum(
        gene_data$median_group == "Low" &
          gene_data$deceased == 1
      ),
      
      Deaths_High = sum(
        gene_data$median_group == "High" &
          gene_data$deceased == 1
      ),
      
      Comparison = "High vs Low",
      HR = HR_med,
      HR_CI_Lower = CI_med[1],
      HR_CI_Upper = CI_med[2],
      Cox_P = p_cox_med,
      Logrank_P = p_logrank_med
    )
    
    # ORIGINAL MEDIAN PLOTTING STYLE
    if (is.finite(p_logrank_med) && p_logrank_med < 0.05) {
      
      p_median <- ggsurvplot(
        km_med,
        data = gene_data,
        pval = TRUE,
        risk.table = TRUE,
        title = paste(
          "TCGA_CHOL", gene, "Median_value",
          sep = "_"
        ),
        surv.median.line = "hv",
        xlab = "Follow up in Months"
      )
      
      median_plots[[gene]] <- p_median
      
      message(
        paste(
          "Significant (p < 0.05) - Generated median plot for:",
          gene
        )
      )
      
      print(p_median)
    }
  }
  
  
  ##########################################################
  # 8B. Q25 CUTOFF: MEDIUM & HIGH VS LOW
  ##########################################################
  
  q25_cutoff <- as.numeric(
    quantile(
      gene_data$TPM,
      probs = 0.25,
      na.rm = TRUE
    )
  )
  
  gene_data$q25_group <- ifelse(
    gene_data$TPM > q25_cutoff,
    "Medium & High",
    "Low"
  )
  
  gene_data$q25_group <- factor(
    gene_data$q25_group,
    levels = c("Low", "Medium & High")
  )
  
  if (length(unique(gene_data$q25_group)) == 2) {
    
    # Kaplan-Meier model
    km_q25 <- survfit(
      Surv(overall_survival, deceased) ~ q25_group,
      data = gene_data
    )
    
    # Log-rank test
    logrank_q25 <- survdiff(
      Surv(overall_survival, deceased) ~ q25_group,
      data = gene_data
    )
    
    p_logrank_q25 <- pchisq(
      logrank_q25$chisq,
      df = length(logrank_q25$n) - 1,
      lower.tail = FALSE
    )
    
    # Cox regression: Medium & High vs Low
    cox_q25 <- coxph(
      Surv(overall_survival, deceased) ~ q25_group,
      data = gene_data
    )
    
    sm_q25 <- summary(cox_q25)
    
    HR_q25 <- sm_q25$coefficients[1, "exp(coef)"]
    p_cox_q25 <- sm_q25$coefficients[1, "Pr(>|z|)"]
    
    CI_q25 <- sm_q25$conf.int[
      1, c("lower .95", "upper .95")
    ]
    
    # Store results
    survival_results[[paste0(gene, "_Q25")]] <- data.frame(
      Gene = gene,
      Cutoff = "Q25",
      Cutoff_Value = q25_cutoff,
      
      Low_N = sum(gene_data$q25_group == "Low"),
      High_N = sum(
        gene_data$q25_group == "Medium & High"
      ),
      
      Deaths_Low = sum(
        gene_data$q25_group == "Low" &
          gene_data$deceased == 1
      ),
      
      Deaths_High = sum(
        gene_data$q25_group == "Medium & High" &
          gene_data$deceased == 1
      ),
      
      Comparison = "Medium & High vs Low",
      HR = HR_q25,
      HR_CI_Lower = CI_q25[1],
      HR_CI_Upper = CI_q25[2],
      Cox_P = p_cox_q25,
      Logrank_P = p_logrank_q25
    )
    
    # ORIGINAL Q25 PLOTTING STYLE
    if (is.finite(p_logrank_q25) && p_logrank_q25 < 0.05) {
      
      p_quant <- ggsurvplot(
        km_q25,
        data = gene_data,
        pval = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        title = paste(
          "TCGA_CHOL", gene, "quantile",
          sep = "_"
        ),
        xlab = "Follow up in Months",
        legend.labs = c("Low", "Medium & High"),
        ncensor.plot = TRUE,
        surv.median.line = "hv",
        ggtheme = theme_bw(),
        palette = c("#2E9FDF", "#e7002f")
      )
      
      quantile_plots[[gene]] <- p_quant
      
      message(
        paste(
          "Significant (p < 0.05) - Generated quantile plot for:",
          gene
        )
      )
      
      print(p_quant)
    }
  }
}


############################################################
# 9. COMBINE RESULTS AND CALCULATE FDR
############################################################

survival_results_df <- bind_rows(survival_results)

if (nrow(survival_results_df) == 0) {
  stop(
    "No valid survival results were generated. ",
    "Check clinical event coding and sample matching."
  )
}

# Adjust Cox and log-rank P-values separately
# within each cutoff family
survival_results_df <- survival_results_df %>%
  group_by(Cutoff) %>%
  mutate(
    Cox_FDR = p.adjust(Cox_P, method = "BH"),
    Logrank_FDR = p.adjust(Logrank_P, method = "BH")
  ) %>%
  ungroup()


############################################################
# 10. CREATE FINAL RESULTS TABLE
############################################################

survival_table <- survival_results_df %>%
  mutate(
    `Expression cutoff` = Cutoff,
    
    `HR (95% CI)` = sprintf(
      "%.3f (%.3f-%.3f)",
      HR,
      HR_CI_Lower,
      HR_CI_Upper
    ),
    
    `Cox P-value` = Cox_P,
    `Cox FDR` = Cox_FDR,
    `Log-rank P-value` = Logrank_P,
    `Log-rank FDR` = Logrank_FDR
  ) %>%
  select(
    Transcript = Gene,
    `Expression cutoff`,
    `Low-expression group (n)` = Low_N,
    `High-expression group (n)` = High_N,
    `Deaths in low group (n)` = Deaths_Low,
    `Deaths in high group (n)` = Deaths_High,
    `HR (95% CI)`,
    `Cox P-value`,
    `Cox FDR`,
    `Log-rank P-value`,
    `Log-rank FDR`
  )


############################################################
# 11. DISPLAY RESULTS
############################################################

print(survival_table, n = Inf)

# Nominally significant log-rank results
print(
  survival_table %>%
    filter(`Log-rank P-value` < 0.05),
  n = Inf
)

# Results with FDR < 0.05
print(
  survival_table %>%
    filter(
      `Cox FDR` < 0.05 |
        `Log-rank FDR` < 0.05
    ),
  n = Inf
)
############################################################
# 12. GENERATE TABLE OF NOMINALLY SIGNIFICANT RESULTS
############################################################

significant_survival_table <- survival_table %>%
  filter(
    `Cox P-value` < 0.05 |
      `Log-rank P-value` < 0.05
  ) %>%
  arrange(`Log-rank P-value`) %>%
  select(
    Transcript,
    `Expression cutoff`,
    `Low-expression group (n)`,
    `High-expression group (n)`,
    `Deaths in low group (n)`,
    `Deaths in high group (n)`,
    `HR (95% CI)`,
    `Cox P-value`,
    `Cox FDR`,
    `Log-rank P-value`,
    `Log-rank FDR`
  )

print(significant_survival_table, n = Inf, width = Inf)
write.csv(significant_survival_table,"CHOL_lncRNA_survival_results/Main_Manuscript_Survival_Table.csv", row.names = TRUE)

