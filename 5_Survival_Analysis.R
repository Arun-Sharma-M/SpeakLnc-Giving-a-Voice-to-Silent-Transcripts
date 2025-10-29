# script to run survival analysis using TCGA data
# setwd("~/Desktop/demo/survivalAnalysis")

#BiocManager::install('TCGAbiolinks')

#library(TCGAbiolinks)
library(survminer)
library(survival)
library(tidyverse)
library(stringr)

# getting clinical data for TCGA-CHOL cohort -------------------
#clinical_CHOL <- GDCquery_clinic("TCGA-CHOL")
clinical_CHOL <- read.csv('clinical_CHOL_forlncRNA_all_35samples_old_new_merged.csv', header = TRUE, sep = ",")
#clinical_CHOL[] <- lapply(clinical_CHOL, gsub, pattern = "'--", replacement = "", fixed = TRUE)

any(colnames(clinical_CHOL) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_CHOL) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_CHOL[,c(18,22,32)]

# looking at some variables associated with survival 
table(clinical_CHOL$vital_status)

# days_to_death, that is the number of days passed from the initial diagnosis to the patient’s death (clearly, this is only relevant for dead patients)
# days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.

# change certain values the way they are encoded
clinical_CHOL$deceased <- ifelse(clinical_CHOL$vital_status == "Alive", FALSE, TRUE)
table(clinical_CHOL$deceased)

# Converting Days to Month by dividing Number of days / 30.417

clinical_CHOL$Months_to_last_follow_up <- clinical_CHOL$days_to_last_follow_up/30.417
clinical_CHOL$Months_to_death <- clinical_CHOL$days_to_death/30.417


# create an "overall survival" variable that is equal to Months_to_death
# for dead patients, and to Months_to_last_follow_up for patients who are still alive
clinical_CHOL$overall_survival <- ifelse(clinical_CHOL$vital_status == "Alive",
                                         clinical_CHOL$Months_to_last_follow_up,
                                         clinical_CHOL$Months_to_death)

#Remove entire rows with NA's in important columns with data for KM plot

#clinical_CHOL <- clinical_CHOL[-c(49:58), ]
#clinical_CHOL$samplename.x <- str_replace_all(clinical_CHOL$samplename.x, '-', '.')

# get mRNA expression data -----------

TCGA_CHOL_TPM_Data_full <- read.csv('TPM_CHOL_lncRNA_Tumour_Deseq.csv', header = TRUE, sep = ",",row.names=1)
TCGA_CHOL_TPM_Data <- TCGA_CHOL_TPM_Data_full[c(7:41)]
head(TCGA_CHOL_TPM_Data)

##########************************************************************##########

# 1. Define the genes to be analyzed in a single list
genes_to_analyze <- c("lnc-CYP2A6-2:1","LINC02362:35", "lnc-ADH4-3:2", "lnc-ADH4-3:3", "CYP4A22-AS1:3", "LINC02362:11", "LINC01093:1", "ALDH1L1-AS2:1","ALDH1L1-AS2:3")

# 2. Get the mRNA expression data and reshape it to a long format
exp_long <- TCGA_CHOL_TPM_Data %>%
  as.data.frame() %>%
  rownames_to_column(var = 'GeneSymbol') %>%
  gather(key = 'samplename', value = 'TPM', -GeneSymbol)

# 3. Create a place to store the plots if you want to view them later
median_plots <- list()
quantile_plots <- list()

# 4. Start the loop for each gene
for (gene in genes_to_analyze) {
  
  # a) Filter data for the current gene and merge with clinical data
  gene_data <- exp_long %>%
    filter(GeneSymbol == gene) %>%
    merge(clinical_CHOL, by.x = 'samplename', by.y = 'samplename')
  
  # b) Handle median-split analysis
  if (nrow(gene_data) > 0) {
    # Denote which cases have higher or lower expression than the median
    gene_data$strata_median <- ifelse(gene_data$TPM >= median(gene_data$TPM), "HIGH", "LOW")
    
    # Fit and plot survival curve for median split
    fit_median <- survfit(Surv(overall_survival, deceased) ~ strata_median, data = gene_data)
    p_median <- ggsurvplot(fit_median,
                           data = gene_data,
                           pval = TRUE,
                           risk.table = TRUE,
                           title = paste("TCGA_CHOL", gene, "Median_value", sep = "_"),
                           surv.median.line = "hv",
                           xlab = "Follow up in Months")
    
    median_plots[[gene]] <- p_median
    message(paste("Generated median plot for:", gene))
    print(p_median)
  }
  
  # c) Handle quantile-split analysis
  if (nrow(gene_data) > 0) {
    # Denote which cases have higher or lower expression than the 75th percentile
    gene_data$strata_quant <- ifelse(gene_data$TPM > quantile(gene_data$TPM, probs = 0.25), "Medium & High", "Low")
    
    # Fit and plot survival curve for quantile split
    fit_quant <- survfit(Surv(overall_survival, deceased) ~ strata_quant, data = gene_data)
    p_quant <- ggsurvplot(fit_quant,
                          data = gene_data,
                          pval = TRUE,
                          conf.int = TRUE,
                          risk.table = TRUE,
                          title = paste("TCGA_CHOL", gene, "quantile", sep = "_"),
                          xlab = "Follow up in Months",
                          legend.labs = c("Medium & High", "Low"),
                          ncensor.plot = TRUE,
                          surv.median.line = "hv",
                          ggtheme = theme_bw(),
                          palette = c("#e7002f", "#2E9FDF"))
    
    quantile_plots[[gene]] <- p_quant
    message(paste("Generated quantile plot for:", gene))
    print(p_quant)
  }
}

# 5. (Optional) Arrange and save all generated plots
# This requires the 'ggpubr' package (which is part of survminer)
# To save all median plots to a PDF:
# ggsave("all_median_plots.pdf", arrange_ggsurvplots(median_plots), width = 10, height = 8)
# To save all quantile plots to a PDF:
# ggsave("all_quantile_plots.pdf", arrange_ggsurvplots(quantile_plots), width = 10, height = 8)
