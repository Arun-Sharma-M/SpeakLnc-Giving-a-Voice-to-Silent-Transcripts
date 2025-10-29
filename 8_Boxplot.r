library(ggplot2)
library(tidyverse)
library(ggpubr) 
library(glue) # Added for string formatting

CHOL_SlimRNA_TPM_in <- read.csv('TPM_CHOL_lncRNA_mRNA_Tumour_Deseq.csv', header = TRUE, sep = ",")

rownames(CHOL_SlimRNA_TPM_in) <- CHOL_SlimRNA_TPM_in$TranscriptSymbol
CHOL_SlimRNA_TPM_in[1:5,1:5]

CHOL_SlimRNA_TPM_in <- CHOL_SlimRNA_TPM_in[ -c(1) ] # Removes entire column be careful
CHOL_SlimRNA_TPM_in[1:5,1:5]

metaData <- read.csv('TCGA_CHOL_Metadata.csv', header = TRUE, sep = ",")
metaData

SlimRNA <- read.csv('Shortlisted_Transcripts.csv', header = FALSE, sep = ",")
SlimRNA<- SlimRNA[[1]]

################################# SlimRNA ##################################

CHOL_SlimRNA_TPM <- CHOL_SlimRNA_TPM_in[rownames(CHOL_SlimRNA_TPM_in) %in% SlimRNA,]      #


#add_column_name_for_first_row
CHOL_SlimRNA_TPM <- CHOL_SlimRNA_TPM %>% 
  data.frame() %>%
  rownames_to_column(var="TranscriptSymbol") %>% 
  as_tibble()

# Rename_transcripts_in order
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "ADH4"] <- "A_ADH4"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "lnc-ADH4-3:2"] <- "B_lnc-ADH4-3:2"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "lnc-ADH4-3:3"] <- "C_lnc-ADH4-3:3"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "ALDH1L1"] <- "D_ALDH1L1"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "ALDH1L1-AS2:1"] <- "E_ALDH1L1-AS2:1"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "ALDH1L1-AS2:3"] <- "F_ALDH1L1-AS2:3"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "CYP2A6"] <- "G_CYP2A6"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "lnc-CYP2A6-2:1"] <- "H_lnc-CYP2A6-2:1"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "CYP4A22"] <- "I_CYP4A22"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "CYP4A22-AS1:3"] <- "J_CYP4A22-AS1:3"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "LINC01093"] <- "K_LINC01093"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "LINC01093:1"] <- "L_LINC01093:1"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "LINC02362"] <- "M_LINC02362"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "LINC02362:11"] <- "N_LINC02362:11"
CHOL_SlimRNA_TPM$TranscriptSymbol[CHOL_SlimRNA_TPM$TranscriptSymbol == "LINC02362:35"] <- "O_LINC02362:35"
# Gathering the columns to have normalized counts to a single column
gathered_SlimRNA_col <- CHOL_SlimRNA_TPM %>%
  gather(colnames(CHOL_SlimRNA_TPM)[2:45], key = "samplename", value = "TPM")

gathered_SlimRNA <- inner_join(metaData, gathered_SlimRNA_col)
View(gathered_SlimRNA)



# --- Load p-values from CSV and prepare for plotting ---
# Assuming your CSV has 'TranscriptSymbol', 'p.value', and 'p.adjusted' columns
pvalues_df <- read.csv('TPM_CHOL_lncRNA_mRNA_Tumour_onlypvalue.csv', header = TRUE, sep = ",")
rownames(pvalues_df) <- pvalues_df$TranscriptSymbol
pvalues_df <- pvalues_df[ -c(1) ] 

pvalues_df_SlimRNA <- pvalues_df[rownames(pvalues_df) %in% SlimRNA,]      
pvalues_df_SlimRNA <- pvalues_df_SlimRNA %>% 
  data.frame() %>%
  rownames_to_column(var="TranscriptSymbol") %>% 
  as_tibble()
# Add a significance column for plotting, if desired
pvalues_df_SlimRNA <- pvalues_df_SlimRNA %>%
  mutate(p.value.signif = case_when(
    pvalue < 0.001 ~ "***",
    pvalue < 0.01 ~ "**",
    pvalue < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# Add a y.position column to place the p-values.
max_tpm_per_transcript <- gathered_SlimRNA %>%
  group_by(TranscriptSymbol) %>%
  summarise(max_tpm = max(TPM, na.rm = TRUE))

pvalues_df_SlimRNA <- inner_join(pvalues_df_SlimRNA, max_tpm_per_transcript, by = "TranscriptSymbol") %>%
  mutate(y.position = log10(max_tpm) + 0.5)

write.csv(pvalues_df_SlimRNA,"Output_TCGA_boxplot\\Pv_SlimRNA_TCGA.csv", row.names = TRUE)

# --- Modified Plotting Code ---

ggplot(gathered_SlimRNA) + 
  geom_boxplot(aes(x = TranscriptSymbol, y = TPM, fill=SampleType), width=0.5, outlier.shape = NA) +
  scale_y_log10() +
  xlab("TranscriptSymbol") +
  ylab("log10(TPM)") +
  scale_fill_manual(values = c("#00FF66", "#FF3300")) +
  ggtitle("TCGA_CHOL_Shortlisted_Transcripts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

#  ,	,	,												








