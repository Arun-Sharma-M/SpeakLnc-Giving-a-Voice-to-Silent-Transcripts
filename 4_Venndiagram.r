#install.packages("ggVennDiagram")

# install.packages("devtools")
#devtools::install_github("gaospecial\\ggVennDiagram")

library(ggVennDiagram)
library(ggplot2)

###############Padj<0.05 and L2FC >2 #############################################
highexp1057_CHOL_SRA_lncRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_lncRNA_CHOL\\Reads_to_TPM_filtering\\\\Input_SRA_CHOL_lncRNA\\Results_SRA_CHOL_lncRNA\\8.pv_l2fc_upr_CHOL_SRA_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp1057_CHOL_SRA_lncRNA)[1]  <- "lncRNASymbol"    # change column name for x column
highexp1057_CHOL_SRA_lncRNA_list <- highexp1057_CHOL_SRA_lncRNA[,-(2:35)]
#rownames(highexp1057_CHOL_SRA_lncRNA_list) <- highexp1057_CHOL_SRA_lncRNA[,1]

highexp123_CHOL_TCGA_lncRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_lncRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_lncRNA\\Results_TCGA_CHOL_lncRNA\\8.pv_l2fc_upr_CHOL_TCGA_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp123_CHOL_TCGA_lncRNA)[1]  <- "lncRNASymbol"    # change column name for x column
highexp123_CHOL_TCGA_lncRNA_list <- highexp123_CHOL_TCGA_lncRNA[,-(2:45)]
#rownames(highexp123_CHOL_TCGA_lncRNA_list) <- highexp123_CHOL_TCGA_lncRNA[,1]

highexp118_CHOL_SRA_mRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_SRA_CHOL_mRNA\\Results_SRA_CHOL_mRNA\\8.pv_l2fc_upr_CHOL_SRA_T_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp118_CHOL_SRA_mRNA)[1]  <- "geneSymbol"    # change column name for x column
highexp118_CHOL_SRA_mRNA_list <- highexp118_CHOL_SRA_mRNA[,-(2:35)]
#rownames(highexp118_CHOL_SRA_mRNA_list) <- highexp118_CHOL_SRA_mRNA[,1]

highexp75_CHOL_TCGA_mRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_mRNA\\Results_TCGA_CHOL_mRNA\\8.pv_l2fc_upr_CHOL_TCGA_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp75_CHOL_TCGA_mRNA)[1]  <- "geneSymbol"    # change column name for x column
highexp75_CHOL_TCGA_mRNA_list <- highexp75_CHOL_TCGA_mRNA[,-(2:46)]
#rownames(highexp75_CHOL_TCGA_mRNA_list) <- highexp75_CHOL_TCGA_mRNA[,1]

x_lncRNA <- list(SRA_lncRNA=sample(highexp1057_CHOL_SRA_lncRNA_list),
          TCGA_lncRNA=sample(highexp123_CHOL_TCGA_lncRNA_list))
ggVennDiagram(x_lncRNA, label_alpha=0) + scale_fill_gradient(low="green",high = "yellow", labs(title = "Four sets Venn Diagram"))      # VennDiagram 1
              
x_mRNA <- list(SRA_mRNA=sample(highexp118_CHOL_SRA_mRNA_list),
          TCGA_mRNA=sample(highexp75_CHOL_TCGA_mRNA_list))
ggVennDiagram(x_mRNA, label_alpha=0) + scale_fill_gradient(low="blue",high = "red")

Merge_SRA_TCGA_lncRNA_Pvalue_L2FC <- intersect(highexp1057_CHOL_SRA_lncRNA_list, highexp123_CHOL_TCGA_lncRNA_list)
Merge_SRA_TCGA_lncRNA_Pvalue_L2FC<-data.frame(Merge_SRA_TCGA_lncRNA_Pvalue_L2FC)
write.csv(Merge_SRA_TCGA_lncRNA_Pvalue_L2FC,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\1.Merge_SRA_TCGA_lncRNA_Pvalue_L2FC.csv")

Merge_SRA_TCGA_mRNA_Pvalue_L2FC <- intersect(highexp118_CHOL_SRA_mRNA_list, highexp75_CHOL_TCGA_mRNA_list)
Merge_SRA_TCGA_mRNA_Pvalue_L2FC<-data.frame(Merge_SRA_TCGA_mRNA_Pvalue_L2FC)
write.csv(Merge_SRA_TCGA_mRNA_Pvalue_L2FC,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\2.Merge_SRA_TCGA_mRNA_Pvalue_L2FC.csv")

#Strip_lncRNA_annotation
highexp1057_CHOL_SRA_lncRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_lncRNA_CHOL\\Reads_to_TPM_filtering\\Input_SRA_CHOL_lncRNA\\Results_SRA_CHOL_lncRNA\\8.pv_l2fc_upr_CHOL_SRA_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp1057_CHOL_SRA_lncRNA_strip)[1]  <- "lncRNASymbol"    # change column name for x column
highexp1057_CHOL_SRA_lncRNA_strip[] <- lapply(highexp1057_CHOL_SRA_lncRNA_strip, gsub, pattern='lnc-', replacement='')
highexp1057_CHOL_SRA_lncRNA_strip[] <- lapply(highexp1057_CHOL_SRA_lncRNA_strip, gsub, pattern='-.*', replacement='')
highexp1057_CHOL_SRA_lncRNA_strip[] <- lapply(highexp1057_CHOL_SRA_lncRNA_strip, gsub, pattern=':.*', replacement='')
highexp1057_CHOL_SRA_lncRNA_striplist <- highexp1057_CHOL_SRA_lncRNA_strip[,-(2:35)]

highexp123_CHOL_TCGA_lncRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_lncRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_lncRNA\\Results_TCGA_CHOL_lncRNA\\8.pv_l2fc_upr_CHOL_TCGA_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp123_CHOL_TCGA_lncRNA_strip)[1]  <- "lncRNASymbol"    # change column name for x column
highexp123_CHOL_TCGA_lncRNA_strip[] <- lapply(highexp123_CHOL_TCGA_lncRNA_strip, gsub, pattern='lnc-', replacement='')
highexp123_CHOL_TCGA_lncRNA_strip[] <- lapply(highexp123_CHOL_TCGA_lncRNA_strip, gsub, pattern='-.*', replacement='')
highexp123_CHOL_TCGA_lncRNA_strip[] <- lapply(highexp123_CHOL_TCGA_lncRNA_strip, gsub, pattern=':.*', replacement='')
highexp123_CHOL_TCGA_lncRNA_striplist <- highexp123_CHOL_TCGA_lncRNA_strip[,-(2:45)]

highexp118_CHOL_SRA_mRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_SRA_CHOL_mRNA\\Results_SRA_CHOL_mRNA\\8.pv_l2fc_upr_CHOL_SRA_T_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp118_CHOL_SRA_mRNA_strip)[1]  <- "geneSymbol"    # change column name for x column
highexp118_CHOL_SRA_mRNA_strip[] <- lapply(highexp118_CHOL_SRA_mRNA_strip, gsub, pattern='-.*', replacement='')
highexp118_CHOL_SRA_mRNA_strip[] <- lapply(highexp118_CHOL_SRA_mRNA_strip, gsub, pattern=':.*', replacement='')
highexp118_CHOL_SRA_mRNA_striplist <- highexp118_CHOL_SRA_mRNA_strip[,-(2:35)]

highexp75_CHOL_TCGA_mRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_mRNA\\Results_TCGA_CHOL_mRNA\\8.pv_l2fc_upr_CHOL_TCGA_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp75_CHOL_TCGA_mRNA_strip)[1]  <- "geneSymbol"    # change column name for x column
highexp75_CHOL_TCGA_mRNA_strip[] <- lapply(highexp75_CHOL_TCGA_mRNA_strip, gsub, pattern='-.*', replacement='')
highexp75_CHOL_TCGA_mRNA_strip[] <- lapply(highexp75_CHOL_TCGA_mRNA_strip, gsub, pattern=':.*', replacement='')
highexp75_CHOL_TCGA_mRNA_striplist <- highexp75_CHOL_TCGA_mRNA_strip[,-(2:46)]

lnc_mRNA_merge_trans <- list(SRA_lncRNA=sample(highexp1057_CHOL_SRA_lncRNA_striplist),
                 TCGA_lncRNA=sample(highexp123_CHOL_TCGA_lncRNA_striplist),
                 SRA_mRNA=sample(highexp118_CHOL_SRA_mRNA_striplist),
                 TCGA_mRNA=sample(highexp75_CHOL_TCGA_mRNA_striplist))
ggVennDiagram(lnc_mRNA_merge_trans, label_alpha=0, edge_lty = "dashed", edge_size ="10") + scale_fill_gradient(low="blue",high = "red")

lnc_mRNA_merge1 <- Venn(lnc_mRNA_merge_trans)
lnc_mRNA_merge <- process_data(lnc_mRNA_merge1)
ggplot() +
  #1.Region count layer
  geom_sf(aes(fill=id), data = venn_region(lnc_mRNA_merge),show.legend = FALSE) +
  #2.Set edge layer, size and colour
  geom_sf(size = 20, lty = "dashed", aes(color = "black"), data = venn_setedge(lnc_mRNA_merge), show.legend = FALSE) +
  #3.set region label layer
  geom_sf_label(aes(label=count), alpha=0, fontface = "bold", data = venn_region(lnc_mRNA_merge), label.size = NA, size = 3) +
  #4.Set title label layer (set_Title)
  geom_sf_text(aes(label = name), colour= "black", data = venn_setlabel(lnc_mRNA_merge), size = 4, fontface = "bold", position = position_dodge(width = 0.2)) +
  theme_void()
 #Adjust edge size
geom_sf()

library(dplyr)
highexp1057_CHOL_SRA_lncRNA_striplist <- data.frame(highexp1057_CHOL_SRA_lncRNA_striplist)
write.csv(highexp1057_CHOL_SRA_lncRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\1.L2FC2_CHOL_SRA_lncRNA_striplist.csv")
highexp123_CHOL_TCGA_lncRNA_striplist <- data.frame(highexp123_CHOL_TCGA_lncRNA_striplist)
write.csv(highexp123_CHOL_TCGA_lncRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\2.L2FC2_CHOL_TCGA_lncRNA_striplist.csv")
highexp118_CHOL_SRA_mRNA_striplist <- data.frame(highexp118_CHOL_SRA_mRNA_striplist)
write.csv(highexp118_CHOL_SRA_mRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\3.L2FC2_CHOL_SRA_mRNA_striplist.csv")
highexp75_CHOL_TCGA_mRNA_striplist <- data.frame(highexp75_CHOL_TCGA_mRNA_striplist)
write.csv(highexp75_CHOL_TCGA_mRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\4.L2FC2_CHOL_TCGA_mRNA_striplist.csv")

###############Padj<0.05, TPM>1 and L2FC >2#############################################

highexp1057_CHOL_SRA_lncRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_lncRNA_CHOL\\Reads_to_TPM_filtering\\Input_SRA_CHOL_lncRNA\\Results_SRA_CHOL_lncRNA\\9.pv_l2fc_TPM_upr_CHOL_SRA_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp1057_CHOL_SRA_lncRNA)[1]  <- "lncRNASymbol"    # change column name for x column
highexp1057_CHOL_SRA_lncRNA_list <- highexp1057_CHOL_SRA_lncRNA[,-(2:35)]
#rownames(highexp1057_CHOL_SRA_lncRNA_list) <- highexp1057_CHOL_SRA_lncRNA[,1]

highexp123_CHOL_TCGA_lncRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_lncRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_lncRNA\\Results_TCGA_CHOL_lncRNA\\9.pv_l2fc_TPM_TPM_CHOL_T_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp123_CHOL_TCGA_lncRNA)[1]  <- "lncRNASymbol"    # change column name for x column
highexp123_CHOL_TCGA_lncRNA_list <- highexp123_CHOL_TCGA_lncRNA[,-(2:45)]
#rownames(highexp123_CHOL_TCGA_lncRNA_list) <- highexp123_CHOL_TCGA_lncRNA[,1]

highexp118_CHOL_SRA_mRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_SRA_CHOL_mRNA\\Results_SRA_CHOL_mRNA\\9.pv_L2FC2_TPM1_CHOL_T_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp118_CHOL_SRA_mRNA)[1]  <- "geneSymbol"    # change column name for x column
highexp118_CHOL_SRA_mRNA_list <- highexp118_CHOL_SRA_mRNA[,-(2:35)]
#rownames(highexp118_CHOL_SRA_mRNA_list) <- highexp118_CHOL_SRA_mRNA[,1]

highexp75_CHOL_TCGA_mRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_mRNA\\Results_TCGA_CHOL_mRNA\\9.pv_l2fc_TPM_CHOL_T_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp75_CHOL_TCGA_mRNA)[1]  <- "geneSymbol"    # change column name for x column
highexp75_CHOL_TCGA_mRNA_list <- highexp75_CHOL_TCGA_mRNA[,-(2:46)]
#rownames(highexp75_CHOL_TCGA_mRNA_list) <- highexp75_CHOL_TCGA_mRNA[,1]

x_lncRNA <- list(SRA_lncRNA=sample(highexp1057_CHOL_SRA_lncRNA_list),
                 TCGA_lncRNA=sample(highexp123_CHOL_TCGA_lncRNA_list))
ggVennDiagram(x_lncRNA, label_alpha=0) + scale_fill_gradient(low="green",high = "yellow")

x_mRNA <- list(SRA_mRNA=sample(highexp118_CHOL_SRA_mRNA_list),
               TCGA_mRNA=sample(highexp75_CHOL_TCGA_mRNA_list))
ggVennDiagram(x_mRNA, label_alpha=0) + scale_fill_gradient(low="blue",high = "red")

Merge_SRA_TCGA_lncRNA_Pvalue_L2FC_TPM <- intersect(highexp1057_CHOL_SRA_lncRNA_list, highexp123_CHOL_TCGA_lncRNA_list)
Merge_SRA_TCGA_lncRNA_Pvalue_L2FC_TPM<-data.frame(Merge_SRA_TCGA_lncRNA_Pvalue_L2FC_TPM)
write.csv(Merge_SRA_TCGA_lncRNA_Pvalue_L2FC_TPM,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\3.Merge_SRA_TCGA_lncRNA_Pvalue_L2FC_TPM.csv")

Merge_SRA_TCGA_mRNA_Pvalue_L2FC_TPM <- intersect(highexp118_CHOL_SRA_mRNA_list, highexp75_CHOL_TCGA_mRNA_list)
Merge_SRA_TCGA_mRNA_Pvalue_L2FC_TPM<-data.frame(Merge_SRA_TCGA_mRNA_Pvalue_L2FC_TPM)
write.csv(Merge_SRA_TCGA_mRNA_Pvalue_L2FC_TPM,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\4.Merge_SRA_TCGA_mRNA_Pvalue_L2FC_TPM.csv")

#Strip_lncRNA_annotation
highexp1057_CHOL_SRA_lncRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_lncRNA_CHOL\\Reads_to_TPM_filtering\\Input_SRA_CHOL_lncRNA\\Results_SRA_CHOL_lncRNA\\9.pv_l2fc_TPM_upr_CHOL_SRA_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp1057_CHOL_SRA_lncRNA_strip)[1]  <- "lncRNASymbol"    # change column name for x column
highexp1057_CHOL_SRA_lncRNA_strip[] <- lapply(highexp1057_CHOL_SRA_lncRNA_strip, gsub, pattern='lnc-', replacement='')
highexp1057_CHOL_SRA_lncRNA_strip[] <- lapply(highexp1057_CHOL_SRA_lncRNA_strip, gsub, pattern='-.*', replacement='')
highexp1057_CHOL_SRA_lncRNA_strip[] <- lapply(highexp1057_CHOL_SRA_lncRNA_strip, gsub, pattern=':.*', replacement='')
highexp1057_CHOL_SRA_lncRNA_striplist <- highexp1057_CHOL_SRA_lncRNA_strip[,-(2:35)]

highexp123_CHOL_TCGA_lncRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_lncRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_lncRNA\\Results_TCGA_CHOL_lncRNA\\9.pv_l2fc_TPM_TPM_CHOL_T_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp123_CHOL_TCGA_lncRNA_strip)[1]  <- "lncRNASymbol"    # change column name for x column
highexp123_CHOL_TCGA_lncRNA_strip[] <- lapply(highexp123_CHOL_TCGA_lncRNA_strip, gsub, pattern='lnc-', replacement='')
highexp123_CHOL_TCGA_lncRNA_strip[] <- lapply(highexp123_CHOL_TCGA_lncRNA_strip, gsub, pattern='-.*', replacement='')
highexp123_CHOL_TCGA_lncRNA_strip[] <- lapply(highexp123_CHOL_TCGA_lncRNA_strip, gsub, pattern=':.*', replacement='')
highexp123_CHOL_TCGA_lncRNA_striplist <- highexp123_CHOL_TCGA_lncRNA_strip[,-(2:45)]

highexp118_CHOL_SRA_mRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_SRA_CHOL_mRNA\\Results_SRA_CHOL_mRNA\\9.pv_L2FC2_TPM1_CHOL_T_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp118_CHOL_SRA_mRNA_strip)[1]  <- "geneSymbol"    # change column name for x column
highexp118_CHOL_SRA_mRNA_strip[] <- lapply(highexp118_CHOL_SRA_mRNA_strip, gsub, pattern='-.*', replacement='')
highexp118_CHOL_SRA_mRNA_strip[] <- lapply(highexp118_CHOL_SRA_mRNA_strip, gsub, pattern=':.*', replacement='')
highexp118_CHOL_SRA_mRNA_striplist <- highexp118_CHOL_SRA_mRNA_strip[,-(2:35)]

highexp75_CHOL_TCGA_mRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_mRNA\\Results_TCGA_CHOL_mRNA\\9.pv_l2fc_TPM_CHOL_T_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp75_CHOL_TCGA_mRNA_strip)[1]  <- "geneSymbol"    # change column name for x column
highexp75_CHOL_TCGA_mRNA_strip[] <- lapply(highexp75_CHOL_TCGA_mRNA_strip, gsub, pattern='-.*', replacement='')
highexp75_CHOL_TCGA_mRNA_strip[] <- lapply(highexp75_CHOL_TCGA_mRNA_strip, gsub, pattern=':.*', replacement='')
highexp75_CHOL_TCGA_mRNA_striplist <- highexp75_CHOL_TCGA_mRNA_strip[,-(2:46)]

lnc_mRNA_merge_trans <- list(SRA_lncRNA=sample(highexp1057_CHOL_SRA_lncRNA_striplist),
                             TCGA_lncRNA=sample(highexp123_CHOL_TCGA_lncRNA_striplist),
                             SRA_mRNA=sample(highexp118_CHOL_SRA_mRNA_striplist),
                             TCGA_mRNA=sample(highexp75_CHOL_TCGA_mRNA_striplist))
ggVennDiagram(lnc_mRNA_merge_trans, label_alpha=0, edge_lty = "dashed", edge_size ="10") + scale_fill_gradient(low="blue",high = "red")

lnc_mRNA_merge1 <- Venn(lnc_mRNA_merge_trans)
lnc_mRNA_merge <- process_data(lnc_mRNA_merge1)
ggplot() +
  #1.Region count layer
  geom_sf(aes(fill=id), data = venn_region(lnc_mRNA_merge),show.legend = FALSE) +
  #2.Set edge layer, size and colour
  geom_sf(size = 20, lty = "dashed", aes(color = "black"), data = venn_setedge(lnc_mRNA_merge), show.legend = FALSE) +
  #3.set region label layer
  geom_sf_label(aes(label=count), alpha=0, fontface = "bold", data = venn_region(lnc_mRNA_merge), label.size = NA, size = 3) +
  #4.Set title label layer (set_Title)
  geom_sf_text(aes(label = name), colour= "black", data = venn_setlabel(lnc_mRNA_merge), size = 4, fontface = "bold", position = position_dodge(width = 0.2)) +
  theme_void()
#Adjust edge size
geom_sf()

library(dplyr)
highexp1057_CHOL_SRA_lncRNA_striplist <- data.frame(highexp1057_CHOL_SRA_lncRNA_striplist)
write.csv(highexp1057_CHOL_SRA_lncRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\1.L2FC2_TPM1_CHOL_SRA_lncRNA_striplist.csv")
highexp123_CHOL_TCGA_lncRNA_striplist <- data.frame(highexp123_CHOL_TCGA_lncRNA_striplist)
write.csv(highexp123_CHOL_TCGA_lncRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\2.L2FC2_TPM1_CHOL_TCGA_lncRNA_striplist.csv")
highexp118_CHOL_SRA_mRNA_striplist <- data.frame(highexp118_CHOL_SRA_mRNA_striplist)
write.csv(highexp118_CHOL_SRA_mRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\3.L2FC2_TPM1_CHOL_SRA_mRNA_striplist.csv")
highexp75_CHOL_TCGA_mRNA_striplist <- data.frame(highexp75_CHOL_TCGA_mRNA_striplist)
write.csv(highexp75_CHOL_TCGA_mRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\4.L2FC2_TPM1_CHOL_TCGA_mRNA_striplist.csv")

#######################################################Padj<0.05 and L2FC >2 TPM>1 COntrol<0.1################################################################

highexp1057_CHOL_SRA_lncRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_lncRNA_CHOL\\Reads_to_TPM_filtering\\Input_SRA_CHOL_lncRNA\\Results_SRA_CHOL_lncRNA\\10.pv_l2fc_CnT_TPM_upr_CHOL_SRA_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp1057_CHOL_SRA_lncRNA)[1]  <- "lncRNASymbol"    # change column name for x column
highexp1057_CHOL_SRA_lncRNA_list <- highexp1057_CHOL_SRA_lncRNA[,-(2:35)]
#rownames(highexp1057_CHOL_SRA_lncRNA_list) <- highexp1057_CHOL_SRA_lncRNA[,1]

highexp123_CHOL_TCGA_lncRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_lncRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_lncRNA\\Results_TCGA_CHOL_lncRNA\\10.pv_l2fc_CnT_TPM_upr_CHOL_TCGA_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp123_CHOL_TCGA_lncRNA)[1]  <- "lncRNASymbol"    # change column name for x column
highexp123_CHOL_TCGA_lncRNA_list <- highexp123_CHOL_TCGA_lncRNA[,-(2:45)]
#rownames(highexp123_CHOL_TCGA_lncRNA_list) <- highexp123_CHOL_TCGA_lncRNA[,1]

highexp118_CHOL_SRA_mRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_SRA_CHOL_mRNA\\Results_SRA_CHOL_mRNA\\10.pv_l2fc_tpmc1_CHOL_T_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp118_CHOL_SRA_mRNA)[1]  <- "geneSymbol"    # change column name for x column
highexp118_CHOL_SRA_mRNA_list <- highexp118_CHOL_SRA_mRNA[,-(2:35)]
#rownames(highexp118_CHOL_SRA_mRNA_list) <- highexp118_CHOL_SRA_mRNA[,1]

highexp75_CHOL_TCGA_mRNA<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_mRNA\\Results_TCGA_CHOL_mRNA\\10.pv_l2fc_CnT_TPM_upr_CHOL_TCGA_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp75_CHOL_TCGA_mRNA)[1]  <- "geneSymbol"    # change column name for x column
highexp75_CHOL_TCGA_mRNA_list <- highexp75_CHOL_TCGA_mRNA[,-(2:46)]
#rownames(highexp75_CHOL_TCGA_mRNA_list) <- highexp75_CHOL_TCGA_mRNA[,1]

x_lncRNA <- list(SRA_lncRNA=sample(highexp1057_CHOL_SRA_lncRNA_list),
                 TCGA_lncRNA=sample(highexp123_CHOL_TCGA_lncRNA_list))
ggVennDiagram(x_lncRNA, label_alpha=0) + scale_fill_gradient(low="green",high = "yellow")

x_mRNA <- list(SRA_mRNA=sample(highexp118_CHOL_SRA_mRNA_list),
               TCGA_mRNA=sample(highexp75_CHOL_TCGA_mRNA_list))
ggVennDiagram(x_mRNA, label_alpha=0) + scale_fill_gradient(low="blue",high = "red")

Merge_SRA_TCGA_lncRNA_Pvalue_L2FC_TPM_cTPM <- intersect(highexp1057_CHOL_SRA_lncRNA_list, highexp123_CHOL_TCGA_lncRNA_list)
Merge_SRA_TCGA_lncRNA_Pvalue_L2FC_TPM_cTPM<-data.frame(Merge_SRA_TCGA_lncRNA_Pvalue_L2FC_TPM_cTPM)
write.csv(Merge_SRA_TCGA_lncRNA_Pvalue_L2FC_TPM_cTPM,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\5.Merge_SRA_TCGA_lncRNA_Pvalue_L2FC_TPM_cTPM.csv")

Merge_SRA_TCGA_mRNA_Pvalue_L2FC_TPM_cTPM <- intersect(highexp118_CHOL_SRA_mRNA_list, highexp75_CHOL_TCGA_mRNA_list)
Merge_SRA_TCGA_mRNA_Pvalue_L2FC_TPM_cTPM<-data.frame(Merge_SRA_TCGA_mRNA_Pvalue_L2FC_TPM_cTPM)
write.csv(Merge_SRA_TCGA_mRNA_Pvalue_L2FC_TPM_cTPM,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\6.Merge_SRA_TCGA_mRNA_Pvalue_L2FC_TPM_cTPM.csv")

#Strip_lncRNA_annotation
highexp1057_CHOL_SRA_lncRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_lncRNA_CHOL\\Reads_to_TPM_filtering\\Input_SRA_CHOL_lncRNA\\Results_SRA_CHOL_lncRNA\\10.pv_l2fc_CnT_TPM_upr_CHOL_SRA_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp1057_CHOL_SRA_lncRNA_strip)[1]  <- "lncRNASymbol"    # change column name for x column
highexp1057_CHOL_SRA_lncRNA_strip[] <- lapply(highexp1057_CHOL_SRA_lncRNA_strip, gsub, pattern='lnc-', replacement='')
highexp1057_CHOL_SRA_lncRNA_strip[] <- lapply(highexp1057_CHOL_SRA_lncRNA_strip, gsub, pattern='-.*', replacement='')
highexp1057_CHOL_SRA_lncRNA_strip[] <- lapply(highexp1057_CHOL_SRA_lncRNA_strip, gsub, pattern=':.*', replacement='')
highexp1057_CHOL_SRA_lncRNA_striplist <- highexp1057_CHOL_SRA_lncRNA_strip[,-(2:35)]

highexp123_CHOL_TCGA_lncRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_lncRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_lncRNA\\Results_TCGA_CHOL_lncRNA\\10.pv_l2fc_CnT_TPM_upr_CHOL_TCGA_lncRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp123_CHOL_TCGA_lncRNA_strip)[1]  <- "lncRNASymbol"    # change column name for x column
highexp123_CHOL_TCGA_lncRNA_strip[] <- lapply(highexp123_CHOL_TCGA_lncRNA_strip, gsub, pattern='lnc-', replacement='')
highexp123_CHOL_TCGA_lncRNA_strip[] <- lapply(highexp123_CHOL_TCGA_lncRNA_strip, gsub, pattern='-.*', replacement='')
highexp123_CHOL_TCGA_lncRNA_strip[] <- lapply(highexp123_CHOL_TCGA_lncRNA_strip, gsub, pattern=':.*', replacement='')
highexp123_CHOL_TCGA_lncRNA_striplist <- highexp123_CHOL_TCGA_lncRNA_strip[,-(2:45)]

highexp118_CHOL_SRA_mRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\SRA_CHOL\\SRA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_SRA_CHOL_mRNA\\Results_SRA_CHOL_mRNA\\10.pv_l2fc_tpmc1_CHOL_T_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp118_CHOL_SRA_mRNA_strip)[1]  <- "geneSymbol"    # change column name for x column
highexp118_CHOL_SRA_mRNA_strip[] <- lapply(highexp118_CHOL_SRA_mRNA_strip, gsub, pattern='-.*', replacement='')
highexp118_CHOL_SRA_mRNA_strip[] <- lapply(highexp118_CHOL_SRA_mRNA_strip, gsub, pattern=':.*', replacement='')
highexp118_CHOL_SRA_mRNA_striplist <- highexp118_CHOL_SRA_mRNA_strip[,-(2:35)]

highexp75_CHOL_TCGA_mRNA_strip<-read.table("E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\TCGA_CHOL\\TCGA_mRNA_CHOL\\Reads_to_TPM_filtering\\Input_TCGA_CHOL_mRNA\\Results_TCGA_CHOL_mRNA\\10.pv_l2fc_CnT_TPM_upr_CHOL_TCGA_mRNA_TPM1ly.csv", header = TRUE, sep = ",")
colnames(highexp75_CHOL_TCGA_mRNA_strip)[1]  <- "geneSymbol"    # change column name for x column
highexp75_CHOL_TCGA_mRNA_strip[] <- lapply(highexp75_CHOL_TCGA_mRNA_strip, gsub, pattern='-.*', replacement='')
highexp75_CHOL_TCGA_mRNA_strip[] <- lapply(highexp75_CHOL_TCGA_mRNA_strip, gsub, pattern=':.*', replacement='')
highexp75_CHOL_TCGA_mRNA_striplist <- highexp75_CHOL_TCGA_mRNA_strip[,-(2:46)]

lnc_mRNA_merge_trans <- list(SRA_lncRNA=sample(highexp1057_CHOL_SRA_lncRNA_striplist),
                             TCGA_lncRNA=sample(highexp123_CHOL_TCGA_lncRNA_striplist),
                             SRA_mRNA=sample(highexp118_CHOL_SRA_mRNA_striplist),
                             TCGA_mRNA=sample(highexp75_CHOL_TCGA_mRNA_striplist))
ggVennDiagram(lnc_mRNA_merge_trans, label_alpha=0, edge_lty = "dashed", edge_size ="10") + scale_fill_gradient(low="blue",high = "red")

lnc_mRNA_merge1 <- Venn(lnc_mRNA_merge_trans)
lnc_mRNA_merge <- process_data(lnc_mRNA_merge1)
ggplot() +
  #1.Region count layer
  geom_sf(aes(fill=id), data = venn_region(lnc_mRNA_merge),show.legend = FALSE) +
  #2.Set edge layer, size and colour
  geom_sf(size = 20, lty = "dashed", aes(color = "black"), data = venn_setedge(lnc_mRNA_merge), show.legend = FALSE) +
  #3.set region label layer
  geom_sf_label(aes(label=count), alpha=0, fontface = "bold", data = venn_region(lnc_mRNA_merge), label.size = NA, size = 3) +
  #4.Set title label layer (set_Title)
  geom_sf_text(aes(label = name), colour= "black", data = venn_setlabel(lnc_mRNA_merge), size = 4, fontface = "bold", position = position_dodge(width = 0.2)) +
  theme_void()
#Adjust edge size
geom_sf()

library(dplyr)
highexp1057_CHOL_SRA_lncRNA_striplist <- data.frame(highexp1057_CHOL_SRA_lncRNA_striplist)
write.csv(highexp1057_CHOL_SRA_lncRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\1.L2FC2_TPM1_C0.1_CHOL_SRA_lncRNA_striplist.csv")
highexp123_CHOL_TCGA_lncRNA_striplist <- data.frame(highexp123_CHOL_TCGA_lncRNA_striplist)
write.csv(highexp123_CHOL_TCGA_lncRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\2.L2FC2_TPM1_C0.1_CHOL_TCGA_lncRNA_striplist.csv")
highexp118_CHOL_SRA_mRNA_striplist <- data.frame(highexp118_CHOL_SRA_mRNA_striplist)
write.csv(highexp118_CHOL_SRA_mRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\3.L2FC2_TPM1_C0.1_CHOL_SRA_mRNA_striplist.csv")
highexp75_CHOL_TCGA_mRNA_striplist <- data.frame(highexp75_CHOL_TCGA_mRNA_striplist)
write.csv(highexp75_CHOL_TCGA_mRNA_striplist,"E:\\Lymphomas\\CHOL\\Cholangiocarcinoma_R_Analysis\\Venn_Diagram_CHOL_dwnr\\4.L2FC2_TPM1_C0.1_CHOL_TCGA_mRNA_striplist.csv")
