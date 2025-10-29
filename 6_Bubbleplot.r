#Data visualisation and network analysis
#The identified gene list was submitted to GeneSetDB "Investigate Human Gene Sets" function with default parameters, and the output was saved as a TSV file (download_GeneSetDB_Reactome_gene_shortlisted.txt) 
# Bubble plot of Reactome Gene Sets
# Load libraries
library(ggplot2)   # For plotting
library(dplyr)     # For data manipulation
library(ggrepel)   # For smart text labels on plots
# 1. Load and preprocess data
# Read the shortlisted gene set data from Reactome
data <- read.csv('download_GeneSetDB_Reactome_gene_shortlisted.txt', header = TRUE, sep = "\t")
# Order by FDR (ascending)
data <- data[order(data$FDR),]
# Add ranking columns
data$Rank_in_order <- seq(1, nrow(data))         
# Scaled rank for plotting
data$Rank <- data$Rank_in_order * 2              
# Keep only relevant columns and rename for clarity
data <- data[, c("Rank_in_order", "geneSet.name", "Genes_in_Overlap_.k.", "FDR", "Rank")]
colnames(data) <- c("Rank_in_order", "process", "genes", "FDR", "Rank")
# Keep only top 10 gene sets for plotting
data <- data[1:10,]
# Save preprocessed data for reference
write.csv(data, "5.GeneSetDB_Reactome_gene_shortlisted.csv", row.names = TRUE)
# 2. Create bubble plot
# Notes:
# - `geom_point` creates the bubbles.
# - `size` aesthetic controls bubble size (genes in overlap).
# - `fill` can be used to color-code each bubble.
# - `geom_hline` adds a significance threshold line (optional).
# - `ggrepel` can be used to avoid overlapping labels if desired.
ggplot(data, 
       aes(x = Rank, 
           # Convert FDR to -log10 scale for better visualization
           y = -log10(FDR),    
           # Bubble size proportional to number of genes in overlap      
           size = genes,       
           # Label each bubble with process name      
           label = process,          
           # Bubble fill aesthetic
           fill = paste0(Rank, ":", process))) +  
  # Bubble outline color
  geom_point(colour = "#2E2E2E", shape = 21) +  
  # Max bubble size
  scale_size_area(max_size = 10) +              
  scale_x_continuous(name = "Rank", limits = c(0, 12.5)) +
  scale_y_continuous(name = "-log10(FDR)") +
  labs(fill = "Rank:Reactome") +
  # Horizontal threshold line
  geom_hline(yintercept = 1.3) +               
  theme(axis.text.x = element_text(family = "Arial", face = "plain", color = "#2E2E2E", size = 10),
        axis.text.y = element_text(family = "Arial", face = "plain", color = "#2E2E2E", size = 10),
        text = element_text(size = 20))
