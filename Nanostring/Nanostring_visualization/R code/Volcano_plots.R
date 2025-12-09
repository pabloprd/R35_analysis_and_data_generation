#Volcano Plots
#By Paul Parodi
#Last modified on October 10th 2025


######### GOAL ###################
#To make nice and modifiable volcano plots to see the spread of gene data

####### Library ##################
library("writexl")
library("readxl")
library('dplyr')
library('tidyverse')
library('tidyr')
library('ggplot2')
library('ggrepel')

###### Loading in the files #############

setwd("/Users/paulparodi/R35_v2/Nanostring/")

#Sex differenitally expressed genes
Sex_DEG <- read_xlsx("Nanostring_analysis/Output/Differential Gene Expression/Male_vs_Female_DGE_table.xlsx")

#Ancestry differentially expressed genes
Ancestry_DEG <- read_xlsx("Nanostring_analysis/Output/Differential Gene Expression/Black_vs_White_DGE_table.xlsx")

###### Regulation #############

regulation_function <- function(dataframe, pval_cutoff = 0.05, logfc_cutoff = 1.0){
  
  dataframe[['Regulation']] <- 'not significant'
  
  dataframe$Regulation[dataframe$adj.P.Val <= pval_cutoff & dataframe$logFC >  logfc_cutoff] <- "upregulated"
  
  dataframe$Regulation[dataframe$adj.P.Val <= pval_cutoff & dataframe$logFC < -logfc_cutoff] <- "downregulated"
  
  return(dataframe)
  
}

####### Volcano Plot Function UNDER CONSTRUCTION ##################
Volcano_plot_function <- function(differentially_expressed_table, pval_cut = 0.05, log_cut = 1.0){
  #Creating the thresholds (p value and fold change)
  vline_threshold <- log_cut
  hline_threshold <- -log10(pval_cut)
  
  genes_of_interest <- differentially_expressed_table %>% #filtering for the genes
    filter(Regulation != 'not significant')
  
  ggplot(data = differentially_expressed_table, 
         aes(x = logFC, y = -log10(adj.P.Val), col = Regulation))+
    geom_vline(xintercept = vline_threshold, col = "black", linetype = 'dashed', linewidth = 0.5) +
    geom_vline(xintercept = -vline_threshold, col = "black", linetype = 'dashed', linewidth = 0.5) +
    geom_hline(yintercept = abs(hline_threshold), col = "black", linetype = 'dashed', linewidth = 0.5) +
    geom_point(size = 2) +
    labs(x = "Log2-Fold Change", y = "Log10(p-value)", 
         color = "Regulation" )+
    geom_text_repel(
      data = genes_of_interest,
      aes(label = Name),
      max.overlaps = Inf) +
    scale_colour_manual(
      values = c("blue", "gray", "red"),
      drop = FALSE
    )+
    theme(
      axis.title = element_text(size =10.0, face = 'bold'),
      axis.text = element_text(size = 10.0, face = 'bold'),
      legend.title = element_text(size = 10.0, face = 'bold'),
      panel.background = element_rect(fill = "white"),             # White background inside plot panel
      plot.background = element_rect(fill = "white"),              # White background around plot
      legend.position = "right",# Minor grid lines even lighter gray
      plot.title = element_text(size = 15.0, hjust = 0.5, face = "bold"))
  
  
}





###### Regulation + Volcano Plotting ###############

Sex_DEG_v2 <- regulation_function(Sex_DEG)
Volcano_plot_function(Sex_DEG_v2)

Ancestry_DEG_v2 <- regulation_function(Ancestry_DEG)
Volcano_plot_function(Ancestry_DEG_v2)


###### Export #############

#Example of how to save as PNG!!
#p <- Single_gene_boxplot(Sex_OXER1, 'Sex','OXER1')   # Assign plot to variable
ggsave("OXER1_boxplot.png", plot = p, width = 6, height = 5, dpi = 300)  # Save as PNG



