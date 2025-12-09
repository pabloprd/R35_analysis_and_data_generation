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
Sex_DEG <- read_xlsx("R35_Nanostring_Data_Analysis/Output/Differential Gene Expression Tables/Male_vs_Female_DGE.xlsx")

#Ancestry differentially expressed genes
Ancestry_DEG <- read_xlsx("R35_Nanostring_Data_Analysis/Output/Differential Gene Expression Tables/Black_vs_White_DGE.xlsx")

####### Volcano Plot Function UNDER CONSTRUCTION ##################
#All I need to do is plot adjusted p value by logFC. 
#And if it is marked as significant make it light up with a different color!
#I also need to set the markers as in 0.05 and 1.0


Volcano_plot_function <- function(differentially_expressed_table){
  #Creating the thresholds (p value and fold change)
  vline_threshold <- 1.0
  hline_threshold <- -log10(0.05)
  
  ggplot(data = differentially_expressed_table, 
         aes(x = logFC, y = -log10(adj.P.Val), col = significant, label = Gene))+
    geom_vline(xintercept = vline_threshold, col = "black", linetype = 'dashed', linewidth = 0.5) +
    geom_vline(xintercept = -vline_threshold, col = "black", linetype = 'dashed', linewidth = 0.5) +
    geom_hline(yintercept = abs(hline_threshold), col = "black", linetype = 'dashed', linewidth = 0.5) +
    geom_point(size = 2) +
    scale_colour_manual(
      values = c("grey", "red"),
    )+
    labs(x = "Log2-Fold Change", y = "Log10(p-value)", 
         col = "Regulation Level" )+
    theme(
      axis.title = element_text(size =10.0, face = 'bold'),
      axis.text = element_text(size = 10.0, face = 'bold'),
      legend.title = element_text(size = 10.0, face = 'bold'),
      panel.background = element_rect(fill = "white"),             # White background inside plot panel
      plot.background = element_rect(fill = "white"),              # White background around plot
      panel.grid.major = element_line(color = "gray80"),           # Major grid lines gray (lighter)
      panel.grid.minor = element_line(color = "gray90"),
      legend.position = "none", # Minor grid lines even lighter gray
      plot.title = element_text(size = 15.0, hjust = 0.5, face = "bold"))
  
  
}


###### Volcano Plotting ###############
Volcano_plot_function(Sex_DEG)
Volcano_plot_function(Ancestry_DEG)


###### Export #############

#Remember to set working directory to the boxplot section
setwd("/Users/paulparodi/R35_v2/Nanostring/R35_Nanostring_Data_Visualization/Output/Volcano Plots/")

#Example of how to save as PNG!!
#p <- Single_gene_boxplot(Sex_OXER1, 'Sex','OXER1')   # Assign plot to variable
ggsave("OXER1_boxplot.png", plot = p, width = 6, height = 5, dpi = 300)  # Save as PNG



