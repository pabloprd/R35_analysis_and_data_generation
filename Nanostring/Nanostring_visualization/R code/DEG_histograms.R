#Histograms
#By Paul Parodi
#Last updated October 14th 2025


######## GOAL ##########
#To make some quality histograms based on differential expressed genes. 

####### LIBRARY ########
library("writexl")
library("readxl")
library('dplyr')
library('tidyverse')
library('tidyr')
library('ggplot2')
library('ggrepel')

####### Loading in the files ##########
setwd("/Users/paulparodi/R35_v2/Nanostring/")

#Sex differenitally expressed genes
Sex_DEG <- read_xlsx("Nanostring_analysis/Output/Differential Gene Expression/Male_vs_Female_DGE_table.xlsx")

#Ancestry differentially expressed genes
Ancestry_DEG <- read_xlsx("Nanostring_analysis/Output/Differential Gene Expression/Black_vs_White_DGE_table.xlsx")




####### Gene selection function #########

Genes_of_interest <- function(dataframe, col_oi, num = 5){
  organized_df <- dataframe %>% arrange(!!sym(col_oi))
  organized_df <- head(organized_df, num)
  output <- organized_df$Name
  return(output)

}

###### Histogram function ###########

#Give people the option of changing the y value. 
Histogram_function <- function(GOI, DEG_table, column, title = NULL){
  #Assuming the column with the genes is called 'Gene'
  DEG_table <- DEG_table %>% filter(Name %in% GOI)
  
  # Plot bar plot of logFC by Gene
  ggplot(DEG_table, aes(x = Name, y = !!sym(column))) +
    geom_col(fill = "steelblue") + 
    labs(title = title,
         x = "Gene", y = paste0(column)) +
    theme_minimal() +
    theme(axis.title = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1),
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))
}

##### Plotting ##############

#!!!!!! SEX !!!!!!!!!!!!
#Sex genes of interest
Sex_GOI <- Genes_of_interest(Sex_DEG, "adj.P.Val")

#Histogram
Histogram_function(Sex_GOI, Sex_DEG, "logFC", "Log Fold Change between Male and Female top 5 genes") 


#!!!! ANCESTRY !!!!!!!!!
#Ancestry genes of interest
Ancestry_GOI <- Genes_of_interest(Ancestry_DEG, "adj.P.Val")

#Histogram
#I can change this to whatever column name that I want!
Histogram_function(Ancestry_GOI, Ancestry_DEG, "logFC", "Log Fold Change between Black and White top 5 genes") 



##### Export ############
#Change to the output folder where I need it!
#Example of how to save as PNG!!
p <- Histogram_function(Sex_GOI, Sex_DEG, "logFC") 

ggsave("Sex_DEG_logFC_histogram", plot = p, width = 6, height = 5, dpi = 300)  # Save as PNG

