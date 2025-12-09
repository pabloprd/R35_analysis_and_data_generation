#Histograms
#By Paul Parodi
#Last updated October 14th 2025


######## GOAL ##########
#To make some quality histograms based on differential expressed genes results
#such as logfc, p value etc. 

####### LIBRARY ########
library("writexl")
library("readxl")
library('dplyr')
library('tidyverse')
library('tidyr')
library('ggplot2')
library('ggrepel')

####### Loading in the files ##########
setwd("/Users/paulparodi/R35_v2/Mass Spec//")

#Sex DGE
Sex_DEG <- read_excel("Mass_Spec_Data_Analysis/Output/Differential Gene Expression Tables/Male_vs_Female_DGE.xlsx")

#Ancestry DGE
Ancestry_DEG <- read_excel("Mass_Spec_Data_Analysis/Output/Differential Gene Expression Tables/EA_vs_AA_DGE.xlsx")


####### Gene selection function #########

Genes_of_interest <- function(dataframe, num = 5, sig = FALSE){
  # Put TRUE rows first
  #Make it so that true is 1 and false is 0 
  if (sig == FALSE){
    organized_df <- dataframe %>% arrange(`adj.P.Val`)
    organized_df <- head(organized_df, num)
    output <- organized_df$GenSymbol
    
  }
  else{
    organized_df <- dataframe[order(-dataframe$significant), ]
    organized_df <- head(organized_df, num)
    output <- organized_df$GenSymbol
    
    
  }
  return(output)
}

###### Histogram function ###########

#Give people the option of changing the y value. 
Histogram_function <- function(GOI, DEG_table, column, title = NULL){
  #Assuming the column with the genes is called 'Gene'
  DEG_table <- DEG_table %>% filter(GenSymbol %in% GOI)
  
  # Plot bar plot of logFC by Gene
  ggplot(DEG_table, aes(x = GenSymbol, y = !!sym(column))) +
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
Sex_GOI <- Genes_of_interest(Sex_DEG)

#Histogram
Histogram_function(Sex_GOI, Sex_DEG, "logFC", "Log Fold Change between Male and Female") 


#!!!! ANCESTRY !!!!!!!!!
#Ancestry genes of interest
Ancestry_GOI <- Genes_of_interest(Ancestry_DEG)

#Histogram
Histogram_function(Ancestry_GOI, Ancestry_DEG, "logFC", "Log Fold Change between Black and White") 



##### Export ############

#Example of how to save as PNG!!
p <- Histogram_function(Sex_GOI, Sex_DEG, "logFC") 

ggsave("Sex_DEG_logFC_histogram", plot = p, width = 6, height = 5, dpi = 300)  # Save as PNG

