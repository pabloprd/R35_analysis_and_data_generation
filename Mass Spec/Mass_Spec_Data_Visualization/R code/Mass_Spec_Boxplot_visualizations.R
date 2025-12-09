#Boxplots [Sex and Ancestry]
#Last updated October 9th 2025
#by Paul Parodi


#Quickest boxplot function for later

######### GOAL #######################
#I can plot nice boxplots from here. 
#All I have to do is modify the function to my liking

#So yeah goal of this file is to make actual nice boxplots.


######### Library ##################
library("writexl")
library("readxl")
library('dplyr')
library('tidyverse')
library('tidyr')
library('ggplot2')
library('ggrepel')

######### Loading in the Files ############
setwd("/Users/paulparodi/R35_v2/Mass Spec/")

normalized_counts <- read_excel("Mass_Spec_Data_Preprocessing/Output/All_confidence_protein_intensities_set1.xlsx")
metadata <- read_excel("Mass_Spec_Data_Preprocessing/Output/R35_Poochon_MS_metadata.xlsx")

########## Box plot Function #########

Single_gene_boxplot <- function(dataframe_long, comparison_group, gene){
  
  #Plotting stuff
  ggplot(dataframe_long, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(color = "black", size = 1, outlier.shape = NA)+
    geom_jitter(width = 0.3, color = "black", alpha = 0.7)+
    scale_fill_brewer(palette = "Set2") + #determining color of inner boxes #What are the different pallets?
    labs(
      title = paste0(comparison_group, " differences in ", gene),
      x = comparison_group ,
      y = paste0("Level of Expression")
    ) +  
    guides(fill = "none")+
    theme_minimal(base_size = 16)+
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))
}




######### Gene of Interest Function ########

#Grab the first word plus the underscore of this filtered_df column name
GOI_function <- function(dataframe, metadata, gene, col, group1, group2){
  #I'm going to need to modify this because it doesn't work with the current output I have
  
  #Dataframe
  filtered_df <- dataframe %>% filter(GenSymbol == gene)
  filtered_df$GenSymbol <- NULL
  
  #Selecting metadata
  table1 <- metadata %>% filter(!!sym(col) == group1)
  table2 <- metadata %>% filter(!!sym(col) == group2)

  #Selecting samples
  table1_samplenames <- table1[['Samples']]
  table2_samplenames <- table2[['Samples']]
  
  
  #selecting the norm counts of the different groups
  table1_selected <- filtered_df[, colnames(filtered_df) %in% table1_samplenames]
  table2_selected <- filtered_df[, colnames(filtered_df) %in% table2_samplenames]
  
  
  #Adding the prefix
  colnames(table1_selected) <- paste0( group1 ,"_", colnames(table1_selected))
  colnames(table2_selected) <- paste0( group2 ,"_", colnames(table2_selected))
  
  composite_df <- cbind(table1_selected, table2_selected)
  

  
  filtered_df_long <- composite_df %>%
    pivot_longer(
      cols = everything(),
      names_to = "Group",
      values_to = "Value"
    ) %>%
    mutate(Group = sub("_.*", "", Group))  # keep only the part before the underscore
  
  return(filtered_df_long)
}


######### Plotting Genes of Interest ########
#To plot the boxplot you have to first select the gene you are interested in
#I wrote a function to change the format of the table so I can easily boxplot it
#Then once it is output into a variable that you name a relevant name, then you can 
#call the boxplot function

#!!!!!! SEX DIFFERENCES !!!!!!!!!!!!!!!!!

#Oxer1
Sex_CRP <- GOI_function(normalized_counts, metadata, 'CRP', 'Sex', 'M', 'F')
Single_gene_boxplot(Sex_OXER1, 'Sex','CRP')


#CCL4
Sex_APOC3 <- GOI_function(normalized_counts, metadata, 'APOC3', 'Sex', 'M', 'F')
Single_gene_boxplot(Sex_CCL4, "Sex", "APOC3")


#!!!! ANCESTRY DIFFERENCES !!!!!!!!!!!!!!!!!!!!!

#Oxer1

Ancestry_CRP <- GOI_function(normalized_counts, metadata, 'CRP', 'Ancestry', 'EA', 'AA')
Single_gene_boxplot(Ancestry_CRP, 'Ancestry','CRP')


#CCL4

Ancestry_CCL4 <- GOI_function(normalized_counts, metadata, 'APOC3', 'Ancestry', 'EA', 'AA')
Single_gene_boxplot(Ancestry_CCL4, "Ancestry", "APOC3")


######### Export #################

#Remember to set working directory to the boxplot section

#Example of how to save as PNG!!!









