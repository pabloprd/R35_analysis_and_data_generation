#Mass Spec Data Preparation
#Paul Parodi
#Last updated: 11/17/2025

################ GOAL ####################
#Modifying the Metadata and the Mass Spec intensity scores to be able
#To easily run analysis on it later.

#Note: Don't have to do any normalization or any other testing with this dataframe
#because it has already been done by Poochon Scientific. 

####################LIBRARIES#########################
library("writexl")
library("readxl")
library('dplyr')
library('tidyverse')
library('tidyr')

############LOADING IN MASS SPEC DATA#####################
#wd C:/Users/pparodi/Documents/R35/R35_MassSpec
setwd('/Users/paulparodi/R35_v2/Mass Spec/')

#European Ancestry vs African Ancestry, c = EA, d = AA. Using this as the selection dataframe.
EA_vs_AA <- read_excel("Mass_Spec_Original_Data/P1841-Appendix3A-32-sera-Normalized-data-statistics.xlsx", sheet = "Proteins-stattistics-c-vs-d")

########### Grabbing the Intensities #######################

selection_dataframe <- EA_vs_AA %>% select(c("Protein FDR Confidence: Combined","GenSymbol",starts_with("0")))

#High 
High_confidence_protein_intensities <- selection_dataframe %>% filter(`Protein FDR Confidence: Combined` == 'High')

#Medium
Medium_confidence_protein_intensities <- selection_dataframe %>% filter(`Protein FDR Confidence: Combined` == 'Medium')

#Low
Low_confidence_protein_intensities <- selection_dataframe %>% filter(`Protein FDR Confidence: Combined` == 'Low')

#All
All_confidence_protein_intensities <- selection_dataframe
All_confidence_protein_intensities$`Protein FDR Confidence: Combined` <- NULL

########### Creating Metadata ############

Samples <- EA_vs_AA %>% select(starts_with("0"))

metadata <- data.frame(Samples = colnames(Samples), Sex = "Unk", Ancestry = "Unk")

metadata <- metadata %>% mutate(
  Sex = case_when(
    grepl("M", Samples) ~ "M",
    grepl("F", Samples) ~ "F",
    TRUE ~ Sex
  ),
  Ancestry = case_when(
    grepl("EA", Samples) ~ "EA",
    grepl("AA", Samples) ~ "AA",
    TRUE ~ Ancestry
  )
)

########### EXPORT ####################
#Export them all into different files
#We will be outputting the different confidences and whole dataframes as well. 

#High Confidence file
write_xlsx(High_confidence_protein_intensities, "Mass_Spec_Data_Preprocessing/Output/High_confidence_protein_intensities_set1.xlsx")

#Medium Confidence file
write_xlsx(Medium_confidence_protein_intensities, "Mass_Spec_Data_Preprocessing/Output/Medium_confidence_protein_intensities_set1.xlsx")

#Low Confidence file
write_xlsx(Low_confidence_protein_intensities, "Mass_Spec_Data_Preprocessing/Output/Low_confidence_protein_intensities_set1.xlsx")

#All Confidence file
write_xlsx(All_confidence_protein_intensities, "Mass_Spec_Data_Preprocessing/Output/All_confidence_protein_intensities_set1.xlsx")

#Whole dataframe
write_xlsx(EA_vs_AA, "Mass_Spec_Data_Preprocessing/Output/Whole_R35_dataframe_set1.xlsx")

#Metadata
write_xlsx(metadata, "Mass_Spec_Data_Preprocessing/Output/R35_Poochon_MS_metadata.xlsx" )








