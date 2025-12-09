#03.5 Gene Tables
#Paul parodi
#October 20th 2025

######### GOAL ##############
#Create a table specifically for each metadata category

##### LIBRARY ###########
library("dplyr")
library('tidyverse')
library('writexl')
library('readxl')

#### LOADING IN THE FILES ###########
setwd("/Users/paulparodi/R35_v2/Nanostring/")

metadata <- read_xlsx("Nanostring_preprocessing/Output/R35_set1_Nanostring_metadata.xlsx")
normalized_counts <- read_xlsx("Nanostring_preprocessing/Output/normalized_counts.xlsx")

#### Gene table function ##########

Gene_table_function <- function(normalized_counts, metadata, col_of_interest){
  
  #Meta table
  meta_table <- metadata %>% dplyr::select(!!sym(col_of_interest), `Sample Name`)
  
  #unique values
  unique_values <- unique(as.list(meta_table[[col_of_interest]]))
  
  #Name
  Name <- normalized_counts[['Name']]
  
  for(val in unique_values){
    #grab a list of all the sample names in metadata that correspond with values
    filtered_samples <- meta_table %>% filter(!!sym(col_of_interest) == val)
    
    #Sample names
    sample_names <- filtered_samples[['Sample Name']]
    
    #Create a new table where you select the columns that are %in% sample_names
    temp_table <- normalized_counts %>% 
      select(all_of(sample_names))
    
    final_table <- cbind(Name, temp_table)
    
    filename <- paste0(val, "_gene_table.xlsx")
    
    write_xlsx(
      final_table,
      path = file.path("Nanostring_analysis/Output/Gene Tables", paste0(filename))
    )
    
    
  }
}



######### Running Gene table function ##############

#Ancestry
Gene_table_function(normalized_counts, metadata, "Ancestry")

#Sex
Gene_table_function(normalized_counts, metadata, "Sex")




