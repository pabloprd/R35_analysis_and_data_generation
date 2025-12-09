#Mass Spec Summary Stats
#Paul Parodi
#October 22nd 2025

######### GOAL ###############
#Goal is to calculate summary stats for Mass Spec data

######## LIBRARY #############

library('readxl')
library('writexl')
library('dplyr')
library('tidyverse')

######## LOADING IN FILES ########
setwd("/Users/paulparodi/R35_v2/Mass Spec/")

#Metadata
metadata <- read_xlsx("Mass_Spec_Data_Preprocessing/Output/R35_Poochon_MS_metadata.xlsx")

#All confidence levels
All_counts <- read_xlsx("Mass_Spec_Data_Preprocessing/Output/All_confidence_protein_intensities_set1.xlsx")

#Other confidence levels:


##### Data prep function ###############

Sum_stats <- function(normalized_counts, metadata, col_of_interest){
  
  meta_table <- metadata %>% dplyr::select(!!sym(col_of_interest), `Samples`)
  
  unique_values <- unique(as.list(meta_table[[col_of_interest]]))
  
  Name <- normalized_counts[['GenSymbol']] #Grabbing just the genes
  #test_dataframe <- test_dataframe[order(test_dataframe$YourColumnName), ] use this for later when you want to order better
  
  
  for(val in unique_values){
    
    #grab a list of all the sample names in metadata that correspond with values
    filtered_samples <- meta_table %>% filter(!!sym(col_of_interest) == val)
    sample_names <- filtered_samples[['Samples']]
    
    #Create a new table where you select the columns that are %in% sample_names
    temp_table <- normalized_counts %>% 
      select(all_of(sample_names))
    
    
    summary_stats_table <- Sum_stats_calculations(temp_table, val)
    
    #I need to bind by index
    Name <- transform(merge(Name, summary_stats_table ,by=0 ,all=TRUE), row.names=Row.names, Row.names=NULL)
    #Changing the column names
    
  }
  Final_results <- Name %>%
    rename("GenSymbol" = "x")
  return(Final_results)
}


##### Summary stat function ############
Sum_stats_calculations <- function(table, symbol){
  
  #input the table as well as the unique symbol
  #calculate these using rowwise functions
  range_name <- paste0(symbol, "_range")
  mean_name  <- paste0(symbol, "_mean")
  min_name   <- paste0(symbol, "_min")
  max_name   <- paste0(symbol, "_max")
  std_name   <- paste0(symbol, "_std")
  
  table[[mean_name]] <- rowMeans(table, na.rm = TRUE)
  table[[range_name]] <- apply(table, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  table[[min_name]]  <- apply(table, 1, min, na.rm = TRUE)
  table[[max_name]]  <- apply(table, 1, max, na.rm = TRUE)
  table[[std_name]] <- apply(table, 1, sd,  na.rm = TRUE)
  
  #selected columns:
  sumstats_table <- table %>%
    dplyr::select(all_of(c(mean_name, range_name, min_name, max_name, std_name)))
  
  return(sumstats_table)
}

##### Plot function ##############

scatterplot_function <- function(dataframe, column1, column2){
  ggplot(dataframe, aes(x = !!sym(column1), y = !!sym(column2))) +
    geom_point()+
    labs(x = paste0(column1), y = paste0(column2), title = paste0("Scatterplot of ", column1, " and ", column2))
  
}


######## Running Summary stats #########


Sex_sumstats <- Sum_stats(All_counts, metadata, 'Sex')

Ancestry_sumstats <- Sum_stats(All_counts, metadata, 'Ancestry')

##### Plotting comparison #############

#Example scatterplots using ancestry
scatterplot_function(Ancestry_sumstats, "EA_std", "AA_std")
scatterplot_function(Ancestry_sumstats, "EA_range", "AA_range")
scatterplot_function(Ancestry_sumstats, "EA_mean", "AA_mean")


####### Exporting #############

#Writing Sex
write_xlsx(Sex_sumstats, "Mass_Spec_Data_Analysis/Output/Sumstats/Sex_summary_stats_set1.xlsx")

#Writing Ancestry
write_xlsx(Ancestry_sumstats, "Mass_Spec_Data_Analysis/Output/Sumstats/Ancestry_summary_stats_set1.xlsx")





