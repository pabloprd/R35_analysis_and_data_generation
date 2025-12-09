#03 Summary Statistics 
#Paul Parodi


######### GOAL ##############
#Calculate the summary statistics the metadata categories of interest

##### LIBRARY ###########
library('dplyr')
library('readxl')
library('writexl')
library('tidyr')

#### LOADING IN THE FILES ###########
#set working directory
setwd("/Users/paulparodi/R35_v2/Nanostring/")

#Load in normalized counts
normalized_counts <- read_xlsx("Nanostring_preprocessing/Output/normalized_counts.xlsx")

#Load in metadata
metadata <- read_xlsx("Nanostring_preprocessing/Output/R35_set1_Nanostring_metadata.xlsx")


Sex_sumstats <- read_xlsx("Nanostring_analysis/Output/Summary Stats/Sex_summary_stats.xlsx")
Ancestry_sumstats <- read_xlsx("Nanostring_analysis/Output/Summary Stats/Ancestry_summary_stats.xlsx")

##### Data prep function ###############

Sum_stats <- function(normalized_counts, metadata, col_of_interest){
  
  meta_table <- metadata %>% dplyr::select(!!sym(col_of_interest), `Sample Name`)
  
  unique_values <- unique(as.list(meta_table[[col_of_interest]]))
  
  Name <- normalized_counts[['Name']] #Grabbing just the genes
  #test_dataframe <- test_dataframe[order(test_dataframe$YourColumnName), ] use this for later when you want to order better
  
  
  for(val in unique_values){
    
    #grab a list of all the sample names in metadata that correspond with values
    filtered_samples <- meta_table %>% filter(!!sym(col_of_interest) == val)
    sample_names <- filtered_samples[['Sample Name']]
    
    #Create a new table where you select the columns that are %in% sample_names
    temp_table <- normalized_counts %>% 
      select(all_of(sample_names))
    
    
    summary_stats_table <- Sum_stats_calculations(temp_table, val)
    
    #I need to bind by index
    Name <- transform(merge(Name, summary_stats_table ,by=0 ,all=TRUE), row.names=Row.names, Row.names=NULL)
    #Changing the column names

  }
  Final_results <- Name %>%
    rename("Name" = "x")
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

##### Plotting function ############
scatterplot_function <- function(dataframe, column1, column2){
  ggplot(dataframe, aes(x = !!sym(column1), y = !!sym(column2))) +
    geom_point()+
    labs(x = paste0(column1), y = paste0(column2), title = paste0("Scatterplot of ", column1, " and ", column2))
  
}
###### Running Summary stats ###########

Sex_sumstats <- Sum_stats(normalized_counts, metadata, 'Sex')
Ancestry_sumstats <- Sum_stats(normalized_counts, metadata, 'Ancestry')


##### Plotting comparison #############

#Example scatterplots
scatterplot_function(Ancestry_sumstats, "black_std", "white_std")
scatterplot_function(Ancestry_sumstats, "black_range", "white_range")
scatterplot_function(Ancestry_sumstats, "black_mean", "white_mean")

##### Export ############

#Exporting Ancestry
write_xlsx(Ancestry_sumstats ,"Nanostring_analysis/Output/Summary Stats/Ancestry_summary_stats.xlsx")

#Exporting Sex
write_xlsx(Sex_sumstats ,"Nanostring_analysis/Output/Summary Stats/Sex_summary_stats.xlsx")




