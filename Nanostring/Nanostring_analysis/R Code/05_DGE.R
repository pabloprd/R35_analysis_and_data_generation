# Differential Gene Expression
# Paul Parodi 
# October 21st 2025



######### GOAL ##############
# Calculate the differential gene expression between two groups in a metadata
# category
##### LIBRARY ###########

library("writexl")
library("readxl")
library("limma")
library("dplyr")

#### LOADING IN THE FILES ###########

setwd("/Users/paulparodi/R35_v2/Nanostring/")

normalized_counts <- read_excel("Nanostring_preprocessing/Output/normalized_counts.xlsx")
metadata <- read_excel("Nanostring_preprocessing/Output/R35_set1_Nanostring_metadata.xlsx")

##### Limma DGE #########

Limma_DGE_test <- function(normalized_counts, metadata_table, testing_column, column1, column2){
  #May have to edit this later to not remove the last three columns as I'll just insert the table with those columns taken out
  Name <- normalized_counts[['Name']]
  normalized_counts$Name <- NULL
  #Creating matrices and groups
  matrix_of_counts <- as.matrix(normalized_counts) # Make the matrix of just the counts. change variable name to be better.
  sample_names <- metadata_table[[testing_column]] # The sample names of the metadata column you are looking at
  group <- factor(sample_names) # Grouping the factors of the sample_names
  design <- model.matrix(~0 + group) # Creating a matrix of the new groups
  colnames(design) <- levels(group) # Make the different levels of the group the column names of the new thing.
  
  #First linear model
  fit <- lmFit(matrix_of_counts, design)
  
  #Contrast matrix
  contrast_name = paste0(column1, 'v', column2) # Make the title
  contrast_formula <- paste0(column1, "-" ,column2) # Putting in the title
  
  contrast.matrix <- limma::makeContrasts(contrasts = setNames(contrast_formula, contrast_name),
                                          levels = design) #Creating the contrast matrix that compares the groups
  
  #More refined linear model
  fit2 <- contrasts.fit(fit, contrast.matrix) #making the linear model based off this contrast matrix
  fit2 <- eBayes(fit2) #running Bayesian correction on the model
  Results <- data.frame(topTable(fit2, adjust.method = 'fdr', number = Inf)) #Outputting the files to results.
  #bind them by the rownames! That's what I should do with the ANOVA!
  
  Results <- transform(merge(Name, Results ,by=0 ,all=TRUE), row.names=Row.names, Row.names=NULL)
  #Changing the column names
  Final_results <- Results %>%
    rename("Name" = "x")
  
  return(Final_results)
}


##### Calling DGE #############

#Remember that with metadata categories with more than two groups you will have to change the names
#To suit that and not confuse yourself

Sex_DGE_table <- Limma_DGE_test(normalized_counts, metadata, 'Sex', "Male", "Female")

Ancestry_DGE_table <- Limma_DGE_test(normalized_counts, metadata, 'Ancestry', "white", "black")

###### EXPORT ###########


#Sex
write_xlsx(Sex_DGE_table,"Nanostring_analysis/Output/Differential Gene Expression/Male_vs_Female_DGE_table.xlsx")

#Ancestry
write_xlsx(Ancestry_DGE_table, "Nanostring_analysis/Output/Differential Gene Expression/Black_vs_White_DGE_table.xlsx")



