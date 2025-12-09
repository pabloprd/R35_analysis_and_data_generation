#04 ANOVA
#Paul Parodi
#October 21st

######### GOAL ##############
# Calculate the ANOVA for each one of the metadata categories and
#output a table comparing all of them directly. 

##### LIBRARY ###########

library("writexl")
library("readxl")
library('dplyr')
library('tidyverse')
library('tidyr')
library('limma')

#### LOADING IN THE FILES ###########
setwd("/Users/paulparodi/R35_v2/Nanostring/")

normalized_counts <- read_xlsx("Nanostring_preprocessing/Output/normalized_counts.xlsx")
metadata <- read_xlsx("Nanostring_preprocessing/Output/R35_set1_Nanostring_metadata.xlsx")

###### ANOVA function ##########

#This is the Anova function that uses the limma package

#I need to check again to see if the order is messed up at all. 
Anova_limma <- function(normalized_counts, metadata, col_of_interest){
  
  Name <- normalized_counts[["Name"]] #Removing and Saving the name of the normalized counts
  normalized_counts$Name <- NULL #Removing the non numeric values of normalized counts
  
  exprs <- as.matrix(normalized_counts) #Creates matrix of normalized counts
  sample_names <- metadata[[col_of_interest]]
  
  group <- factor(sample_names)
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  groups <- colnames(design) #stores the names of the columns in the design matrix (which are group names like control and treatment) into a variable called groups
  contrast_pairs <- combn(groups, 2, simplify = FALSE)
  contrast_strings <- sapply(contrast_pairs, function(pair) paste(pair, collapse = "-"))
  
  contrast.matrix <- do.call(makeContrasts, 
                             c(lapply(contrast_strings, function(c) 
                               parse(text = c)[[1]]), list(levels = design)))
  
  fit <- lmFit(exprs, design) #Making linear model
  fit2 <- contrasts.fit(fit, contrast.matrix) #
  fit2 <- eBayes(fit2)  #A new fitted model that contains results for my specified contrasts, rather than just the original group coefficients. 
  Results <- topTable(fit2, adjust.method = 'fdr', number = Inf) #Outputting the files to results.
  
  Results <- transform(merge(Name, Results ,by=0 ,all=TRUE), row.names=Row.names, Row.names=NULL)
  #Changing the column names
  Final_results <- Results %>%
    rename("Name" = "x")
  
  return(Final_results)
}


##### Calling ANOVA ##########
#Calling the Anova from the function

#Income broad
Income_broad_anova <- Anova_limma(normalized_counts, metadata, 'income_broad')

#Sex
Sex_anova <- Anova_limma(normalized_counts, metadata, 'Sex')

#Ancestry
Ancestry_anova <- Anova_limma(normalized_counts, metadata, "Ancestry")

###### EXPORT #############

#Income broad
write_xlsx(Income_broad_anova, "Nanostring_analysis/Output/ANOVA/Income_broad_ANOVA_table.xlsx")

#Sex broad
write_xlsx(Sex_anova, "Nanostring_analysis/Output/ANOVA/Sex_ANOVA_table.xlsx")


#Ancestry broad
write_xlsx(Ancestry_anova, "Nanostring_analysis/Output/ANOVA/Ancestry_ANOVA_table.xlsx")




















##### Export ###############