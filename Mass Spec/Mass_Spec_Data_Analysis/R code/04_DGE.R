#Differential Gene Expression
#By Paul Parodi
#October 22nd 2025



######### GOAL ##############
#Poochon Scientific already ran their tests on the data
#However, we want to run our own tests just for better validity
#We are going to be using limma which is a great package for just this thing. 


######## Library ############

library("writexl")
library("readxl")
library("limma")
library("dplyr")

####### Loading files ########
setwd("/Users/paulparodi/R35_v2/Mass Spec/")



#All counts
All_counts <- read_excel("Mass_Spec_Data_Preprocessing/Output/All_confidence_protein_intensities_set1.xlsx")

#Can test this on other counts as well, high, medium, low

Metadata <- read_excel("Mass_Spec_Data_Preprocessing/Output/R35_Poochon_MS_metadata.xlsx")

###### Limma Differential Gene Expression function #############

Limma_DGE_test <- function(normalized_counts, metadata_table, testing_column, column1, column2){
#May have to edit this later to not remove the last three columns as I'll just insert the table with those columns taken out
Name <- normalized_counts[['GenSymbol']]
normalized_counts$GenSymbol <- NULL
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
  rename("GenSymbol" = "x")

return(Final_results)
}


###### Running DGE ###########

Sex_DGE <- Limma_DGE_test(All_counts, Metadata, 'Sex', 'M', 'F')

Ancestry_DGE <- Limma_DGE_test(All_counts, Metadata, 'Ancestry', 'EA', 'AA')

############### Export ##################

#Sex
write_xlsx(Sex_DGE, "Mass_Spec_Data_Analysis/Output/Differential Gene Expression Tables/Male_vs_Female_DGE.xlsx")

#Ancestry
write_xlsx(Ancestry_DGE, "Mass_Spec_Data_Analysis/Output/Differential Gene Expression Tables/EA_vs_AA_DGE.xlsx")



