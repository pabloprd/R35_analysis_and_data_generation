#01 Quality Control and normalization
#Paul Parodi


########### GOAL ################
#Run the quality control and normalization of the data

######### LIBRARY ##############
library('readxl')
library('writexl')
library('readr')
library('tidyr')
library('dplyr')
library("purrr") #functional programming tools for lists and vectors
library('nanostringr') #Tools for analyzing Nanostring ncounter data

######### LOADING IN THE FILES #############

setwd("/Users/paulparodi/R35_v2/Nanostring/Nanostring_original_data/")

rcc_files <- list.files("R35_RCCs", pattern = "\\.RCC$", full.names = TRUE)

all_counts <- lapply(rcc_files, parse_counts)
all_counts_table <- reduce(all_counts, full_join, by = c("Code.Class", "Name", "Accession"))

#Master metadata file
mmetadata <- read_csv('/Users/paulparodi/R35_v2/Nanostring/Nanostring_original_data/R35_complete_metadata.csv')


######## Quality Control #############
#Making table of positive controls
pos_controls <- all_counts_table %>% filter(Code.Class == 'Positive') %>% group_by(Name)
#Making table of negative controls
neg_controls <- all_counts_table %>% filter(Code.Class == 'Negative') %>% group_by(Name)

#Calculating the negative geometric means.
#Why do I need to calculate this?
'Geometric means quantify the noise which is used to set thresholds for detecting true signal'
'Perform background correction'

neg_geom_means <- sapply(neg_controls[, grep('^NA', names(neg_controls))], function(x) exp(mean(log(x+1))))

#Making table of housekeeping controls
hk_controls <- all_counts_table[all_counts_table$Code.Class == 'Housekeeping', ]

#Calculating the housekeeping means
hk_means <- sapply(hk_controls[, grep("^NA", names(hk_controls))], mean)





############VISUALIZING QC METRICS##################
#I may just call the boxplot function that I wrote. 
#I just loaded in the MooreNanostring package I made!

#NEGATIVE CONTROLS
boxplot(as.matrix(neg_controls[ , grep("^NA", names(neg_controls))]),
        main = "Negative Control Counts Across Samples",
        ylab = "Counts")
abline(h = neg_geom_means, col = "red", lty = 2)


#POSITIVE CONTROLS
boxplot(as.matrix(pos_controls[ , grep("^NA", names(pos_controls))]),
        main = "Positive Control Counts Across Samples",
        ylab = "Counts")
###############FLAG OUTLIERS##################

#Reinvestigate this. 

#Flagging samples where negative control are more than they should
threshold <- 2 * median(neg_geom_means)
flagged_samples <- names(neg_geom_means)[neg_geom_means > threshold]


#Flagging samples where housekeeping genes are more or less than they should
hk_threshold <- median(hk_means) -2 * mad(hk_means)
flagged_hk <- names(hk_means)[hk_means < hk_threshold]

flagged_samples
#One sample so far was flagged for having a wack housekeeping gene count compared to the others
#that one is NA22344-001 which is a 21 year old white male earning 100-150,000 dollars a year with some 
#college. 

#His positive control count and negative control count is also not bad. For now I am deciding to keep him.




#########NORMALIZATION#########################
#Normalization
#Normalize using housekeeping genes (log2 gene count - mean(log2(housekeeping gene counts))
normalized_counts_OG <-HKnorm(all_counts_table) 

#Getting rid of the other columns in normalized counts and returning a table with the gene index

#converting to R dataframe:
normalized_counts <- as.data.frame(normalized_counts_OG)
normalized_counts$Code.Class <- NULL
normalized_counts$Accession <- NULL

####### EXPORT ############

#Exporting normalized counts 
write_xlsx(normalized_counts, "/Users/paulparodi/R35_v2/Nanostring/Nanostring_preprocessing/Output/normalized_counts.xlsx")