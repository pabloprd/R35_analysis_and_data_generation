

#I will have to split this into a few different functions. #First is creating the new group names. 

column_renaming <- function(df_of_samples, metadata, col, sample_col){

#Making dictionary from metadata
sample_to_featuredict <- setNames(metadata[[paste0(col)]], metadata[[paste0(sample_col)]])


#Making a copy of the dataframe to edit
df_copy <- df_of_samples
namelist <- list()

#For each sample name in the lupus samples data frame find the key and add the  _group to it
for(i in names(sample_to_featuredict)){
  newcol <- paste0(i, "_", sample_to_featuredict[[i]], "_group")
  namelist <- append(namelist, newcol)}

#New dictionary mapping it
sample_to_newsample_name_dict <- setNames(namelist, metadata[[paste0(sample_col)]])


#Renaming the column names so they have group in it. 
df_copy <- df_copy %>% rename(!!!setNames(names(sample_to_newsample_name_dict), unname(sample_to_newsample_name_dict)))
df_copy <- data.frame(df_copy) #Making dataframe

return(df_copy)
}

protein_stats_runner <- function(df_of_samples, metadata, col, sample_col){

  
#######################PAIRWISE TESTS AND RENAMING####################################    
#Renaming the columns by groups
df_copy <- column_renaming(df_of_samples, metadata, col, sample_col)

#turning all 0s into NAs so I can calculate correct averages. Important!!
df_copy[df_copy == 0] <- NA

#Making the genes the rows
rownames(df_copy) <- df_copy[['GenSymbol']]
df_copy$GenSymbol <- NULL
protein_matrix <- as.data.frame(df_copy[ , !(names(df_copy) == "GenSymbol") ])

#Making a data frame of the different pairwise tests we can run for each group
combinations <- t(combn(unique(metadata[[paste0(col)]]), 2))
pairwise_tests <- as.data.frame(combinations)
colnames(pairwise_tests) <- c("Group1", "Group2")




###############CALCULATING THE AVERAGE###########################
#for each group in column calculate the average for each row
#Append the average to a new column

#This one is tricky, may come back and edit
for(item in unique(metadata[[paste0(col)]]) ){ #this line works
  group_name <- (paste0(item, "_group")) #This works
  #selecting column names that have the rght group
  df_temp <- df_copy %>% select(which(grepl(paste0(group_name), names(df_copy))))
  df_temp_copy <- df_temp
  df_temp_copy$average <- rowMeans(df_temp, na.rm = TRUE)
  df_copy[[paste0(group_name,".avg")]] <-  df_temp_copy$average
}


###############CALCULATING THE T-TEST###########################################

#REMEMEBER TO REDO THIS WITH A NEW DATAFRAME BECAUSE IT WILL COUNT THE AVERAGES AND THE LOG2 STUFF AS VALUES
group_list = list()

#Iterate through each pairwise test and then run 
for(i in 1:nrow(pairwise_tests)){
  group1 <- (pairwise_tests[i, 1])
  group2 <- (pairwise_tests[i, 2])
  group1_avg_name <- (paste0(group1, "_group.avg"))
  group2_avg_name <- (paste0(group2, "_group.avg"))
  group_name <- paste0("Ratio(", group1, "_vs_", group2, ")")
  
  #Getting the ratio
  df_copy[[paste0(group_name)]] <-  df_copy[[paste0(group1_avg_name)]]/ df_copy[[paste0(group2_avg_name)]]
  #MS_Lupus_samples_copy[[paste0("Log2(", group_name, ")")]] <-  log2(MS_Lupus_samples_copy[[paste0(group1_avg_name)]]/ MS_Lupus_samples_copy[[paste0(group2_avg_name)]])
  
  #I want to select every group1 thing, and select every group2 thing, and then do the test test
  group1_df <- df_copy %>% select(which(grepl(paste0(group1), names(df_copy))))
  group2_df <- df_copy %>% select(which(grepl(paste0(group2), names(df_copy))))
  
  #Not just do the t-test 
  #I want to select every group1 thing, and select every group2 thing, then call it 
  #MS_Lupus_samples_copy[[paste0("Ttest(", group_name, ")")]] <- t.test(group1_df, group2_df, var.equal = FALSE )
  
  #Then do the ratio, send message to Emily talking about it tho
  
  #group_list <- append(group_list, group_name)

  ########################THE END###################################  
  
}
return(df_copy)
}
