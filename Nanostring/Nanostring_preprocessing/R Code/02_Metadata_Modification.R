#02 Metadata Modification
#Paul Parodi



####### GOAL #############
#To modify the metadata so that it can more easily handled and worked with down the line

##### LIBRARY ###########
library('dplyr')
library('tidyverse')
library('readxl')
library('writexl')

#### LOADING IN THE FILES ###########
setwd("/Users/paulparodi/R35_v2/Nanostring/")

mmetadata <- read_csv("Nanostring_original_data/R35_complete_metadata.csv")

##### MODIFYING THE METADATA ########
#1. Renaming
colnames(mmetadata)[colnames(mmetadata) == "Master Sample"] <- "Sample Name"
colnames(mmetadata)[colnames(mmetadata) == "Income Level"] <- "income_group"


#2. Eliminating the NAs from the metadata
mmetadata <- mmetadata %>% filter(!is.na(`Sample Name`))


#3. CREATING NEW EARNING COLUMN
mmetadata <- mmetadata %>% mutate(income_num = as.numeric(str_replace_all(str_extract(income_group, 
                                                                                      "\\d{1,3}(,\\d{3})*"), ",", "")))

#4. INCOME
#above and below 80,0000 just to make sure we have enough in each sample
#Tweak it. 
mmetadata  <- mmetadata  %>%
  mutate(
    income_broad = case_when(
      income_num >= 100000  ~ "HI",
      income_num >= 41000 & income_num <= 99999 ~ "MI",
      income_num <= 40000 ~ "LI",
      TRUE ~ NA_character_
    )
  )


#5. DISCRIMINATION SECTION!
#Making Discrimination categorical HD = High discrimination, LD = Low Discrimination 
mmetadata  <- mmetadata  %>%
  mutate(
    # Self discrimination
    Discrimination = case_when(
      `Cumulative Score Yourself` >= 4 ~ "HD",
      `Cumulative Score Yourself` >= 0 & `Cumulative Score Yourself` <= 3 ~ "LD",
      TRUE ~ NA_character_
    ),
    # Proxy discrimination
    Proxy_discrimination = case_when(
      `Cumulative Score Some Close to you` >= 4 ~ "HD",
      `Cumulative Score Some Close to you` >= 0 & `Cumulative Score Some Close to you` <= 3 ~ "LD",
      TRUE ~ NA_character_
    ),
    # Total discrimination
    Total_Discrimination = case_when(
      `Cumulative Score Discrimination` >= 7 ~ "HD",
      `Cumulative Score Discrimination` >= 0 & `Cumulative Score Discrimination` <= 6 ~ "LD",
      TRUE ~ NA_character_
    )
  )



#6. EDUCATION
#Making Education shorter and easier to read
#high school = HS, some college = SC, college degree or higher = CDOH
mmetadata <-mmetadata %>%
  mutate(Education = case_when(
    Education == "high school" ~ "HS",
    Education == "some college" ~ "SC",
    Education == "college degree or higher" ~ "CDOH",
    TRUE ~ NA_character_  # or another label for unmatched cases
  ))

##############  INDEX SCORING  ##########################################


#8. Scoring income numbers
mmetadata  <- mmetadata  %>%
  mutate(
    income_score = case_when(
      income_num >= 100000  ~ 0,
      income_num >= 41000 & income_num <= 99999 ~ 1,
      income_num <= 40999 ~ 2,
      TRUE ~ NA_integer_
    )
  )

#9. Scoring education
mmetadata <-mmetadata %>%
  mutate(Education_score = case_when(
    Education == "HS" ~ 2,
    Education == "SC" ~ 1,
    Education == "CDOH" ~ 0,
    TRUE ~ NA_integer_  # or another label for unmatched cases
  ))

#10 Scoring Discrimination
mmetadata  <- mmetadata  %>%
  mutate(
    # Self discrimination
    Discrimination_score = case_when(
      `Discrimination` == "HD" ~ 1,
      `Discrimination` == "LD" ~ 0,
      TRUE ~ NA_integer_
    ),
    # Proxy discrimination
    Proxy_discrimination_score = case_when(
      `Proxy_discrimination` == "HD" ~ 1,
      `Proxy_discrimination` == "LD" ~ 0,
      TRUE ~ NA_integer_
    ),
    # Total discrimination
    Total_Discrimination_score = case_when(
      `Total_Discrimination` == "HD" ~ 1,
      `Total_Discrimination` == "LD" ~ 0,
      TRUE ~ NA_integer_
    )
  )



#11. Scoring both
#D score = Discrimination score
mmetadata$D_score = mmetadata$Discrimination_score + mmetadata$Education_score + mmetadata$income_score

#P score = Proxy score  
mmetadata$PD_score = mmetadata$Proxy_discrimination_score + mmetadata$Education_score + mmetadata$income_score

#TD score = Total Discrimination score  
mmetadata$TD_score = mmetadata$Total_Discrimination_score + mmetadata$Education_score + mmetadata$income_score



#12.
#Turning it into score Low_SES, Mid_SES, High_SES:
mmetadata <-mmetadata %>%
  mutate(TD_score_label = case_when(
    TD_score >= 4  ~ "High_SES",
    TD_score <= 3 & TD_score >= 2 ~ "Mid_SES",
    TD_score <= 1 ~ "Low_SES",
    TRUE ~ NA_character_  #or another label for unmatched cases
  ))


##### EXPORTING THE METADATA #########

#Exporting the metadata
write_xlsx(mmetadata ,"Nanostring_preprocessing/Output/R35_set1_Nanostring_metadata.xlsx")