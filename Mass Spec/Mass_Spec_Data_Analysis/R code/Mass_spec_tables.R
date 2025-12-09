#Mass Spec data exploratory analysis
#Paul Parodi
#Last updated: 10/14/2025

####################LIBRARIES#########################


library("writexl")
library("readxl")
library('dplyr')
library('tidyverse')
library('tidyr')
library('ggplot2')
############LOADING IN MASS SPEC DATA#####################
#wd C:/Users/pparodi/Documents/R35/R35_MassSpec
setwd('/Users/paulparodi/Documents/UMD_projects/R35/R35_MassSpec')
#P1841-Appendix1A-32-Sera-MS-data-Proteins

Mass_spec_dataframe <- read_excel("P1841-Appendix3A-32-sera-Normalized-data-statistics.xlsx")

#Read only the desired sheet
#a = Female, b = Male
Female_vs_male <- read_excel("P1841-Appendix3A-32-sera-Normalized-data-statistics.xlsx", sheet = "Proteins-stattistics-a-vs-b")

#European Ancestry vs African Ancestry, c = EA, d = AA
EA_vs_AA <- read_excel("P1841-Appendix3A-32-sera-Normalized-data-statistics.xlsx", sheet = "Proteins-stattistics-c-vs-d")


#############MAKING THE METADATA#################
Samples <- Female_vs_male %>% select(starts_with("0"))
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

#############TOP AND BOTTOM SELECTION FUNCTIONS###########

#I want to select only the top 5 genes of interest
top_selection <- function(table, column_name, number){
  table1 <- table %>% slice_max(order_by = .data[[column_name]], n = number)
  return(table1)
}

bottom_selection <- function(table, column_name, number){
  table1 <- table %>% slice_min(order_by = .data[[column_name]], n = number)
  return(table1)
}


############BOXPLOT FUNCTION########################
Single_gene_boxplot_exp <- function(metadataframe, x, y,  plot_title = NA, x_axis_title = NA, y_axis_title = NA, y_limit = -5, line_size = 1.0, point_size = 2, palette_set  = "Set1",
                                    x_axis_title_size = 16, y_axis_title_size = 16, plot_title_size = 16){
 
  #Merging the metadata 

  
  #Plotting stuff
  ggplot(metadataframe, aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])) +
    geom_boxplot(color = "black", size = line_size, outlier.shape = NA) +
    geom_jitter(width = 0.3, color = "black", alpha = 0.7, size = point_size) +#adding in the dots
    scale_fill_brewer(palette = paste0(palette_set)) + #determining color of inner boxes #What are the different pallets?
    labs(
      title = paste(plot_title),
      x = x_axis_title,
      y = y_axis_title
    ) +  
    theme_minimal(base_size = 16) +
    theme(
      axis.title = element_text(size = x_axis_title_size),
      axis.text = element_text(size = y_axis_title_size),
      plot.title = element_text(size = plot_title_size, hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14))
}
############CREATING JUST THE TOTAL COUNTS DATAFRAME###############
Sex_counts <- Female_vs_male
Ancestry_counts <- EA_vs_AA

#ANCESTRY
#Only selecting the column names that start with '0'
Sex_counts <- Sex_counts %>% select(GenSymbol, starts_with("0"))
rownames(Sex_counts) <- Sex_counts$GenSymbol
Ancestry_counts <- Ancestry_counts %>% select(GenSymbol, starts_with("0"))
rownames(Ancestry_counts) <- Ancestry_counts$GenSymbol
############FEMALE VS MALE############################
colnames(Female_vs_male)
#Columns of interest
#target protein
#core pathway
#description
#Average
#call (Ratioa/b) > 1/25 or <0.8 p <0.05
#Ratio(a/b)
#subcellular location

#female vs male differentially expressed
FvsM_DE <- top_selection(Female_vs_male,"call (Ratio(a/b)>1.25 or <0.8, p<0.05)", 3)
curated_FvsM_DE <- select(FvsM_DE, GenSymbol, "Protein names", "Core pathway", "Pathway", "Sub-pathway1","Subcellular location [CC]", "Function [CC]",
                          "Gene ontology (biological process)", "Gene ontology (cellular component)", "Gene ontology (molecular function)")
barplot(table(FvsM_DE$`Core pathway`), cex.names = 0.6, cex.axis = 0.5) #adding in the do

############FEMALE VS MALE BOXPLOTTING###########################
protein_boxplot_sel <- function(counts, gene, metadata){
  
  #Making dataframe of just that gene index
  dataframe <- counts[gene,]
  
  #Pivoting the table so that it can be used for boxplot
  dataframe <- pivot_longer(dataframe, cols = -GenSymbol, names_to = "Samples", values_to = "Expression")
  
  #Merging with dataframe
  dataframe %>% inner_join(metadata, by = "Samples")
}


#DOPO Male vs Female
Dopo_counts <- protein_boxplot_sel(Sex_counts, "DOPO", metadata)
Single_gene_boxplot_exp(Dopo_counts, 'Sex', "Expression", "DOPO Expression between Females and Males", "Sex", "Expression")


#FCG3A
FCG3A_counts <- protein_boxplot_sel(Sex_counts, "FCG3A", metadata)
Single_gene_boxplot_exp(FCG3A_counts, 'Sex', "Expression", "FCG3A Expression between Females and Males", "Sex", "Expression")


#	CSF1R
CSF1R_counts <- protein_boxplot_sel(Sex_counts, "CSF1R", metadata)
Single_gene_boxplot_exp(CSF1R_counts, 'Sex', "Expression", "CSF1R Expression between Females and Males", "Sex", "Expression")


##############AA VS EA###################################

#BARPLOT


colnames(EA_vs_AA)

#European ancestry vs african ancestry differential expression
EAvsAA_DE <- top_selection(EA_vs_AA,"call (Ratio(C/D)>1.25 or <0.8, p<0.05)", 96)
EAvsAA_DE

#Curated version
curated_EAvsAA_DE <- select(EAvsAA_DE, GenSymbol, "Protein names", "Core pathway", "Pathway", "Sub-pathway1", "Subcellular location [CC]" , "Function [CC]",
                            "Gene ontology (biological process)", "Gene ontology (cellular component)", "Gene ontology (molecular function)")

#fillna with unknown. 
curated_EAvsAA_DE$`Core pathway`[is.na(curated_EAvsAA_DE$`Core pathway`)] <- "Undetermined"

#I need to group the curated version into a graph I'm pretty sure. 
#What variable should I use for that?
#There are just so many NAs
#For the ones that have multiple categories it marks them as NA in the core pathway column. I should just plot the core pathways right now

#Let's make a barplot of core pathways. 
#X is the core pathway labels, y is the counts
#par(mar = c(5.1, 4.1, 4.1, 2.1)) normal margins
par(mar = c(8, 4, 4, 2) + 1.0)
barplot(table(curated_EAvsAA_DE$`Core pathway`), cex.names = 0.6, cex.axis = 0.5,las = 2) #adding in the do


#Making the barplot legible (need to take out the brackets. )
par(mar = c(8, 4, 4, 2) + 15.0)
barplot(table(curated_EAvsAA_DE$Pathway), cex.names = 0.6, cex.axis = 0.5, las = 2) #adding in the do


#BOXPLOT
ggplot(merged_gene_dataframe, aes(x = .data[[column]], y = Expression, fill = .data[[column]])) +
  geom_boxplot(color = "black", size = line_size, outlier.shape = NA) +  # Thicker box
  geom_jitter(width = 0.3, color = "black", alpha = 0.7, size = point_size) + #adding in the dots
  scale_fill_brewer(palette = paste0(palette_set)) + #determining color of inner boxes #What are the different pallets?
  labs(
    title = paste(plot_title, gene),
    x = column,
    y = y_axis_title
  )




############AA VS EA BOXPLOTTING###########################

#Targets we were looking for:

#CRP Male vs Female
CRP_counts <- protein_boxplot_sel(Sex_counts, "CRP", metadata)
Single_gene_boxplot_exp(CRP_counts, 'Ancestry', "Intensity", "CRP Expression between African Ancestry and European Ancestry", "Ancestry", "Intensity")


#TGB1 (NO SIGNIFICANT DIFFERENCE)
TGFB1_counts <- protein_boxplot_sel(Sex_counts, "TGFB1", metadata)
Single_gene_boxplot_exp(TGFB1_counts, 'Ancestry', "Intensity", "TGFB1 Expression between African Ancestry and European Ancestry", "Ancestry", "Intensity")



#IL10
IL10_counts <- protein_boxplot_sel(Sex_counts, "IL10", metadata)
Single_gene_boxplot_exp(IL10_counts, 'Ancestry', "Expression", "IL10 Expression between African Ancestry and European Ancestry", "Ancestry", "Intensity")

#######EXPORTING TO Excel FOR NANA###########################


write_xlsx(FvsM_DE, "Sex_DPE_Massspec.xlsx")
write_xlsx(EAvsAA_DE, "Ancestry_DPE_Massspec.xlsx")

#######EXPORTING IMPORTANT TABLES TO EXCEL#####################
#Metadata table
write_xlsx(metadata, "R35_massspec_metadatav2.xlsx")

#Curated tables (Differential expressed)
write_xlsx(curated_EAvsAA_DE, "DE_ancestry_proteins.xlsx")
write_xlsx(curated_FvsM_DE, "DE_sex_proteins.xlsx")

#Male and Female tables
write_xlsx(Female_vs_male, "Sex_table.xlsx")
write_xlsx(EA_vs_AA, "Ancestry_table.xlsx")

