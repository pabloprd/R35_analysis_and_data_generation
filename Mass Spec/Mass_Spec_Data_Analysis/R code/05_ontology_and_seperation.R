#Ontology and seperation
#October 23rd 2025
#Paul Parodi


########## GOAL ###############
# To organize the Poochon Ontology dataframes that we have. 


######## Library #############
library("readxl")
library("writexl")
library("dplyr")
library("tidyverse")
library("org.Hs.eg.db")
library('GO.db')
library("topGO")
library("BiocGenerics")
library("parallel")
library("AnnotationDbi")
library('ggkegg')
library('tidygraph')
library("KEGGREST")
library("GOxploreR")


####### Loading in files ###########

setwd("/Users/paulparodi/R35_v2/Mass Spec/")

#Loading in the whole Poochon dataframe
Poochon_dataframe <- read_excel("Mass_Spec_Data_Preprocessing/Output/Whole_R35_dataframe_set1.xlsx")

Poochon_dataframe$`Gene ontology IDs`
###### All the GO columns of interest to view ################

test <- c("sortID", "GenSymbol", "Protein names", "Core pathway", 
          "Pathway", "Sub-pathway1", "Cross-reference (KEGG)", 
          "Cross-reference (Reactome)", "Cross-reference (UniPathway)",
          "Cross-reference (PANTHER)", "Function [CC]", 
          "Involvement in disease", "Gene ontology (GO)", 
          "Gene ontology (cellular component)", 
          "Gene ontology (molecular function)", 
          "Gene ontology IDs", "Subcellular location [CC]")


############## Separation function ###########################

#This function is mean to create a dataframe of just two columns that seperate each
# column of interest according to how you need it. 

seperation_function <- function(dataframe, gene_column, col_of_interest, new_colname){
  
  #Making a smaller selective dataframe
  col_selection_df <- setNames(data.frame(dataframe[[gene_column]]), gene_column)
  
  col_selection_df[[new_colname]] <- dataframe[[col_of_interest]]
  
  #Seperation
  col_selection_df <- col_selection_df %>%
    separate_rows(!!sym(new_colname), sep = ";") %>%
    mutate(!!sym(new_colname) := str_trim(!!sym(new_colname)))
  
  
  #removing the GO boxes
  col_selection_df <- col_selection_df %>%
    mutate(!!sym(new_colname) := str_remove_all(!!sym(new_colname), "\\[GO:[0-9]+\\]"))
  
  #removing hsa prefixes
  col_selection_df <- col_selection_df %>%
    mutate(!!sym(new_colname) := str_remove_all(!!sym(new_colname), "hsa:"))
  

  
  #Removing any blank rows
  col_selection_df <- col_selection_df %>% filter(!!sym(new_colname) != "" )
  
  return(col_selection_df)
  
}

######## GENE ONTOLOGY  RANKING FUNCTION ################

Gene_ontology_function <- function(GOI){
  
  GOI_entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, keys=GOI, keytype="SYMBOL", columns="ENTREZID")
  
  GOI_entrez_ids[['ENTREZID']] <- as.integer(GOI_entrez_ids[['ENTREZID']])
  
  GO_table <- Gene2GOTermAndLevel(genes = GOI_entrez_ids, organism = "Homo sapiens", domain = "BP")
  
  GO_table[['GOID']] <- GO_table[['GO ID']]
  
  GO_table[['GO ID']] <- NULL
  
  GO_table[['ENTREZID']] <- GO_table[['Entrezgene ID']]
  
  GO_table[['Entrezgene ID']] <- NULL
  
  GO_definition_table <- data.frame(AnnotationDbi::select(GO.db, keys = GO_table[['GOID']], 
                                                          columns = c("TERM"), keytype = "GOID"))
  
  
  #Now I just need to append the terms to the GO table by the GO ID
  Final_table <- inner_join(GO_table, GO_definition_table, by = "GOID", relationship = "many-to-many")
  
  Final_table <- inner_join(GOI_entrez_ids, Final_table, by = "ENTREZID", relationship = "many-to-many")
  
  Final_table <- distinct(Final_table)
  
  return(Final_table)
  
}



############# Running separations ######################

GO_dataframe <- seperation_function(Poochon_dataframe, "GenSymbol", "Gene ontology (biological process)", "GO")

Kegg_dataframe <- seperation_function(Poochon_dataframe, "GenSymbol", "Cross-reference (KEGG)", "KEGG IDs")

Reactome_dataframe <- seperation_function(Poochon_dataframe, "GenSymbol", "Cross-reference (Reactome)", "Reactome IDs")

Core_pathway_dataframe <- seperation_function(Poochon_dataframe, "GenSymbol", "Core pathway", "Core pathway")

Pathway_dataframe <- seperation_function(Poochon_dataframe, "GenSymbol", "Pathway", "Pathway")

Sub_pathway_dataframe <- seperation_function(Poochon_dataframe, "GenSymbol", "Sub-pathway1", "Pathway")

GO_ids_dataframe <- seperation_function(Poochon_dataframe, "GenSymbol", "Gene ontology IDs","GOIDs")

Accession_ids_dataframe <- seperation_function(Poochon_dataframe, "GenSymbol", "Accession","Accession")


########### RUNNING RANKING ##################
#here is where we run our own ontology with ranking
GO_levels_dataframe <- Gene_ontology_function(Poochon_dataframe$GenSymbol)

############ Export ####################

#Exporting all gene ontology definitions
write_xlsx(GO_dataframe, "Mass_Spec_Data_Analysis/Output/Gene Ontology tables/Gene_ontology_dataframe.xlsx")

#Exporting all Kegg IDs
write_xlsx(Kegg_dataframe, "Mass_Spec_Data_Analysis/Output/Gene Ontology tables/Kegg_id_dataframe.xlsx" )

#Exporting Reactome
write_xlsx(Reactome_dataframe, "Mass_Spec_Data_Analysis/Output/Gene Ontology tables/Reactome_dataframe.xlsx" )

#Exporting pathways
write_xlsx(Core_pathway_dataframe, "Mass_Spec_Data_Analysis/Output/Gene Ontology tables/Core_pathway_dataframe.xlsx")

write_xlsx(Pathway_dataframe, "Mass_Spec_Data_Analysis/Output/Gene Ontology tables/Pathway_dataframe.xlsx")

write_xlsx(Sub_pathway_dataframe, "Mass_Spec_Data_Analysis/Output/Gene Ontology tables/Sub_pathway_dataframe.xlsx")

#Exporting GO ids
write_xlsx(GO_ids_dataframe, "Mass_Spec_Data_Analysis/Output/Gene Ontology tables/GO_IDs_dataframe.xlsx")

#Exporting accession ids
write_xlsx(Accession_ids_dataframe, "Mass_Spec_Data_Analysis/Output/Gene Ontology tables/Accession_IDs_dataframe.xlsx")

#Exporting GO ranking
write_xlsx(GO_levels_dataframe, "Mass_Spec_Data_Analysis/Output/Gene Ontology tables/GO_levels_dataframe.xlsx")






