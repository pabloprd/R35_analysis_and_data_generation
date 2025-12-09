#Gene Ontology Directed Acylic Graph
#Last updated: October 15th 2025
#By: Paul Parodi 

#DECOMMISSIONED AS I JUST TRANSFERRED EVERTYIHING TO Gene Ontology and Gene enrichment

######### GOAL ##############
#Detailing my exact plan:
# Grab my ontology graphs
# Take the genes from them
# Get the entrez ids
# Get the levels of each specific definition
# Export the tables so that they can be graphed. Maybe a directed acylic graph but
# at the very least be plotted at the same level.




######## Library ############
library("writexl")
library("readxl")
library('dplyr')
library('tidyverse')
library('tidyr')
library('ggplot2')
library('ggrepel')
library("org.Hs.eg.db")
library('GO.db')
library("topGO")
library("BiocGenerics")
library("parallel")
library("AnnotationDbi")
library("GOxploreR")

######## Loading in the files ##########

setwd("/Users/paulparodi/R35_v2/Nanostring/R35_Nanostring_Data_Analysis/Output/Gene Ontology tables/")

#Loading in the Gene ontology to quickly get genes of interest
Ancestry_ontology_table <- read_excel("Ancestry_GO_definition_table.xlsx")


############ ONTOLOGY FUNCTION #####################
#I wouldn't have to call this Gene ontology function if I could just save it correctly. 
#I could also just transfer this to the Gene Ontology and Enrichment file. 


Gene_ontology_function <- function(genes, col_1 = 'GOALL', col_2 = 'ONTOLOGYALL'){
  #If you want the specific definitions just change col_1 to GO and change col_2 to ONTOLOGY
  
  #Mapping GO to the Annotation database (How we get the GO IDs!)
  go_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = genes,
    columns = c(paste0(col_1), paste0(col_2)),
    keytype = "SYMBOL")
  
  #Selecting the columns I want and filtering them down to only biological processes
  go_map <- go_map %>%   dplyr::select(all_of(c("SYMBOL", col_1, col_2)))
  go_map <- go_map %>% filter(.data[[col_2]] == "BP")
  
  #Getting the verbal definitions from the GO database
  GO_terms <- data.frame(AnnotationDbi::select(GO.db, keys = go_map[[paste0(col_1)]], 
                                               columns = c("TERM"), keytype = "GOID"))
  #Binding the dataframes together
  go_map <- cbind(go_map, GO_terms)
  
  return(data.frame(go_map))
  
}

############ ENTREZ ID FUNCTION ####################


GO_function_and_rank <- function(table){
  
  table_GOI <- unlist(as.list(unique(table %>% dplyr::select("SYMBOL"))))
  table_entrez_ids <- select(org.Hs.eg.db, keys=table_GOI, keytype="SYMBOL", columns="ENTREZID")
  table_ids <- as.list(table_entrez_ids[['ENTREZID']])
  GOID_rank_table <- Gene2GOTermAndLevel(genes = table_ids, organism = "Homo sapiens", domain = "BP")
  Gene_ontology_table <- Gene_ontology_function(as.matrix(table_ids$SYMBOL), 'GO', 'ONTOLOGY')
  level_vector <- setNames(GOID_rank_table$Level, GOID_rank_table$`GO ID` )
  Updated_GO_table <-   Gene_ontology_table %>% mutate(Level  = level_vector[GOID])
  Updated_GO_table$GO <- NULL
  
  return( Updated_GO_table)
}


####### Getting the ENTREZ IDS ####################



#Maybe I move this to the 05_Gene_Ontology_enrichment

Ancestry_GOI <- unlist(as.list(unique(Ancestry_ontology_table %>% dplyr::select("SYMBOL"))))
Ancestry_entrez_ids <- select(org.Hs.eg.db, keys=Ancestry_GOI, keytype="SYMBOL", columns="ENTREZID")
Ancestry_ids <- as.list(Ancestry_entrez_ids[['ENTREZID']])
Gene2GOTermAndLevel(genes = Ancestry_ids, organism = "Homo sapiens", domain = "BP")

test_entrez <- Entrez_ids_and_gene2GO(Ancestry_ontology_table)
test_GO <- Gene_ontology_function(as.matrix(entrez_ids$SYMBOL), 'GO', 'ONTOLOGY')
#I should create a new dataframe with just the symbol, the Ontology, the GOID, the term, and the label
#It may take a while though. What is the most efficient way to do this??

#maybe an efficient way would be to make a dictionary from test_GO levels and GOID

#This vector solution may be important for later. 
level_vector <- setNames(test_entrez$Level, test_entrez$`GO ID` )
level_vector
test_GO <- test_GO %>% mutate(Level  = level_vector[GOID])
test_GO
########## Mapping gene definitions ##################



###### 




#### DAG function ##########


##### Plotting #########


###### EXPORT ############

