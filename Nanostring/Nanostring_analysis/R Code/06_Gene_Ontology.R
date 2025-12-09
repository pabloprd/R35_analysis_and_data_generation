# 06 Gene Ontology
# Paul Parodi 
# October 21st 2025


######### GOAL ##############
# Get gene ontology for R35 Data
# Use the GO explore R package to organize the data in a way
# that can easily be divded into categories. 

#SOMETHING TO WATCH OUT FOR: The GO.db package I call is not super comprehensive. It's just a good way to organize stuff
#If you truly want a comprehensive set of definition call on Gene ontology. And also, if you have a ton of genes
#you'll just have to call in gene enrichment. 

##### LIBRARY ###########
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
#### LOADING IN THE FILES ###########
setwd("/Users/paulparodi/R35_v2/Nanostring/")

#Load in Differential Gene Expression tables
Ancestry_DGE_table <- read_excel("Nanostring_analysis/Output/Differential Gene Expression/Black_vs_White_DGE_table.xlsx")
Sex_DGE_table <- read_excel("Nanostring_analysis/Output/Differential Gene Expression/Male_vs_Female_DGE_table.xlsx")


####### Genes of interest function ##############

Genes_of_interest <- function(DGE_table, pvalcut = 0.05, logcut = 1){
  #logFC filter
  biglog <- DGE_table %>% filter( logFC >= logcut )
  smalllog <-  DGE_table %>% filter(logFC <= -logcut )

  New_table <- rbind(biglog, smalllog)
  #Stacking them
  New_table <- New_table %>% filter(`adj.P.Val` <= pvalcut)
  
  print(New_table[['Name']])
  return(New_table[['Name']])
  
}

######## GENE ONTOLOGY FUNCTION ################

Gene_ontology_function <- function(GOI){
  
  GOI_entrez_ids <- select(org.Hs.eg.db, keys=GOI, keytype="SYMBOL", columns="ENTREZID")
  
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



########## Running Gene Ontology ################
Ancestry_GOI <- Genes_of_interest(Ancestry_DGE_table, pvalcut = 0.05, logcut = 0.5)
Ancestry_GO_table <- Gene_ontology_function(Ancestry_GOI )

#These are not of interset. This is just an example:

Sex_GOI <- Genes_of_interest(Sex_DGE_table, pvalcut = 1.0, logcut = 0.5)
Sex_GO_table <- Gene_ontology_function(Sex_GOI)

##### EXPORT ###########

#Ancestry GO
##GOI means gene of interest 
write_xlsx(Ancestry_GO_table ,"Nanostring_analysis/Output/Ontology/Ancestry_GOI_GO_table.xlsx")

#Sex GO
write_xlsx(Sex_GO_table ,"Nanostring_analysis/Output/Ontology/Sex_GOI_GO_table.xlsx")













