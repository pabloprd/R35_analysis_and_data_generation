#Gene Enrichment and Trimming for Plotting
#October 23rd 2025
#Paul Parodi

############# GOAL #################
#Calculate the gene enrichment for the genes of interest


########### Library ################
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
library('pheatmap')
library("enrichplot")
library("clusterProfiler")
library("UniProt.ws")

#Do I want to add a significance function in DEG or in here? I think it would possibly be better in DEG.
#Nvm I may want to change this again and again depending on what we're doing. I can just call it from another file. 

########## Loading in files ############
setwd("/Users/paulparodi/R35_v2/Mass Spec/")

#Ancestry
Ancestry_DEG <- read_excel("Mass_Spec_Data_Analysis/Output/Differential Gene Expression Tables/EA_vs_AA_DGE.xlsx")

#Accession dataframe
Accession_dataframe <- read_excel("Mass_Spec_Data_Analysis/Output/Gene Ontology tables/Accession_IDs_dataframe.xlsx")


#Kegg dataframe
Kegg_df <- read_excel("Mass_Spec_Data_Analysis/Output/Gene Ontology tables/Kegg_id_dataframe.xlsx")

#GO ids
GO_id_df <- read_excel("Mass_Spec_Data_Analysis/Output/Gene Ontology tables/GO_IDs_dataframe.xlsx")

#Gene ontology dataframe
Gene_ontology_df <- read_excel("Mass_Spec_Data_Analysis/Output/Gene Ontology tables/Gene_ontology_dataframe.xlsx")

#Gene levels ontology dataframe 
Gene_levels_ontology_df <- read_excel("Mass_Spec_Data_Analysis/Output/Gene Ontology tables/GO_levels_dataframe.xlsx")

####### Marking Significance Function ##############

significance_function <- function(dataframe, pval_cutoff = 0.05, logfc_cutoff = 1.0){
  
  dataframe[['significance']] <- NA
  
  dataframe$significance[dataframe$adj.P.Val <= pval_cutoff & dataframe$logFC >  logfc_cutoff] <- "significant"
  
  dataframe$significance[dataframe$adj.P.Val <= pval_cutoff & dataframe$logFC < -logfc_cutoff] <- "significant"
  
  return(dataframe)
  
}

########## Gene enrichment function ############
#Should I modify this so it either takes the gene symbols or just the gene names?
gene_enrichment_function <- function(dataframe, goidtable, pvalcutoff = 0.05, qvalcutoff = 0.2){

  dataframe <- dataframe %>% filter(significance == "significant")
  new_protein_table <- goidtable %>% filter(goidtable$GenSymbol %in% dataframe$GenSymbol)
  genes <- new_protein_table$GenSymbol
  genes <- unique(genes)
  print(genes)
  
  go2entrez <- AnnotationDbi::select(org.Hs.eg.db,
                      keys = new_protein_table$GOIDs,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "GO")  # or "GO" for direct (non-inferred) terms
  
  go2entrez_filtered <- go2entrez %>% filter( ONTOLOGY == 'BP' )

  results <- enrichGO(
    gene =   genes,
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = pvalcutoff,
    pAdjustMethod = "BH",
    qvalueCutoff = qvalcutoff,
    minGSSize = 3,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE)
  
return(data.frame(results))
}
######## KEGG enrichment function ###############

Kegg_enrichment_function <- function(dataframe, kegg_table, pvalcutoff = 0.05, qvalcutoff = 0.2){
  dataframe <- dataframe %>% filter(significance == "significant")
  new_protein_table <- kegg_table %>% filter(kegg_table$GenSymbol %in% dataframe$GenSymbol)
  proteins <- unique(new_protein_table$`KEGG IDs`)
  
  results <- enrichKEGG(
    gene = new_protein_table$`KEGG IDs`,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = pvalcutoff,
    pAdjustMethod = "BH",
    minGSSize = 2,
    maxGSSize = 500,
    qvalueCutoff = qvalcutoff,
    use_internal_data = FALSE)
  
  return(data.frame(results))
}

########## Running: GENE AND PROTEIN ENRICHMENT #################

#ANCESTRY::::
# marking which genes are significant for the ancestry category
Ancestry_DEG <- significance_function(Ancestry_DEG)


# GENE ENRICHMENT (should be 0.05 and 0.2 respectively as they are the pvalue and the qvalue)
Ancestry_gene_enrichment <- gene_enrichment_function(Ancestry_DEG, GO_id_df, 0.5, 0.5)

# KEGG (change the last numbers to fit your needs better. 0.05 and 0.2 respectively) I just have them like this 
# So they can produce nonzero results here for example. 
Ancestry_kegg <- Kegg_enrichment_function(Ancestry_DEG, Kegg_df, 0.5, 0.5)




########## Trimming dataframes ###################

#Goal of this dataframe is to trim them down to only the genes of interest
Dataframe_trim <- function(dataframe, definition_table, gene_column){
  dataframe <- dataframe %>% filter(significance == "significant")
  new_table <- definition_table %>% filter(definition_table[[paste0(gene_column)]] %in% dataframe$GenSymbol)
  return(new_table)
}

########## Running: Trimming ###################

#Trimming the genes that are not in the Ancestry differentially expressed genes or definitions
Ancestry_ranked <- Dataframe_trim(Ancestry_DEG, Gene_levels_ontology_df, "SYMBOL")
Ancestry_definitions <- Dataframe_trim(Ancestry_DEG, Gene_ontology_df, "GenSymbol")


########## Export ##############

#Exporting enrichment
write_xlsx(Ancestry_gene_enrichment, "Mass_Spec_Data_Analysis/Output/Gene Enrichment Tables/Ancestry_GO_gene_enrichment.xlsx")
write_xlsx(Ancestry_kegg, "Mass_Spec_Data_Analysis/Output/Gene Enrichment Tables/Ancestry_kegg_protein_enrichment.xlsx")



#Exporting trimmed

write_xlsx(Ancestry_ranked, "Mass_Spec_Data_Analysis/Output/Seperated GO tables/Ancestry_ranked_genes.xlsx")
write_xlsx(Ancestry_definitions, "Mass_Spec_Data_Analysis/Output/Seperated GO tables/Ancestry_ontology.xlsx")






