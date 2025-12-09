#Chord Plot Ontology Visualizations
#By Paul Parodi
#Last updated October 14th 2025

###### GOAL ########
#To plot the ONTOLOGY (not enrichment) of the genes of interest per metdata group

##### LIBRARY ########
library("readxl")
library("writexl")
library("dplyr")
library("stringr")
library("tidyverse")
library('tidyr')
library('ggplot2')
library("DEP")
library("circlize")
library("GOplot")


##### LOADING IN THE FILES ##########
setwd("/Users/paulparodi/R35_v2/Nanostring/")
Ancestry_GO_table <- read_xlsx("Nanostring_analysis/Output/Ontology/Ancestry_GOI_GO_table.xlsx")
Sex_GO_table <- read_xlsx("Nanostring_analysis/Output/Ontology/Sex_GOI_GO_table.xlsx")

###### Chord Plot function ##########

chord_diagram_plot <- function(dataframe, col1, col2, rank = NA){
  
  #Assuming that the table you loaded in has rank in it
  if(!is.na(rank)){
    dataframe <- dataframe %>% filter(Level == rank)
  }
  else{
    dataframe = dataframe
  }

    
  chordDiagram(dataframe[, c(col1, col2)], 
               annotationTrack = "grid", 
               preAllocateTracks = 1)
  # Add vertical labels (rotated 90 degrees) on the allocated track
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    circos.text(x = mean(xlim), y = ylim[1] + mm_y(2), labels = sector.name,
                facing = "downward",   # 'downward' or 'reverse.clockwise' rotates text vertically
                niceFacing = TRUE,     # automatic rotation adjustment for readability
                adj = c(0.5, 0),      # center horizontally, align bottom vertically
                cex = 0.7)             # label font size
  }, bg.border = NA)
  
  
}

##### Plotting ##############

#ANCESTRY
#I could really dive into this later and make it so that you 
#can select certain genes associated with whatever definition you're interested in
#Chord diagram plot of only a certain ranking!!!!
chord_diagram_plot(Ancestry_GO_table, "SYMBOL", "TERM", 7)

#SEX
chord_diagram_plot(Sex_GO_table, "SYMBOL", "TERM", 2)


##### Export ############


#Export PNGS as needed.
ggsave("Ontology_plot.png", plot = p, width = 6, height = 5, dpi = 300)  # Save as PNG


