#Chord Plot Ontology Visualizations
#By Paul Parodi
#Last updated October 27th 2025
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
setwd("/Users/paulparodi/R35_v2/Mass Spec/")

Ancestry_ranked_plot <- read_xlsx("Mass_Spec_Data_Analysis/Output/Seperated GO tables/Ancestry_ranked_genes.xlsx")

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

#Example
#This is the command to filter to only one gene. 
only_IL10 <- Ancestry_ranked_plot %>% filter(`SYMBOL` == "IL10")

#Chord diagram plot of only a certain ranking!!!!
chord_diagram_plot(Ancestry_ranked_plot, "SYMBOL", "TERM", 2)


#Later try loading in the whole definition system to try and categorize genes of interest


##### Export ############

#setwd as needed. 

#Export PNGS as needed.
ggsave("Ontology_plot.png", plot = p, width = 6, height = 5, dpi = 300)  # Save as PNG


