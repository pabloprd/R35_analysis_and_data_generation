#Volcano plot and boxplot modification for Alison
#By Paul Parodi
#Last updated: August 1st
###########GOALS#########################
#Goal: Remake some of the Mass Spec figures
#Enlarge the volcano plot, reduce the range of y axis and x axis
#Make clearer
#Only mark the significant proteins in the middle
#Make dots color blind friendly
#Fonts of significant proteins bigger
#More descriptive x-axis
#better graph title


###########WORKFLOW######################
#Load in the neccessary tables from the folder
#Mark up regulated and down regulated for volcano plot
#volcano plot it
#Change boxplots colors and titles
#Export images into files. 

###########LIBRARY#####################
library("writexl")
library("readxl")
library('dplyr')
library('tidyverse')
library('tidyr')
library('ggplot2')
library('ggrepel')
###########LOADING IN####################
setwd('/Users/paulparodi/Documents/UMD_projects/R35/R35_MassSpec')

#Female and male table
Female_vs_male <- read_excel("Sex_table.xlsx")

#European and African Ancestry
EA_vs_AA <- read_excel("Ancestry_table.xlsx")

#Metadata
metadata <- read_excel("R35_massspec_metadatav2.xlsx")
###########METADATA TWEAKING###########
#changing metadata labels for better plotting
metadata$Sex[metadata$Sex == 'F'] <- 'Female'
metadata$Sex[metadata$Sex == 'M'] <- 'Male'

##########VOLCANO PLOT PREP##############

#I'm going to have to mark as upregulated, downregulated, and None
#Will also have to take out all the ones that are not marked as high confidence. 
#Then volcano plot them over. 

#This will be based on whether it is marked with a two and whether it is negative or a postive
#Make a new column called: regulation_level
sigfig <- -log10(0.05)
Female_vs_male <- data.frame(Female_vs_male)


Female_vs_male_regulation_df <- Female_vs_male[, c('log10.ttest.Ratio.a.b..', 'Protein.FDR.Confidence..Combined',
                                              'log2.Ratio.a.b.', 'call..Ratio.a.b..1.25.or..0.8..p.0.05.', 'GenSymbol' )]
Female_vs_male_regulation_df$regulation.level <- 'None'
# First, initialize the 'regulation.level' column with a default value, e.g., NA
Female_vs_male_regulation_df[['regulation.level']] <- NA


#I need an extra line here removing all the no non high confidence values
Female_vs_male_regulation_df <- Female_vs_male_regulation_df[
  Female_vs_male_regulation_df$Protein.FDR.Confidence..Combined == "High", 
]

# Then, use logical indexing to assign 'UP' and 'DOWN'
# The '&' operator performs element-wise logical AND
Female_vs_male_regulation_df[['regulation.level']][
  Female_vs_male_regulation_df[['log10.ttest.Ratio.a.b..']] > sigfig & 
    Female_vs_male_regulation_df[['log2.Ratio.a.b.']] > log2(1.25)
] <- 'UP'

Female_vs_male_regulation_df[['regulation.level']][
  Female_vs_male_regulation_df[['log10.ttest.Ratio.a.b..']] > sigfig & 
    Female_vs_male_regulation_df[['log2.Ratio.a.b.']] < log2(0.8)
] <- 'DOWN'

Female_vs_male_regulation_df$regulation.level[
  is.na(Female_vs_male_regulation_df$`call..Ratio.a.b..1.25.or..0.8..p.0.05.`) 
] <- "NO"



#########VOLCANO PLOTTING#################
top_genes <- Female_vs_male_regulation_df %>%
  arrange(call..Ratio.a.b..1.25.or..0.8..p.0.05.) %>%   # adjust column as needed
  slice(1:3)

vline_thresholds <- log2(c(1.25, 0.8))
hline_threshold <- -log10(0.05) # Should this be +1.3?

title <- "Protein Abundance Profiling: Female versus Male" # define title

Female_vs_male_volcanoplot <- ggplot(data = Female_vs_male_regulation_df, 
       aes(x = log2.Ratio.a.b., y = log10.ttest.Ratio.a.b.., 
           col = regulation.level, label = GenSymbol)) +
  geom_vline(xintercept = 0, col = "black", size = 0.2) +    # y axis line
#  geom_hline(yintercept = 0, col = "black", size = 1) +    # x axis line
  geom_vline(xintercept = vline_thresholds, col = "black", linetype = 'dashed', size = 0.5) +
  geom_hline(yintercept = abs(hline_threshold), col = "black", linetype = 'dashed', size = 0.5) +
  geom_point(size = 2) +
  scale_colour_manual(
    values = c("#8C8B8B", "red"),
    labels = c( "Not significant", "Upregulated"),
    drop = FALSE
  )+
  geom_text_repel(
    data = top_genes,
    aes(label = GenSymbol),
    max.overlaps = Inf
  )+
  ggtitle(title)+
  labs(x = "Log2-Fold Change(Female/Male)", y = "Log10(p-value)", 
       col = "Regulation Level" )+
  theme(
  axis.title = element_text(size =10.0, face = 'bold'),
  axis.text = element_text(size = 10.0, face = 'bold'),
  legend.title = element_text(size = 10.0, face = 'bold'),
  panel.background = element_rect(fill = "white"),             # White background inside plot panel
    plot.background = element_rect(fill = "white"),              # White background around plot
    panel.grid.major = element_line(color = "gray80"),           # Major grid lines gray (lighter)
    panel.grid.minor = element_line(color = "gray90"),
  legend.position = "none", # Minor grid lines even lighter gray
  plot.title = element_text(size = 15.0, hjust = 0.5, face = "bold"))



#########BOXPLOT FUNCTIONS####################
#Function that does the boxplots:
Single_gene_boxplot_exp <- function(metadataframe, x, y,  plot_title = NA, x_axis_title = NA, y_axis_title = NA, y_limit = -5, line_size = 1.0, point_size = 2, palette_set  = "Set1",
                                    x_axis_title_size = 16, y_axis_title_size = 16, plot_title_size = 20){
  
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
    guides(fill = "none")+
    theme_minimal(base_size = 16) +
    theme(
      axis.title = element_text(size = x_axis_title_size),
      axis.text = element_text(size = y_axis_title_size),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(size = plot_title_size, hjust = 0.5, face = "bold"))
  
}



#Function that selects the genes
protein_boxplot_sel <- function(counts, gene, metadata){
  
  #Making dataframe of just that gene index
  dataframe <- counts[gene, ]
  
  #Pivoting the table so that it can be used for boxplot
  dataframe <- pivot_longer(dataframe, cols = -GenSymbol, names_to = "Samples", values_to = "Expression")
  
  #Merging with dataframe
  dataframe %>% inner_join(metadata, by = "Samples")
}

##########CREATING TOTAL COUNTS########################
Sex_counts <- Female_vs_male
Ancestry_counts <- EA_vs_AA

Sex_counts <- Sex_counts %>% select(GenSymbol, starts_with("0"))
Ancestry_counts <- Ancestry_counts %>% select(GenSymbol, starts_with("0"))
rownames(Sex_counts) <- Sex_counts$GenSymbol
rownames(Ancestry_counts) <- Ancestry_counts$GenSymbol

########SEX BOXPLOT###########################

#DOPO Male vs Female
Dopo_counts <- protein_boxplot_sel(Sex_counts, "DOPO", metadata)
Single_gene_boxplot_exp(Dopo_counts, 'Sex', "Expression", "DOPO Relative Abundance", "Sex", "Abundance", palette_set  = "Set2")
Dopo_boxplot <- Single_gene_boxplot_exp(Dopo_counts, 'Sex', "Expression", "DOPO Relative Abundance", "Sex", "Abundance", palette_set  = "Set2")


#FCG3A
FCG3A_counts <- protein_boxplot_sel(Sex_counts, "FCG3A", metadata)
Single_gene_boxplot_exp(FCG3A_counts, 'Sex', "Expression", "FCG3A Relative Abundance", "Sex", "Abundance", palette_set  = "Set2")
FCG3A_boxplot <- Single_gene_boxplot_exp(FCG3A_counts, 'Sex', "Expression", "FCG3A Relative Abundance", "Sex", "Abundance", palette_set  = "Set2")


#	CSF1R
CSF1R_counts <- protein_boxplot_sel(Sex_counts, "CSF1R", metadata)
Single_gene_boxplot_exp(CSF1R_counts, 'Sex', "Expression", "CSF1R Relative Abundance", "Sex", "Abundance", palette_set  = "Set2")
CSFR1_boxplot <- Single_gene_boxplot_exp(CSF1R_counts, 'Sex', "Expression", "CSF1R Relative Abundance", "Sex", "Abundance", palette_set  = "Set2")
 

######EXPORTING BOXPLOTS AS PNGS#####################
ggsave("DOPO_boxplot_MvF.png", plot = Dopo_boxplot, width = 6, height = 4, dpi = 300)
ggsave("FCG3A_boxplot_MvF.png", plot = FCG3A_boxplot, width = 6, height = 4, dpi = 300)
ggsave("CSFR1_boxplot_MvF.png", plot = CSFR1_boxplot, width = 6, height = 4, dpi = 300)



######EXPORTING VOLCANO PLOTS AS PNGS##############
ggsave("Female_vs_Male_volcano_plot_v1.png", plot = Female_vs_male_volcanoplot, width = 6, height = 4, dpi = 300)

