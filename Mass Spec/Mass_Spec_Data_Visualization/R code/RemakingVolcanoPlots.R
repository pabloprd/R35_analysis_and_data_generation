#REMEMBER TO SAVE TABLES/IMAGES/RESULTS AT CERTAIN CHECKPOINTS SO I DON'T HAVE TO REMAKE STUFF
#Try to make each R script focus on one specific job. 

###########INSTALL AND LOAD IN LIBRARIES##########################

#BiocManager::install('writexl')
#BiocManager::install('readxl')
#BiocManager::install('dplyr')
#BiocManager::install('tidyverse')
#BiocManager::install('tidyr')
#BiocManager::install('ggplot2')


library("writexl")
library("readxl")
library('dplyr')
library('tidyverse')
library('tidyr')
library('ggplot2')
library('ggrepel')


#####SET WORKING DIRECTORY###############################
setwd('/Users/paulparodi/Documents/UMD_projects/R35/R35_MassSpec')

#############LOAD IN THE CORRECT DATA TABLE#######################
read_excel("P1841-Appendix3A-32-sera-Normalized-data-statistics.xlsx")

#a = Female, b = Male
Female_vs_male <- read_excel("P1841-Appendix3A-32-sera-Normalized-data-statistics.xlsx", sheet = "Proteins-stattistics-a-vs-b")

#European Ancestry vs African Ancestry, c = EA, d = AA
EA_vs_AA <- read_excel("P1841-Appendix3A-32-sera-Normalized-data-statistics.xlsx", sheet = "Proteins-stattistics-c-vs-d")

#Make dataframe
EA_vs_AA <- data.frame(EA_vs_AA)

#Potentially save these files seperately so they are easier to load. 


########ADDING IN A COLUMN SAYING UP OR DOWNREGULATED#############
#This will be based on whether it is marked with a two and whether it is negative or a postive
#Make a new column called: regulation_level
sigfig <- -log10(0.05)
EA_vs_AA_regulation_df <- EA_vs_AA[, c('log10.ttest.Ratio.C.D..', 'log2.Ratio.C.D.', 'Protein.FDR.Confidence..Combined', 'call..Ratio.C.D..1.25.or..0.8..p.0.05.', 'GenSymbol' )]
EA_vs_AA_regulation_df$regulation.level <- 'None'
# First, initialize the 'regulation.level' column with a default value, e.g., NA
EA_vs_AA_regulation_df[['regulation.level']] <- NA



EA_vs_AA_regulation_df <- EA_vs_AA_regulation_df[
  EA_vs_AA_regulation_df$Protein.FDR.Confidence..Combined == "High", 
]
# Then, use logical indexing to assign 'UP' and 'DOWN'
# The '&' operator performs element-wise logical AND
EA_vs_AA_regulation_df[['regulation.level']][
  EA_vs_AA_regulation_df[['log10.ttest.Ratio.C.D..']] > sigfig & 
  EA_vs_AA_regulation_df[['log2.Ratio.C.D.']] > log2(1.25)
] <- 'UP'

EA_vs_AA_regulation_df[['regulation.level']][
  EA_vs_AA_regulation_df[['log10.ttest.Ratio.C.D..']] > sigfig & 
  EA_vs_AA_regulation_df[['log2.Ratio.C.D.']] < log2(0.8)
] <- 'DOWN'

EA_vs_AA_regulation_df$regulation.level[
  is.na(EA_vs_AA_regulation_df$`call..Ratio.C.D..1.25.or..0.8..p.0.05.`) 
] <- "NO"





########VOLCANO PLOTTING RACE#################################
#c('log10.ttest.Ratio.C.D..', 'log2.Ratio.C.D.', 'call..Ratio.C.D..1.25.or..0.8..p.0.05.', 'GenSymbol' )]
# For log2 fold-change cutoffs at 1.25x and 0.8x:

top_genes <- EA_vs_AA_regulation_df %>%
  arrange(call..Ratio.C.D..1.25.or..0.8..p.0.05.) %>%   # adjust column as needed
  slice(1:50)

vline_thresholds <- log2(c(1.25, 0.8))
hline_threshold <- -log10(0.05) # Should this be +1.3?

title <- "European Ancestry vs African Ancestry Plot" # define title

EA_vs_AA_vplot <- ggplot(data = EA_vs_AA_regulation_df, 
       aes(x = log2.Ratio.C.D., y = log10.ttest.Ratio.C.D.., 
           col = regulation.level, label = GenSymbol)) +
  geom_vline(xintercept = vline_thresholds, col = "black", linetype = 'dashed') +
  geom_hline(yintercept = abs(hline_threshold), col = "black", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_colour_manual(
    values = c("blue", "gray", "red"),
    labels = c("Downregulated", "Not significant", "Upregulated"),
    drop = FALSE
  )+
  geom_text_repel(
    data = top_genes,
    aes(label = GenSymbol),
    max.overlaps = Inf
  )+
  ggtitle(title) +
labs(x = "Log2-Fold Change(European Ancestry/ African Ancestry)", y = "Log10(p-value)", 
     col = "Regulation Level" )+
  theme(
    axis.title = element_text(size =10.0, face = 'bold'),
    axis.text = element_text(size = 10.0, face = 'bold'),
    legend.title = element_text(size = 10.0, face = 'bold'),
    panel.background = element_rect(fill = "white"),             # White background inside plot panel
    plot.background = element_rect(fill = "white"),              # White background around plot
    panel.grid.major = element_line(color = "gray80"),           # Major grid lines gray (lighter)
    panel.grid.minor = element_line(color = "gray90"),
    plot.title = element_text(size = 15.0, hjust = 0.5, face = "bold"))


EA_vs_AA_vplot
  

#######EXPERIMENTAL HEATMAP############################
#Keep an eye out if these are actually significant. 
top_genes <- EA_vs_AA_regulation_df %>%
  arrange(regulation.level) %>%   # adjust column as needed
  slice(1:50)
#Remember to step back and plan. 



ggplot(top_genes, aes(x = log2.Ratio.C.D. , y = GenSymbol, fill = log10.ttest.Ratio.C.D..)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Heatmap using ggplot2",
       x = "Fold Change",
       y = "Protein",
       fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############HEATMAPS##################################
#Selecting the dataframe of just high confidence proteins
EA_vs_AA_high_conf_df <- EA_vs_AA[EA_vs_AA$Protein.FDR.Confidence..Combined == "High", ]

#Select just the samples
EA_vs_AA_samples <- EA_vs_AA_high_conf_df %>% select(GenSymbol, starts_with("X0"))

#Make long format for easier plotting
long_data <- EA_vs_AA_samples %>%
  pivot_longer(
    cols = c(starts_with("X0")),  # replace with column names you want
    names_to = "Metric",
    values_to = "Value"
  )

#Changing long data
long_data <- long_data[long_data$GenSymbol %in% c("NUCB1", "CRP", "TGFB1"), ]


#Plot their abundance levels
ggplot(long_data, aes(x = Metric, y = GenSymbol, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "gray", mid = "yellow", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Heatmap using ggplot2",
       x = "Sample",
       y = "Abundance",
       fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Why is my heatmap showing up differently than the MAss spec peoples?
#They have relative abundance between 0 and 1. My value goes from 0 to 4. A bit confusing
#Should ask about that. 



########EXPORTING THE IMAGES INTO A PICTURES FOLDER##############
ggsave("EA_vs_AA_vplot_v1.png", plot = EA_vs_AA_vplot, width = 6, height = 4, dpi = 300)




