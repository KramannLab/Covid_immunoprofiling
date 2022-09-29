# Copyright (c) [2022] [Felix Schreibing]
# feschreibing@ukaachen.de

# Use CrossTalkeR results to investigate and plot NK cell-related cell-cell interactions ----

library(ggplot2)
library(dplyr)
library(tidyverse)
library(readxl)
library(cowplot)
library(viridis)
library(ggpubr)

indir = '~/sciebo/COVID_project/CrossTalkR/'
outdir = '~/sciebo/COVID_project/01 Results/WP4/4.4 NK-like TEMRA related cell-cell interactions/'


# Load 'LR analysis data' output tables from CrossTalkR ----
CCI_mild_vs_healthy = read_excel(paste0(indir, 'LR_analysis_mild_vs_healthy.xlsx'), sheet = 1)
CCI_severe_vs_healthy = read_excel(paste0(indir, 'LR_analysis_severe_vs_healthy.xlsx'), sheet = 1)
CCI_severe_vs_mild = read_excel(paste0(indir, 'LR_analysis_severe_vs_mild.xlsx'), sheet = 1)

CCI_mild_vs_healthy = CCI_mild_vs_healthy %>%
  mutate(condition = 'mild_vs_healthy')

CCI_severe_vs_healthy = CCI_severe_vs_healthy %>%
  mutate(condition = 'severe_vs_healthy')

CCI_severe_vs_mild = CCI_severe_vs_mild %>%
  mutate(condition = 'severe_vs_mild')

CCI = rbind(CCI_mild_vs_healthy, CCI_severe_vs_healthy, CCI_severe_vs_mild)


# Plot interesting differential CCIs ----
NK_related_interactions = c('HLA-E-KLRC1', 'HLA-C-KIR2DL1', 'HLA-C-KIR2DL2', 'HLA-C-KIR2DL3', 'HLA-B-KIR3DL1', 'HLA-A-KIR3DL2', 'KITLG-KIT', 'FLT3LG-FLT3',
                 'IL2-IL2RA', 'IL2-IL2RB', 'IL2-IL2RG', 'IL15-IL15RA', 'IL15-IL2RA', 'IL15-IL2RB', 'IL15-IL2RG', 'KLRB1-CLEC2D', 'CLEC2D-KLRB1',
                 'LTA-LTBR', 'LTA-TNFRSF1A', 'LTA-TNFRSF1B' ,'LTA-TNFRSF14', 'LTB-LTBR', 'FASLG-FAS', 'IL21-IL21R', 'IL21-IL2RG', 'CXCL8-CXCR1',
                 'CXCL8-CXCR2', 'CCL19-CCR7', 'SELL-SELPLG', 'SELPLG-SELL', 'CCL3-CCR1', 'CCL3-CCR4', 'CCL3-CCR5')

# define ligand cluster cell types
cell_type = c('Plasma', 'NK', 'Neu', 'Mono', 'Mast', 'Macro', 'Epi', 'DC', 'CD4', 'B')

# define comparisons
comparisons = c('mild_vs_healthy', 'severe_vs_healthy', 'severe_vs_mild')


# Plot CrossTalkeR results for NK cell-related interactions
for(comparison in comparisons){
  
  df = CCI %>%
    filter(Receptor.Cluster == 'CD8+ NK-like TEMRA cells') %>%
    mutate(interaction = paste0(Ligand, '-', Receptor)) %>%
    filter(interaction %in% NK_related_interactions) %>%
    filter(Ligand.Cluster %in% cell_type) %>%
    filter(condition == comparison)
  
  pdf(file = paste0(outdir, 'NK_related_CCI_NK-like_TEMRA_dotplot_', comparison, '_.pdf'), width = 8, height = 5)
  print(ggplot(df, mapping = aes(x = interaction, y = Ligand.Cluster, size = abs(LRScore)))  +
          geom_point(color = ifelse(df$LRScore < 0, '#330A5F', '#ED6925')) +
          theme_cowplot() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
                axis.text.y = element_text(size = 8, hjust = 1),
                axis.ticks = element_blank(),
                strip.text.x = element_text(size = 6),
                panel.grid.major = element_line(colour = 'gray'),
                panel.border = element_rect(colour = 'black')) +
          scale_size(trans = 'reverse')
  )
  dev.off()
  
}


# Comparing LR scores for selected interactions between mild and severe conditions ----
# Generating a data frame with single interaction values for mild and severe conditions
CCI_mild = read_excel(paste0(indir, 'LR_analysis_single_mild.xlsx'), sheet = 1)
CCI_severe = read_excel(paste0(indir, 'LR_analysis_single_severe.xlsx'), sheet = 1)

CCI_mild = CCI_mild %>%
  mutate(condition = 'mild')

CCI_severe = CCI_severe %>%
  mutate(condition = 'severe')

single_CCI = rbind(CCI_mild, CCI_severe)

single_CCI = single_CCI %>%
  mutate(interaction = paste0(Ligand, '-', Receptor))

# Define interesting interactions
interactions = c('IL15-IL2RB', 'IFNG-IFNGR1', 'IL18-CD48', 'IL27-IL27RA', 'CCL2-CCR5', 'CCL3-CCR4', 'CCL3-CCR5', 'CCL4-CCR5',
                 'HLA-A-KIR3DL2', 'HLA-B-KLRD1', 'HLA-C-KIR2DL3', 'HLA-E-KLRC1', 'HLA-E-KLRC2', 'HLA-E-KLRD1', 'HLA-E-KLRK1', 'HLA-F-KIR3DL2')

for(i in interactions){
  
  # Subset for interesting interactions with only CD8+ T cells as Receptor Cluster
  test = single_CCI %>%
    filter(interaction == i) %>%
    filter(Receptor.Cluster != 'CD4' & Receptor.Cluster != 'NK' & Receptor.Cluster != 'Plasma',
           Receptor.Cluster != 'B', Receptor.Cluster != 'DC', Receptor.Cluster != 'Epi',
           Receptor.Cluster != 'Macro', Receptor.Cluster != 'Mast', Receptor.Cluster != 'Mega',
           Receptor.Cluster != 'Mono', Receptor.Cluster != 'Neu')
  
  
  # Plot comparison
  pdf(file = paste0(outdir, 'CCI_severe_vs_mild_', i, '_.pdf'), width = 4, height = 5)
  print(ggplot(test, aes(x = condition, y = LRScore)) +
          geom_boxplot(width = 0.5, color = '#00204D', lwd = 2) +
          geom_jitter(color = '#808080') +
          stat_compare_means() +
          theme_cowplot() +
          theme(text = element_text(size = 8),
                axis.text.x = element_text(angle = 90, size = 8),
                axis.title.x = element_blank()) +
          xlab(i)
  )
  dev.off()
  
}