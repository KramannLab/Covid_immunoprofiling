# Copyright (c) [2022] [Felix Schreibing]
# feschreibing@ukaachen.de


# Identification of CD8+ NK-like TEMRA cell marker genes ----

library(Seurat)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(tidyverse)


# Load the sc dataset ----
indir = '~/sciebo/COVID_project/Robjects/'
outdir = '~/sciebo/COVID_project/01 Results/WP2/2.2 NK-like TEMRA cell marker genes/'

sc = readRDS(paste0(indir, 'integrated.RNA.Tcells.annotated.revision.filtered.rds'))
Idents(sc) = 'integrated_annotations'


# Calculate CD8+ NK-like TEMRA cell marker genes ----
NK_TEMRA_markers <- FindMarkers(sc, ident.1 = 'CD8+ NK-like TEMRA cells', min.pct = 0.25, logfc.threshold = 0.25)

write.csv(NK_TEMRA_markers, paste0(outdir, 'NK_like_TEMRA_marker_genes.txt'))


# plot CD8+ NK-like TEMRA cell marker genes ----
NK_TEMRA_markers = rownames_to_column(NK_TEMRA_markers, 'gene')

colors <- ifelse(
  NK_TEMRA_markers$avg_log2FC < -0.5 & NK_TEMRA_markers$p_val_adj < 0.05, '#440154',
  ifelse(NK_TEMRA_markers$avg_log2FC > 0.5 & NK_TEMRA_markers$p_val_adj < 0.05, '#FDE725',
         'black'))
colors[is.na(colors)] <- 'black'
names(colors)[colors == '#FDE725'] <- 'up'
names(colors)[colors == 'black'] <- 'not significant'
names(colors)[colors == '#440154'] <- 'down'

pdf(file = paste0(outdir, 'NK_like_TEMRA_marker_genes_volcano_plot.pdf'), width = 11, height = 9)
EnhancedVolcano(
  NK_TEMRA_markers,
  lab = NK_TEMRA_markers$gene,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  xlim = c(-3, 3),
  selectLab = c('KLRF1', 'KLRC2', 'KIR2DL3', 'KIR3DL1', 'KIR3DL2', 'FCGR3A', 'KLRD1', 'IL7R'),
  pCutoff = 0.05,
  FCcutoff = 0.5,
  colCustom = colors,
  pointSize = 3,
  colAlpha = 1,
  cutoffLineType = "solid",
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  labSize = 3.0,
  labCol = 'black',
  labFace = 'bold',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  legendPosition = 'right',)
dev.off()