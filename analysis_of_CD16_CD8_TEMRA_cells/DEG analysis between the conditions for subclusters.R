# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de
# edited by Felix Schreibing

# Differential gene expression analysis for CD16+ CD8+ TEMRA subtypes between the conditions----

library(Seurat)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)


# Load the subclustered sc datasets ----
indir = '~/sciebo/COVID_project/Robjects/'
outdir = '~/sciebo/COVID_project/01 Results/WP5/5.4 DEG analysis between the conditions for subclusters/'

NK_TEMRA = readRDS(file = paste0(indir, 'NK-like_TEMRA_subclusters.rds'))


# Redefine the differential gene expression function ----
get_markers = function(sc, condition, control, cell.type, only.sig = TRUE, adj.pval.cutoff = 0.05, logfc.threshold = 0.25, out.dir = 'diff_genes'){
  # Take care if no cells
  markers = data.frame()
  condition.subset = tryCatch(subset(sc, NK_TEMRA_subtype.condition == paste(cell.type, condition)), 
                              error = function(e) {data.frame()})
  control.subset = tryCatch(subset(sc, NK_TEMRA_subtype.condition == paste(cell.type, control)), 
                            error = function(e) {data.frame()})
  
  if (ncol(condition.subset) > 20 & ncol(control.subset) > 20){
    markers = FindMarkers(sc, ident.1 = paste(cell.type, condition), 
                          ident.2 = paste(cell.type, control), min.pct = 0.25, logfc.threshold = logfc.threshold)
    
    if (only.sig){
      markers = markers[markers$p_val_adj < adj.pval.cutoff,]
    }
    
    
    markers$gene = rownames(markers)
    cols = colnames(markers)
    cell.type.name = gsub('/', '_',cell.type)
    cell.type.name = gsub(' ', '_',cell.type.name)
    
    if (nrow(markers) > 0){
      write.table(markers[,c('gene',cols[-6])], 
                  file = paste0(out.dir, '/integrated.diff.genes.', cell.type.name,
                                '.', condition, '.vs.', control, '.txt'),
                  quote = FALSE, sep = '\t', row.names = FALSE)
    }
  }
  return(markers)
}


# Run DEG analysis ----
cell.types = levels(Idents(NK_TEMRA))
NK_TEMRA$NK_TEMRA_subtype.condition = paste(NK_TEMRA$NK_TEMRA_subset_annotations, NK_TEMRA$condition_collapsed)
Idents(NK_TEMRA) = 'NK_TEMRA_subtype.condition'


# DEGanalysis subtype vs subtype across conditions
for (cell.type in cell.types){
  # mild vs. healthy
  condition = 'mild'
  control = 'healthy'
  
  markers = get_markers(NK_TEMRA, condition = condition, 
                        control = control, 
                        cell.type = cell.type, only.sig = FALSE,
                        out.dir = outdir)
  
  # severe vs. healthy
  condition = 'severe'
  control = 'healthy'
  
  markers = get_markers(NK_TEMRA, condition = condition, 
                        control = control, 
                        cell.type = cell.type, only.sig = FALSE,
                        out.dir = outdir)
  
  # severe vs. mild
  condition = 'severe'
  control = 'mild'
  
  markers = get_markers(NK_TEMRA, condition = condition, 
                        control = control, 
                        cell.type = cell.type, only.sig = FALSE,
                        out.dir = outdir)
}


# Plot DGE results for subcluster 2 between severe and mild condition ----
genes = read.delim(paste0(outdir, 'integrated.diff.genes.CD16+_CD8+_TEMRA_2.severe.vs.mild.txt'))

colors <- ifelse(
  genes$avg_log2FC < -0.25 & genes$p_val_adj < 0.05, '#8AFF6B',
  ifelse(genes$avg_log2FC > 0.25 & genes$p_val_adj < 0.05, '#FF5D5D',
         'black'))
colors[is.na(colors)] <- 'black'
names(colors)[colors == '#FF5D5D'] <- 'up'
names(colors)[colors == 'black'] <- 'not significant'
names(colors)[colors == '#8AFF6B'] <- 'down'

pdf(file = paste0(outdir, 'DEGs_NK_TEMRA_subset_2_severe_vs_mild.pdf'), width = 7, height = 7)
EnhancedVolcano(
  genes,
  lab = genes$gene,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  pCutoff = 0.05,
  FCcutoff = 0.25,
  selectLab = c('HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DPA1', 'IFITM1', 'IFITM3', 'KLRC2', 'KIR3DL2'),
  xlim = c(-4.5, 3),
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