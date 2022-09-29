# Copyright (c) [2022] [Felix Schreibing]
# feschreibing@ukaachen.de

# Calculate different functional scores to characterize subclusters of CD16+ CD8+ TEMRA cells----

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)


# Load the subclustered sc datasets ----
indir = '~/sciebo/COVID_project/Robjects/'
outdir = '~/sciebo/COVID_project/01 Results/WP5/5.2 characterization of subclusters/'

NK_TEMRA = readRDS(file = paste0(indir, 'NK-like_TEMRA_subclusters.rds'))


# Load defined color schemes ----
NK.TEMRA.subset.colors = c('CD16+_CD8+_TEMRA_1' = '#8FD744',
                           'CD16+_CD8+_TEMRA_2' = '#9700AD',
                           'CD16+_CD8+_TEMRA_3' = '#00F1FF',
                           'CD16+_CD8+_TEMRA_4' = '#FFFA00',
                           'CD16+_CD8+_TEMRA_5' = '#FF5D5D',
                           'CD16+_CD8+_TEMRA_6' = '#000000')


# Calculate marker genes for every subcluster ----
# Plot marker genes as Heatmap
NK_TEMRA_markers <- FindAllMarkers(NK_TEMRA, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)

write.csv(NK_TEMRA_markers, paste0(outdir, 'NK_TEMRA_marker_genes.txt'))

top10 = NK_TEMRA_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

pdf(file = paste0(outdir, 'NK_TEMRA_subsets_marker_genes.pdf'), width = 17, height = 15)
DoHeatmap(NK_TEMRA, features = top10$gene, 
          group.colors = NK.TEMRA.subset.colors, group.bar.height = 0.05,
          angle = 90, size = 3) + 
  scale_fill_gradientn(colors = c("#00204D", "#7C7B78", "#FFEA46"))
dev.off()

# Plot top 5 marker genes for CD8+ NK-like TEMRA subtypes as dotplot
top5 = NK_TEMRA_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

feature = unique(top5$gene)

pdf(file = paste0(outdir, 'NK_TEMRA_subsets_marker_genes_dotplot.pdf'), width = 6, height = 9)
DotPlot(object = NK_TEMRA, features = feature) +
  coord_flip() +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8, hjust = 1),
        axis.ticks = element_blank()) +
  scale_color_viridis(option = 'E')
dev.off()


# Characterize the CD8+ NK-like TEMRA subsets by calculating different scores ----
# Load HAY_BONE_MARROW_NK_CELLS gene set (M39204)
# https://www.gsea-msigdb.org/gsea/msigdb/cards/HAY_BONE_MARROW_NK_CELLS
NK_cell_BM = read.csv(paste0(outdir, 'HAY_BONE_MARROW_NK_CELLS .txt'))


# Load REACTOME_APOPTOSIS gene set (M15303)
# https://www.gsea-msigdb.org/gsea/msigdb/cards/REACTOME_APOPTOSIS.html
apoptosis = read.csv(paste0(outdir, 'REACTOME_APOPTOSIS.txt'))


# Define the necessary scores
# exhaustion and cytotoxicity scores adapted from https://www.biorxiv.org/content/10.1101/2020.03.19.998658v2
exhaustion_signature = c("HAVCR2", "PDCD1", "ENTPD1", "TNFRSF9", "CCL3", "PHLDA1", "SIRPG", "CTLA4", "TIGIT",
                         "SNAP47", "WARS", "CD27", "RGS1", "TNFRSF1B", "CD27-AS1", "ACP5", "AFAP1L2", "MYO7A",
                         "AKAP5", "ADGRG1", "TOX", "LAYN", "CD38", "ITGAE", "RGS2", "CCND2", "FKBP1A", "MTHFD1",
                         "GALM", "CSF1", "SNX9", "ICOS", "LYST", "TPI1", "SARDH", "GZMB", "CREM", "HLA-DMA", "PRDM1",
                         "PKM", "MYO1E", "FUT8", "DUSP4", "RAB27A", "HLA-DRA", "UBE2F-SCLY", "IGFLR1", "CXCL13", "IFNG",
                         "UBE2F", "ITM2A", "ID3", "CD2BP2", "CHST12", "CTSD", "STAT3", "BST2", "CXCR6", "RALGDS", "VCAM1",
                         "TRAFD1", "SYNGR2", "VAPA", "IFI35", "CD63", "NAB1", "PARK7", "MIR4632", "YARS", "PRKAR1A", "IL2RB",
                         "WHRN", "HMGN3", "MTHFD2", "NDFIP2", "PRF1", "PRDX5", "LAG3", "MS4A6A", "LINC00299")

cytotoxicity_signature = c("GZMA", "GZMB", "GZMM", "GZMK", "GZMH", "PRF1", "CD8A", "CD8B")

NK_cell_signature = as.character(NK_cell_BM$HAY_BONE_MARROW_NK_CELLS)

apoptosis_signature = as.character(apoptosis$REACTOME_APOPTOSIS)


# Assign the scores to each cell 
NK_TEMRA <- AddModuleScore(NK_TEMRA, list(exhaustion_signature), 
                        ctrl = length(exhaustion_signature), name = "exhaustion_score")
NK_TEMRA <- AddModuleScore(NK_TEMRA, list(cytotoxicity_signature), 
                        ctrl = length(cytotoxicity_signature), name = "cytotoxicity_score")
NK_TEMRA <- AddModuleScore(NK_TEMRA, list(NK_cell_signature), 
                        ctrl = length(NK_cell_signature), name = "NK_cell_signature_score")
NK_TEMRA <- AddModuleScore(NK_TEMRA, list(apoptosis_signature), 
                        ctrl = length(apoptosis_signature), name = "apoptosis_signature")


# Plot the exhaustion score per cell subtype as violin plot
pdf(file = paste0(outdir, 'NK_TEMRA_subsets_exhaustion_score_violin_plot.pdf'), width = 4, height = 4)
VlnPlot(NK_TEMRA, features = 'exhaustion_score1', pt.size = 0, cols = NK.TEMRA.subset.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  xlab('exhaustion score')
dev.off()


# Plot the cytotoxicity score per cell subtype as violin plot ----
pdf(file = paste0(outdir, 'NK_TEMRA_subsets_cytotoxicity_score_violin_plot.pdf'), width = 4, height = 4)
VlnPlot(NK_TEMRA, features = 'cytotoxicity_score1', pt.size = 0, cols = NK.TEMRA.subset.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  xlab('cytotoxicity score')
dev.off()


# Plot the NK cell signature score per cell subtype as violin plot ----
pdf(file = paste0(outdir, 'NK_TEMRA_subsets_NK_cell_signature_score_violin_plot.pdf'), width = 4, height = 4)
VlnPlot(NK_TEMRA, features = 'NK_cell_signature_score1', pt.size = 0, cols = NK.TEMRA.subset.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  xlab('NK cell signature score')
dev.off()


# Plot the HIV-induced T cell apoptosis signature score per cell subtype as violin plot ----
pdf(file = paste0(outdir, 'NK_TEMRA_subsets_apoptosis_signature_score_violin_plot.pdf'), width = 4, height = 4)
VlnPlot(NK_TEMRA, features = 'apoptosis_signature1', pt.size = 0, cols = NK.TEMRA.subset.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  xlab('apoptosis signature score')
dev.off()


# Plot dotplot with certain NK-cell related genes to characterize the subclusters further ----
NK_features = c('FCGR3A', 'NCAM1', 'KLRD1', 'KLRC1', 'KLRC2', 'KLRF1', 'KLRK1', 'KIR2DL1', 'KIR2DL3', 'KIR3DL1', 'KIR3DL2', 'ITGAE', 'GZMA',
             'GZMB', 'GZMK', 'GNLY', 'PFN1', 'LAMP1', 'CD38', 'CD69', 'HLA-DRA', 'IL7R', 'IKZF2', 'TCF7', 'TBX21', 'EOMES')

pdf(file = paste0(outdir, 'NK_TEMRA_subsets_NK_markers_dotplot.pdf'), width = 6, height = 9)
DotPlot(object = NK_TEMRA, features = NK_features) +
  coord_flip() +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8, hjust = 1),
        axis.ticks = element_blank()) +
  scale_colour_gradient2(low = '#00204D', mid = '#7C7B78', high = '#FFEA46')
dev.off()


# Plot selected gene and ADT marker expression ----
DefaultAssay(NK_TEMRA) = 'ADT'

pdf(file = paste0(outdir, 'NK_TEMRA_CD279_violin_plot.pdf'), width = 4, height = 3)
VlnPlot(NK_TEMRA, features = 'CD279', pt.size = 0, cols = NK.TEMRA.subset.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  ylab('CD279 expression level')
dev.off()

pdf(file = paste0(outdir, 'NK_TEMRA_CD127_violin_plot.pdf'), width = 4, height = 3)
VlnPlot(NK_TEMRA, features = 'CD127', pt.size = 0, cols = NK.TEMRA.subset.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  ylab('CD127 expression level')
dev.off()

pdf(file = paste0(outdir, 'NK_TEMRA_HLA-DR_violin_plot.pdf'), width = 4, height = 3)
VlnPlot(NK_TEMRA, features = 'HLA-DR', pt.size = 0, cols = NK.TEMRA.subset.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  ylab('HLA-DR expression level')
dev.off()

pdf(file = paste0(outdir, 'NK_TEMRA_CD161_violin_plot.pdf'), width = 4, height = 3)
VlnPlot(NK_TEMRA, features = 'CD161', pt.size = 0, cols = NK.TEMRA.subset.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  ylab('CD161 expression level')
dev.off()


# Compare functional scores in CD16+ CD8+ TEMRA-1 cells between mild and severe conditions ----
# Isolate cells from mild and severe groups for each subcluster separately and calculate the scores in these cells per condition
celltypes = as.character(levels(Idents(NK_TEMRA)))

for(i in celltypes){
  
  outdir = paste0('~/sciebo/COVID_project/01 Results/WP5/5.2 characterization of subclusters/scores_', i, '_severe_vs_mild/')
  
  subset_1 = subset(NK_TEMRA, subset = NK_TEMRA_subset_annotations == i)
  Idents(subset_1) = 'condition_collapsed'
  
  mild_vs_severe_subset_1 = subset(subset_1, idents = c('mild', 'severe'))
  
  
  # Extract score data to plot with ggplot
  df = mild_vs_severe_subset_1@meta.data
  
  df = df %>%
    select('condition_collapsed', 'exhaustion_score1', 'cytotoxicity_score1', 'NK_cell_signature_score1', 'apoptosis_signature1')
  
  
  # Plot the exhaustion score for subset one cells as violin plot and compare between mild and severe groups
  pdf(file = paste0(outdir, 'exhaustion_score_violin_plot.pdf'), width = 5, height = 5)
  print(ggplot(df, aes(x = condition_collapsed, y = exhaustion_score1, fill = condition_collapsed)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = 'gray') +
    stat_compare_means() +
    theme_cowplot() +
    theme(text = element_text(size = 8),
          axis.text.x = element_text(angle = 90, size = 8),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = c('#8AFF6B', '#FF5D5D')) +
    xlab('exhaustion score')
  )
  dev.off()
  
  
  # Plot the cytotoxicity score for subset one cells as violin plot and compare between mild and severe groups
  pdf(file = paste0(outdir, 'cytotoxicity_score_violin_plot.pdf'), width = 5, height = 5)
  print(ggplot(df, aes(x = condition_collapsed, y = cytotoxicity_score1, fill = condition_collapsed)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = 'gray') +
    stat_compare_means() +
    theme_cowplot() +
    theme(text = element_text(size = 8),
          axis.text.x = element_text(angle = 90, size = 8),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = c('#8AFF6B', '#FF5D5D')) +
    xlab('cytotoxicity score')
  )
  dev.off()
  
  
  # Plot the NK cell signature score for subset one cells as violin plot and compare between mild and severe groups
  pdf(file = paste0(outdir, 'NK_cell_signature_score_violin_plot.pdf'), width = 5, height = 5)
  print(ggplot(df, aes(x = condition_collapsed, y = NK_cell_signature_score1, fill = condition_collapsed)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = 'gray') +
    stat_compare_means() +
    theme_cowplot() +
    theme(text = element_text(size = 8),
          axis.text.x = element_text(angle = 90, size = 8),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = c('#8AFF6B', '#FF5D5D')) +
    xlab('NK cell signature score')
  )
  dev.off()
  
  
  # Plot the NK cell signature score for subset one cells as violin plot and compare between mild and severe groups
  pdf(file = paste0(outdir, 'apoptosis_score_violin_plot.pdf'), width = 5, height = 5)
  print(ggplot(df, aes(x = condition_collapsed, y = apoptosis_signature1, fill = condition_collapsed)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = 'gray') +
    stat_compare_means() +
    theme_cowplot() +
    theme(text = element_text(size = 8),
          axis.text.x = element_text(angle = 90, size = 8),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = c('#8AFF6B', '#FF5D5D')) +
    xlab('apoptosis score')
  )
  dev.off()
  
}