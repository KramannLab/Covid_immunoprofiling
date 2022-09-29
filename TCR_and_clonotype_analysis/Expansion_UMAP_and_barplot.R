# Clonal Expansion UMAP projection and barplot
# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de
# edited by Felix Schreibing

# T cell receptor clonality UMAP projection and barplot ----

library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(viridis)
library(ggplot2)


# Load the intgrated sc dataset ----
indir = '~/sciebo/COVID_project/Robjects/'
outdir = '~/sciebo/COVID_project/01 Results/WP3/3.1 Clonality UMAP and barplot/'

sc = readRDS(paste0(indir, 'integrated.RNA.Tcells.annotated.revision.filtered.rds'))

DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'


# Plot clonotype expansion groups on UMAP ----
# Define the clonotype groups according to their size
sc$clonotype_cut = cut(sc$clonotype_size, breaks = c(0, 1, 5, 20, 100, Inf), 
                       labels = rev(c('Hyperexpanded (100 < x <= Inf)', 
                                      'Large (20 < x <= 100)', 
                                      'Medium (5 < x <= 20)', 
                                      'Small (1 < x <= 5)', 
                                      'Single (0 < x <= 1)')))

# Order of groups for plot
group.order = c(NA, rev(c('Hyperexpanded (100 < x <= Inf)', 
                          'Large (20 < x <= 100)', 
                          'Medium (5 < x <= 20)', 
                          'Small (1 < x <= 5)', 
                          'Single (0 < x <= 1)')))

pdf(file = paste0(outdir, 'clonal_expansion_UMAP.pdf'), width = 9.8)
DimPlot(sc, group.by = 'clonotype_cut', 
        order = TRUE,
        cols = viridis(5))
dev.off()


# Bar charts of clonotype expansion groups separated ----
sc$celltype.condition = paste0(sc$integrated_annotations, '.', sc$condition_collapsed)

df = sc@meta.data
df[is.na(df$Clonotype),'clonotype_cut'] = NA


# Plot clonotype groups per cell type
pdf(file = paste0(outdir, 'Barplot_expansion_per_cell_type_all_cells.pdf'), width = 6, height = 4)
ggplot(df, aes(x = integrated_annotations, 
               fill = factor(clonotype_cut, levels = group.order, exclude = NULL))) + 
  geom_bar(position = 'fill') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 1, hjust = 1, color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = 'black'),
        axis.ticks = element_blank()) +
  scale_fill_viridis(discrete = TRUE, na.value = 'lightgrey') +
  labs(fill = 'Clonal expansion',
       y =  'Relative abundance')
dev.off()


# Plot clonotype groups per condition
pdf(file = paste0(outdir, 'Barplot_expansion_per_condition_all_cells.pdf'), width = 5, height = 4)
ggplot(df, aes(x = condition_collapsed,
               fill = factor(clonotype_cut, levels = group.order, exclude = NULL))) + 
  geom_bar(position = 'fill') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 1, hjust = 1, color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = 'black'),
        axis.ticks = element_blank()) +
  scale_fill_viridis(discrete = TRUE, na.value = 'lightgrey') +
  labs(fill = 'Clonal expansion',
       y =  'Relative abundance')
dev.off()


# Plot clonotype groups per celltype separately for each condition
pdf(file = paste0(outdir, 'Barplot_expansion_per_cell_type_separated_by_condition.pdf'), width = 9, height = 4.5)
ggplot(df, aes(x = celltype.condition,
               fill = factor(clonotype_cut, levels = group.order, exclude = NULL))) + 
  geom_bar(position = 'fill') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 1, hjust = 1, color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = 'black'),
        axis.ticks = element_blank()) +
  scale_fill_viridis(discrete = TRUE, na.value = 'lightgrey') +
  labs(fill = 'Clonal expansion',
       y =  'Relative abundance')
dev.off()