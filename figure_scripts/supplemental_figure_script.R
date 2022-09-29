# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de
# edited by Felix Schreibing

# Supplemental figure scripts ----

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(viridis)


# Load the data (whole dataset with invariant-like T cell subsets) ----
indir = '~/sciebo/COVID_project/Robjects/'
outdir = '~/sciebo/COVID_project/01 Results/WP_Supp/'

sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))

DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'


cell.type.colors = c('CD8+ naive T cells' = '#f58231',
                     'CD8+ central memory T cells' = '#ffe119',
                     'CD8+ CD73+ regulatory T cells' = '#fabebe',
                     'CD8+ TEMRA cells' = '#000075',
                     'CD8+ NK-like TEMRA cells' = '#568198',
                     'CD8+ effector memory T cells 1' = '#70036a',
                     'CD8+ effector memory T cells 2' = '#911eb4',
                     'CD8+ cycling effector T cells' = '#f032e6',
                     'CD8+ NK-like early effector T cells' = '#ff8ff9',
                     'Atypical NKT cells' = '#aaffc3',
                     'CD8+ exhausted T cells' = '#e6194B',
                     'MAIT cells' = '#999999',
                     'Gamma Delta T cells' = '#000000')

# General QC plots ----
patients = c('31', '32',
             '1', '2', '15',
             '5', '6', '7',
             '3', '11', '19',
             '12', '13', '18')


# Plot cell count per sample in bar chart
data = as.data.frame(table(sc$orig.ident))
colnames(data) = c('sample', 'count')
data$sample = as.character(sub('_GEX_SURF', '', data$sample))
data$sample = factor(data$sample, 
                     levels = patients)

p1 = ggplot(data, aes(x = sample, y = log10(count))) +
  geom_bar(aes(fill = sample), 
           alpha = 1, stat = 'identity') +
  scale_fill_viridis(discrete = TRUE) +
  xlab('') + 
  ylab('Number of valid cells (log10)') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 7),
        legend.position = 'none')


# Plot UMI count per sample in violin plot
data = data.frame(cell = names(Idents(sc)), 
                  UMI = as.numeric(sc$nCount_RNA), 
                  gene = as.numeric(sc$nFeature_RNA),
                  sample = sc$orig.ident)
data$sample = as.character(sub('_GEX_SURF', '', data$sample))
data$sample = factor(data$sample, 
                     levels = patients)

p2 = ggplot(data, aes(x = sample, y = log10(UMI))) +
  geom_boxplot(aes(color = sample),
               outlier.size = 0.5) +
  geom_violin(aes(color = sample, fill = sample), 
              alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  xlab('') + 
  ylab('Number of UMIs (log10)') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 7),
        axis.ticks = element_blank(),
        legend.position = 'none')

# Plot gene count per sample in violin plot
p3 = ggplot(data, aes(x = sample, y = log10(gene))) +
  geom_boxplot(aes(color = sample),
               outlier.size = 0.5) +
  geom_violin(aes(color = sample, fill = sample), 
              alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  xlab('') + 
  ylab('Number of genes (log10)') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 7),
        axis.ticks = element_blank(),
        legend.position = 'none')


pdf(file = paste0(outdir, 'per_patient_scRNA_QC.pdf'), width = 4, height = 6)
plot_grid(p1, p2, p3, 
          ncol = 1, 
          align = 'v', 
          rel_heights = c(1, 1, 1.3))
dev.off()


# Plot percentage of mitochondrial reads per condition
Idents(sc) = 'condition_collapsed'

pdf(file = paste0(outdir, 'Percentage_mitochondrial_reads.pdf'), width = 5, height = 5)
VlnPlot(sc, features = 'percent.mt', pt.size = 0) +
  scale_fill_viridis(discrete = TRUE, option = 'D')
dev.off()


# Bar charts with average cell type distribution ----
cell.table = data.frame(cell = colnames(sc), condition = sc$condition_collapsed,
                        cluster = sc$integrated_annotations, 
                        patient = sc$patient)

df = cell.table %>%
  group_by(patient) %>%
  count(cluster, name = 'count', .drop = FALSE) %>%
  add_count(wt = count, name = 'total') %>%
  mutate(fraction = count/total) %>%
  as.data.frame

condition = c(rep('healthy', 2*13),
              rep(c('mild',
                    'severe', 
                    'mild',
                    'severe'), each = 3*13))

df = cbind(df, condition)

df = df %>% 
  group_by(condition, cluster) %>%
  summarize(mean = mean(fraction)) %>%
  as.data.frame

pdf(file = paste0(outdir, 'integrated_Tcells_barchart_celltype_condition_average.pdf'), width = 6, height = 5)
ggplot(df, aes(fill = cluster, y = mean, x = condition)) + 
  geom_bar(stat = 'identity', position = position_fill(reverse = TRUE)) + 
  labs(y = 'Average proportion', x = element_blank(), fill = 'Cell type') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1,
                                   color = 'black'),
        axis.text.y = element_text(color = 'black'),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = cell.type.colors)
dev.off()


# Bar charts with cell type distribution per patient ----
# Change order of patients on X axis according to conditions
cell.table$patient = factor(cell.table$patient, levels = c('31', '32', '1', '2', '15', '5', '6', '7', '3', '11', '19', '12', '13', '18'))

pdf(file = paste0(outdir, 'integrated_Tcells_barchart_celltype_patient.pdf'), width = 7, height = 3)
ggplot(cell.table, aes(fill = cluster, x = patient)) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  labs(y = 'Proportion', x = element_blank(), fill = 'Cell type') +
  theme_classic() +
  theme(axis.text.x = element_text(color = 'black'),
        axis.text.y = element_text(color = 'black'),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = cell.type.colors)
dev.off()


# Plot UMAPs for supplemental materials ----
# Patient UMAP projection
sc$patient = factor(sc$patient, levels = c('31', '32', '1', '2', '15', '5', '6', '7', '3', '11', '19', '12', '13', '18'))

pdf(file = paste0(outdir, 'UMAP_patient_projection.pdf'), width = 5.5, height = 5)
DimPlot(sc, group.by = 'patient') +
  scale_color_viridis(discrete = TRUE, option = 'B')
dev.off()


# Condition UMAP projection
pdf(file = paste0(outdir, 'UMAP_condition_projection.pdf'), width = 6, height = 5)
DimPlot(sc, group.by = 'condition_collapsed') +
  scale_color_viridis(discrete = TRUE, option = 'D')
dev.off()


# Clonotype availability UMAP projection
sc$clonotype_NA = ifelse(is.na(sc$patient_clonotype) == TRUE, 'clonotype not available', 'clonotype available')

Idents(sc) = 'clonotype_NA'

pdf(file = paste0(outdir, 'UMAP_clonotype_availability.pdf'), width = 7, height = 5)
DimPlot(sc) +
  scale_color_viridis(discrete = TRUE, option = 'D')
dev.off()


# Bar chart with TCR V genes in each cluster (show top 20 per cluster) ----
Idents(sc) = 'integrated_annotations'
TCR_outdir = paste0(outdir, 'TRA and TRB gene usage plots/')

pdf(file = paste0(TCR_outdir, 'integrated_Tcells_TCRV_genes.pdf'), width = 12)
for (cell.type in names(cell.type.colors)){
  subset = subset(sc, idents = cell.type)
  
  # Get all TCR V genes
  cell.table = data.frame(cell = colnames(subset), cluster = Idents(subset),
                          tcrv = subset$TCR_V_GENE)
  
  # Get top 20 TCR V genes
  plot.table = cell.table %>% 
    group_by(tcrv) %>% 
    count() %>% 
    arrange(desc(n)) %>% 
    as.data.frame %>% head(20)
  
  print(ggplot(cell.table %>% filter(tcrv %in% plot.table$tcrv), 
               aes(x = tcrv, fill = as.factor(cluster))) + 
          geom_bar() + 
          labs(y = 'Count', x = element_blank(), fill = 'Cell type') +
          ggtitle(paste0(cell.type, ' (n = ', nrow(cell.table), ')')) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_fill_manual(values = cell.type.colors[cell.type]))
}
dev.off()


# Cell-cycle analysis ----
pdf(file = paste0(outdir, 'integrated_cells_Cell_cycle_score.pdf'), width = 4, height = 4)
VlnPlot(sc, features = 'S.Score', 
        pt.size = 0, 
        group.by = 'integrated_annotations', 
        cols = cell.type.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank()) +
  labs(y = 'S cell cycle score') +
  NoLegend() +
  ggtitle('')

VlnPlot(sc, features = 'G2M.Score', 
        pt.size = 0, 
        group.by = 'integrated_annotations', 
        cols = cell.type.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank()) +
  labs(y = 'G2M cell cycle score') +
  NoLegend() +
  ggtitle('')
dev.off()


# ATTENTION: CHANGE DATASET TO FILTERED DATASET
# Load the data (filtered dataset) ----
indir = '~/sciebo/COVID_project/Robjects/'
sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.revision.filtered.rds'))
DefaultAssay(sc) = 'ADT'
Idents(sc) = 'integrated_annotations'


# Plot ADT expression of IL7R (CD127) ----
pdf(file = paste0(outdir, 'CD127_ADT_expression.pdf'), width = 10, height = 7)
VlnPlot(sc, features = 'CD127', pt.size = 0) +
  coord_flip() +
  scale_fill_manual(values = cell.type.colors)
dev.off()
