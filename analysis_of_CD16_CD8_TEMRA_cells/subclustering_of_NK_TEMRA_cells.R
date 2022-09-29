# Copyright (c) [2022] [Felix Schreibing]
# feschreibing@ukaachen.de

# Subclustering of CD16+ CD8+ TEMRA cells (previously termed CD8+ NK-like TEMRA cells)

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)


# Load the intgrated sc dataset ----
indir = '~/sciebo/COVID_project/Robjects/'
outdir = '~/sciebo/COVID_project/01 Results/WP5/5.1 subclustering of TEMRA and NK-like TEMRA cells/CD8+ NK-like TEMRA cells/'

sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.revision.filtered.rds'))


# Subset the CD8+ NK-like TEMRA cells ----
subset_1 = subset(sc, subset = integrated_annotations == 'CD8+ NK-like TEMRA cells')

# Run the standard Seurat workflow ----
# Normalizing the data
subset_1 = NormalizeData(subset_1)

# Find Variable Features
subset_1 = FindVariableFeatures(subset_1)

# Scaling the data
subset_1 <- ScaleData(subset_1)

# Perform linear dimensional reduction
subset_1 <- RunPCA(subset_1)

# Determine dimensionality of the dataset
ElbowPlot(subset_1)


# Run non-linear dimensional reduction to check for batch effects ----
subset_1 <- RunUMAP(subset_1, dims = 1:20, reduction = 'pca')
DimPlot(subset_1, reduction = 'umap', group.by = 'patient')


# Run Harmony to re-integrate subset ----
subset1.harmony = subset_1 %>%
  RunHarmony(group.by.vars = 'patient', plot_convergence = FALSE)

subset1.harmony.embed = Embeddings(subset1.harmony, 'harmony')
subset1.harmony.embed


# Run UMAP and clustering using harmony embeddings ----
NK_like_TEMRA = subset1.harmony%>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(resolution = 1)

DimPlot(NK_like_TEMRA)


# Determine the number of cells per cluster to potentially change the resolution
df = NK_like_TEMRA@meta.data

df %>%
  group_by(seurat_clusters) %>%
  summarize(count = n())


# Change subcluster names and define a color scheme ----
# Change subcluster names
NK_like_TEMRA = RenameIdents(NK_like_TEMRA, `0` = 'CD16+_CD8+_TEMRA_1', `1` = 'CD16+_CD8+_TEMRA_2',
                             `2` = 'CD16+_CD8+_TEMRA_3', `3` = 'CD16+_CD8+_TEMRA_4',
                             `4`  = 'CD16+_CD8+_TEMRA_5', `5` = 'CD16+_CD8+_TEMRA_6')

NK_like_TEMRA$NK_TEMRA_subset_annotations = Idents(NK_like_TEMRA)

# Define colors for CD8+ NK-like TEMRA cell subsets
NK.TEMRA.subset.colors = c('CD16+_CD8+_TEMRA_1' = '#8FD744',
                           'CD16+_CD8+_TEMRA_2' = '#9700AD',
                           'CD16+_CD8+_TEMRA_3' = '#00F1FF',
                           'CD16+_CD8+_TEMRA_4' = '#FFFA00',
                           'CD16+_CD8+_TEMRA_5' = '#FF5D5D',
                           'CD16+_CD8+_TEMRA_6' = '#000000')

NK_like_TEMRA$NK_TEMRA_subset_annotations = factor(NK_like_TEMRA$NK_TEMRA_subset_annotations, 
                                                   levels = names(NK.TEMRA.subset.colors))

DimPlot(NK_like_TEMRA, cols = NK.TEMRA.subset.colors)

saveRDS(NK_like_TEMRA, file = paste0(outdir, 'NK-like_TEMRA_subclusters.rds'))


# Visualize CD8+ NK-like TEMRA cell subsets as UMAP ----
pdf(file = paste0(outdir, 'UMAP_NK_TEMRA_subsets.pdf'), width = 9)
DimPlot(NK_like_TEMRA, reduction = "umap", pt.size = 2, cols = NK.TEMRA.subset.colors)
dev.off()

DimPlot(NK_like_TEMRA, group.by = 'clonetype_cut')


# Visualize cell subset distribution per condition and per patient ----
cell.table = data.frame(cell = colnames(NK_like_TEMRA), condition = NK_like_TEMRA$condition_collapsed,
                        cluster = Idents(NK_like_TEMRA), 
                        patient = NK_like_TEMRA$patient)

df = cell.table %>%
  group_by(patient) %>%
  count(cluster, name = 'count', .drop = FALSE) %>%
  add_count(wt = count, name = 'total') %>%
  mutate(fraction = count/total) %>%
  as.data.frame

condition = c(rep('healthy', 2*6),
              rep(c('mild',
                    'severe', 
                    'mild',
                    'severe'), each = 3*6))

df = cbind(df, condition)

df = df %>% 
  group_by(condition, cluster) %>%
  summarize(mean = mean(fraction)) %>%
  as.data.frame


# Plot subtype distribution per condition
pdf(file = paste0(outdir, 'NK_TEMRA_barchart_per_condition.pdf'), width = 5, height = 5)
ggplot(df, aes(fill = cluster, y = mean, x = condition)) + 
  geom_bar(stat = 'identity', position = position_fill(reverse = TRUE)) + 
  labs(y = 'Average proportion', x = element_blank(), fill = 'Cell type') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1,
                                   color = 'black'),
        axis.text.y = element_text(color = 'black'),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = NK.TEMRA.subset.colors)
dev.off()


# Plot subtype distribution per patient
cell.table$patient = factor(cell.table$patient, levels = c('31', '32', '1', '2', '15', '5', '6', '7', '3', '11', '19', '12', '13', '18'))

pdf(file = paste0(outdir, 'NK_TEMRA_barchart_per_patient.pdf'), width = 9, height = 5)
ggplot(cell.table, aes(fill = cluster, x = patient)) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  labs(y = 'Proportion', x = element_blank(), fill = 'Cell type') +
  theme_classic() +
  theme(axis.text.x = element_text(color = 'black'),
        axis.text.y = element_text(color = 'black'),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = NK.TEMRA.subset.colors)
dev.off()


# Visualize TCR clonality of CD8+ NK-like TEMRA cell subsets as barplot ----
group.order = c(NA, rev(c('Hyperexpanded (100 < x <= Inf)', 
                          'Large (20 < x <= 100)', 
                          'Medium (5 < x <= 20)', 
                          'Small (1 < x <= 5)', 
                          'Single (0 < x <= 1)')))

df2 = NK_like_TEMRA@meta.data

pdf(file = paste0(outdir, 'barplot_expansion_per_NK_TEMRA_subtype.pdf'), width = 6, height = 4)
ggplot(df2, aes(x = NK_TEMRA_subset_annotations, 
                fill = factor(clonetype_cut, levels = group.order, exclude = NULL))) + 
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

pdf(file = paste0(outdir, 'barplot_expansion_per_condition_in_NK_TEMRA.pdf'), width = 5, height = 3)
ggplot(df2, aes(x = condition_collapsed,
                fill = factor(clonetype_cut, levels = group.order, exclude = NULL))) + 
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