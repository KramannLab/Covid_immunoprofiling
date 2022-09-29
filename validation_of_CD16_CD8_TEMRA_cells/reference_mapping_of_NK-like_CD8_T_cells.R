# Copyright (c) [2022] [Felix Schreibing]
# feschreibing@ukaachen.de

# Integrating our CD8+ T cell dataset (query) with the reference CD8+ T cell dataset (reference) using Seurat ----
# to verify the existence of a CD16+ CD8+ T cell population in the reference dataset ----

# Source 
# https://satijalab.org/seurat/articles/integration_mapping.html

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(viridis)

memory.limit(17590000000000)

ref.indir = '~/sciebo/COVID_project/Robjects/'
que.indir = '~/sciebo/COVID_project/01 Results/WP5/5.10 validation of subpopulations/large_dataset_CD8_subsets/'
outdir = '~/sciebo/COVID_project/01 Results/WP5/5.10 validation of subpopulations/'


# Load the Seurat objects ----
query = readRDS(paste0(ref.indir, 'integrated.RNA.Tcells.annotated.revision.filtered.rds'))
ref = readRDS(paste0(que.indir, 'CD8_large_dataset_clustered.rds'))

# Add information 'query' and 'reference'
query$dataset = 'query'
ref$dataset = 'reference'

# Adjust the metadata of the reference dataset to match the query dataset
ref$patient = ref$PatientID
ref$condition_collapsed = ref$CoVID.19.severity

ref$PatientID = NULL
ref$CoVID.19.severity = NULL

# Rename identities of reference dataset
ref = RenameIdents(ref, `0` = 'CD8+_1', `1` = 'CD8+_2',
                    `2` = 'CD8+_3', `3` = 'CD8+_4',
                    `4`  = 'CD8+_5', `5` = 'CD8+_6',
                    `6` = 'CD8+_7', `7` = 'CD8+_8',
                    `8` = 'CD8+_9', `9` = 'CD8+_10',
                    `10` = 'CD8+_11', `11` = 'CD8+_12',
                    `12` = 'CD8+_13', `13` = 'CD8+_14',
                    `14` = 'CD8+_15', `15` = 'CD8+_16',
                    `16` = 'CD8+_17', `17` = 'CD8+_18',
                    `18` = 'CD8+_19')

ref$CD8_annotations = Idents(ref)
DimPlot(ref)

# Define colors for CD8+ NK-like TEMRA cell subsets
CD8.colors = c('CD8+_1' = '#3DE356',
               'CD8+_2' = '#57FA37',
               'CD8+_3' = '#8FD744',
               'CD8+_4' = '#D6ED40',
               'CD8+_5' = '#0CE8B6',
               'CD8+_6' = '#00F1FF',
               'CD8+_7' = '#006BFF',
               'CD8+_8' = '#0026F5',
               'CD8+_9' = '#000000',
               'CD8+_10' = '#7B00D1',
               'CD8+_11' = '#AB00C7',
               'CD8+_12' = '#B54FFF',
               'CD8+_13' = '#FF42D2',
               'CD8+_14' = '#FF428B',
               'CD8+_15' = '#F54640',
               'CD8+_16' = '#E85F48',
               'CD8+_17' = '#FF824F',
               'CD8+_18' = '#FFCA0D',
               'CD8+_19' = '#FFFB00')

ref$CD8_annotations = factor(ref$CD8_annotations, 
                             levels = names(CD8.colors))


# Visualize reference dataset as UMAP ----
pdf(file = paste0(outdir, 'UMAP_CD8_reference_dataset.pdf'), width = 7.5)
DimPlot(ref, reduction = "umap", cols = CD8.colors)
dev.off()

# grouped by condition
pdf(file = paste0(outdir, 'UMAP_CD8_reference_dataset_condition.pdf'), width = 8)
DimPlot(ref, reduction = "umap", group.by = 'condition_collapsed') +
  scale_color_viridis(discrete = TRUE, option = 'D')
dev.off()

# grouped by patient
pdf(file = paste0(outdir, 'UMAP_CD8_reference_dataset_patient.pdf'), width = 9.5)
DimPlot(ref, reduction = "umap", group.by = 'patient') +
  scale_color_viridis(discrete = TRUE, option = 'B')
dev.off()


# Integrate the two datasets ----
# Create list of the two Seurat objects to integrate
data.list = list(query, ref)

for (i in 1:length(data.list)) {
  data.list[[i]] <- NormalizeData(data.list[[i]])
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = 2000)
}

# Identify anchors for integration
data.list = list(query, ref)

CD8.anchors = FindIntegrationAnchors(object.list = data.list, dims = 1:30)

# Integrate the data
CD8.integrated <- IntegrateData(anchorset = CD8.anchors, dims = 1:30)

# Save integrated dataset before running harmony ----
saveRDS(CD8.integrated, file = paste0(outdir, 'CD8_integrated_before_harmony.rds'))

# Switch to integrated assay
DefaultAssay(CD8.integrated) <- "integrated"


# Run the standard workflow ----
# Scaling the data
CD8.integrated <- ScaleData(CD8.integrated)

# Perform linear dimensional reduction
CD8.integrated <- RunPCA(CD8.integrated)

# Determine dimensionality of the dataset
ElbowPlot(CD8.integrated, ndims = 50)

# Run non-linear dimensional reduction to check for batch effects
CD8.integrated <- RunUMAP(CD8.integrated, dims = 1:20, reduction = 'pca')
DimPlot(CD8.integrated, reduction = 'umap', group.by = 'patient')


# Run Harmony ----
CD8.integrated.harmony = CD8.integrated %>%
  RunHarmony(group.by.vars = 'patient', plot_convergence = FALSE, assay.use = 'integrated')

CD8.integrated.harmony.embed = Embeddings(CD8.integrated.harmony, 'harmony')
CD8.integrated.harmony.embed


# Run UMAP and clustering using harmony embeddings ----
integrated_CD8_dataset = CD8.integrated.harmony%>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(resolution = 0.4)


# Change subcluster names and define a color scheme ----
# Change subcluster names
integrated_CD8_dataset = RenameIdents(integrated_CD8_dataset, `0` = 'CD8+_T_1', `1` = 'CD8+_T_2',
                             `2` = 'CD8+_T_3', `3` = 'CD8+_T_4',
                             `4`  = 'CD8+_T_5', `5` = 'CD8+_T_6',
                             `6` = 'CD8+_T_7', `7` = 'CD8+_T_8',
                             `8` = 'CD8+_T_9', `9` = 'CD8+_T_10',
                             `10` = 'CD8+_T_11')

integrated_CD8_dataset$CD8_subset_annotations = Idents(integrated_CD8_dataset)

# Define colors for CD8+ NK-like TEMRA cell subsets
CD8.subset.colors = c('CD8+_T_1' = '#3DE356',
                           'CD8+_T_2' = '#8FD744',
                           'CD8+_T_3' = '#00F1FF',
                           'CD8+_T_4' = '#0CA3E8',
                           'CD8+_T_5' = '#0D73FF',
                           'CD8+_T_6' = '#4209BA',
                           'CD8+_T_7' = '#B54FFF',
                           'CD8+_T_8' = '#E848D6',
                           'CD8+_T_9' = '#FF5D5D',
                           'CD8+_T_10' = '#FFCA0D',
                           'CD8+_T_11' = '#FFFA00')

integrated_CD8_dataset$CD8_subset_annotations = factor(integrated_CD8_dataset$CD8_subset_annotations, 
                                                   levels = names(CD8.subset.colors))

DimPlot(integrated_CD8_dataset, cols = CD8.subset.colors)

saveRDS(integrated_CD8_dataset, file = paste0(outdir, 'CD8_integrated_with_harmony.rds'))


# Visualize integrated dataset as UMAP ----
pdf(file = paste0(outdir, 'UMAP_CD8_integrated.pdf'), height = 7, width = 8)
DimPlot(integrated_CD8_dataset, reduction = "umap", cols = CD8.subset.colors, raster = FALSE)
dev.off()

# grouped by condition
pdf(file = paste0(outdir, 'UMAP_CD8_integrated_condition.pdf'), height = 7, width = 8)
DimPlot(integrated_CD8_dataset, reduction = "umap", group.by = 'condition_collapsed') +
  scale_color_viridis(discrete = TRUE, option = 'D')
dev.off()

# grouped by patient
pdf(file = paste0(outdir, 'UMAP_CD8_integrated_patient.pdf'), width = 11)
DimPlot(integrated_CD8_dataset, reduction = "umap", group.by = 'patient') +
  scale_color_viridis(discrete = TRUE, option = 'B')
dev.off()

# grouped by dataset
pdf(file = paste0(outdir, 'UMAP_CD8_integrated_dataset.pdf'), height = 7, width = 8)
DimPlot(integrated_CD8_dataset, reduction = "umap", group.by = 'dataset') +
  scale_color_viridis(discrete = TRUE, option = 'D')
dev.off()


# Subset cells from the query and the reference dataset to plot as separate UMAPs ----
ref_subset = subset(integrated_CD8_dataset, subset = dataset == 'reference')
query_subset = subset(integrated_CD8_dataset, subset = dataset == 'query')

# Plot only cells from reference dataset
pdf(file = paste0(outdir, 'UMAP_CD8_integrated_only_reference_cells.pdf'), height = 7, width = 8)
DimPlot(ref_subset, reduction = "umap", cols = CD8.subset.colors, raster = FALSE)
dev.off()

# Plot only cells from query dataset
pdf(file = paste0(outdir, 'UMAP_CD8_integrated_only_query_cells.pdf'), height = 7, width = 8)
DimPlot(query_subset, reduction = "umap", cols = CD8.subset.colors, raster = FALSE)
dev.off()


# Analyze the proportion of query clusters in the integrated clusters ----
# Create dataframe with proportions
Idents(integrated_CD8_dataset) = 'integrated_annotations'

DimPlot(integrated_CD8_dataset)

test = subset(integrated_CD8_dataset, idents = c('CD8+ effector memory T cells 2', 'CD8+ TEMRA cells', 'CD8+ effector memory T cells 1',
                                          'CD8+ NK-like TEMRA cells', 'CD8+ central memory T cells', 'CD8+ NK-like early effector T cells',
                                          'CD8+ naive T cells', 'CD8+ CD73+ regulatory T cells', 'CD8+ cycling effector T cells',
                                          'CD8+ exhausted T cells'))


df = test@meta.data

df = df %>%
  select('integrated_annotations', 'CD8_subset_annotations') %>%
  group_by(integrated_annotations) %>%
  add_count('integrated_annotations', name = 'total_cells') %>%
  group_by(integrated_annotations, CD8_subset_annotations) %>%
  add_count('CD8_subset_annotations', name = 'count_cells_cluster') %>%
  mutate(frequency = count_cells_cluster/total_cells) %>%
  select('integrated_annotations', 'CD8_subset_annotations', 'frequency') %>%
  unique()

# Create new dataframe for plotting as heatmap
freq = data.frame(integrated_cluster = df$CD8_subset_annotations, original_cluster = df$integrated_annotations, proportion = df$frequency)

freq.1 = data.frame(original_cluster = rep(c('CD8+ effector memory T cells 2', 'CD8+ TEMRA cells', 'CD8+ effector memory T cells 1',
                 'CD8+ NK-like TEMRA cells', 'CD8+ central memory T cells', 'CD8+ NK-like early effector T cells',
                 'CD8+ naive T cells', 'CD8+ CD73+ regulatory T cells', 'CD8+ cycling effector T cells',
                 'CD8+ exhausted T cells'), each = 1*11),
                 integrated_cluster = rep(c('CD8+_T_1', 'CD8+_T_2', 'CD8+_T_3', 'CD8+_T_4', 'CD8+_T_5', 'CD8+_T_6', 'CD8+_T_7',
                                        'CD8+_T_8', 'CD8+_T_9', 'CD8+_T_10', 'CD8+_T_11')))

data = merge(freq, freq.1, by = c('integrated_cluster', 'original_cluster'), all = TRUE)
data[is.na(data)] = 0
data = dcast(data, original_cluster ~ integrated_cluster)
data = column_to_rownames(data, var = 'original_cluster')

# Plot heatmap
pdf(file = paste0(outdir, 'heatmap_ref_vs_query_CD8_populations.pdf'), height = 6, width = 9)
pheatmap(data, border_color = "black")
dev.off()




# # NOW WE KNOW WHICH CELLS IN THE REFERENCE DATASET RESEMBLE OUR CD16+ CD8+ TEMRA POPULATION ----
# NOW WE WANT TO IDENTIFY WHERE IN THE UMAP OF THE REFERENCE CD8 T CELL DATASET THESE NK-LIKE CD8 T CELLS ARE LOCATED ----
# Isolate Cell IDs of all cells in integrated dataset, that cluster with CD8+ NK-like TEMRA cells (cluster name: CD8+_T_1)
outdir = '~/sciebo/COVID_project/01 Results/WP5/5.10 validation of subpopulations/validation_of_CD16_CD8_subpopulations/'

Idents(integrated_CD8_dataset) = 'CD8_subset_annotations'
NK_like_subset = WhichCells(integrated_CD8_dataset, idents = 'CD8+_T_1')

# Identify these cells in the reference dataset based on cell ID
ref$cell.names = WhichCells(ref)
ref$NK_like_subset = ifelse(ref$cell.names %in% NK_like_subset, 'NK-like subset', 'others')

# Plot NK-like CD8+ T cells
Idents(ref) = 'NK_like_subset'

pdf(file = paste0(outdir, 'UMAP_NK-like_CD8_in_reference_dataset.pdf'), width = 8)
DimPlot(ref, cols = c('darkgray','#3DE356'), raster = FALSE)
dev.off()

# Plot FCGR3A expression in NK-like CD8+ T cells compared to all other populations
pdf(file = paste0(outdir, 'FCGR3A_in_reference_NK-like_CD8_cells.pdf'), width = 5, height = 5)
VlnPlot(ref, features = 'FCGR3A', cols = c('darkgray','#3DE356'), pt.size = 0)
dev.off()


# Save the NK-like CD8+ T cells from the reference dataset as a separate Seurat object ----
# Run the standard workflow
# Normalizing the data
large_NK_like_subset = NormalizeData(large_NK_like_subset)

# Find Variable Features
large_NK_like_subset = FindVariableFeatures(large_NK_like_subset)

# Scaling the data
large_NK_like_subset <- ScaleData(large_NK_like_subset)

# Perform linear dimensional reduction
large_NK_like_subset <- RunPCA(large_NK_like_subset)

# Determine dimensionality of the dataset
ElbowPlot(large_NK_like_subset)

# Run non-linear dimensional reduction to observe batch effects
large_NK_like_subset <- RunUMAP(large_NK_like_subset, dims = 1:20, reduction = 'pca')
DimPlot(large_NK_like_subset, reduction = 'umap', group.by = 'PatientID')


# Run Harmony
large_NK_like_subset.harmony = large_NK_like_subset %>%
  RunHarmony(group.by.vars = 'PatientID', plot_convergence = FALSE)

large_NK_like_subset.harmony.embed = Embeddings(large_NK_like_subset.harmony, 'harmony')
large_NK_like_subset.harmony.embed


# Run UMAP and clustering using harmony embeddings
ref_NK_like_TEMRA = large_NK_like_subset.harmony%>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(resolution = 1)

DimPlot(ref_NK_like_TEMRA)

saveRDS(ref_NK_like_TEMRA, paste0(outdir, 'large_NK_like_dataset_clustered.rds'))


# Rename identities of reference dataset
ref_NK_like_TEMRA = RenameIdents(ref_NK_like_TEMRA, `0` = 'CD16+_CD8+_1', `1` = 'CD16+_CD8+_2',
                   `2` = 'CD16+_CD8+_3', `3` = 'CD16+_CD8+_4',
                   `4`  = 'CD16+_CD8+_5', `5` = 'CD16+_CD8+_6',
                   `6` = 'CD16+_CD8+_7', `7` = 'CD16+_CD8+_8',
                   `8` = 'CD16+_CD8+_9', `9` = 'CD16+_CD8+_10',
                   `10` = 'CD16+_CD8+_11', `11` = 'CD16+_CD8+_12',
                   `12` = 'CD16+_CD8+_13', `13` = 'CD16+_CD8+_14',
                   `14` = 'CD16+_CD8+_15')

ref_NK_like_TEMRA$CD16_CD8_annotations = Idents(ref_NK_like_TEMRA)
DimPlot(ref_NK_like_TEMRA)

# Define colors for CD8+ NK-like TEMRA cell subsets
CD16.colors = c('CD16+_CD8+_1' = '#3DE356',
               'CD16+_CD8+_2' = '#57FA37',
               'CD16+_CD8+_3' = '#D6ED40',
               'CD16+_CD8+_4' = '#0CE8B6',
               'CD16+_CD8+_5' = '#00F1FF',
               'CD16+_CD8+_6' = '#0026F5',
               'CD16+_CD8+_7' = '#000000',
               'CD16+_CD8+_8' = '#7B00D1',
               'CD16+_CD8+_9' = '#AB00C7',
               'CD16+_CD8+_10' = '#B54FFF',
               'CD16+_CD8+_11' = '#FF42D2',
               'CD16+_CD8+_12' = '#F54640',
               'CD16+_CD8+_13' = '#FF824F',
               'CD16+_CD8+_14' = '#FFCA0D',
               'CD16+_CD8+_15' = '#FFFB00')

ref_NK_like_TEMRA$CD16_CD8_annotations = factor(ref_NK_like_TEMRA$CD16_CD8_annotations, 
                             levels = names(CD16.colors))


# Visualize reference dataset as UMAP ----
pdf(file = paste0(outdir, 'UMAP_NK_like_reference_subset_dataset.pdf'), width = 7.5)
DimPlot(ref_NK_like_TEMRA, reduction = "umap", cols = CD16.colors)
dev.off()

# grouped by condition
pdf(file = paste0(outdir, 'UMAP_NK_like_reference_subset_dataset_condition.pdf'), width = 8)
DimPlot(ref_NK_like_TEMRA, reduction = "umap", group.by = 'CoVID.19.severity') +
  scale_color_viridis(discrete = TRUE, option = 'D')
dev.off()

# grouped by patient
pdf(file = paste0(outdir, 'UMAP_NK_like_reference_subset_dataset_patient.pdf'), width = 9.5)
DimPlot(ref_NK_like_TEMRA, reduction = "umap", group.by = 'PatientID') +
  scale_color_viridis(discrete = TRUE, option = 'B')
dev.off()