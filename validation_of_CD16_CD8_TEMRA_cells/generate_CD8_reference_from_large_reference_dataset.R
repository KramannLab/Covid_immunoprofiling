# Copyright (c) [2020] [Felix Schreibing]
# feschreibing@ukaachen.de

# Use large reference COVID-19 dataset
# from https://www.cell.com/cell/fulltext/S0092-8674(21)00148-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421001483%3Fshowall%3Dtrue ----
# Subset all CD8+ T cells in samples from frozen PBMCs and 5´-seuquencing to generate a reference CD8+ T cell dataset ----
# to validate the existence of a CD16+ CD8+ TEMRA population ----

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(tidyverse)
library(ggpubr)
library(readxl)


# Load the data ----
indir = '~/sciebo/COVID_project/01 Results/WP5/5.7 validation of signature score/large_COVID_dataset_CELL/'
outdir = '~/sciebo/COVID_project/01 Results/WP5/5.10 validation of subpopulations/'

memory.limit(17592100000000)


# Generate a large CD8+ T cell reference dataset ----
# Load the large reference dataset files and subset for CD8+ T cells
files = c('selected030.h', 'selected3060.h', 'selected6090.h', 'selected90120.h',
          'selected120150.h', 'selected150180.h', 'selected180196.h')


for(file in files){
  
  name = as.character(file)
  
  file = readRDS(paste0(indir, file, '.RDS'))
  
  DefaultAssay(file) = 'RNA'
  
  Idents(file) = 'majorType'
  
  CD8 = subset(file, idents = 'CD8')
  
  saveRDS(CD8, paste0(outdir, 'large_dataset_CD8_subsets/', name, '.CD8.rds') )
  
}

# Combine the split CD8 datasets into larger CD8 datasets
CD8_1 = readRDS(paste0(outdir, 'large_dataset_CD8_subsets/', 'selected030.h.CD8.rds'))
CD8_2 = readRDS(paste0(outdir, 'large_dataset_CD8_subsets/', 'selected3060.h.CD8.rds'))

data = merge(CD8_1, y = CD8_2, project = "large_CD8_dataset")
saveRDS(data, paste0(outdir, 'large_dataset_CD8_subsets/CD8_1_2.rds'))

CD8_3 = readRDS(paste0(outdir, 'large_dataset_CD8_subsets/', 'selected6090.h.CD8.rds'))
CD8_4 = readRDS(paste0(outdir, 'large_dataset_CD8_subsets/', 'selected90120.h.CD8.rds'))

data = merge(CD8_3, y = CD8_4, project = "large_CD8_dataset")
saveRDS(data, paste0(outdir, 'large_dataset_CD8_subsets/CD8_3_4.rds'))

CD8_5 = readRDS(paste0(outdir, 'large_dataset_CD8_subsets/', 'selected120150.h.CD8.rds'))
CD8_6 = readRDS(paste0(outdir, 'large_dataset_CD8_subsets/', 'selected150180.h.CD8.rds'))
CD8_7 = readRDS(paste0(outdir, 'large_dataset_CD8_subsets/', 'selected180196.h.CD8.rds'))

data = merge(CD8_5, y = c(CD8_6, CD8_7), project = "large_CD8_dataset")
saveRDS(data, paste0(outdir, 'large_dataset_CD8_subsets/CD8_5_6_7.rds'))


# Combine the large CD8 datasets into one large CD8 reference dataset
CD8_1_2 = readRDS(paste0(outdir, 'large_dataset_CD8_subsets/', 'CD8_1_2.rds'))
CD8_3_4 = readRDS(paste0(outdir, 'large_dataset_CD8_subsets/', 'CD8_3_4.rds'))

data = merge(CD8_1_2, y = CD8_3_4, project = "large_CD8_dataset")
saveRDS(data, paste0(outdir, 'large_dataset_CD8_subsets/CD8_1_2_3_4.rds'))

CD8_1_2_3_4 = readRDS(paste0(outdir, 'large_dataset_CD8_subsets/', 'CD8_1_2_3_4.rds'))
CD8_5_6_7 = readRDS(paste0(outdir, 'large_dataset_CD8_subsets/', 'CD8_5_6_7.rds'))

sc = merge(CD8_1_2_3_4, y = CD8_5_6_7, project = "large_CD8_dataset")
saveRDS(sc, paste0(outdir, 'large_dataset_CD8_subsets/CD8_large_dataset.rds'))


# Filter for only 5´libraries and for frozen PBMC samples ----
data.dir = '~/sciebo/COVID_project/01 Results/WP5/5.10 validation of subpopulations/large_dataset_CD8_subsets/'
sc = readRDS(paste0(data.dir, 'CD8_large_dataset.rds'))

subset = subset(sc, subset = Single.cell.sequencing.platform == "10X 5'")
subset = subset(subset, subset = Sample.type == 'frozen PBMC')


# Re-integrate the CD8+ reference dataset and perform clustering ----
# Run the standard Seurat workflow
# Normalizing the data
subset = NormalizeData(subset)

# Find Variable Features
subset = FindVariableFeatures(subset)

# Scaling the data
subset <- ScaleData(subset)

# Perform linear dimensional reduction
subset <- RunPCA(subset)

# Determine dimensionality of the dataset
ElbowPlot(subset)

# Run non-linear dimensional reduction and check for batch effects
subset <- RunUMAP(subset, dims = 1:20, reduction = 'pca')
DimPlot(subset, reduction = 'umap', group.by = 'PatientID')


# Run Harmony to re-integrate the data
subset.harmony = subset %>%
  RunHarmony(group.by.vars = 'PatientID', plot_convergence = FALSE)

subset.harmony.embed = Embeddings(subset.harmony, 'harmony')
subset.harmony.embed


# Run UMAP and clustering using harmony embeddings
CD8 = subset.harmony%>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(resolution = 1)

DimPlot(CD8, group.by = 'PatientID')
DimPlot(CD8, group.by = 'CoVID.19.severity')

saveRDS(CD8, paste0(outdir, 'large_dataset_CD8_subsets/CD8_large_dataset_clustered.rds'))