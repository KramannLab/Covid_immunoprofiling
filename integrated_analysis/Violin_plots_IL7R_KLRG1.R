# Copyright (c) [2022] [Felix Schreibing]
# feschreibing@ukaachen.de

# IL7R and KLRG1 expression in CD8+ T cell populations ----

library(Seurat)


# Load the intgrated sc dataset ----
indir = '~/sciebo/COVID_project/Robjects/'
outdir = '~/sciebo/COVID_project/01 Results/WP1/1.2 cell types and annotation of CD8 T cells/'

sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.revision.filtered.rds'))
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


# Generate a VlnPlot showing expression of IL7R and KLRG1 ----
pdf(paste0(outdir, 'IL7R_VlnPlot.pdf'), width = 9, height = 5)
VlnPlot(sc, features = c('IL7R'), cols = cell.type.colors, pt.size = 0)
dev.off()

pdf(paste0(outdir, 'KLRG1_VlnPlot.pdf'), width = 9, height = 5)
VlnPlot(sc, features = c('KLRG1'), cols = cell.type.colors, pt.size = 0)
dev.off()
