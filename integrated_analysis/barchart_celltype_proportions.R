# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de
# adapted by Felix Schreibing
# feschreibing@ukaachen.de


# Barplots cell type distribution per condition and patient ----

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

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

# Load the data ----
indir = '~/sciebo/COVID_project/Robjects/'
outdir = '~/sciebo/COVID_project/Robjects/01 Results/WP1/1.2 cell types and annotation of CD8 T cells/'

sc = readRDS(paste0(indir, 'integrated.RNA.Tcells.annotated.revision.filtered.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'

# Bar charts with average cell type distribution per condition ----
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