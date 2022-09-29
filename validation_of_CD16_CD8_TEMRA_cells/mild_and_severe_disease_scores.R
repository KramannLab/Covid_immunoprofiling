# Copyright (c) [2022] [Felix Schreibing]
# feschreibing@ukaachen.de

# Define disease severity signatures in CD16+ CD8+ T cells ----

library(Seurat)
library(harmony)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(viridis)

memory.limit(17590000000000)


# Load the data ----
que.indir = '~/sciebo/COVID_project/Robjects/'
ref.indir = '~/sciebo/COVID_project/01 Results/WP5/5.10 validation of subpopulations/validation_of_CD16_CD8_subpopulations/'
marker.dir = '~/sciebo/COVID_project/01 Results/WP2/2.2 NK-like TEMRA cell marker genes/'
outdir = '~/sciebo/COVID_project/01 Results/WP5/5.10 validation of subpopulations/'

query = readRDS(paste0(que.indir, 'NK-like_TEMRA_subclusters.rds'))
ref = readRDS(paste0(ref.indir, 'large_NK_like_dataset_clustered.rds'))


# Perform DEG analysis in CD16+ CD8+ TEMRA cells from our dataset ----
# between mild and severe conditions ----
Idents(ref) = 'CoVID.19.severity'
Idents(query) = 'condition_collapsed'

DEGs_quer = FindMarkers(query, ident.1 = 'mild', ident.2 = 'severe', min.pct = 0.25)
DEGs_quer_export = rownames_to_column(DEGs_quer, var = 'gene')
write_xlsx(DEGs_quer_export, paste0(outdir, 'DEGs_mild_vs_severe_all_CD16_CD8_TEMRA_cells.xlsx'))

# Load genes differentially expressed in NK-like TEMRA cells compared to all other cell types in our dataset
NK_TEMRA_markers = read.csv(paste0(marker.dir, 'NK_like_TEMRA_marker_genes.txt'))


# Generate a score for each of the two disease conditions including ----
# all significant DEGs between mild vs. severe that intersect with highly significant NK-TEMRA markers ----
NK_TEMRA_sig = NK_TEMRA_markers %>%
  filter(p_val_adj < 0.00000000000000001) %>%
  filter(avg_log2FC > 0.25)

DEGs_quer_sig = DEGs_quer %>%
  filter(p_val_adj < 0.05) %>%
  filter(avg_log2FC > 0.25 | avg_log2FC < -0.25) %>%
  rownames_to_column('gene')

intersecting_genes = DEGs_quer_sig %>%
  filter(gene %in% NK_TEMRA_sig$X)

# generate a mild signature score and a severe signature score
mild_genes = intersecting_genes %>%
  filter(avg_log2FC > 0) %>%
  filter(p_val_adj < 0.05)

severe_genes = intersecting_genes %>%
  filter(avg_log2FC < 0) %>%
  filter(p_val_adj < 0.05)

mild_score = mild_genes$gene
severe_score = severe_genes$gene


# Calculate mild and severe disease scores for every cell in the query and the reference dataset ----
query = AddModuleScore(query, list(mild_score), ctrl = length(mild_score), name = 'mild_score')
query = AddModuleScore(query, list(severe_score), ctrl = length(severe_score), name = 'severe_score')
ref = AddModuleScore(ref, list(mild_score), ctrl = length(mild_score), name = 'mild_score')
ref = AddModuleScore(ref, list(severe_score), ctrl = length(severe_score), name = 'severe_score')


# Plots for query dataset ----
# Plot mild and severe disease scores for query dataset 
pdf(paste0(outdir, 'mild_disease_score_Feature_plot_query.pdf'), width = 4, height = 4)
FeaturePlot(query, features = 'mild_score1', pt.size = 2) +
  scale_color_viridis(option = 'B')
dev.off()

pdf(paste0(outdir, 'severe_disease_score_Feature_plot_query.pdf'), width = 4, height = 4)
FeaturePlot(query, features = 'severe_score1', pt.size = 2) +
  scale_color_viridis(option = 'B')
dev.off()


# Plot score values per condition as boxplot for query dataset...
df_que = query@meta.data %>%
  select(patient, condition_collapsed, mild_score1, severe_score1)

pairwise_comparisons1 = list(c('healthy', 'mild'), c('mild', 'severe'), c('healthy', 'severe'))

# ...for mild score
pdf(paste0(outdir, 'mild_disease_score_boxplot_query.pdf'), width = 7, height = 7)
ggplot(df_que, aes(x = condition_collapsed, y = mild_score1, fill = condition_collapsed)) +
  geom_boxplot(width = 0.7, lwd = 1) +
  stat_compare_means(comparisons = pairwise_comparisons1, label.y = c(1.6, 1.7, 1.8)) +
  stat_compare_means(label.y = 2.1) + # for global Anova p-value
  theme_cowplot() +
  ggtitle('CD16+ CD8+ mild score') +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  scale_fill_viridis(discrete = TRUE, option = 'B')
dev.off()

# ...for severe score
pdf(paste0(outdir, 'severe_disease_score_boxplot_query.pdf'), width = 7, height = 7)
ggplot(df_que, aes(x = condition_collapsed, y = severe_score1, fill = condition_collapsed)) +
  geom_boxplot(width = 0.7, lwd = 1) +
  stat_compare_means(comparisons = pairwise_comparisons1, label.y = c(1.8, 2, 2.2)) +
  stat_compare_means(label.y = 2.5) + # for global Anova p-value
  theme_cowplot() +
  ggtitle('CD16+ CD8+ severe score') +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  scale_fill_viridis(discrete = TRUE, option = 'B')
dev.off()


# Plots for reference dataset ----
# Plot mild and severe disease scores for reference dataset
pdf(paste0(outdir, 'mild_disease_score_Feature_plot_ref.pdf'), width = 5, height = 5)
FeaturePlot(ref, features = 'mild_score1', pt.size = 0.5) +
  scale_color_viridis(option = 'B')
dev.off()

pdf(paste0(outdir, 'severe_disease_score_Feature_plot_ref.pdf'), width = 5, height = 5)
FeaturePlot(ref, features = 'severe_score1', pt.size = 0.5) +
  scale_color_viridis(option = 'B')
dev.off()


# Plot score values per condition as boxplot for reference dataset...
# Order conditions of reference dataset in the way they are ordered in query dataset
Idents(ref) = 'CoVID.19.severity'
levels = c('control', 'mild/moderate', 'severe/critical')
Idents(ref) = factor(Idents(ref), levels = levels)

df_ref = ref@meta.data %>%
  select(PatientID, CoVID.19.severity, Outcome, mild_score1, severe_score1, medication = COVID.19.related.medication.and.anti.microbials)

# GROUPED BY CONDITON
# ...for mild score
pairwise_comparisons2 = list(c('control', 'mild/moderate'), c('mild/moderate', 'severe/critical'), c('control', 'severe/critical'))

pdf(paste0(outdir, 'mild_disease_score_boxplot_reference_condition.pdf'), width = 7, height = 7)
ggplot(df_ref, aes(x = CoVID.19.severity, y = mild_score1, fill = CoVID.19.severity)) +
  geom_boxplot(width = 0.7, lwd = 1) +
  stat_compare_means(comparisons = pairwise_comparisons2, label.y = c(1.6, 1.7, 1.8)) +
  stat_compare_means(label.y = 2.1) + # for global Anova p-value
  theme_cowplot() +
  ggtitle('CD16+ CD8+ mild score') +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  scale_fill_viridis(discrete = TRUE, option = 'B')
dev.off()

# ...for severe score
pdf(paste0(outdir, 'severe_disease_score_boxplot_reference_condition.pdf'), width = 7, height = 7)
ggplot(df_ref, aes(x = CoVID.19.severity, y = severe_score1, fill = CoVID.19.severity)) +
  geom_boxplot(width = 0.7, lwd = 1) +
  stat_compare_means(comparisons = pairwise_comparisons2, label.y = c(1.8, 2.0, 2.2)) +
  stat_compare_means(label.y = 2.6) + # for global Anova p-value
  theme_cowplot() +
  ggtitle('CD16+ CD8+ severe score') +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  scale_fill_viridis(discrete = TRUE, option = 'B')
dev.off()



# GROUPED BY OUTCOME
# ...for mild score
df_ref$Outcome = factor(df_ref$Outcome, levels = c('control', 'discharged', 'deceased'))
pairwise_comparisons3 = list(c('control', 'discharged'), c('discharged', 'deceased'), c('control', 'deceased'))

pdf(paste0(outdir, 'mild_disease_score_boxplot_reference_outcome.pdf'), width = 7, height = 7)
ggplot(df_ref, aes(x = Outcome, y = mild_score1, fill = Outcome)) +
  geom_boxplot(width = 0.7, lwd = 1) +
  stat_compare_means(comparisons = pairwise_comparisons3, label.y = c(1.6, 1.7, 1.8)) +
  stat_compare_means(label.y = 2.1) + # for global Anova p-value
  theme_cowplot() +
  ggtitle('CD16+ CD8+ mild score') +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  scale_fill_viridis(discrete = TRUE, option = 'B')
dev.off()

# ...for severe score
pdf(paste0(outdir, 'severe_disease_score_boxplot_reference_outcome.pdf'), width = 7, height = 7)
ggplot(df_ref, aes(x = Outcome, y = severe_score1, fill = Outcome)) +
  geom_boxplot(width = 0.7, lwd = 1) +
  stat_compare_means(comparisons = pairwise_comparisons3, label.y = c(1.8, 2.0, 2.2)) +
  stat_compare_means(label.y = 2.6) + # for global Anova p-value
  theme_cowplot() +
  ggtitle('CD16+ CD8+ severe score') +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, size = 8),
        axis.title.x = element_blank()) +
  scale_fill_viridis(discrete = TRUE, option = 'B')
dev.off()


# Plot scores for all CD8+ T cell types in the whole query and reference dataset ----
# Whole query
sc_query = readRDS(paste0(que.indir, 'integrated.RNA.Tcells.annotated.revision.filtered.rds'))

sc_query = AddModuleScore(sc_query, list(mild_score), ctrl = length(mild_score), name = 'mild_score')
sc_query = AddModuleScore(sc_query, list(severe_score), ctrl = length(severe_score), name = 'severe_score')

pdf(paste0(outdir, 'mild_disease_score_whole_query.pdf'), width = 5, height = 5)
FeaturePlot(sc_query, features = 'mild_score1') +
  scale_color_viridis(option = 'B')
dev.off()

pdf(paste0(outdir, 'severe_disease_score_whole_query.pdf'), width = 5, height = 5)
FeaturePlot(sc_query, features = 'severe_score1') +
  scale_color_viridis(option = 'B')
dev.off()

# Plot score values in query dataset per celltype as VlnPlot
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

Idents(sc_query) = 'integrated_annotations'

pdf(paste0(outdir, 'mild_disease_score_whole_query_violin_plot.pdf'), width = 8, height = 7)
VlnPlot(sc_query, features = 'mild_score1', cols = cell.type.colors, pt.size = 0)
dev.off()

pdf(paste0(outdir, 'severe_disease_score_whole_query_violin_plot.pdf'), width = 8, height = 7)
VlnPlot(sc_query, features = 'severe_score1', cols = cell.type.colors, pt.size = 0)
dev.off()


# Whole reference
sc_ref = readRDS('~/sciebo/COVID_project/01 Results/WP5/5.10 validation of subpopulations/large_dataset_CD8_subsets/CD8_large_dataset_clustered.rds')

sc_ref = AddModuleScore(sc_ref, list(mild_score), ctrl = length(mild_score), name = 'mild_score')
sc_ref = AddModuleScore(sc_ref, list(severe_score), ctrl = length(severe_score), name = 'severe_score')

pdf(paste0(outdir, 'mild_disease_score_whole_reference.pdf'), width = 5, height = 5)
FeaturePlot(sc_ref, features = 'mild_score1') +
  scale_color_viridis(option = 'B')
dev.off()

pdf(paste0(outdir, 'severe_disease_score_whole_reference.pdf'), width = 5, height = 5)
FeaturePlot(sc_ref, features = 'severe_score1') +
  scale_color_viridis(option = 'B')
dev.off()

# Plot score values in reference dataset per celltype as VlnPlot
NK_like = WhichCells(ref)

sc_ref$cell_name = WhichCells(sc_ref)
sc_ref$NK_like = ifelse(sc_ref$cell_name %in% NK_like, 'NK-like cells', 'others')

Idents(sc_ref) = 'NK_like'

pdf(paste0(outdir, 'mild_disease_score_whole_reference_violin_plot.pdf'), width = 5, height = 7)
VlnPlot(sc_ref, features = 'mild_score1', cols = c('darkgray','#3DE356'), pt.size = 0)
dev.off()

pdf(paste0(outdir, 'severe_disease_score_whole_reference_violin_plot.pdf'), width = 5, height = 7)
VlnPlot(sc_ref, features = 'severe_score1', cols = c('darkgray','#3DE356'), pt.size = 0)
dev.off()