# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de
# edited by Felix Schreibing

# TF and pathway footprints for CD16+ CD8+ TEMRA subtypes----

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(progeny)
library(dorothea)
library(tibble)
library(pheatmap)
library(tidyverse)
library(tidyr)
library(viper)
library(fs)
source('sc_source/sc_source.R')


# Load the subclustered sc datasets ----
indir = '~/sciebo/COVID_project/Robjects/'
outdir = '~/sciebo/COVID_project/01 Results/WP5/5.3 TF and pathway footprints of subclusters/'
setwd(outdir)

NK_TEMRA = readRDS(file = paste0(indir, 'NK-like_TEMRA_subclusters.rds'))

# Load defined color schemes ----
NK.TEMRA.subset.colors = c('CD16+_CD8+_TEMRA_1' = '#8FD744',
                           'CD16+_CD8+_TEMRA_2' = '#9700AD',
                           'CD16+_CD8+_TEMRA_3' = '#00F1FF',
                           'CD16+_CD8+_TEMRA_4' = '#FFFA00',
                           'CD16+_CD8+_TEMRA_5' = '#FF5D5D',
                           'CD16+_CD8+_TEMRA_6' = '#000000')


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


# Run DGEA where we keep all genes and no padj cutoff ----
NK_TEMRA$NK_TEMRA_subtype.condition = paste(NK_TEMRA$NK_TEMRA_subset_annotations, NK_TEMRA$condition_collapsed)
Idents(NK_TEMRA) = 'NK_TEMRA_subtype.condition'


cell.types = names(NK.TEMRA.subset.colors)

for (cell.type in cell.types){
  # mild vs healthy
  case = 'mild'
  control = 'healthy'
  
  markers = get_markers(sc = NK_TEMRA, condition = case, control = control, 
                        cell.type = cell.type, only.sig = FALSE, 
                        logfc.threshold = 0, out.dir = paste0(outdir, 'diff_genes_for_dorothea_NK_TEMRA_subclusters'))
  
  
  # severe vs healthy
  case = 'severe'
  control = 'healthy'
  
  markers = get_markers(sc = NK_TEMRA, condition = case, control = control, 
                        cell.type = cell.type, only.sig = FALSE, 
                        logfc.threshold = 0, out.dir = paste0(outdir, 'diff_genes_for_dorothea_NK_TEMRA_subclusters'))
  
  
  # severe vs mild
  case = 'severe'
  control = 'mild'
  
  markers = get_markers(sc = NK_TEMRA, condition = case, control = control, 
                        cell.type = cell.type, only.sig = FALSE, 
                        logfc.threshold = 0, out.dir = paste0(outdir, 'diff_genes_for_dorothea_NK_TEMRA_subclusters'))
}


# Run DoRothEA on contrasts ----
# Prepare human DoRothEA regulons
dorothea.path = 'https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv'
dorothea_regulon_human = read.csv(dorothea.path)


# Group regulons
regulon = dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c('A','B','C')) %>% 
  split(.$tf) %>%
  map(function(dat) {
    tf = dat %>% distinct(tf) %>% pull()
    targets = setNames(dat$mor, dat$target)
    likelihood = dat$likelihood
    list(tfmode = targets, likelihood = likelihood)
  })

# Run dorothea
# mild vs. healthy
diff.indir = paste0(outdir, 'diff_genes_for_dorothea_NK_TEMRA_subclusters/')
out.dir =  paste0(outdir, 'TF_activity_NK_TEMRA_subclusters/')

case = 'mild'
control = 'healthy'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, out.dir = out.dir, 
                  dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, case, '_vs_', control, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()


# severe vs. healthy
case = 'severe'
control = 'healthy'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, out.dir = out.dir, 
                  dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, case, '_vs_', control, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()


# severe vs. mild
case = 'severe'
control = 'mild'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, out.dir = out.dir, 
                  dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, case, '_vs_', control, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()


# Plot significant TFs ----
# TF data for all three comparisons is only available for CD16+ CD8+ T cell subsets 1; 3 and 4
celltypes = c('CD16+_CD8+_TEMRA_1', 'CD16+_CD8+_TEMRA_3', 'CD16+_CD8+_TEMRA_4')
# Load the TF statistics

for(i in celltypes){
  
  df_mild_healthy = read.delim(paste0(out.dir, 'mild_vs_healthy_tf_activity_', i, '.txt'))
  df_severe_healthy = read.delim(paste0(out.dir, 'severe_vs_healthy_tf_activity_', i, '.txt'))
  df_severe_mild = read.delim(paste0(out.dir, 'severe_vs_mild_tf_activity_', i, '.txt'))
  
  df_mild_healthy$comparison = 'mild_vs_healthy'
  df_severe_healthy$comparison = 'severe_vs_healthy'
  df_severe_mild$comparison = 'severe_vs_mild'
  
  tf.tables = rbind(df_mild_healthy, df_severe_healthy, df_severe_mild)
  rownames(tf.tables) = NULL
  
  # Highlight the significant TFs
  sig.regulons = tf.tables %>% 
    filter(FDR < 0.05) %>% 
    select(Regulon) %>% 
    unique()
  
  sig.regulons = sig.regulons$Regulon
  
  # Filter all significant TFs
  tf.plot = tf.tables %>%
    filter(Regulon %in% sig.regulons) %>%
    mutate(significance = ifelse(FDR < 0.001, '***',
                                 ifelse(FDR < 0.01, '**',
                                        ifelse(FDR < 0.05, '*', 'ns'))))
  
  # Plot NES as heatmap for the three comparisons
  pdf(file = paste0(out.dir, i, '_sig_TF_heatmap.pdf'), width = 7, height = 5)
  print(ggplot(tf.plot, aes(x = Regulon, y = comparison, fill = NES)) +
    geom_tile(color = 'black', lwd = 1) +
    geom_text(aes(label = significance), color = "black", size = 4) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, size = 14)) +
    scale_fill_gradient2(low = '#8AFF6B', mid = 'grey', high = '#FF5D5D', limits = c(-7, 7))
  )
  dev.off()
  
}


# Run PROGENy ----
# Re-define function to compute p-values and make violin plots for Progeny pathway inference results
compute_stats = function(df, celltype, pathways, conditions, plot = FALSE){
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(viridis))
  
  p = NULL
  df_sub = df %>% filter(condition %in% conditions,
                         Pathway %in% pathways,
                         cell.type %in% celltype)
  
  # p.adjust for number of cells * number of pathways tested
  num_test = nrow(df_sub) * length(pathways)
  
  stats = df_sub %>% group_by(Pathway) %>% 
    nest() %>%
    mutate(wilcox = map(data, function(df){
      stest = wilcox.test(Activity ~ condition, data = df, alternative = 'two.sided')
      broom::tidy(stest) %>%
        mutate(corr_pvalue = p.adjust(p.value, method = 'BH', n = num_test))})) %>%
    select(-data) %>%
    unnest(wilcox) %>%
    ungroup() %>%
    arrange(corr_pvalue) %>% 
    as.data.frame
  
  if(plot){
    p = ggplot(df_sub, aes(x = condition, y = Activity, fill = condition)) +
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), colour = 'lightgrey') +
      theme_classic() + 
      scale_fill_manual(values = c('#8AFF6B', '#FF5D5D')) +
      theme(strip.background = element_rect(fill = 'lightgrey'),
            axis.ticks = element_blank()) +
      xlab('') +
      ylab('Pathway activity') +
      facet_wrap( ~ Pathway)
  }
  return(list('stats' = stats, 'p' = p))
}

# Prepare the data
NK_TEMRA$condition_NK_TEMRA_subtype = paste(NK_TEMRA$condition_collapsed, NK_TEMRA$NK_TEMRA_subset_annotations, sep = '.')
Idents(NK_TEMRA) = 'condition_NK_TEMRA_subtype'

NK_TEMRA = progeny(NK_TEMRA, scale = FALSE, organism = 'Human', top = 500, 
             perm = 1, return_assay = TRUE)
NK_TEMRA = ScaleData(NK_TEMRA, assay = 'progeny')

# Create dataframe of clusters
CellsClusters = data.frame(Cell = names(Idents(NK_TEMRA)),
                           CellType = as.character(Idents(NK_TEMRA)),
                           stringsAsFactors = FALSE)

# Transform to data frame
progeny_scores_df = as.data.frame(t(GetAssayData(NK_TEMRA, slot = 'scale.data', assay = 'progeny'))) %>%
  rownames_to_column('Cell') %>%
  gather(Pathway, Activity, -Cell)

# Match Progeny scores with the clusters
progeny_scores_df = inner_join(progeny_scores_df, CellsClusters)

# Summarize Progeny scores 
summarized_progeny_scores = progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

# Create dataframe for plotting
summarized_progeny_scores_df = summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# Plot heatmap
outdir = '~/sciebo/COVID_project/01 Results/WP5/5.3 TF and pathway footprints of subclusters/Pathway_activity_NK_TEMRA_subclusters/'

annotation = data.frame(condition = sapply(strsplit(row.names(summarized_progeny_scores_df),'\\.'), `[`, 1),
                        row.names = row.names(summarized_progeny_scores_df))
celltype = sapply(strsplit(row.names(summarized_progeny_scores_df),'\\.'), `[`, 2)

condition = c('healthy', 'mild', 'severe')
annotation$condition = factor(annotation$condition, levels = condition)
condition.colors = viridis(3)
names(condition.colors) = condition

paletteLength = 100
progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                      length.out = ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength,
                      max(summarized_progeny_scores_df),
                      length.out = floor(paletteLength/2)))

pdf(file = paste0(outdir, 'progeny_heatmap.pdf'), width = 10, height = 5)
pheatmap(t(summarized_progeny_scores_df[,-1]),
         fontsize = 9,
         fontsize_row = 9,
         color = colorRampPalette(c('#00204D', '#FFFFFF', '#FFEA46'))(paletteLength), 
         breaks = progenyBreaks,
         main = '', 
         angle_col = 90,
         treeheight_col = 0,  
         border_color = NA,
         annotation_col = annotation,
         annotation_colors = list(condition = condition.colors),
         labels_col = celltype,
         cluster_cols = FALSE,
         annotation_names_col = FALSE)
dev.off()


# Compute p_values for selected PROGENy pathways in selected NK-like TEMRA cell subsets and plot them ----
progeny_scores_df$condition = sapply(strsplit(progeny_scores_df$CellType,'\\.'), `[`, 1)
progeny_scores_df$cell.type = sapply(strsplit(progeny_scores_df$CellType,'\\.'), `[`, 2)


# PI3K signaling in CD16+_CD8+_TEMRA_1 cells
# severe vs. mild
out.dir = paste0(outdir, 'PI3K_severe_vs_mild/')

conditions = c('mild', 'severe')
celltype = 'CD16+_CD8+_TEMRA_1'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# severe vs. healthy
out.dir = paste0(outdir, 'PI3K_severe_vs_healthy/')

conditions = c('healthy', 'severe')
celltype = 'CD16+_CD8+_TEMRA_1'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# mild vs. healthy
out.dir = paste0(outdir, 'PI3K_mild_vs_healthy/')

conditions = c('healthy', 'mild')
celltype = 'CD16+_CD8+_TEMRA_1'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# PI3K signaling in CD16+_CD8+_TEMRA_2 cells
# severe vs. mild
out.dir = paste0(outdir, 'PI3K_severe_vs_mild/')

conditions = c('mild', 'severe')
celltype = 'CD16+_CD8+_TEMRA_2'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# severe vs. healthy
out.dir = paste0(outdir, 'PI3K_severe_vs_healthy/')

conditions = c('healthy', 'severe')
celltype = 'CD16+_CD8+_TEMRA_2'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# mild vs. healthy
out.dir = paste0(outdir, 'PI3K_mild_vs_healthy/')

conditions = c('healthy', 'mild')
celltype = 'CD16+_CD8+_TEMRA_2'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# PI3K signaling in CD16+_CD8+_TEMRA_3 cells
# severe vs. mild
out.dir = paste0(outdir, 'PI3K_severe_vs_mild/')

conditions = c('mild', 'severe')
celltype = 'CD16+_CD8+_TEMRA_3'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# severe vs. healthy
out.dir = paste0(outdir, 'PI3K_severe_vs_healthy/')

conditions = c('healthy', 'severe')
celltype = 'CD16+_CD8+_TEMRA_3'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# mild vs. healthy
out.dir = paste0(outdir, 'PI3K_mild_vs_healthy/')

conditions = c('healthy', 'mild')
celltype = 'CD16+_CD8+_TEMRA_3'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# PI3K signaling in CD16+_CD8+_TEMRA_4 cells
# severe vs. mild
out.dir = paste0(outdir, 'PI3K_severe_vs_mild/')

conditions = c('mild', 'severe')
celltype = 'CD16+_CD8+_TEMRA_4'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# severe vs. healthy
out.dir = paste0(outdir, 'PI3K_severe_vs_healthy/')

conditions = c('healthy', 'severe')
celltype = 'CD16+_CD8+_TEMRA_4'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# mild vs. healthy
out.dir = paste0(outdir, 'PI3K_mild_vs_healthy/')

conditions = c('healthy', 'mild')
celltype = 'CD16+_CD8+_TEMRA_4'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# PI3K signaling in CD16+_CD8+_TEMRA_5 cells
# severe vs. mild
out.dir = paste0(outdir, 'PI3K_severe_vs_mild/')

conditions = c('mild', 'severe')
celltype = 'CD16+_CD8+_TEMRA_5'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# severe vs. healthy
out.dir = paste0(outdir, 'PI3K_severe_vs_healthy/')

conditions = c('healthy', 'severe')
celltype = 'CD16+_CD8+_TEMRA_5'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# mild vs. healthy
out.dir = paste0(outdir, 'PI3K_mild_vs_healthy/')

conditions = c('healthy', 'mild')
celltype = 'CD16+_CD8+_TEMRA_5'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# PI3K signaling in CD16+_CD8+_TEMRA_6 cells
# severe vs. mild
out.dir = paste0(outdir, 'PI3K_severe_vs_mild/')

conditions = c('mild', 'severe')
celltype = 'CD16+_CD8+_TEMRA_6'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# severe vs. healthy
out.dir = paste0(outdir, 'PI3K_severe_vs_healthy/')

conditions = c('healthy', 'severe')
celltype = 'CD16+_CD8+_TEMRA_6'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# mild vs. healthy
out.dir = paste0(outdir, 'PI3K_mild_vs_healthy/')

conditions = c('healthy', 'mild')
celltype = 'CD16+_CD8+_TEMRA_6'
pathways = c('PI3K')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()