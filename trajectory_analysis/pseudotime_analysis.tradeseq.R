# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de
# edited by Hyojin Kim [2022]

#---- Pseudotime analysis on T cell sub-populations with Slingshot and tradeSeq

library(tradeSeq)
library(Seurat)
library(slingshot)
library(ggplot2)
library(viridis)
library(ggbeeswarm)
library(ggthemes)
library(patchwork)
library(SingleCellExperiment)
library(cowplot)
library(gridExtra)
library(dplyr)
library(scales)
library(PseudotimeDE)
set.seed(42)
'%ni%' = Negate('%in%')


indir = '~/Felix_revision/input/'
outdir = '~/Felix_revision_2021/trajectory_analysis/'
source('~/Felix_revision_2021/CovidEpiMap/sc_source/sc_source.R')
if(!dir.exists(outdir)) dir.create(outdir);
'%ni%' = Negate('%in%')



sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.revision.filtered.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'



# Exclude the below cell types as these have a separate 
excluded.cell.types = 'CD8+ CD73+ regulatory T cells'
cell.types = unique(sc$integrated_annotations)
cell.types = cell.types[cell.types %ni% excluded.cell.types]
sc.subset = subset(sc, integrated_annotations %in% cell.types)
Idents(sc.subset) = factor(Idents(sc.subset), levels = unique(Idents(sc.subset)))
cell.type.colors <- cell.type.colors[cell.types]


# Run slingshot on UMAP embedding
start.clus = 'CD8+ naive T cells'
reduction = 'umap'
sds = slingshot(Embeddings(sc.subset, reduction)[,1:2], clusterLabels = Idents(sc.subset), start.clus = start.clus)
sc.subset@tools[['slingshot']] = SlingshotDataSet(sds)
pseudotime = slingPseudotime(sds)


# Plot slingshot curves
pseudotime = slingPseudotime(sds)
pseudotime %>% as.data.frame() %>%
                write.table(paste0(outdir, "slingshot.pseudotime.txt"), sep="\t", row.names = TRUE, col.names = TRUE)
curves = colnames(pseudotime)
palette = viridis(100, end = 0.95)



# add reduceDim
rd <- Embeddings(sc.subset, reduction)[,1:2] %>% as.matrix()
rd <- rd[colnames(sc.subset),]


pdf(file = paste0(outdir, 'integrated_Tcells_slingshot_curves.pdf'), width = 7.2)
plot(sds$reducedDim, col = cell.type.colors[as.vector(sc.subset$integrated_annotations)], pch = 16, cex = 0.5)
lines(SlingshotDataSet(sds), lwd = 2, col = 'black')
# Show which curve cells belong to
par(mfrow = c(2, 3))
for (curve in curves) {
	colors = palette[cut(pseudotime[,curve], breaks = 100)]
	print(plot(sds$reducedDim, col = colors, pch = 16, cex = 0.5, main = curve) +
	lines(SlingshotDataSet(sds), lwd = 2, col = 'black'))
}
dev.off()




# Add pseudotimes (arclength) to meta data for visualisation
pseudotime2 = slingPseudotime(sds, na = FALSE)
sc.subset$slingshot_pseudotime = pseudotime2[,2]
sc.subset$slingshot_pseudotime_curve1 = pseudotime[,1]
sc.subset$slingshot_pseudotime_curve2 = pseudotime[,2]
saveRDS(sc.subset, file = paste0(indir, 'integrated.RNA.Tcells.pseudotime.subset.rds'))


pdf(file = paste0(outdir, 'integrated_Tcells_slingshot_pseudotime.pdf'), width = 7.8)
FeaturePlot(sc.subset, 
            feature = 'slingshot_pseudotime', 
            cols = c('red', 'yellow')) + 
  ggtitle('') +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(colour = 'Pseudotime')
dev.off()



# Condition density along pseudotime 
df = sc.subset@meta.data
df$integrated_annotations = factor(df$integrated_annotations, 
                                    levels = names(cell.type.colors))

# Active mild and severe conditions
p1 = pseudo_density(data = df %>% 
                   filter(condition_collapsed %in% c('mild', 'severe')), 
                 x = 'slingshot_pseudotime_curve1')
pdf(file = paste0(outdir, 'integrated_Tcells_density_mild_severe.curve1.pdf'), width = 5, height = 5)
p1
dev.off()


p2 = pseudo_density(data = df %>%
                   filter(condition_collapsed %in% c('mild', 'severe')),
                 x = 'slingshot_pseudotime_curve2')
pdf(file = paste0(outdir, 'integrated_Tcells_density_mild_severe.curve2.pdf'), width = 5, height = 5)
p2
dev.off()


# Healthy, recovered mild and recovered severe
p3 = pseudo_density(data = df %>% 
                 filter(condition_collapsed %in% c('healthy', 'mild', 'severe')), 
               x = 'slingshot_pseudotime_curve1')
pdf(file = paste0(outdir, 'integrated_Tcells_density_healthy_mild_severe.curve1.pdf'), width = 6, height = 5)
p3
dev.off()


p4 = pseudo_density(data = df %>%
                 filter(condition_collapsed %in% c('healthy', 'mild', 'severe')),
               x = 'slingshot_pseudotime_curve2')
pdf(file = paste0(outdir, 'integrated_Tcells_density_healthy_mild_severe.curve2.pdf'), width = 6, height = 5)
p4
dev.off()


# Test difference in condition densitiy along pseudotime (active mild vs severe)
df = df %>% 
  filter(condition_collapsed %in% c('mild', 'severe'))
df$condition_collapsed = droplevels(df$condition_collapsed)

# Kolmogorov-Smirnof method to test difference in distribution 
# Use this method over Wilcoxon rank test as it is more sensitive to shape changes in distribution
# and ours appear to be bimodal.

# Lineage 1 (p < 2.2e-16, D = 0.31004)
ks.test(df[df$condition_collapsed == 'mild','slingshot_pseudotime_curve1'],
        df[df$condition_collapsed == 'severe','slingshot_pseudotime_curve1'])


# Lineage 2 (p < 2.2e-16, D = 0.40237)
ks.test(df[df$condition_collapsed == 'mild','slingshot_pseudotime_curve2'],
        df[df$condition_collapsed == 'severe','slingshot_pseudotime_curve2'])



#---- Fit NB-GAM on gene expression data using tradeSeq
# Fit negative binomial generalized additive model (NB-GAM) on active infection subset

subset = subset(sc.subset, condition_collapsed %in% c('mild', 'severe'))
subset = FindVariableFeatures(subset, nfeatures = 10000, verbose = FALSE)
var.genes = VariableFeatures(subset)
subset = subset(subset, features = var.genes)
subset$condition_collapsed = droplevels(subset$condition_collapsed)



# As tradeSeq relies on a negative binomial count distribution, use raw counts
counts = as.matrix(subset[['RNA']]@counts)
pseudotime = slingPseudotime(sds, na = FALSE)[colnames(counts),]
cell.weights = slingCurveWeights(sds)[colnames(counts),]



# Fit smoothed average gene expression profile along pseudotime using a NB-GAM
# with a condition-specific smoother for each lineage (run on cluster, took 2.5 days)
BPPARAM = BiocParallel::bpparam()
BPPARAM$workers = 4

sce = fitGAM(counts = counts, pseudotime = pseudotime,
            cellWeights = cell.weights, nknots = 6,
            conditions = as.factor(unname(subset$condition_collapsed)),
            parallel = TRUE, 
            BPPARAM = BPPARAM)

saveRDS(sce, file = 'sce.parallel.active.subset.rds')


