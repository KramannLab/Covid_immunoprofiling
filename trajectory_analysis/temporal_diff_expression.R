# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de
# editted by Hyojin Kim [2022] 

#---- Temporal differential gene expression analysis with tradeSeq
library(tradeSeq)
library(ggplot2)
library(viridis)
library(Seurat)
library(pheatmap)
library(SingleCellExperiment)
library(gridExtra)
library(slingshot)
library(fgsea)
library(msigdbr)
set.seed(42)


## final data
indir = '~/Felix_revision/input/'
outdir = '~/Felix_revision_2021/trajectory_analysis/'
source('~/Felix_revision_2021/CovidEpiMap/sc_source/sc_source.R')
if(!dir.exists(outdir)) dir.create(outdir);
'%ni%' = Negate('%in%')


# Get slingshot data
sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.pseudotime.subset.rds'))
sds = sc@tools$slingshot

# NB-GAM model based on active mild/severe subset
sce = readRDS(file = paste0(indir, 'sce.parallel.active.subset.rds'))
sce$crv <- sce$slingshot


#---- Visualize knots on UMAP

pdf(file = paste0(outdir, 'tradeseq_knots.pdf'))
plotGeneCount(sds, models = sce, clusters = apply(slingClusterLabels(sds), 1, which.max))
dev.off()


#---- Test for DE between start- and end-point of lineages (progenitor markers)
# Temporal DE
start.res = startVsEndTest(sce, l2fc = log2(2), lineages = TRUE)
start.res$padj_lineage1 = p.adjust(start.res$pvalue_lineage1, 'fdr')
start.res$padj_lineage2 = p.adjust(start.res$pvalue_lineage2, 'fdr')

order = order(start.res$waldStat, decreasing = TRUE)
start.order = start.res[order,]
start.order$gene = rownames(start.order)

# Write to file
write.table(start.order[,c(14,1:13)], 
			file = paste0(outdir, 'start.vs.end.test.active.subset.txt'),
			sep = '\t', 
			row.names = FALSE,
			quote = FALSE)

# Plot top ranking genes
sig.subset.1 = start.res[which(start.res$padj_lineage1 < 0.05),]
sig.subset.2 = start.res[which(start.res$padj_lineage2 < 0.05),]
order.1 = order(sig.subset.1$waldStat_lineage1, decreasing = TRUE)
order.2 = order(sig.subset.2$waldStat_lineage2, decreasing = TRUE)
start.genes.1 = rownames(sig.subset.1[order.1,])
start.genes.2 = rownames(sig.subset.2[order.2,])

file.prefix.1 = 'start_vs_end_test_lineage1'
file.prefix.2 = 'start_vs_end_test_lineage2'

plot_top_genes(gene.set = start.genes.1, 
               model = sce,
               n = 150, 
               out.dir = outdir, 
               file.name = paste0(file.prefix.1, '_top_genes'))

plot_top_genes(gene.set = start.genes.2, 
               model = sce, 
               n = 150,
               out.dir = outdir, 
               file.name = paste0(file.prefix.2, '_top_genes'))


# Heatmap of top genes 
genes = start.genes.1[1:100]

n = 50
smooth = predictSmooth(sce, gene = genes, nPoints = n, tidy = FALSE)
smooth = smooth[,c((3*n+1):(n*4),(n+1):(2*n),(2*n+1):(3*n),1:n)]

# Annotations for columns
col.data = data.frame(data = colnames(smooth))
col.data$lineage = c(rep('lineage 2', 2*n), rep('lineage 1', 2*n))
col.data$condition = c(rep(c(rep('severe', n), rep('mild', n)),2))
rownames(col.data) = col.data$data
col.data$data = NULL

# Annotation colours
annotation.colors = list(lineage = c('lineage 1' = viridis(20)[10],
                                     'lineage 2' = viridis(20)[17]),
                         condition = c('mild' = viridis(2)[1],
                                       'severe' = viridis(2)[2]))

# Define colour palette
paletteLength = 100
myColor = colorRampPalette(c(viridis(paletteLength)[1], 
                             'white', 
                             viridis(paletteLength)[paletteLength]))(paletteLength)

smooth = scale(t(smooth))
myBreaks = c(seq(min(smooth), 0, 
                 length.out = ceiling(paletteLength/2) + 1),
             seq(max(smooth)/paletteLength, 
                 max(smooth), 
                 length.out = floor(paletteLength/2)))

pdf(file = paste0(outdir, 'start_vs_end_test_DEG_heatmap.top100.pdf'), width = 5, height = 8.0)
pheatmap(t(smooth),
         cluster_cols = FALSE,
         show_rownames = TRUE, 
         show_colnames = FALSE,
         breaks = myBreaks,
         fontsize_row = 6,
         treeheight_row = 0,
         annotation_col = col.data[c('condition', 'lineage')],
	 color = myColor,
         annotation_colors = annotation.colors,
         border = NA)
dev.off()



# GSEA
bg.genes = names(sce)
order.1 = order(start.res$waldStat_lineage1, decreasing = TRUE)
order.2 = order(start.res$waldStat_lineage2, decreasing = TRUE)

stats.1 = start.res[order.1,]$waldStat_lineage1
stats.2 = start.res[order.2,]$waldStat_lineage2
names(stats.1) = rownames(start.res[order.1,])
names(stats.2) = rownames(start.res[order.2,])
stats.1 = stats.1[!is.na(stats.1)]
stats.2 = stats.2[!is.na(stats.2)]

## GO
try(run_gsea(bg.genes = bg.genes, stats = stats.1, 
		category = 'C5', subcategory = 'BP',
		out.dir = outdir, plot.title = 'GO',
		file.prefix = file.prefix.1, n = 150)
		)

try(run_gsea(bg.genes = bg.genes, stats = stats.2, 
		category = 'C5', subcategory = 'BP',
		out.dir = outdir, plot.title = 'GO',
		file.prefix = file.prefix.2, n = 150)
		)

## PID
try(run_gsea(bg.genes = bg.genes, stats = stats.1, 
		category = 'C2', subcategory = 'PID',
		out.dir = outdir, plot.title = 'PID',
		file.prefix = file.prefix.1, n = 150)
		)

try(run_gsea(bg.genes = bg.genes, stats = stats.2, 
		category = 'C2', subcategory = 'PID',
		out.dir = outdir, plot.title = 'PID',
		file.prefix = file.prefix.2, n = 150)
		)

## Immunological signature
try(run_gsea(bg.genes = bg.genes, stats = stats.1, 
		category = 'C7',
		out.dir = outdir, plot.title = 'Immunological Signature',
		file.prefix = file.prefix.1, n = 150)
		)

try(run_gsea(bg.genes = bg.genes, stats = stats.2, 
		category = 'C7',
		out.dir = outdir, plot.title = 'Immunological Signature',
		file.prefix = file.prefix.2, n = 150)
		)


#---- Test for DE between end-points of lineages (differentiated markers)

# Temporal DE
end.res = diffEndTest(sce, l2fc = log2(2))
end.res$padj = p.adjust(end.res$pvalue, 'fdr')

order = order(end.res$waldStat, decreasing = TRUE)
end.order = end.res[order,]
end.order$gene = rownames(end.order)
end.genes = rownames(end.order)[which(end.order$padj < 0.05)]

# Write to file
write.table(end.order[,c(6,1:5)], 
			file = paste0(outdir, 'diff.end.test.active.subset.txt'),
			sep = '\t', 
			row.names = FALSE,
			quote = FALSE)

# Plot top ranking genes
file.prefix = 'diff_end_test'

plot_top_genes(gene.set = end.genes, 
               model = sce, 
               n = length(end.genes),
               out.dir = outdir, 
               file.name = paste0(file.prefix, '_top_genes'))


## GSEA

stats = end.order$waldStat
names(stats) = rownames(end.order)
stats = stats[!is.na(stats)]


# GO
try(run_gsea(bg.genes = bg.genes, stats = stats, 
		category = 'C5', subcategory = 'BP',
		out.dir = outdir, plot.title = 'GO',
		file.prefix = file.prefix, n = 150)
		)

# PID
try(run_gsea(bg.genes = bg.genes, stats = stats, 
		category = 'C2', subcategory = 'PID',
		out.dir = outdir, plot.title = 'PID',
		file.prefix = file.prefix, n = 150)
		)

## Immunological signature
try(run_gsea(bg.genes = bg.genes, stats = stats, 
		category = 'C7',
		out.dir = outdir, plot.title = 'Immunological Signature',
		file.prefix = file.prefix, n = 150)
		)


#---- Test DE where the lineages split

# Temporal DE
# Use lower l2fc threshold to capture more subtle changes
split.res = earlyDETest(sce, knots = c(2, 3), l2fc = log2(1.5))
split.res$padj = p.adjust(split.res$pvalue, 'fdr')

order = order(split.res$waldStat, decreasing = TRUE)
split.order = split.res[order,]
split.order$gene = rownames(split.order)
split.genes = rownames(split.order)[which(split.order$padj < 0.05)]

# Write to file
write.table(split.order[,c(6,1:5)], 
			file = paste0(outdir, 'lineage.split.test.active.subset.txt'),
			sep = '\t', 
			row.names = FALSE,
			quote = FALSE)

# Plot top ranking genes
file.prefix = 'lineage_split_test'
plot_top_genes(gene.set = split.genes, 
               model = sce,
               n = length(split.genes),
               out.dir = outdir, 
               file.name = paste0(file.prefix, '_top_genes'))


## GSEA

stats = split.order$waldStat
names(stats) = rownames(split.order)
stats = stats[!is.na(stats)]

## GO
try(run_gsea(bg.genes = bg.genes, stats = stats, 
         category = 'C5', subcategory = 'BP',
         out.dir = outdir, plot.title = 'GO',
         file.prefix = file.prefix, n = 150)
	 )


## PID
try(run_gsea(bg.genes = bg.genes, stats = stats, 
         category = 'C2', subcategory = 'PID',
         out.dir = outdir, plot.title = 'PID',
         file.prefix = file.prefix, n = 150)
	)



# Immunological signature
try(run_gsea(bg.genes = bg.genes, stats = stats, 
         category = 'C7',
         out.dir = outdir, plot.title = 'Immunological Signature',
         file.prefix = file.prefix, n = 150)
	)


#---- Test for DE genes between conditions

# Temporal DE
condition.res = conditionTest(sce, l2fc = log2(2),  global = TRUE, lineages = TRUE) #pairwise = TRUE)
condition.res$padj_lineage1 = p.adjust(condition.res$pvalue_lineage1, 'fdr')
condition.res$padj_lineage2 = p.adjust(condition.res$pvalue_lineage2, 'fdr')

order1 = order(condition.res$waldStat_lineage1, decreasing = TRUE)
condition.order_lineage1 = condition.res[order1,]
condition.order_lineage1$gene = rownames(condition.order_lineage1)

order2 = order(condition.res$waldStat_lineage2, decreasing = TRUE)
condition.order_lineage2 = condition.res[order2,]
condition.order_lineage2$gene = rownames(condition.order_lineage2)


# Write to file
write.table(condition.order_lineage1[,c(12,4,5,6,10)] , 
            file = paste0(outdir, 'condition.lineage1.active.subset.txt'),
            sep = '\t', 
            row.names = FALSE,
            quote = FALSE)

write.table(condition.order_lineage2[,c(12,7,8,9,11)] ,
            file = paste0(outdir, 'condition.lineage2.active.subset.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

# Plot top ranking genes
sig.subset.1 = condition.res[which(condition.res$padj_lineage1 < 0.05),]
sig.subset.2 = condition.res[which(condition.res$padj_lineage2 < 0.05),]

order.1 = order(sig.subset.1$waldStat_lineage1, decreasing = TRUE)
order.2 = order(sig.subset.2$waldStat_lineage2, decreasing = TRUE)
genes.1 = rownames(sig.subset.1[order.1,])
genes.2 = rownames(sig.subset.2[order.2,])

file.prefix.1 = 'condition_test_lineage1'
file.prefix.2 = 'condition_test_lineage2'

plot_top_genes(gene.set = genes.1, model = sce, n = length(genes.1), 
				out.dir = outdir, file.name = paste0(file.prefix.1, '_top_genes'))

plot_top_genes(gene.set = genes.2, model = sce, n = length(genes.2), 
				out.dir = outdir, file.name = paste0(file.prefix.2, '_top_genes'))



# GSEA
order.1 = order(condition.res$waldStat, decreasing = TRUE)
stats.1 = condition.res[order.1,]$waldStat
names(stats.1) = rownames(condition.res[order.1,])
stats.1 = stats.1[!is.na(stats.1)]

# GO
try(run_gsea(bg.genes = bg.genes, stats = stats.1, 
         category = 'C5', subcategory = 'BP',
         out.dir = outdir, plot.title = 'GO',
         file.prefix = file.prefix.1, n = 150)
	)


# PID
try(run_gsea(bg.genes = bg.genes, stats = stats.1, 
         category = 'C2', subcategory = 'PID',
         out.dir = outdir, plot.title = 'PID',
         file.prefix = file.prefix.1, n = 150)
	)


# Immunological signature
try(run_gsea(bg.genes = bg.genes, stats = stats.1, 
         category = 'C7',
         out.dir = outdir, plot.title = 'Immunological Signature',
         file.prefix = file.prefix.1, n = 150)
	)




