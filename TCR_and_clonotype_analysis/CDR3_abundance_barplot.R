#!/hpcwork/noco0007/tool/run_anaconda3/anaconda3/envs/salmon/bin/Rscript
# source /hpcwork/noco0007/tool/run_anaconda3/anaconda3/etc/profile.d/conda.sh
# conda activate salmon

# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de
# edited by Hyojin Kim and Felix Schreibing

# Analysis of CDR3-sequence abundance per condition ----

library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(tidytext)
library(viridis)
library(cowplot)
library(textshape)
library(tibble)


# Load the data ----
indir = '~/sciebo/COVID_project/Robjects/'
outdir = '~/sciebo/COVID_project/01 Results/WP3/3.4 CDR3 abundance/'

sc = readRDS(paste0(indir, 'integrated.RNA.Tcells.annotated.revision.filtered.rds'))

DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'


# Analysis of CDR3 sequence abundance ----
# Arrange the TCR data ----
# For all cells extract meta data as data frame
df = sc@meta.data %>% rownames_to_column()

# Count CDR3 sequences of all cells (per condition collapsed)
data = df %>% 
  group_by(condition_collapsed) %>% 
  dplyr::count(TCR_cdr3s_aa, name = 'TCR_cdr3s_aa_count') %>% 
  arrange(desc(TCR_cdr3s_aa_count)) %>% 
  as.data.frame()

df$clonotype_cut = cut(df$clonotype_size, breaks = c(0, 1, 5, 20, 100, Inf),
                       labels = rev(c('Hyperexpanded (100 < x <= Inf)',
                                      'Large (20 < x <= 100)',
                                      'Medium (5 < x <= 20)',
                                      'Small (1 < x <= 5)',
                                      'Single (0 < x <= 1)')))

df[is.na(df$Clonotype),'clonotype_cut'] = NA

# Separate CDR3 sequences for TRA and TRB and separate multiple TCR_V_GENES or TCR_J_GENES 
df.orig = df
df = df %>% 
  separate_rows(TCR_cdr3s_aa, sep = ';') %>%
  separate_rows(TCR_V_GENE, sep = ';') %>% ## added by Hyojin Kim
  separate_rows(TCR_J_GENE, sep = ';') %>% ## added by Hyojin Kim
  as.data.frame

# The data frame now has the following format:
#             TCR_cdr3s_aa TCR_V_GENE TCR_D_GENE TCR_J_GENE TCR_C_GENE
# 1       TRA:CAVTSGANSKLTF     TRAV21       None     TRAJ56 TRAC;TRBC1
# 2       TRA:CAVTSGANSKLTF     TRAV21       None    TRBJ1-1 TRAC;TRBC1
# 3       TRA:CAVTSGANSKLTF   TRBV10-3       None     TRAJ56 TRAC;TRBC1
# 4       TRA:CAVTSGANSKLTF   TRBV10-3       None    TRBJ1-1 TRAC;TRBC1
# 5       TRB:CAIRGHWNTEAFF     TRAV21       None     TRAJ56 TRAC;TRBC1
# 6       TRB:CAIRGHWNTEAFF     TRAV21       None    TRBJ1-1 TRAC;TRBC1
# 7       TRB:CAIRGHWNTEAFF   TRBV10-3       None     TRAJ56 TRAC;TRBC1
# 8       TRB:CAIRGHWNTEAFF   TRBV10-3       None    TRBJ1-1 TRAC;TRBC1
# 9  TRA:CAMARGRNSGGSNYKLTF TRAV14/DV4       None     TRAJ53 TRAC;TRBC2
# 10 TRA:CAMARGRNSGGSNYKLTF TRAV14/DV4       None    TRBJ2-5 TRAC;TRBC2


# remove all TCR_cdr3s_aa-TCR_V_GENE-TCR_J_GENE combinations that do not have the same TCR chain (TRA or TRB)
df <- df %>% 
  mutate(TCR_cdr3s_chain = strsplit(TCR_cdr3s_aa, '[:]') %>% lapply(., function(x) x[1]) ) %>%
  mutate(TCR_V_GENE_chain = strsplit(TCR_V_GENE, '[V]') %>% lapply(., function(x) x[1]) ) %>%
  mutate(TCR_J_GENE_chain = strsplit(TCR_J_GENE, '[J]') %>% lapply(., function(x) x[1]) ) %>% 
  mutate(TEST = ifelse(as.character(TCR_cdr3s_chain) == as.character(TCR_V_GENE_chain),1,0)) %>% 
  filter(TEST==1) %>% 
  mutate(TEST.II = ifelse(as.character(TCR_cdr3s_chain) == as.character(TCR_J_GENE_chain),1,0)) %>%
  filter(TEST.II==1)

# The data frame now has the following format:
#                         rowname orig.ident TCR_cdr3s_chain TCR_V_GENE_chain TCR_J_GENE_chain TEST TEST.II
# 1 3_GEX_SURF_AAACGGGGTTGAGTTC-1 3_GEX_SURF             TRA              TRA              TRA    1       1
# 2 3_GEX_SURF_AAACGGGGTTGAGTTC-1 3_GEX_SURF             TRB              TRB              TRB    1       1
# 3 3_GEX_SURF_AACCGCGGTCCTCCAT-1 3_GEX_SURF             TRA              TRA              TRA    1       1


# Add new column with TCR chain
df$TCR_chain = 'TRA'
df[grep('TRB:', df$TCR_cdr3s_aa),'TCR_chain'] = 'TRB'

# Remove "TRA:" and "TRB:" from the CDR3 sequence name
df = df %>% 
  separate_rows(TCR_cdr3s_aa, sep = ':') %>%
  filter(TCR_cdr3s_aa != 'TRA' & TCR_cdr3s_aa != 'TRB')

# Arrange CDR3 sequence in the following format 'TCR_V_GENE_CDR3_TCR_J_GENE'
df$TCR_cdr3s_aa = gsub('TR(A|B):', '', df$TCR_cdr3s_aa)
df$TCR_cdr3s_aa_origin = df$TCR_cdr3s_aa
df$TCR_cdr3s_aa = paste0(df$TCR_V_GENE, "_", df$TCR_cdr3s_aa, "_", df$TCR_J_GENE)
df_clone_cut <- df %>% select("condition_collapsed", "TCR_chain", "TCR_cdr3s_aa", "clonotype_cut") %>%
  unique()

# Summarize data
tcr.data = df %>% 
  group_by(condition_collapsed, TCR_chain) %>%
  dplyr::count(TCR_cdr3s_aa) %>%
  arrange(desc(n)) %>%
  group_by(condition_collapsed, TCR_chain) %>%
  mutate(total_n = sum(n)) %>% 
  mutate(percent = (n / total_n)*100)

  
# Plot the CDR3 sequence abundance per condition for top 15 clones ----
chain_set = c('TRB', 'TRA')
for ( chain in c(chain_set) ) { # added by Hyojin 
  pdf(file = paste0(outdir, 'CDR3_abundance_all_cells_', chain, '_chain.pdf'), height = 5)
  print(tcr.data %>% 
          filter(TCR_chain == chain) %>% # edited by Hyoijn 
          group_by(condition_collapsed) %>% 
          top_n(n = 15, wt = percent) %>%
          ggplot(aes(x = reorder_within(TCR_cdr3s_aa, -percent, condition_collapsed), 
                     y = percent, 
                     fill = condition_collapsed)) +
          geom_bar(stat = 'identity', show.legend = FALSE) +
          scale_x_reordered() +
          cowplot::theme_cowplot() +
          theme(axis.ticks = element_blank(),
                axis.text.x = element_text(size = 8, 
                                           angle = 90, 
                                           hjust = 1, 
                                           vjust = 0.5),
                axis.text.y = element_text(size = 8),
                axis.title.x = element_blank()) +
          facet_grid(~ condition_collapsed, 
                     scales = 'free_x', 
                     space = 'free') +
          scale_fill_viridis(discrete = TRUE, option = 'viridis') +
          labs(y = paste0(chain, ' (%)')) # edited by Hyojin 
  )
  dev.off()
}