#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(here)

here::i_am("README.md")

# Set working directory
setwd(here())

tf_target_mapping = fread("data/regulons/TF_target_mapping_filtered_merged_K562_with_ppi_with_dnase_with_atac_with_motifs_with_ccres.tsv", nThread=10)
coexpr_mat = fread("data/coexpression/Bicor_ENCODE_bulk_K562.tsv")
coexpr_mat$V1 = colnames(coexpr_mat)
coexpr_adj_list = melt(coexpr_mat, id.vars=c("V1"))
colnames(coexpr_adj_list) = c("tf", "gene_symbol", "corr")
coexpr_adj_list %>% filter((tf != gene_symbol) & 
				(tf %in% unique(tf_target_mapping$tf))) -> coexpr_adj_list
left_join(tf_target_mapping, coexpr_adj_list, by=c("tf", "gene_symbol")) %>% 
	mutate(is_coexpressed = (!is.na(corr) & (abs(corr) > 0.6))) -> tf_target_mapping
fwrite(tf_target_mapping, "data/regulons/TF_target_mapping_filtered_merged_K562_with_ppi_with_dnase_with_atac_with_motifs_with_ccres.tsv", sep="\t", col.names=T, row.names=F, quote=F)

