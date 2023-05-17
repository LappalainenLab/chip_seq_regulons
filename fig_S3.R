library(data.table)
library(here)
library(jtools)
library(stringr)
library(biomaRt)
library(ggplot2)
library(fastglm)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(doParallel)
 
registerDoParallel(cores=10)
here::i_am("extract_tss.R")
 
 
# Plot style
source("fig_style.R")

temp = fread("data/enrich_analysis/per_gene_motif_stats.tsv")

plot = ggplot(temp, aes(x = n_tfs)) +
        geom_histogram() +
        geom_vline(xintercept = c(mean(temp$n_tfs)), linetype="dashed", color = c("red")) +
        gtex_v8_figure_theme() + ylab("TFs") + xlab("Number of TFs per gene") +
        annotate("text", x = 140, y=1000, label=paste0("mean: ", round(mean(temp$n_tfs), 0)))

ggsave(filename = "plots/N_tfs_with_motif_per_gene.png", plot = plot, device = "png", dpi = 600, width=4, height=4)


temp = fread("data/enrich_analysis/per_tf_motif_stats.tsv")

plot = ggplot(temp, aes(x = n_genes)) +
        geom_histogram() +
        geom_vline(xintercept = c(mean(temp$n_genes)), linetype="dashed", color = c("red")) +
        gtex_v8_figure_theme() + ylab("Genes") + xlab("Number of genes per TF") +
        annotate("text", x = 14000, y=10, label=paste0("mean: ", round(mean(temp$n_genes), 0)))

ggsave(filename = "plots/N_genes_with_motif_per_tf.png", plot = plot, device = "png", dpi = 600, width=4, height=4)
