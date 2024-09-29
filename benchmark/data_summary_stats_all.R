#!/usr/bin/env Rscript

# Load required packages
library(data.table)
library(dplyr)
library(doParallel)
library(here)

# Register parallel processing
registerDoParallel(cores=20)

# Set working directory
setwd(here())


# Define data directories
data_dir = "data/regulons/"
processed_data_dir = "data/1-dataset_stats/"

# List of cell lines
cells = c("K562", "HepG2", "MCF7", "GM12878", "A549", "HeLa-S3", "HCT116", "SK-N-SH", "MCF10A", "IMR90", "PC3", "H9", "A375", "OCI-Ly7", "PANC-1", "BE2C", "GM23338", "Caco2", "BJ", "GM12891", "HT29", "Jurkat", "A673", "G401", "HT1080", "PC9", "GM12892", "PFSK1", "U-87MG", "Calu3", "HUES64", "H1", "HFF", "Karpas422", "NCI-H460", "MG-63-3", "SK-MEL-5", "GM23248", "DAOY", "SJSA1")


# Iterate through each cell line
ds_stats <- rbindlist(foreach (ct = cells) %do% {
	# Define the input file
	file1 = paste0(processed_data_dir, "per_tf_stats_S2Mb_", ct, ".tsv")
    file2 = paste0(processed_data_dir, "per_tf_stats_M2Kb_", ct, ".tsv")
    file3 = paste0(processed_data_dir, "per_tf_stats_S2Kb_", ct, ".tsv")
    file4 = paste0(processed_data_dir, "per_tf_stats_S100Kb_", ct, ".tsv")
    file5 = paste0(processed_data_dir, "per_tf_stats_M100Kb_", ct, ".tsv")

	# Read the data
	temp1 = fread(file1, nThread=10)
    temp2 = fread(file2, nThread=10)
    temp3 = fread(file3, nThread=10)
    temp4 = fread(file4, nThread=10)
    temp5 = fread(file5, nThread=10)
	data.table(ct, mean(temp1$n_genes),  mean(temp2$n_genes),  mean(temp3$n_genes), mean(temp4$n_genes), mean(temp5$n_genes))
})

colnames(ds_stats) = c("cell_line", "S2Mb", "M2Kb", "S2Kb", "S100Kb", "M100Kb")
fwrite(ds_stats, paste0(processed_data_dir, "dataset_summary_stats_all.tsv"), sep="\t", col.names=T, row.names=F, quote=F)

