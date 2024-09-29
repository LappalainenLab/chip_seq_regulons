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
foreach (ct = cells) %do% {
	# Define the input file
	file = paste0(data_dir, "TF_target_mapping_filtered_merged_", ct, "_with_ppi_with_dnase_with_atac_with_motifs_with_ccres_cleaned.tsv")	
	# Read the data
	data = fread(file, nThread=10)

	data %>% select(c(n_motifs,
                          n_motifs_hoco,
                          is_dnase,
                          is_S2Mb,
                          is_S2Kb,
                          is_M2Kb,
                          is_S100Kb,
                          is_M100Kb,
                          is_atac,
                          ensembl_gene_id,
                          tpm_total,
                          is_ppi,
                          tf,
                          gene_symbol,
                          accession)) %>%
                  distinct() -> data

	# Filter the data based on method conditions
	data %>% filter(is_S2Mb | is_M2Kb | is_S2Kb | is_M100Kb | is_S100Kb) -> data_full

        print(colnames(data))
	# Iterate through each method
	foreach(m = c("S2Mb", "M2Kb", "S2Kb", "M100Kb", "S100Kb")) %do% {
		# Filter interactions for the current regulon
		data_full %>% filter(!!as.name(paste0("is_", m)) == T) -> data
		
		# Calculate number of TFs per gene
		data %>% 
			group_by(ensembl_gene_id) %>% 
			summarise(n_tfs = n_distinct(tf)) -> temp
		
		# Write the result to a file
		fwrite(temp, paste0(processed_data_dir, "per_gene_stats_", m, "_", ct, ".tsv"), sep="\t", col.names=T, row.names=F, quote=F)
		
		# Calculate number of genes per TF
		data %>% 
			group_by(tf) %>% 
			summarise(n_genes = n_distinct(ensembl_gene_id)) -> temp
		
		# Write the result to a file
		fwrite(temp, paste0(processed_data_dir, "per_tf_stats_", m, "_", ct, ".tsv"), sep="\t", col.names=T, row.names=F, quote=F)
		}
}
