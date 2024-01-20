# Load required packages
library(data.table)
library(dplyr)
library(doParallel)

# Register parallel processing
registerDoParallel(cores=20)

# Set working directory
setwd("/proj/lappalainen_lab1/users/marii/chip_seq_ann")


# Define data directories
data_dir = "data/regulons/"
processed_data_dir = "data/1-dataset_stats/"

# List of cell lines
cells = c("K562", "HepG2", "MCF7", "GM12878")

# Iterate through each cell line
foreach (ct = cells) %do% {
	# Define the input file
	file = paste0(data_dir, ct, "_all_regulons.tsv")	
	# Read the data
	data = fread(file, nThread=10)

	# Filter the data based on method conditions
	data %>% filter(is_method_1 | is_method_2 | is_method_3) -> data_full
	data_full <- rename(data, all_of(c(is_S2Mb = "is_method_1",
					is_M2Kb = "is_method_2",
					is_S2Kb = "is_method_3")))
        print(colnames(data))
	# Iterate through each method
	foreach(m = c("S2Mb", "M2Kb", "S2Kb")) %do% {
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
