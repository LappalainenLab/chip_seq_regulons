# Suppress package startup messages to prevent clutter
suppressPackageStartupMessages({
library(dplyr)
library(tidyverse)
library(tidyr)
library(biomaRt)
library(data.table)
library(stringr)
library(parallel)
library(foreach)
library(here)
library(argparser)
})

# Set working directory and define data directories
# The working directory is set to a specific path.
setwd("/proj/lappalainen_lab1/users/marii/chip_seq_ann/")
data_dir = "data/encode_bulk_rna/"
processed_data_dir = "data/tss_mapping/"

# Create a connection to the Ensembl database
ensembl <- biomaRt::useEnsembl(
        biomart = "ensembl", 
        dataset = "hsapiens_gene_ensembl",
        host = "http://www.ensembl.org",
	mirror = "useast"
)

# Parse command-line arguments using the argparser package
# This script expects 'cell_line' and 'reps' as input arguments.
# 'cell_line' is the name of a cell line, and 'reps' are replicate filenames separated by commas.
p <- arg_parser("Get input file")
p <- add_argument(p, "cell_line", help="cell_line_name", type="character")
p <- add_argument(p, "reps", help="replicate filenames comma separated", type="character") 
argv <- parse_args(p)

# Extract cell_line and replicate filenames from parsed arguments
cell_line = argv$cell_line
reps = unlist(str_split(argv$reps, pattern=","))

# Define a function to remove version from IDs
remove_version <- function(x) {
        unlist(str_split(x, pattern="[.]"))[1]
}

# Check if the input data files have 5 columns
if (ncol(fread(paste0(data_dir, reps[1], ".tsv"), nThread=10)) == 5){
        # Data processing for 5-column input files
        
        # Define columns to separate target_ids
        target_id_columns <- c("ensembl_transcript", "ensembl_gene", "havana_gene", "havana_transcript", "transcript_symbol", "gene_symbol", "gene_id", "gene_type")
        
        # Read and process data for each replicate
	count_mats = foreach(rep = reps) %do% {
		replicate_1 = fread(paste0(data_dir, rep, ".tsv"), nThread=10)
		# Separate target_ids into multiple columns
		replicate_1 %>% separate(target_id, target_id_columns, sep="[|]")
	}

	# Calculate average expression between replicates
	count_matrix <- count_mats %>% reduce(inner_join, by=target_id_columns)
	count_matrix$tpm_total <- count_matrix %>% dplyr::select(starts_with("tpm")) %>% rowMeans()
	count_matrix = count_matrix[order(count_matrix$ensembl_gene),]
	
	
	# Calculate maximum gene expression among transcripts for each gene
	count_matrix %>% 
		group_by(ensembl_gene) %>%
		mutate(max_tpm = max(tpm_total)) -> count_matrix
	
	# Filter for transcripts with the highest expression
	count_matrix %>%
		filter(tpm_total == max_tpm & tpm_total > 0) -> expressed
	
	# Filter non-expressed genes
	count_matrix %>%
		filter(! ensembl_gene %in% unique(expressed$ensembl_gene)) -> non_expressed
	
	# Remove version from ensembl_gene_id and ensembl_transcript_id
	expressed$ensembl_gene_id <- sapply(expressed$ensembl_gene, remove_version)
	expressed$ensembl_transcript_id <- sapply(expressed$ensembl_transcript, remove_version)
	
	# Retrieve gene ids, strand, TSS, and exons coordinates for expressed genes
	symbols_ensembl_expressed <- as.data.table(biomaRt::getBM(
	        attributes = c("external_gene_name", 
	                       "ensembl_gene_id", 
	                       "chromosome_name", 
	                       "strand", 
	                       "ensembl_transcript_id_version", 
	                       "transcription_start_site", 
	                       "transcript_biotype", 
	                       "ensembl_transcript_id"),
	        filters = c("ensembl_gene_id"),
	        values = expressed$ensembl_gene_id,
	        mart = ensembl
	))
	
	symbols_ensembl_expressed[, strand_num := ifelse(strand == 1, "+", "-")]
	symbols_ensembl_expressed[, chr := paste0("chr", chromosome_name)]
	
	# Join data for expressed genes
	expressed %>% left_join(symbols_ensembl_expressed, by=c("ensembl_transcript_id", "ensembl_gene_id")) -> expressed_joint
	
	# Remove version from ensembl_gene_id and ensembl_transcript_id
	non_expressed$ensembl_gene_id <- sapply(non_expressed$ensembl_gene, remove_version)
	non_expressed$ensembl_transcript_id <- sapply(non_expressed$ensembl_transcript, remove_version)
	
	
	# Retrieve gene ids, strand, TSS, and exons coordinates for non-expressed genes
	symbols_ensembl_non_expressed <- as.data.table(biomaRt::getBM(
	        attributes = c("external_gene_name", 
	                       "ensembl_gene_id",
	                       "chromosome_name",
	                       "strand", 
	                       "ensembl_transcript_id_version", 
	                       "transcription_start_site",
	                       "transcript_biotype",
				"ensembl_transcript_id"),
	        filters = c("ensembl_gene_id"),
        	values = non_expressed$ensembl_gene_id,
        	mart = ensembl
	)) # TSS matches exon_chrom_start if strand == 1, and matches exon_chrom_end if strand == -1
	symbols_ensembl_non_expressed[, strand_num := ifelse(strand == 1, "+", "-")]
	symbols_ensembl_non_expressed[, chr := paste0("chr", chromosome_name)]
	
	# Join data for non-expressed genes
	non_expressed %>% left_join(symbols_ensembl_non_expressed, by=c("ensembl_transcript_id", "ensembl_gene_id")) -> non_expressed_joint
	
	# Select most abundant TSS fpr non-expressed genes
	non_expressed_joint %>% 
		group_by(ensembl_gene_id, transcription_start_site) %>%
		mutate(tss_count=n()) %>% 
		ungroup() %>%
		group_by(ensembl_gene) %>%
		mutate(max_tss_count = max(tss_count)) %>% 
		filter(tss_count == max_tss_count) -> non_expressed_joint
	
	# Join data for expressed and non-expressed genes
	symbols_ensembl = bind_rows(expressed_joint, non_expressed_joint)
	
	# Convert data to BED file format
	bed = unique(symbols_ensembl[c("chr", "transcription_start_site", "transcription_start_site",
        	                         "ensembl_gene_id", "ensembl_transcript", 
        	                         "strand_num", "tpm_total", "transcript_biotype",
        	                         "gene_symbol", "gene_type")])
	bed %>% drop_na() -> bed
	
	# Write data to file
	write.table(bed, paste0(processed_data_dir, cell_line, "_gene_TSS_highest_expressed.bed"), quote=F, sep="\t", row.names = F, col.names = F)
} else {
        # Data processing for 19-column input files
        
        # Read and process data for each replicate
        count_mats = foreach(rep = reps) %do% {
                replicate_1 = fread(paste0(data_dir, rep, ".tsv"), nThread=10)
		replicate_1 %>% dplyr::select(colnames(replicate_1)[1:10])
        }
 
        # Calculate average expression between replicates
        count_matrix <- count_mats %>% reduce(inner_join, by=c("transcript_id", "gene_id"))
        count_matrix$tpm_total <- count_matrix %>% dplyr::select(starts_with("TPM")) %>% rowMeans()
	count_matrix = count_matrix[order(count_matrix$gene_id),]
	
	# Calculate maximum gene expression among transcripts for each gene
	count_matrix %>%
	        group_by(gene_id) %>%
        	mutate(max_tpm = max(tpm_total)) -> count_matrix
	
	# Filter for transcripts with the highest expression
	count_matrix %>%
        	filter(tpm_total == max_tpm & tpm_total > 0) -> expressed
	
	# Filter non-expressed genes
	count_matrix %>%
        	filter(! gene_id %in% unique(expressed$gene_id)) -> non_expressed
	# Remove version from ensembl_gene_id and ensembl_transcript_id
	expressed$ensembl_gene_id <- sapply(expressed$gene_id, remove_version)
	expressed$ensembl_transcript_id <- sapply(expressed$transcript_id, remove_version)
	
	# Retrieve gene ids, strand, TSS, and exons coordinates for expressed genes
	symbols_ensembl_expressed <- as.data.table(biomaRt::getBM(
        	attributes = c("external_gene_name",
                	       "ensembl_gene_id",
                	       "chromosome_name",
                	       "strand",
                	       "gene_biotype",
                	       "ensembl_transcript_id_version",
                	       "transcription_start_site",
                	       "transcript_biotype",
                	        "ensembl_transcript_id"),
        	filters = c("ensembl_gene_id"),
        	values = as.vector(expressed$ensembl_gene_id),
        	mart = ensembl
	)) # TSS matches exon_chrom_start if strand == 1, and matches exon_chrom_end if strand == -1
	symbols_ensembl_expressed[, strand_num := ifelse(strand == 1, "+", "-")]
	symbols_ensembl_expressed[, chr := paste0("chr", chromosome_name)]
	
	# Join data for expressed genes
	expressed %>% left_join(symbols_ensembl_expressed, by=c("ensembl_transcript_id", "ensembl_gene_id")) -> expressed_joint
	
	# Remove version from ensembl_gene_id and ensembl_transcript_id
	non_expressed$ensembl_gene_id <- sapply(non_expressed$gene_id, remove_version)
	non_expressed$ensembl_transcript_id <- sapply(non_expressed$transcript_id, remove_version)
	
	
	# Retrieve gene ids, strand, TSS, and exons coordinates for non-expressed genes
	symbols_ensembl_non_expressed <- as.data.table(biomaRt::getBM(
	        attributes = c("external_gene_name",
        	               "ensembl_gene_id",
        	               "chromosome_name",
        	               "gene_biotype",   
        	               "strand",
        	               "ensembl_transcript_id_version",
        	               "transcription_start_site",
        	               "transcript_biotype",
        	                "ensembl_transcript_id"),
        	filters = c("ensembl_gene_id"),
        	values = non_expressed$ensembl_gene_id,
        	mart = ensembl
	)) # TSS matches exon_chrom_start if strand == 1, and matches exon_chrom_end if strand == -1
	symbols_ensembl_non_expressed[, strand_num := ifelse(strand == 1, "+", "-")]
	symbols_ensembl_non_expressed[, chr := paste0("chr", chromosome_name)]
	
	# Join data for non-expressed genes
	non_expressed %>% left_join(symbols_ensembl_non_expressed, by=c("ensembl_transcript_id", "ensembl_gene_id")) -> non_expressed_joint
	
	# Select most abundant TSS fpr non-expressed genes
	non_expressed_joint %>%
        	group_by(ensembl_gene_id, transcription_start_site) %>%
        	mutate(tss_count=n()) %>%
        	ungroup() %>%
        	group_by(gene_id) %>%
        	mutate(max_tss_count = max(tss_count)) %>%
        	filter(tss_count == max_tss_count) -> non_expressed_joint
	
	# Join data for expressed and non-expressed genes
	symbols_ensembl = bind_rows(expressed_joint, non_expressed_joint)
	
	# Convert data to BED file format
	bed = unique(symbols_ensembl[c("chr", "transcription_start_site", "transcription_start_site",
	                                 "ensembl_gene_id", "ensembl_transcript_id",
	                                 "strand_num", "tpm_total", "transcript_biotype",
	                                 "external_gene_name", "gene_biotype")])
 
	bed %>% drop_na() -> bed
	
	# Write data to file
	write.table(bed, paste0(processed_data_dir, cell_line, "_gene_TSS_highest_expressed.bed"), quote=F, sep="\t", row.names = F, col.names = F)
}
