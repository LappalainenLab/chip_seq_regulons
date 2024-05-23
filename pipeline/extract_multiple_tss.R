# Load required packages
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
library(httr)
})


# Set working directory
setwd("/proj/lappalainen_lab1/users/marii/chip_seq_ann/")

# Create a connection to the Ensembl database
ensembl <- biomaRt::useEnsembl(
        biomart = "ensembl",
        dataset = "hsapiens_gene_ensembl",
        #host = "https://feb2023.archive.ensembl.org", # release 109
	#host = "https://apr2020.archive.ensembl.org", # release 100
        host = "https://jul2023.archive.ensembl.org", # release 110
        mirror = "uswest"
)

# Parse command-line arguments using the argparser package
# This script expects 'cell_line' and 'reps' as input arguments.
# 'cell_line' is the name of a cell line, and 'reps' are replicate filenames separated by commas.
p <- arg_parser("Get input file")
p <- add_argument(p, "cell_line", help="cell_line_name", type="character")
p <- add_argument(p, "reps", help="replicate filenames comma separated", type="character")
p <- add_argument(p, "path", help="path to working directory", type="character")
argv <- parse_args(p)

# Define data directories
data_dir = paste0(argv$path, "/encode_bulk_rna/")
processed_data_dir = paste0(argv$path, "/tss_mapping/")

# Extract cell_line and replicate filenames
cell_line = argv$cell_line
reps = unlist(str_split(argv$reps, pattern=","))

# Define a function to remove version from IDs
remove_version <- function(x) {
        unlist(str_split(x, pattern="[.]"))[1]
}

if (ncol(fread(paste0(data_dir, reps[1], ".tsv"), nThread=10)) == 5){
        # Data processing for 5-column input files
        
        # Read and parse data for each replicate
        count_mats = foreach(rep = reps) %do% {
                replicate_1 = fread(paste0(data_dir, rep, ".tsv"), nThread=10)

                # Separate target_ids into multiple columns
                replicate_1 %>% separate(target_id, c("ensembl_transcript", "ensembl_gene", "havana_gene","havana_transcript","transcript_symbol", "gene_symbol", "gene_id", "gene_type"), sep="[|]")
        }

        # Calculate average expression between replicates
        count_matrix <- count_mats %>% reduce(inner_join, by=c("ensembl_transcript", "ensembl_gene", "havana_gene", "havana_transcript","transcript_symbol", "gene_symbol", "gene_id", "gene_type"))
        count_matrix$tpm_total <- count_matrix %>% dplyr::select(starts_with("tpm")) %>% rowMeans()
        count_matrix$est_counts <- count_matrix %>% dplyr::select(starts_with("est_counts")) %>% rowMeans()
        count_matrix = count_matrix[order(count_matrix$ensembl_gene),]

        # Calculate median expression of a gene across isoforms
	count_matrix %>% 
		group_by(ensembl_gene) %>%
		mutate(median_tpm = median(tpm_total)) -> count_matrix
	
	# Select top 50% expressed isoforms of each gene
	count_matrix %>%
		filter(tpm_total > median_tpm & est_counts > 5.0) -> expressed
	
	
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
	)) # TSS matches exon_chrom_start if strand == 1, and matches exon_chrom_end if strand == -1
	symbols_ensembl_expressed[, strand_num := ifelse(strand == 1, "+", "-")]
	symbols_ensembl_expressed[, chr := paste0("chr", chromosome_name)]
	
	# Join data for expressed genes
	expressed %>% left_join(symbols_ensembl_expressed, by=c("ensembl_transcript_id", "ensembl_gene_id")) -> expressed_joint
	
		
	# Join data for expressed and non-expressed genes
	symbols_ensembl = expressed_joint
	
	# Convert data to BED file format
	bed = unique(symbols_ensembl[c("chr", "transcription_start_site", "transcription_start_site",
        	                         "ensembl_gene_id", "ensembl_transcript", 
                	                 "strand_num", "tpm_total", "transcript_biotype",
                        	         "gene_symbol", "gene_type")])
	bed %>% drop_na() -> bed
        
	# Write data to file
	write.table(bed, paste0(processed_data_dir, cell_line, "_gene_TSS_50_expressed_50_occur_enc_v110.bed"), quote=F, sep="\t", row.names = F, col.names = F)
} else {
        # Data processing for 19-column input files
        
        # Read and parse data for each replicate
        count_mats = foreach(rep = reps) %do% {
                replicate_1 = fread(paste0(data_dir, rep, ".tsv"), nThread=10)
		replicate_1 %>% dplyr::select(colnames(replicate_1)[1:10])
        }
   
        # Calculate average expression between replicates
        count_matrix <- count_mats %>% reduce(inner_join, by=c("transcript_id", "gene_id"))
        count_matrix$tpm_total <- count_matrix %>% dplyr::select(starts_with("TPM")) %>% rowMeans()
        count_matrix$expected_count <- count_matrix %>% dplyr::select(starts_with("expected_count")) %>% rowMeans()
        count_matrix = count_matrix[order(count_matrix$gene_id),]
	
        # Calculate median expression of a gene across isoforms
	count_matrix %>%
	        group_by(gene_id) %>%
        	mutate(median_tpm = median(tpm_total)) -> count_matrix
	
	# Select top 50% expressed isoforms of each gene
	count_matrix %>%
        	filter(tpm_total > median_tpm & expected_count > 5.0) -> expressed
	
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
	
	# Join data for non-expressed genes
	expressed %>% left_join(symbols_ensembl_expressed, by=c("ensembl_transcript_id", "ensembl_gene_id")) -> expressed_joint
	
	symbols_ensembl = expressed_joint

	# Convert data to BED file format
	bed = unique(symbols_ensembl[, c("chr", "transcription_start_site", "transcription_start_site",
        	                         "gene_id", "transcript_id",
                	                 "strand_num", "tpm_total", "transcript_biotype",
                        	         "external_gene_name", "gene_biotype")])
	
	bed %>% drop_na() -> bed
	
	# Write data to file
	write.table(bed, paste0(processed_data_dir, cell_line, "_gene_TSS_50_expressed_50_occur_enc_v110.bed"), quote=F, sep="\t", row.names = F, col.names = F)
}
