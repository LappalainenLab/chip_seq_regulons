#!/usr/bin/env Rscript

# Load required packages while suppressing startup messages
suppressPackageStartupMessages({
library(GenomicRanges)
library(tidyr)
library(stringr)
library(data.table)
library(dplyr)
library(argparser)
library(here)
})

here::i_am("README.md")

# Set working directory
setwd(here())

# Parse command-line arguments
# Create a parser
p <- arg_parser("Get input file")

# Add command line arguments
p <- add_argument(p, "cell_line", help="cell_line_name", type="character")
p <- add_argument(p, "atac", help="ATAC-Seq accession", type="character")
p <- add_argument(p, "dnase", help="DNAse-Seq accession", type="character")
p <- add_argument(p, "path", help="data path", type="character")
argv <- parse_args(p)

processed_data_dir = paste0(argv$path, "/regulons/")

# Read TF target mapping data
data_tf = fread(paste0(processed_data_dir, "TF_target_mapping_filtered_merged_", argv$cell_line, ".tsv"), nThread=10)

# Read protein-protein interaction (PPI) data
data_ppi = fread(paste0(argv$path, "/9606.protein.links.full.v11.5.with.names.txt"), nThread=10)

# Join TF target mapping with PPI data
left_join(data_tf, data_ppi, by=c("tf" = "protein1_name", "gene_symbol" = "protein2_name")) %>% 
    select(-c("protein1", "protein2")) %>%
    mutate(is_ppi = !is.na(combined_score)) -> temp

# Filter out rows with 'gene_symbol' equal to "."
temp %>% filter(!gene_symbol == ".") -> tf_target_mapping

# Create binary annotations
tf_target_mapping$tpm_total = as.numeric(tf_target_mapping$tpm_total)
replace_na(tf_target_mapping, list(tpm_total=0)) -> tf_target_mapping

# DNAse annotations 

if (file.exists(paste0(argv$path, "/encode_dnase/", argv$dnase, ".bed"))){
	# Read DNAse data
	dnase_data = fread(paste0(argv$path, "/encode_dnase/", argv$dnase, ".bed"), header=F)
	
	# Create GenomicRanges objects for TF target mapping and DNAse data
	gr_1 = GRanges(seqnames = tf_target_mapping$chr,
	                        ranges = IRanges(tf_target_mapping$start_tss, tf_target_mapping$end_tss),
	                        peak_name = tf_target_mapping$peak_name,
	                        ensembl_transcript = tf_target_mapping$ensembl_transcript)
	gr_2 = GRanges(seqnames = dnase_data$V1,
	                                ranges = IRanges(dnase_data$V2, dnase_data$V3))

	# Count overlaps between TF target mapping and DNAse data
	cre = as.data.table(gr_1)
	cre$n_dnase_peaks = countOverlaps(gr_1, gr_2, type = "any", maxgap = 2000)
	cre %>% select(start, peak_name, ensembl_transcript, n_dnase_peaks) -> cre
	
	# Join TF target mapping with DNAse count data
	tf_target_mapping <- left_join(tf_target_mapping, cre, by = c("start_tss" = "start", "peak_name" = "peak_name", "ensembl_transcript" = "ensembl_transcript"))
	
	# Create a binary variable 'is_dnase' based on the DNAse count
tf_target_mapping$is_dnase = as.logical(tf_target_mapping$n_dnase_peaks)
} else {
	tf_target_mapping$is_dnase = NA
}


# ATAC annotations 

if (file.exists(paste0(argv$path, "/encode_atac/", argv$atac, ".bed"))){
	# Read ATAC data
	atac_data = fread(paste0(argv$path, "/encode_atac/", argv$atac, ".bed"), header=F)
	
	# Create GenomicRanges objects for TF target mapping and ATAC data
	gr_1 = GRanges(seqnames = tf_target_mapping$chr,
	                        ranges = IRanges(tf_target_mapping$start_tss, tf_target_mapping$end_tss),
	                        peak_name = tf_target_mapping$peak_name,
	                        ensembl_transcript = tf_target_mapping$ensembl_transcript)
	gr_2 = GRanges(seqnames = atac_data$V1,
	                                ranges = IRanges(atac_data$V2, atac_data$V3))	

	# Count overlaps between TF target mapping and ATAC data
	cre = as.data.table(gr_1)
	cre$n_atac_peaks = countOverlaps(gr_1, gr_2, type = "any", maxgap = 2000)
	cre %>% select(start, peak_name, ensembl_transcript, n_atac_peaks) -> cre
	
	# Join TF target mapping with ATAC count data
	tf_target_mapping <- left_join(tf_target_mapping, cre, by = c("start_tss" = "start", "peak_name" = "peak_name", "ensembl_transcript" = "ensembl_transcript"))
	
	# Create a binary variable 'is_atac' based on the ATAC count
	tf_target_mapping$is_atac = as.logical(tf_target_mapping$n_atac_peaks)
} else {
	tf_target_mapping$is_atac = NA
}

# Check if 'n_motifs' columns exist and create binary variables based on motif methods
if (as.logical(sum(str_locate(colnames(data_tf), "n_motifs"), na.rm=T))){
        tf_target_mapping$is_motif_method1 = as.logical(tf_target_mapping$n_motifs_method1)
        tf_target_mapping$is_motif_method2 = as.logical(tf_target_mapping$n_motifs_method2)
        tf_target_mapping$is_motif_method3 = as.logical(tf_target_mapping$n_motifs_method3)
        
}

# Write the result to a file
fwrite(tf_target_mapping, paste0(processed_data_dir, "TF_target_mapping_filtered_merged_", argv$cell_line, "_with_ppi_with_dnase_with_atac.tsv"), col.names=T, row.names=F, quote=F, sep="\t")

