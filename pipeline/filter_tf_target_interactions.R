# Load required packages while suppressing startup messages
suppressPackageStartupMessages({
library(tidyr)
library(data.table)
library(stringr)
library(argparser)
library(dplyr)
})


# Set working directory
setwd("/proj/lappalainen_lab1/users/marii/chip_seq_ann")

# Parse command-line arguments
# Create a parser
p <- arg_parser("Get input file")

# Add command-line argument for cell line name
p <- add_argument(p, "cell_line", help="cell_line_name", type="character")
p <- add_argument(p, "path", help="data path", type="character")
argv <- parse_args(p)

# Function to standardize column names in data frames
get_col_annot = function(data){
        col_names = c(
                "chr_peak", "start_peak", "end_peak", "peak_name",
                "score", "peak_strand", "color_start", "color_end",
                "color", "chr", "start_tss", "end_tss",
                "ensembl_gene_id", "ensembl_transcript",
                "strand_num", "tpm_total", "transcript_biotype",
                "gene_symbol", "gene_type", "distance"
        )
        
	if (ncol(data) == 20){
	
	        colnames(data) = col_names
	}
	if (ncol(data) == 16) {
	        colnames(data) = c(col_names[1:5], col_names[10:20])
	}
	if (ncol(data) < 16) {
	        colnames(data) = c(col_names[1:3], col_names[10:20])
		print(head(data))
	}
	if (ncol(data) == 17) {
	        colnames(data) = c(col_names[1:6], col_names[10:20])
	}
	return(data)
}


# S2 methods: highest expressed isoform + most frequent non-expressed

# Read data from a file
tf_target_mapping_highest_meth_1 = fread(paste0(argv$path, "/mappings/remap2022_", argv$cell_line, "_macs2_hg38_v1_0_no_black_sorted_mapped_highest_expressed_enc_v110.bed"), header=F)

# Standardize column names
tf_target_mapping_highest_meth_1 = get_col_annot(tf_target_mapping_highest_meth_1)

# Filter interactions based on the distance filters: 
# +/- 1Mb around TSS for the S2Mb method and +/- 1Kb around TSS for the S2Kb method
tf_target_mapping_highest_meth_1 %>%
	mutate(is_S2Mb = abs(distance) < 1000000,
	       is_S2Kb = abs(distance) < 1000) -> tf_target_mapping_highest


# M2Kb method: top 50% expressed isoforms + top 50% most abundant isoforms for non-expressed genes

# Read data from a file
tf_target_mapping_50 = fread(paste0(argv$path, "/mappings/remap2022_", argv$cell_line, "_macs2_hg38_v1_0_no_black_sorted_mapped_top50_expressed_enc_v110.bed"), header=F)

# Standardize column names
tf_target_mapping_50 = get_col_annot(tf_target_mapping_50)

# Filter interactions based on the distance filters: +/- 1Kb around TSS
tf_target_mapping_50 %>% 
	mutate(is_M2Kb = abs(distance) < 1000) -> tf_target_mapping_50

# Join S2 and M2 regulons and filter out negative interactions
full_join(tf_target_mapping_highest, tf_target_mapping_50,
		by = join_by(chr_peak, 
		             start_peak, 
		             end_peak, 
		             peak_name, 
		             score, 
		             peak_strand, 
		             chr, start_tss, 
		             end_tss, 
		             ensembl_gene_id, 
		             ensembl_transcript, 
		             strand_num, tpm_total, 
		             transcript_biotype, 
		             gene_symbol, 
		             gene_type))  %>% 
        mutate(is_S2Mb = replace_na(is_S2Mb, FALSE),
               is_S2Kb = replace_na(is_S2Kb, FALSE), 
               is_M2Kb = replace_na(is_M2Kb, FALSE)) -> joint_mapping


needed_cols = c("chr", "start_tss", "end_tss",
                "ensembl_gene_id", "ensembl_transcript",
                "strand_num", "tpm_total", "transcript_biotype",
                "gene_symbol", "gene_type",
                "chr_peak", "start_peak", "end_peak", "peak_name",
                "score", "peak_strand", "distance.x", "distance.y")

# Format data frame
if (ncol(joint_mapping) < length(needed_cols)){
	differ = setdiff(needed_cols, colnames(joint_mapping))
	for (col in differ){
		joint_mapping[col] = NA
	}
}

# Calculate number of mapped non-redundant ChIP-Seq peaks per transcript
joint_mapping %>%
	group_by(peak_name, ensembl_transcript, start_tss, transcript_biotype) %>% 
	mutate(n_peaks_M2Kb = sum(is_M2Kb), 
		n_peaks_S2Kb = sum(is_S2Kb), 
		n_peaks_S2Mb = sum(is_S2Mb)) %>%
	ungroup() %>%
	distinct(peak_name, ensembl_transcript, start_tss, transcript_biotype, .keep_all = TRUE) -> tf_target_mapping

# Split peak_name into 'accession' and 'tf'
tf_target_mapping %>%
	rowwise() %>%
	mutate(accession = unlist(str_split(peak_name, pattern="[.]"))[1],
		tf = unlist(str_split(peak_name, pattern="[.]"))[2]) -> tf_target_mapping

# Write the data frame to a TSV file
write.table(
        tf_target_mapping,
        paste0(argv$path, "/regulons/TF_target_mapping_filtered_merged_", argv$cell_line, "_enc_v110.tsv"),
        quote=F, sep="\t", row.names = F, col.names = T
        )

