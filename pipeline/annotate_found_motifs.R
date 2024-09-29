#!/usr/bin/env Rscript

# Load required packages while suppressing startup messages
suppressPackageStartupMessages({
library(GenomicRanges)
library(tidyr)
library(stringr)
library(data.table)
library(dplyr)
library(argparser)
library(foreach)
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
p <- add_argument(p, "path", help="data path", type="character")
argv <- parse_args(p)


# Setting data paths
processed_data_dir = paste0(argv$path, "/regulons/")
motif_data_dir = paste0(argv$path, "/motifs_ann/", argv$cell_line)

# Read regulons 
data_tf = fread(paste0(processed_data_dir, "TF_target_mapping_filtered_merged_", argv$cell_line, "_with_ppi_with_dnase_with_atac.tsv"), nThread=10)

# Substring to filter out
substring_to_exclude <- "hoco"

# List files in the directory
files <- list.files(paste0(motif_data_dir))

# ------------------------------------------------------------------------------------------------------

# Motif annotation with HOMER db 

# Filter out files that do not contain the substring
filtered_files <- files[!grepl(substring_to_exclude, files)]

data_tf_temp_ann <- rbindlist(foreach (motif_file = filtered_files) %dopar% {
	print(motif_file)
	motif_data = fread(paste0(motif_data_dir, "/", motif_file), header=T)
	## Filter out peaks with no annotated motifs
	
	motif_data %>% filter(nzchar(.data[[colnames(motif_data)[ncol(motif_data)]]])) -> motif_data
	temp_tf = toupper(unlist(str_split(motif_file, pattern="_"))[1])
	data_tf %>% filter(tf == temp_tf) -> data_tf_temp
	print(temp_tf)
	gr_1 = GRanges(seqnames = data_tf_temp$chr, ranges = IRanges(data_tf_temp$start_peak, data_tf_temp$end_peak), peak_name = data_tf_temp$peak_name, ensembl_transcript = data_tf_temp$ensembl_transcript)
	
	gr_2 = GRanges(seqnames = motif_data$Chr, ranges = IRanges(motif_data$Start, motif_data$End), peak_name = motif_data$PeakID, motif_name = motif_data %>% select(colnames(motif_data)[ncol(motif_data)]))
	
	
	overlap <- findOverlaps(gr_1, gr_2, type="equal", maxgap=1)
	
	data_tf_temp_ann = data_tf_temp[queryHits(overlap)]
	motif_data[subjectHits(overlap)] %>% select(colnames(motif_data)[ncol(motif_data)]) -> data_tf_temp_ann$motifs
	
	data_tf_temp_ann %>% distinct() -> data_tf_temp_ann
	data_tf_temp_ann[, n_motifs := str_count(motifs, "\\(")]
	data_tf_temp_ann
})

if (nrow(data_tf_temp_ann) == 0){
	temp = data_tf
	temp$n_motifs = NA
	temp$motifs = NA
} else {
	left_join(data_tf, data_tf_temp_ann) -> temp
}
#almost there, just need file writing

# ------------------------------------------------------------------------------------------------------

# Motif annotation with HOCOMOCO db
 
# Filter out files that do not contain the substring
filtered_files <- files[grepl(substring_to_exclude, files)]

data_tf_temp_ann <- rbindlist(foreach (motif_file = filtered_files) %dopar% {
        print(motif_file)
        motif_data = fread(paste0(motif_data_dir, "/", motif_file), header=T)
        ## Filter out peaks with no annotated motifs
        motif_data %>% filter(nzchar(.data[[colnames(motif_data)[ncol(motif_data)]]])) -> motif_data
        temp_tf = toupper(unlist(str_split(motif_file, pattern="_"))[1])
        data_tf %>% filter(tf == temp_tf) -> data_tf_temp
        print(temp_tf)   
        gr_1 = GRanges(seqnames = data_tf_temp$chr, ranges = IRanges(data_tf_temp$start_peak, data_tf_temp$end_peak), peak_name = data_tf_temp$peak_name, ensembl_transcript = data_tf_temp$ensembl_transcript)

        gr_2 = GRanges(seqnames = motif_data$Chr, 
			ranges = IRanges(motif_data$Start, motif_data$End), 
			peak_name = motif_data$PeakID, 
			motif_name_hoco = motif_data %>% select(colnames(motif_data)[ncol(motif_data)]))


        overlap <- findOverlaps(gr_1, gr_2, type="equal", maxgap=1)

        data_tf_temp_ann = data_tf_temp[queryHits(overlap)]
        motif_data[subjectHits(overlap)] %>% select(colnames(motif_data)[ncol(motif_data)]) -> data_tf_temp_ann$motifs_hoco

        data_tf_temp_ann %>% distinct() -> data_tf_temp_ann
        data_tf_temp_ann[, n_motifs_hoco := str_count(motifs_hoco, "\\(")]
	data_tf_temp_ann
})

if (nrow(data_tf_temp_ann) == 0){
        temp = temp
	temp$n_motifs_hoco = NA
        temp$motifs_hoco = NA
} else {
        left_join(temp, data_tf_temp_ann) -> temp
}

temp %>% filter(!is.na(n_motifs) | !is.na(n_motifs_hoco)) %>% nrow
print(nrow(temp))

# Write the result to a file
fwrite(temp, paste0(processed_data_dir, "TF_target_mapping_filtered_merged_", argv$cell_line, "_with_ppi_with_dnase_with_atac_with_motifs.tsv"), col.names=T, row.names=F, quote=F, sep="\t")

