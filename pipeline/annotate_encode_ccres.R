#!/usr/bin/Rscript

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
p <- add_argument(p, "path", help="data path", type="character")
argv <- parse_args(p)

processed_data_dir = paste0(argv$path, "/regulons/")

# Read TF target mapping data
tf_target_mapping = fread(paste0(processed_data_dir, "TF_target_mapping_filtered_merged_", argv$cell_line, "_with_ppi_with_dnase_with_atac_with_motifs.tsv"), nThread=10)

temp = tf_target_mapping
# PLS annotations 

if (file.exists(paste0(argv$path, "/encode_ccre/", argv$cell_line, "_PLS.bed"))){
	# Read PLS data
	pls_data = fread(paste0(argv$path, "/encode_ccre/", argv$cell_line, "_PLS.bed"), header=F)
	
	# Create GenomicRanges objects for TF target mapping and PLS data
        gr_1 = GRanges(seqnames = temp$chr, 
			ranges = IRanges(temp$start_peak, temp$end_peak), 
			peak_name = temp$peak_name, 
			ensembl_transcript = temp$ensembl_transcript)


	gr_2 = GRanges(seqnames = pls_data$V1,
	                                ranges = IRanges(pls_data$V2, pls_data$V3))

	overlap <- findOverlaps(gr_1, gr_2, type="any", maxgap=1)
 
        data_tf_temp_ann = temp[queryHits(overlap)]
 
        data_tf_temp_ann %>% distinct() -> data_tf_temp_ann

        data_tf_temp_ann$is_pls = T

        left_join(temp, data_tf_temp_ann) -> temp       
	temp[is.na(temp$is_pls)]$is_pls <- F 
} else {
        temp = temp
        temp$is_pls = NA
} 

# pELS annotations

if (file.exists(paste0(argv$path, "/encode_ccre/", argv$cell_line, "_pELS.bed"))){
        # Read pELS data
        pELS_data = fread(paste0(argv$path, "/encode_ccre/", argv$cell_line, "_pELS.bed"), header=F)

        # Create GenomicRanges objects for TF target mapping and pELS data
        gr_1 = GRanges(seqnames = temp$chr,
                        ranges = IRanges(temp$start_peak, temp$end_peak),
                        peak_name = temp$peak_name,
                        ensembl_transcript = temp$ensembl_transcript)


        gr_2 = GRanges(seqnames = pELS_data$V1,
                                        ranges = IRanges(pELS_data$V2, pELS_data$V3))

        overlap <- findOverlaps(gr_1, gr_2, type="any", maxgap=1)
 
        data_tf_temp_ann = temp[queryHits(overlap)]
 
        data_tf_temp_ann %>% distinct() -> data_tf_temp_ann

        data_tf_temp_ann$is_pels = T

        left_join(temp, data_tf_temp_ann) -> temp
temp[is.na(temp$is_pels)]$is_pels <- F 
} else {
        temp = temp
        temp$is_pels = NA
}


# dELS annotations

if (file.exists(paste0(argv$path, "/encode_ccre/", argv$cell_line, "_dELS.bed"))){
        # Read dELS data
        dELS_data = fread(paste0(argv$path, "/encode_ccre/", argv$cell_line, "_dELS.bed"), header=F)

        # Create GenomicRanges objects for TF target mapping and dELS data
        gr_1 = GRanges(seqnames = temp$chr,
                        ranges = IRanges(temp$start_peak, temp$end_peak),
                        peak_name = temp$peak_name,
                        ensembl_transcript = temp$ensembl_transcript)


        gr_2 = GRanges(seqnames = dELS_data$V1,
                                        ranges = IRanges(dELS_data$V2, dELS_data$V3))

        overlap <- findOverlaps(gr_1, gr_2, type="any", maxgap=1)
 
        data_tf_temp_ann = temp[queryHits(overlap)]
 
        data_tf_temp_ann %>% distinct() -> data_tf_temp_ann

        data_tf_temp_ann$is_dels = T

        left_join(temp, data_tf_temp_ann) -> temp
        temp[is.na(temp$is_dels)]$is_dels <- F 
} else {
        temp = temp
        temp$is_dels = NA
}

# Write the result to a file
write.table(temp, paste0(processed_data_dir, "TF_target_mapping_filtered_merged_", argv$cell_line, "_with_ppi_with_dnase_with_atac_with_motifs_with_ccres.tsv"), col.names=T, row.names=F, quote=F, sep="\t")

