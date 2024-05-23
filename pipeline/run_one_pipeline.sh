#!/bin/bash

# Configure
# CELLS="sample_cell" # Name of the cell type consistent with ReMap2022 file name
# RNA_FILES="sample_RNA_1,sample_RNA_2" # Comma-separated RNA-Seq replicate data accessions 
# DNASE_FILE="sample_DNAse" # DNAse-Seq data accession
# ATAC_FILE="sample_ATAC" # ATAC-Seq data accession

CELLS="K562" # Name of the cell type consistent with ReMap2022 file name
RNA_FILES="ENCFF823JHX,ENCFF940GYO" # Comma-separated RNA-Seq replicate data accessions 
DNASE_FILE="ENCFF274YGF" # DNAse-Seq data accession
ATAC_FILE="ENCFF558BLC" # ATAC-Seq data accession

DATA_PATH="./data"
CHIP_FILE_NAME="remap2022_${CELLS}_macs2_hg38_v1_0"

echo "Starting the pipeline"

if [ ! -d "${DATA_PATH}/tss_mapping" ] 
then
    	mkdir "${DATA_PATH}/tss_mapping"
fi

if [ ! -d "${DATA_PATH}/mappings" ] 
then 
     	mkdir "${DATA_PATH}/mappings"
fi 

# Extract highest expressed TSS
Rscript ./pipeline/extract_tss.R $CELLS $RNA_FILES $DATA_PATH
echo "Extracted TSS 1"


# Extract top 50% expressed TSS 
Rscript ./pipeline/extract_multiple_tss.R $CELLS $RNA_FILES $DATA_PATH
echo "Extracted TSS 2"


# Intersect ChIP-Seq with Blacklist regions
bedtools intersect -wa \
			-a "${DATA_PATH}/chipseq_data/${CELLS}/${CHIP_FILE_NAME}.bed" \
			-b "${DATA_PATH}/blacklist/ENCFF023CZC.bed" \
			-v > "${DATA_PATH}/chipseq_data/${CELLS}/${CHIP_FILE_NAME}_no_black.bed"
echo "Removed blacklist regions"

# Sort BED files
sort -k1,1 -k2,2n "${DATA_PATH}/chipseq_data/${CELLS}/${CHIP_FILE_NAME}_no_black.bed" > \
		"${DATA_PATH}/chipseq_data/${CELLS}/${CHIP_FILE_NAME}_no_black_sorted.bed"

sort -k1,1 -k2,2n "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_highest_expressed_enc_v110.bed" > \
                "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_highest_expressed_enc_v110_sorted.bed"

sort -k1,1 -k2,2n "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_50_expressed_50_occur_enc_v110.bed" > \
                "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_50_expressed_50_occur_enc_v110_sorted.bed"
echo "Sorted files"


# Method 1 and 3
bedtools closest -d -a "${DATA_PATH}/chipseq_data/${CELLS}/${CHIP_FILE_NAME}_no_black_sorted.bed" \
                    -b "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_highest_expressed_enc_v110_sorted.bed" > "${DATA_PATH}/mappings/${CHIP_FILE_NAME}_no_black_sorted_mapped_highest_expressed_enc_v110.bed"
echo "Mapped S2"


# Method 2
bedtools closest -d -a "${DATA_PATH}/chipseq_data/${CELLS}/${CHIP_FILE_NAME}_no_black_sorted.bed" \
-b "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_50_expressed_50_occur_enc_v110_sorted.bed" > "${DATA_PATH}/mappings/${CHIP_FILE_NAME}_no_black_sorted_mapped_top50_expressed_enc_v110.bed"
echo "Mapped M2"


# Distance filtering of interactions
Rscript ./pipeline/filter_tf_target_interactions.R  $CELLS $DATA_PATH
echo "Filtered"


# Annotation with PPI, ATAC and DNAse
Rscript ./pipeline/annotate.R $CELLS $ATAC_FILE $DNASE_FILE $DATA_PATH
echo "Annotated"
