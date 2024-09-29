#!/bin/bash

# Configure
CELLS=$1 # Name of the cell type consistent with ReMap2022 file name
RNA_FILES=$2 # Comma-separated RNA-Seq replicate data accessions 
DNASE_FILE=$3 # DNAse-Seq data accession
ATAC_FILE=$4 # ATAC-Seq data accession

echo $CELL $RNA_FILES $DNASE_FILE $ATAC_FILE

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

if [ ! -d "${DATA_PATH}/regulons" ] 
then 
     	mkdir "${DATA_PATH}/regulons"
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

sort -k1,1 -k2,2n "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_highest_expressed.bed" > \
                "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_highest_expressed_sorted.bed"

sort -k1,1 -k2,2n "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_50_expressed_50_occur.bed" > \
                "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_50_expressed_50_occur_sorted.bed"
echo "Sorted files"


# Method 1 and 3
bedtools closest -d -a "${DATA_PATH}/chipseq_data/${CELLS}/${CHIP_FILE_NAME}_no_black_sorted.bed" \
                    -b "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_highest_expressed_sorted.bed" > "${DATA_PATH}/mappings/${CHIP_FILE_NAME}_no_black_sorted_mapped_highest_expressed.bed"
echo "Mapped S2"


# Method 2
bedtools closest -d -a "${DATA_PATH}/chipseq_data/${CELLS}/${CHIP_FILE_NAME}_no_black_sorted.bed" \
-b "${DATA_PATH}/tss_mapping/${CELLS}_gene_TSS_50_expressed_50_occur_sorted.bed" > "${DATA_PATH}/mappings/${CHIP_FILE_NAME}_no_black_sorted_mapped_top50_expressed.bed"
echo "Mapped M2"


# Distance filtering of interactions
Rscript ./pipeline/filter_tf_target_interactions.R  $CELLS $DATA_PATH
echo "Filtered"


# Annotation with PPI, ATAC and DNAse
Rscript ./pipeline/annotate.R $CELLS $ATAC_FILE $DNASE_FILE $DATA_PATH
echo "Annotated"

# Annotation with TFBS
Rscript ./pipeline/annotate_found_motifs.R $CELLS $DATA_PATH
echo "Annotated motifs"

# Annotation with ENCODE cCREs"
Rscript ./pipeline/annotate_encode_ccres.R $CELLS $DATA_PATH
echo "Annotated cCREs"

# Running cleanup
filename=data/regulons/*$CELLS*ccres.tsv
Rscript ./pipeline/dataset_cleanup.R $filename

# Final file cleaning up
cd data/regulons
ls TF_target_mapping_filtered_merged_*$CELLS*.tsv | grep -v cleaned | xargs rm
