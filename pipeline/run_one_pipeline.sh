#!/bin/bash
#SBATCH -A snic2022-5-560
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH -J TF_target_mapping

cd /proj/lappalainen_lab1/users/marii/chip*/

# Loading modules
module load R_packages/4.1.1


# Configure
CELLS="try_K562"
RNA_FILES="ENCFF823JHX,ENCFF940GYO"
DNASE_FILE="ENCFF274YGF"
ATAC_FILE="ENCFF558BLC"

echo "Starting the pipeline"

# Extract highest expressed TSS
Rscript ./pipeline/extract_tss.R $CELLS $RNA_FILES
echo "Extracted TSS 1"


# Extract top 50% expressed TSS 
Rscript ./pipeline/extract_multiple_tss.R $CELLS $RNA_FILES
echo "Extracted TSS 2"


# Intersect ChIP-Seq with Blacklist regions
bedtools intersect -wa \
			-a "data/chipseq_data/${CELLS}/remap2022_${CELLS}_macs2_hg38_v1_0.bed" \
			-b "data/blacklist/ENCFF023CZC.bed" \
			-v > "data/chipseq_data/${CELLS}/remap2022_${CELLS}_macs2_hg38_v1_0_no_black.bed"
echo "Removed blacklist regions"

if [ ! -d "data/tss_mapping" ] 
then
	mkdir data/tss_mapping
fi

if [ ! -d "data/mappings" ] 
then 
        mkdir data/mappings
fi 

# Sort BED files
sort -k1,1 -k2,2n "data/chipseq_data/${CELLS}/remap2022_${CELLS}_macs2_hg38_v1_0_no_black.bed" > \
		"data/chipseq_data/${CELLS}/remap2022_${CELLS}_macs2_hg38_v1_0_no_black_sorted.bed"

sort -k1,1 -k2,2n "data/tss_mapping/${CELLS}_gene_TSS_highest_expressed.bed" > \
                "data/tss_mapping/${CELLS}_gene_TSS_highest_expressed_sorted.bed"

sort -k1,1 -k2,2n "data/tss_mapping/${CELLS}_gene_TSS_50_expressed_50_occur.bed" > \
                "data/tss_mapping/${CELLS}_gene_TSS_50_expressed_50_occur_sorted.bed"
echo "Sorted files"


# Method 1 and 3
bedtools closest -d -a "data/chipseq_data/${CELLS}/remap2022_${CELLS}_macs2_hg38_v1_0_no_black_sorted.bed" \
                    -b "data/tss_mapping/${CELLS}_gene_TSS_highest_expressed_sorted.bed" > "data/mappings/remap2022_${CELLS}_macs2_hg38_v1_0_no_black_sorted_mapped_highest_expressed.bed"
echo "Mapped S2"


# Method 2
bedtools closest -d -a "data/chipseq_data/${CELLS}/remap2022_${CELLS}_macs2_hg38_v1_0_no_black_sorted.bed" \
-b "data/tss_mapping/${CELLS}_gene_TSS_50_expressed_50_occur_sorted.bed" > "data/mappings/remap2022_${CELLS}_macs2_hg38_v1_0_no_black_sorted_mapped_top50_expressed.bed"
echo "Mapped M2"


# Distance filtering of interactions
Rscript ./pipeline/filter_tf_target_interactions.R  $CELLS
echo "Filtered"


# Annotation with PPI, ATAC and DNAse
Rscript ./pipeline/annotate.R $CELLS $ATAC_FILE $DNASE_FILE
echo "Annotated"
