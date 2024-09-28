#!/bin/bash

source ~/.bashrc
conda activate chip_seq_reg


CELLS="GM23248"
DATA_PATH="./data"

chip_file=$DATA_PATH"/regulons/TF_target_mapping_filtered_merged_"$CELLS"_with_ppi_with_dnase_with_atac.tsv"

p=$DATA_PATH"/motifs/"
for i in $p/*motif;
do
	file=${i##*/}  
	base=${file%.*}

	echo $base

	# Get respective TF peak files
	Rscript ./pipeline/select_tf_peaks.R $base $chip_file $DATA_PATH"/motifs_ann/"$CELLS
        
        # Count the number of lines in the file
        line_count=$(wc -l < $DATA_PATH"/motifs_ann/"$CELLS"/"$base"_temp_peak_file.bed")

        # Check if the line count is equal to 1
        if [ "$line_count" -eq 0 ]; then
		rm $DATA_PATH"/motifs_ann/"$CELLS"/"$base"_temp_peak_file.bed"
		continue
        fi

	# Annotate motifs
	annotatePeaks.pl $DATA_PATH"/motifs_ann/"$CELLS"/"$base"_temp_peak_file.bed" data/hg38/ \
	            -m $i > $DATA_PATH"/motifs_ann/"$CELLS"/"$base"_out_hoco.txt"

	# Count the number of lines in the file
	line_count=$(wc -l < $DATA_PATH"/motifs_ann/"$CELLS"/"$base"_out_hoco.txt")

	# Check if the line count is equal to 1
	if [ "$line_count" -eq 1 ]; then
		rm $DATA_PATH"/motifs_ann/"$CELLS"/"$base"_out_hoco.txt"
	fi
	rm $DATA_PATH"/motifs_ann/"$CELLS"/"$base"_temp_peak_file.bed"
done
