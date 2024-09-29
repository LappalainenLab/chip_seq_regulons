# chip_seq_regulons

Specifying cellular context of transcription factor regulons for exploring context-specific gene regulation programs

### Installation

The installation of the dependencies required for the repository, requires `conda`. Further instructions on installing Anaconda could be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

```
git clone https://github.com/LappalainenLab/chip_seq_regulons.git
cd chip_seq_regulons/
conda env create --file environment.yaml
conda activate chip_seq_reg
```

### Annotation pipeline:

Contains scripts to run the S2Mb, S100Kb, S2Kb, M100Kb, and M2Kb TF-target gene annotation pipelines

1. Put all necessary data files, i.e. ReMap and RNA-Seq data, to respective data directories (see data directory structure below)
2. Run annotation pipeline with 
```
bash pipeline/run_one_pipeline.sh sample_cell sample_RNA_1,sample_RNA_2 sample_DNAse sample_ATAC
```
where sample_cell is the name of the cell line; sample_RNA_1,sample_RNA_2 are comma-separated RNA-Seq replicate data accessions (i.e. file names); sample_DNAse is the DNAse-Seq data accession (i.e. file name); sample_ATAC is the ATAC-Seq data accession (i.e. file name).


### Benchmarking:

To replicate manuscript results:

1. Collect regulons from the respective Zenodo repository.
2. Run benchmarking scripts with
```
# Compare regulons
python benchmark/benchmark_with_decoupler_{cell_line}_comparison.py

# Compare filtering strategies
python benchmark/benchmark_with_decoupler_{cell_line}_filtering.py
```
3. To extract regulons' statistics
```
Rscript benchmark/data_summary_stats.R
```

4. To perform enrichment analysis of TF-target gene pairs in PPI interactions
```
Rscript benchmark/test_ppi_enrich.R
```

5. To perform enrichment analysis of TF-target gene pairs in biological networks for the K562 cell line (see Methods section of the manuscript for deatils)
```
Rscript benchmark/test_net_enrich.R
```

6. To perform enrichment analysis of TF-target gene pairs in coexpression networks for the K562 cell line
```
Rscript benchmark/test_coexp_enrich.R
```

### Case studies:

Notebooks used to conduct cancer case studies are located at `case_studies`

### Figures:

Scripts to reproduce figures from the manuscripts can be run by
```
Rscript figures/manuscript/{script_name}.R
```

## Citation
> Minaeva, M., Domingo, J., Rentzsch, P., & Lappalainen, T. (2024). 
> Specifying cellular context of transcription factor regulons for exploring context-specific gene regulation programs.
> bioRxiv, 2023-12. https://doi.org/10.1101/2023.12.31.573765
