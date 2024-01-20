# chip_seq_regulons

This repo contains code used in preprint Specifying cellular context of transcription factor regulons for exploring context-specific gene regulation prog$

### Installation

```
git clone https://github.com/LappalainenLab/chip_seq_regulons.git
conda create --name chip_seq_reg --file environment.txt
conda activate chip_seq_reg
```

### Annotation pipeline:

Contains scripts to run the S2Mb, S2Kb and M2Kb TF-target gene annotation pipelines

1. Put all necessary data files, i.e. ReMap and RNA-Seq data, to respective data directories (see data directory structure below)
2. Run annotation pipeline with 
```
bash pipeline/run_one_pipeline.sh
```

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
