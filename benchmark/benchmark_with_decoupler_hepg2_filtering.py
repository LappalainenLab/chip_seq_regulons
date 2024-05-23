import pandas as pd
import numpy as np
import decoupler as dc
import os, sys
sys.path.append("/proj/lappalainen_lab1/users/marii/chip_seq_ann")

# Read the knockout expression matrix and metadata
mat = pd.read_csv('data/KnockTF2/knockTF_expr.csv', index_col=0)
obs = pd.read_csv('data/KnockTF2/knockTF_meta.csv', index_col=0)

# Read TF-target mapping data
# data = pd.read_table('data/regulons/HepG2_all_regulons.tsv', sep="\t")
data = pd.read_table("data/regulons/TF_target_mapping_filtered_merged_HepG2_with_ppi_with_dnase_with_atac_with_motifs_with_ccres_cleaned.tsv", sep="\t")

# Filter the knockout experiments based on logFC and cell line
msk = obs['logFC'] < -1
hepg2_mask = obs["Biosample.Name"] == "HepG2"
mat = mat[msk & hepg2_mask]
obs = obs[msk & hepg2_mask]

# Extract TF-target interactions for S2Mb, M2Kb and S2Kb methods without applied filters
s2mb = data.loc[data.is_S2Mb, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
m2kb = data.loc[data.is_M2Kb, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
s2kb = data.loc[data.is_S2Kb, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])

# Extract TF-target interactions for S2Mb, M2Kb and S2Kb methods with applied TF binding motif filter
s2mb_motif = data.loc[data.is_S2Mb & ((data.n_motifs > 0) | (data.n_motifs_hoco > 0)), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
m2kb_motif = data.loc[data.is_M2Kb & ((data.n_motifs > 0) | (data.n_motifs_hoco > 0)), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
s2kb_motif = data.loc[data.is_S2Kb & ((data.n_motifs > 0) | (data.n_motifs_hoco > 0)), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])

# Extract TF-target interactions for S2Mb, M2Kb and S2Kb methods with applied open chromatin filter
s2mb_open = data.loc[data.is_S2Mb & (data.is_dnase | data.is_atac), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
m2kb_open = data.loc[data.is_M2Kb & (data.is_dnase | data.is_atac), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
s2kb_open = data.loc[data.is_S2Kb & (data.is_dnase | data.is_atac), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])

# Extract TF-target interactions for S2Mb and M2Kb methods with applied cCRE filter
s2mb_ccre = data.loc[data.is_S2Mb & (data.is_pls | data.is_pels | data.is_dels), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
m2kb_ccre = data.loc[data.is_M2Kb & (data.is_pls | data.is_pels | data.is_dels), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
s2kb_ccre = data.loc[data.is_S2Kb & (data.is_pls | data.is_pels | data.is_dels), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])

# Rename columns
for dataset in [s2mb, m2kb, s2kb, s2mb_motif, m2kb_motif, s2kb_motif, s2mb_open, \
		m2kb_open, s2kb_open, s2mb_ccre, m2kb_ccre, s2kb_ccre]:
    dataset.columns = ["source", "target"]

# Create randomized networks by shuffling TF<->target gene associations in the S2Mb, M2Kb and S2Kb networks
s2mb_rand = dc.shuffle_net(s2mb, target='target', weight=None).drop_duplicates(['source', 'target'])
m2kb_rand =  dc.shuffle_net(m2kb, target='target', weight=None).drop_duplicates(['source', 'target'])
s2kb_rand =  dc.shuffle_net(s2kb, target='target', weight=None).drop_duplicates(['source', 'target'])


# Build dictionary of networks to test
nets = {
    'S2Mb': s2mb,
    'M2Kb': m2kb,
    'S2Kb': s2kb,
    'S2Mb_motif': s2mb_motif, 
    'M2Kb_motif': m2kb_motif, 
    'S2Kb_motif': s2kb_motif, 
    'S2Mb_open': s2mb_open,
    'M2Kb_open': m2kb_open,
    'S2Kb_open': s2kb_open,
    'S2Mb_ccre': s2mb_ccre, 
    'M2Kb_ccre': m2kb_ccre, 
    'S2Kb_ccre': s2kb_ccre,
    'S2Mb_random': s2mb_rand,
    'M2Kb_random': m2kb_rand,
    'S2Kb_random': s2kb_rand}

# Define decoupler parameters for each network
decouple_kws = {
    'S2Mb': {'min_n': 5, 'weight': None},
    'M2Kb': {'min_n': 5, 'weight': None},
    'S2Kb': {'min_n': 5, 'weight': None},
    'S2Mb_motif': {'min_n': 5, 'weight': None},
    'M2Kb_motif': {'min_n': 5, 'weight': None},
    'S2Kb_motif': {'min_n': 5, 'weight': None},
    'S2Mb_open': {'min_n': 5, 'weight': None},
    'M2Kb_open': {'min_n': 5, 'weight': None},
    'S2Kb_open': {'min_n': 5, 'weight': None},
    'S2Mb_ccre': {'min_n': 5, 'weight': None},
    'M2Kb_ccre': {'min_n': 5, 'weight': None},
    'S2Kb_ccre': {'min_n': 5, 'weight': None},
    'S2Mb_random': {'min_n': 5, 'weight': None},
    'M2Kb_random': {'min_n': 5, 'weight': None},
    'S2Kb_random': {'min_n': 5, 'weight': None}
}

# Run benchmark pipeline
df = dc.benchmark(mat, obs, nets, perturb='TF', sign=-1, verbose=True, decouple_kws=decouple_kws, min_exp=2)

# Save the results to a TSV file
df.to_csv("data/2-plot_decoupler_filter_benchmark_across_methods/hepg2_filtering_benchmark.tsv", sep="\t", index=False)
