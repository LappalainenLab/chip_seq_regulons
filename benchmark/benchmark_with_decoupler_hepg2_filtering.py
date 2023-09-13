import pandas as pd
import numpy as np
import decoupler as dc
import os, sys
sys.path.append("/proj/lappalainen_lab1/users/marii/chip_seq_ann")

# Read the knockout expression matrix and metadata
mat = pd.read_csv('data/KnockTF2/knockTF_expr.csv', index_col=0)
obs = pd.read_csv('data/KnockTF2/knockTF_meta.csv', index_col=0)

# Read TF-target mapping data
data = pd.read_table('data/regulons/TF_target_mapping_filtered_merged_HepG2_with_motifs_with_ppi_with_dnase_with_atac_with_dist_score.tsv', sep="\t")

# Filter the knockout experiments based on logFC and cell line
msk = obs['logFC'] < -1
hepg2_mask = obs["Biosample.Name"] == "HepG2"
mat = mat[msk & hepg2_mask]
obs = obs[msk & hepg2_mask]

# Extract TF-target interactions for S2Mb, M2Kb and S2Kb methods without applied filters
s2mb = data.loc[data.is_method_1, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
m2kb = data.loc[data.is_method_2, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
s2kb = data.loc[data.is_method_3, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])

# Extract TF-target interactions for S2Mb, M2Kb and S2Kb methods with applied TF binding motif filter
s2mb_motif = data.loc[data.is_method_1 & (data.n_motifs_method1 > 0), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
m2kb_motif = data.loc[data.is_method_2 & (data.n_motifs_method2 > 0), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
s2kb_motif = data.loc[data.is_method_3 & (data.n_motifs_method3 > 0), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])

# Extract TF-target interactions for S2Mb, M2Kb and S2Kb methods with applied expression filter
s2mb_exp = data.loc[data.is_method_1 & (data.tpm_total > 0), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
m2kb_exp = data.loc[data.is_method_2 & (data.tpm_total > 0), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
s2kb_exp = data.loc[data.is_method_3 & (data.tpm_total > 0), ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])

# Extract TF-target interactions for S2Mb, M2Kb and S2Kb methods with applied DNAse filter
s2mb_dnase = data.loc[data.is_method_1 & data.is_dnase, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
m2kb_dnase = data.loc[data.is_method_2 & data.is_dnase, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
s2kb_dnase = data.loc[data.is_method_3 & data.is_dnase, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"]) 

# Extract TF-target interactions for S2Mb and M2Kb methods with applied ATAC filter
s2mb_atac = data.loc[data.is_method_1 & data.is_atac, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
m2kb_atac = data.loc[data.is_method_2 & data.is_atac, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
s2kb_atac = data.loc[data.is_method_3 & data.is_atac, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])

# Rename columns
for dataset in [s2mb, m2kb, s2kb, s2mb_motif, m2kb_motif, s2kb_motif, s2mb_dnase, \
		m2kb_dnase, s2kb_dnase, s2mb_atac, m2kb_atac, s2kb_atac, s2mb_exp, \
		m2kb_exp, s2kb_exp]:
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
    'S2Mb_dnase': s2mb_dnase,
    'M2Kb_dnase': m2kb_dnase,
    'S2Kb_dnase': s2kb_dnase,
    'S2Mb_atac': s2mb_atac, 
    'M2Kb_atac': m2kb_atac, 
    'S2Kb_atac': s2kb_atac,
    'S2Mb_exp': s2mb_exp,
    'M2Kb_exp': m2kb_exp,
    'S2Kb_exp': s2kb_exp,
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
    'S2Mb_exp': {'min_n': 5, 'weight': None},
    'M2Kb_exp': {'min_n': 5, 'weight': None},
    'S2Kb_exp': {'min_n': 5, 'weight': None},
    'S2Mb_dnase': {'min_n': 5, 'weight': None},
    'M2Kb_dnase': {'min_n': 5, 'weight': None},
    'S2Kb_dnase': {'min_n': 5, 'weight': None},
    'S2Mb_atac': {'min_n': 5, 'weight': None},
    'M2Kb_atac': {'min_n': 5, 'weight': None},
    'S2Kb_atac': {'min_n': 5, 'weight': None},
    'S2Mb_random': {'min_n': 5, 'weight': None},
    'M2Kb_random': {'min_n': 5, 'weight': None},
    'S2Kb_random': {'min_n': 5, 'weight': None}
}

# Run benchmark pipeline
df = dc.benchmark(mat, obs, nets, perturb='TF', sign=-1, verbose=True, decouple_kws=decouple_kws)

# Save the results to a TSV file
df.to_csv("data/enrich_analysis/hepg2_filtering_benchmark.tsv", sep="\t", index=False)
