import pandas as pd
import numpy as np
import decoupler as dc
import os, sys
sys.path.append("/proj/lappalainen_lab1/users/marii/chip_seq_ann")

# Read the knockout expression matrix and metadata
mat = pd.read_csv('data/KnockTF2/knockTF_expr.csv', index_col=0)
obs = pd.read_csv('data/KnockTF2/knockTF_meta.csv', index_col=0)

# Read TF-target mapping data
data = pd.read_table('data/regulons/K562_regulon.tsv', sep="\t")

# Filter the knockout experiments based on logFC and cell line
msk = obs['logFC'] < -1
k562_mask = obs["Biosample.Name"] == "K562"
mat = mat[msk & k562_mask]
obs = obs[msk & k562_mask]

# Extract TF-target interactions for S2Mb and M2Kb methods with applied ATAC filter
s2mb_atac = data.loc[data.is_method_1 & data.is_atac, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])
m2kb_atac = data.loc[data.is_method_2 & data.is_atac, ["tf", "gene_symbol"]].drop_duplicates(["tf", "gene_symbol"])

# Rename columns
s2mb_atac.columns = ["source", "target"]
m2kb_atac.columns = ["source", "target"]

# Read and preprocess other regulons
regnet = pd.read_table("data/trans_networks/regnet_human.tsv", sep="\t", header=None).iloc[:,[0, 2]].drop_duplicates()
regnet.columns = ["source", "target"]

trrust = pd.read_table("data/trans_networks/trrust_rawdata.human.tsv", sep="\t", header=None).iloc[:,0:2].drop_duplicates()
trrust.columns = ["source", "target"]

chip_atlas = pd.read_table("data/trans_networks/ChIP-Atlas_target_genes_K562.tsv", sep="\t").drop_duplicates()
chip_atlas.columns = ["source", "target"]

collectri = dc.get_collectri(organism='human', split_complexes=False)
dorothea = dc.get_dorothea(organism='human')

tfs = set.intersection(
    set(collectri.source), set(m2kb_atac.source), set(regnet.source),
    set(trrust.source), set(chip_atlas.source), set(dorothea.source)
)

# Filter regulons to shared TFs
s2mb_atac = s2mb_atac[s2mb_atac.source.isin(tfs)]
m2kb_atac = m2kb_atac[m2kb_atac.source.isin(tfs)]
chip_atlas = chip_atlas[chip_atlas.source.isin(tfs)]
regnet = regnet[regnet.source.isin(tfs)]
trrust = trrust[trrust.source.isin(tfs)]
collectri = collectri[collectri.source.isin(tfs)]
dorothea = dorothea[dorothea.source.isin(tfs)]

# Build dictionary of networks to test
nets = {
    'S2Mb': s2mb_atac, 
    'M2Kb': m2kb_atac, 
    "ChIP-Atlas": chip_atlas,
    'CollecTri': collectri,
    'Dorothea': dorothea,
    'RegNet': regnet,
    'TRRUST': trrust,
}

# Define decoupler parameters for each network
decouple_kws = {
    'S2Mb': {'min_n': 5, 'weight': None},
    'M2Kb': {'min_n': 5, 'weight': None},
    'ChIP-Atlas': {'min_n': 5, 'weight': None},
    'Dorothea': {'min_n': 5, 'weight': None},
    'RegNet': {'min_n': 5, 'weight': None},
    'TRRUST': {'min_n': 5, 'weight': None},
    'CollecTri': {'min_n': 5, 'weight': None}
}

# Run benchmark pipeline
df = dc.benchmark(mat, obs, nets, perturb='TF', sign=-1, verbose=True, decouple_kws=decouple_kws)

# Get some statistics
collectri = df.loc[df.net == "CollecTri"]
print(collectri.loc[(collectri.method == "consensus_estimate") & (collectri.metric == "mcauprc")].describe())
# CollecTri mean MCAUPRC 0.627082

m2kb = df.loc[df.net == "M2Kb"] 
print(m2kb.loc[(m2kb.method == "consensus_estimate") & (m2kb.metric == "mcauprc")].describe())
# M2Kb mean MCAUPRC 0.587219

# Save the results to a TSV file
df.to_csv("data/enrich_analysis/k562_comparison_benchmark.tsv", sep="\t", index=False)
