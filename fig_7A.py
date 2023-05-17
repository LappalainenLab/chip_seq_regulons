import json
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns
import matplotlib.pyplot as plt
import gseapy as gp

raw_data = pd.read_csv("data/network_analysis/centrality_metrics_m2_full_network.tsv", sep="\t")
raw_data.index = raw_data.name
with_motifs = raw_data.loc[:, raw_data.columns.str.endswith("motifs")]
with_ppi = raw_data.loc[:, raw_data.columns.str.endswith("ppi")]
init_cols = list(set(raw_data).difference(set(with_motifs)))
init_cols = list(set(init_cols).difference(set(with_ppi)))
init_net = raw_data.loc[:, init_cols]

init_net = init_net.drop(["name", "n_components"], axis=1)
with_motifs = with_motifs.drop(["n_components_with_motifs"], axis=1)
with_ppi = with_ppi.drop(["n_components_with_ppi"], axis=1)

cor_mat = pd.DataFrame(sp.stats.spearmanr(init_net).correlation)
cor_mat.columns = ["Eigen", "AuthS", "OutD", "Betw", "InD", "Page", "HubS", "Closs"] 
cor_mat.index = ["Eigen", "AuthS", "OutD", "Betw", "InD", "Page", "HubS", "Closs"]

sns.clustermap(cor_mat, 
	annot=True, 
	annot_kws={"size": 7}, 
	cmap = sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True))
plt.gcf().set_size_inches(5, 5)
plt.savefig("plots/centrality_corr_init_fig.png",  dpi=600)
plt.close()
