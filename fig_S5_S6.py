import pandas as pd
from tqdm import tqdm
import numpy as np
import igraph as ig
from create_ground_networks import *
import matplotlib.pyplot as plt
from collections import defaultdict
import gseapy as gp
plt.rcParams["figure.figsize"] = 5, 7

metrics_df_comb = pd.read_csv("data/network_analysis/centrality_metrics_m2_full_network.tsv", header=0, sep="\t")

mets = ["authority_score",
	"closeness",
	"eigenvector_centrality",
	"hub_score",
	"in_degree",
	"out_degree",
	"pagerank",
	"vertex_betweenness"]

for met in mets:
    print(met)
    fig, axs = plt.subplots()
    check_cor(metrics_df_comb.dropna(subset=[met, met + "_with_motifs"])[met].to_numpy(), metrics_df_comb.dropna(subset=[met, met + "_with_motifs"])[met + "_with_motifs"].to_numpy(), 
	x_lab="no_extra", y_lab="with_motifs", title=met, ax=axs)
    plt.gcf().set_size_inches(7, 5)
    plt.savefig("plots/compare_networks_with_motifs_" + met + ".png", dpi=600)
    plt.close()


for met in mets:
    print(met)
    fig, axs = plt.subplots()
    check_cor(metrics_df_comb.dropna(subset=[met, met + "_with_ppi"])[met].to_numpy(), metrics_df_comb.dropna(subset=[met, met + "_with_ppi"])[met + "_with_ppi"].to_numpy(),
        x_lab="no_extra", y_lab="with_ppi", title=met, ax=axs)
    plt.gcf().set_size_inches(7, 5)
    plt.savefig("plots/compare_networks_with_ppi_" + met + ".png", dpi=600)
    plt.close()
