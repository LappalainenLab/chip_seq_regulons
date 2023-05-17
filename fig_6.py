import seaborn as sns
import igraph as ig
from matplotlib import cm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable 
from matplotlib.colors import Normalize

d2n = ig.Graph.Read(f="data/enrich_analysis/d2n_fig_6_graph.gml", format="gml")

norm1 = ScalarMappable(norm=Normalize(min(d2n.vs["color"]), max(d2n.vs["color"])), cmap=sns.color_palette("icefire", as_cmap=True, n_colors=len(d2n.vs["color"])))
palette = ig.drawing.colors.PrecalculatedPalette(sns.color_palette("icefire", as_cmap=False, n_colors=len(d2n.vs["color"])))


fig, ax = plt.subplots()
ig.plot(d2n, target=ax, layout=d2n.layout_random(),
    palette = palette,
    vertex_frame_width=0.2,
    vertex_label_dist = 1,
    vertex_label_size=4.5,
    vertex_size = d2n.vs["degree"],
    vertex_label=d2n.vs["name"],
    vertex_color=d2n.vs["color"],
    edge_width=0.05,
    bbox = (5000, 4000),
    margin = 0,
    keep_aspect_ratio=True,
    edge_arrow_size=0.001)
fig.colorbar(norm1, orientation="vertical", label='Degree difference')
plt.gcf().set_size_inches(5, 5)
plt.savefig("plots/d2n_network_degree_plus_priority_random.png", dpi=600)
plt.close()
