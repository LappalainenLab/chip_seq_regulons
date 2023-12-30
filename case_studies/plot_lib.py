import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.patches as mpatches
import matplotlib
from plotnine import (ggplot, aes, geom_point, 
                      position_dodge2, facet_wrap, 
                      ylab, ggtitle, scale_size, 
                      theme_bw, ggsave, theme,
                      element_text, element_blank, xlab,
                      guides, guide_legend, scale_color_manual,
                      geom_vline, element_line, scale_y_discrete, coord_cartesian
                     )
import textwrap



def plot_enrich(df, palette, size="Combined Score"):
    def wraping_func(text):
        #return [textwrap.wrap(wraped_text, 30) for wraped_text in text]
        return ["\n".join(textwrap.wrap(wraped_text, 25)) for wraped_text in text]

    p = (ggplot(df, aes(x="-log10(Adj. P-value)", y="Term", color="Network", size=size)) + 
        geom_point() +
        geom_vline(xintercept=1.3, linetype="dashed") +
        scale_y_discrete(breaks=df['Term'].unique().tolist(), 
                           labels=wraping_func) +
        ylab("") + 
        theme_bw() +
        scale_color_manual(palette) +
            theme(axis_text = element_text(size=9),
                    axis_text_x = element_text(size=9),
                    legend_background = element_blank(),
                    axis_title = element_text(size=9),
                    legend_text = element_text(size=9),
                    legend_title = element_text(size=9),
                    legend_position = "bottom",
                    rect=element_blank(),
                    axis_line=element_line(),
                    axis_ticks=element_line()) +
            xlab("-log10(Adj. P-value)") +
            guides(color=guide_legend(ncol=1,
                                      byrow=True,
                               title="Regulon",
                               override_aes = {"size": 6},
                                    title_position = "top",
                            title_hjust =0.5),
                  size=guide_legend(ncol=1,
                                      byrow=True,
                               title="Combined Score",
                                    title_position = "top",
                            title_hjust =0.5)) + 
            scale_size(range=[3, 8]))
    
    return p


def plot_bubbles(megaheat, scale=0.35):  
    megaheat = megaheat.sort_values(["score", "Network"])
    fig, ax = plt.subplots(1, 1,  dpi=720, #figsize=(5.924, 6.693)) #poster figsize
                           figsize=(2.55906, 4.92126), )
    ns = (megaheat["log_p"].values * scale * plt.rcParams["lines.markersize"]) ** 2
    ax.grid(axis='x')
    nc = colors.TwoSlopeNorm(vcenter=0)
    scatter = ax.scatter(x=megaheat["Network"], y=megaheat["TF"], s=ns, c=megaheat["score"], 
                cmap="RdBu_r", norm=nc,  clip_on=False)
    ax.set_axisbelow(True)
    ax.set_xlabel("Regulon", fontsize=7)
    yticks = ax.get_yticklabels(minor=False, which=None)
    ax.set_yticklabels(yticks, fontdict={'fontsize': 9})
    
    xticks = ax.get_xticklabels(minor=False, which=None)
    ax.set_xticklabels(xticks, rotation=30, ha='right', rotation_mode='anchor', fontdict={'fontsize': 9})
    
    handles, labels = scatter.legend_elements(
            prop="sizes",
            num=3,
            fmt="{x:.2f}",
            func=lambda s: np.sqrt(s) / plt.rcParams["lines.markersize"] / scale)

    ax.legend(
            handles,
            [x[:-1] for x in labels],
            title="-log10(p_value)",
            frameon=False,           
            bbox_to_anchor=(-0.2, -0.475),
            loc="lower left",
            labelspacing= 1.0,
            columnspacing = 0.05,
            fontsize='9',
            title_fontsize='9',
            ncols = len(labels),
            handletextpad=0.5
        )


     # Add colorbar
    clb = fig.colorbar(
            scatter,
            shrink=0.25,
            aspect=10,
            orientation='horizontal',
            anchor=(1.0, 0.25),
            location="bottom"
        )
    clb.ax.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
    fig.axes[1].tick_params(axis="x", labelsize=9)
    tick_locator = matplotlib.ticker.MaxNLocator(nbins=5)
    clb.locator = tick_locator
    clb.update_ticks()
    clb.update_ticks()
    clb.ax.set_title("Score", loc="center", size='9')
    ax.margins(x=0.25, y=0.05)
    sns.despine(fig=fig, ax=ax)
    return ax