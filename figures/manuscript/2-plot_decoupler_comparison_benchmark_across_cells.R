library(ggplot2); theme_set(theme_bw())
library(data.table)
library(dplyr)
library(ggpubr)
library(ggpattern)
library(here)
library(see)
library(grid)
library(gridExtra)
library(tidyr)
library(see)
library(cowplot)
palette = c()

mcf7 = fread("data/enrich_analysis/mcf7_comparison_benchmark.tsv")
k562 = fread("data/enrich_analysis/k562_comparison_benchmark.tsv")
hepg2 = fread("data/enrich_analysis/hepg2_comparison_benchmark.tsv")
mcf7[, "cells"] = "MCF-7"
k562[, "cells"] = "K-562"
hepg2[, "cells"] = "HepG2"

bind_rows(k562, mcf7, hepg2) -> data

data %>% pivot_wider(names_from = metric, values_from = score, values_fn=mean) %>%
	filter(method == "consensus_estimate")-> to_plot
to_plot$net = factor(to_plot$net, levels=c("S2Mb", "M2Kb", "CollecTri", "Dorothea", "ChIP-Atlas", "RegNet", "TRRUST"))
to_plot$cells = factor(to_plot$cells, levels=c("K-562", "MCF-7", "HepG2"))

p1 = ggplot(to_plot , aes(x = net, y = mcauprc, group = net, fill=cells)) +
	geom_col(position = position_dodge2(), color = "black") +
        scale_fill_okabeito(palette = "black_original", order=c(3, 4, 5)) +
        theme_pubr(legend = "right") +
        theme(axis.text = element_text(size=9),
		axis.text.x = element_text(size=9, angle = 30, vjust=0.75),
                legend.background = element_blank(),
                axis.title = element_text(size=9),
                legend.text = element_text(size=9),
                legend.title = element_text(size=9),
                legend.spacing.y = unit(0.1, 'pt')) +
	xlab("Regulon") +
	ylab("MCAUPRC") +
        guides(fill=guide_legend(ncol=1,byrow=T,
                           keywidth=3,
                           keyheight=3,
                           default.unit="mm",
			   title="Cell line"))


ggsave(plot = p1, filename = "plots/decoupler_bench_scatter_comparison_all_lines.svg", width=100, height=70, units="mm", dpi=720, device = "svg")
