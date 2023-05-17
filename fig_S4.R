library(foreach)
library(gridExtra)
library(stringr)
library(data.table)
library(ggstatsplot)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

here::i_am("fig_style.R")
source("fig_style.R")



to_plot = fread("data/enrich_analysis/target_dist_for_dbs.tsv")

plot = ggplot(to_plot, aes(x = N)) +
	geom_histogram() +
	xlab("Number of target genes per TF") +
	ylab("Genes") +
	facet_wrap(~ db, scales = "free") +
	gtex_v8_figure_theme()

ggsave(plot=plot, filename = "plots/target_dist_for_dbs.png", device = "png", dpi = 600, width = 5, height = 4, limitsize = FALSE)

