library(tidyr)
library(foreach)
library(psych)
library(foreach)
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
source("fig_style.R")


to_plot = fread("data/network_analysis/direct_target_overlap_in_pert.tsv")

print(quantile(to_plot$fraction))

plot = ggplot(to_plot, aes(x = fraction)) +
        geom_histogram() +
	ylab("TFs") +
        xlab("% overlap") +
        gtex_v8_figure_theme()
ggsave(plot=plot, filename = "plots/TF_enrich_overlaps.png", device = "png", dpi = 600, width = 4, height = 4)
