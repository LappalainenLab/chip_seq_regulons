library(tidyr)
library(psych)
library(foreach)
library(ggplot2); theme_set(theme_bw())
library(data.table)
library(dplyr)

source("fig_style.R")

Cor_dt = fread("data/network_analysis/summary_correlation_with_centrality.tsv")

p <- ggplot(Cor_dt, aes(x= metric, y = cent_met)) +
  geom_point(aes(color=r, size=log10_padj)) +
  scale_color_gradient2(low = "#377EB8", high = "#FF7F00", midpoint = 0) +
  theme(panel.grid = element_blank()) +
  ylab("") +
  xlab("") +
  scale_y_discrete(labels = c("BetW", "HubS", "Eigen")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = p, filename = "plots/summary_correlation_with_centrality.png", device = "png", dpi = 600, width = 4, height = 4)
