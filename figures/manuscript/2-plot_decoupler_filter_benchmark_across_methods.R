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


data = fread("data/enrich_analysis/mcf7_filtering_benchmark.tsv")

data %>% separate(net, c("group", "filter"), "_") %>%
        mutate(filter = replace_na(filter, "no")) %>%
        filter(!(filter %in% c("exp+atac", "exp")))-> data

data %>% pivot_wider(names_from = metric, values_from = score, values_fn=mean) %>%
        filter(method == "consensus_estimate")-> to_plot

to_plot$group = factor(to_plot$group, levels=c("S2Mb", "M2Kb", "S2Kb"))
to_plot$filter =  factor(to_plot$filter, levels=c("random", "no", "exp", "dnase", "atac"))

p1 = ggplot(to_plot , aes(x = filter, y = mcauprc, group = filter, fill=group)) +
	geom_col(position = position_dodge2(), color = "black") +
        scale_fill_okabeito(palette = "black_original", order=c(6, 4, 7)) +
        theme_pubr(legend = "right") +
        theme(axis.text = element_text(size=9),
		axis.text.x = element_text(size=9, angle = 30, vjust=0.75),
                legend.background = element_blank(),
                axis.title = element_text(size=9),
                legend.text = element_text(size=9),
                legend.title = element_text(size=9),
                legend.spacing.y = unit(0.1, 'pt')) +
	xlab("Filter") +
	ylab("MCAUPRC") +
        guides(fill=guide_legend(ncol=1,byrow=T,
                           keywidth=3,
                           keyheight=3,
                           default.unit="mm",
			   title="Cell line"))


ggsave(plot = p1, filename = "plots/decoupler_bench_filtering_mcf7_across_methods.svg", width=100, height=70, units="mm", dpi=720, device = "svg")
