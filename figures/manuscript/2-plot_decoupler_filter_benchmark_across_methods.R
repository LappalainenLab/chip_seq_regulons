library(ggplot2); theme_set(theme_bw())
library(data.table)
library(dplyr)
library(ggpubr)
library(here)
library(see)
library(grid)
library(gridExtra)
library(tidyr)
library(see)
library(cowplot)
library(foreach)



cell_line = c("k562", "hepg2", "mcf7")

foreach(cl = cell_line) %do% {

	input_file <- sprintf("data/2-plot_decoupler_filter_benchmark_across_methods/%s_filtering_benchmark.tsv", cl)
	data = fread(input_file)
	
	data %>% separate(net, c("group", "filter"), "_") %>%
        	mutate(filter = replace_na(filter, "no")) %>%
        	filter(!(filter %in% c("exp+atac", "exp", "polr2a")))-> data

	data %>% pivot_wider(names_from = metric, values_from = score, values_fn=mean) %>%
        	filter(method == "consensus_estimate")-> to_plot

	data %>% filter((method == "consensus_estimate") & (metric == "mcauprc")) -> data

	data$group = factor(data$group, levels=c("S2Mb", "M2Kb", "S2Kb"))
	data$filter =  factor(data$filter, levels=c("random", "no", "motif", "dnase", "atac"), labels=c("shuffled \nnetwork", "no filter", "motif", "dnase", "atac"))



	p1 = ggbarplot(data, x="filter", y="score", fill="group", group="filter", add = c("mean_ci"), position = position_dodge()) +
        	scale_fill_okabeito(palette = "black_original", order=c(6, 4, 7)) +
        	theme_pubr(legend = "right") +
        	theme(axis.text = element_text(size=9),
			axis.text.x = element_text(size=9, angle = 30, vjust=0.75, hjust=1),
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
				   title="Method"))

	plot_file = sprintf("plots/2-plot_decoupler_filter_benchmark_across_methods/decoupler_bench_filtering_across_methods_%s.svg", cl)
	ggsave(plot = p1, filename = plot_file, width=95, height=65, units="mm", dpi=720, device = "svg")
}
