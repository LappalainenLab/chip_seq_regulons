#------------------------------------------------------------------------------
# This script visualizes results of benchmarking of different filtering 
# strategies for our approaches
#------------------------------------------------------------------------------

library(ggplot2); theme_set(theme_bw())
library(data.table)
library(foreach)
library(dplyr)
library(ggpubr)
library(here)
library(see)
library(grid)
library(gridExtra)
library(tidyr)
library(cowplot)

#---------------fixing-paths------------------------
here::i_am("README.md")
data_dir = "data/2-plot_decoupler_filter_benchmark_across_methods/"
plot_dir = "plots/s5_s7-plot_decoupler_filtering_benchmark/"
#---------------------------------------------------

cells = c("mcf7", "hepg2", "k562")

foreach(ct = cells) %do% {
	data = fread(paste0(data_dir, ct, "_filtering_benchmark.tsv"))

	data %>% separate(net, c("group", "filter"), "_") %>%
		mutate(filter = replace_na(filter, "no")) %>%
		filter(filter != "exp+atac")-> data
	
	data %>% pivot_wider(names_from = metric, values_from = score, values_fn=mean) %>%
		filter(method == "consensus_estimate")-> to_plot
	
	to_plot$group = factor(to_plot$group, levels=c("S2Mb", "M2Kb", "S2Kb"))
	to_plot$filter =  factor(to_plot$filter, levels=c("random", "no", "motif", "open", "ccre"))
	p1 = ggplot(data = to_plot, aes(x = mcauroc, y = mcauprc, shape = filter)) +
		geom_point() + 
		facet_wrap(~group) +
		theme_pubr() +
		theme(axis.text = element_text(size=9),
			legend.position = c(0.5, 0.1),
			legend.background = element_blank(),
			legend.direction="horizontal",
                	axis.title = element_text(size=9),
			legend.text = element_text(size=9),
                	legend.title = element_text(size=9))
	
	ggsave(plot = p1, filename = paste0(plot_dir, "decoupler_bench_scatter_", ct, "_filtering.svg"), width=210, height=70, units="mm", dpi=720, device = "svg")
	
	data$group = factor(data$group, levels=c("S2Mb", "M2Kb", "S2Kb"))
	
	p2 = ggviolin(data = data %>% filter(metric == "mcauprc") %>% mutate(filter = replace_na(filter, "no"),
										group = factor(group, levels=c("S2Mb", "M2Kb", "S2Kb"))),
					x = "filter",
					y = "score",
					add = "boxplot", 
					add.params = list(fill = "white"))+
		geom_pwc(method = "wilcox_test", label = "p.adj.format", ref.group="random",
				p.adjust.method = "fdr", vjust = 0.5, label.size = 2, bracket.nudge.y = 0.1, tip.length = 0, hide.ns=T) +
		geom_pwc(method = "wilcox_test", label = "p.adj.format", ref.group="no",
        	                p.adjust.method = "fdr", vjust = 0.5, label.size = 2, bracket.nudge.y = 0.3, tip.length = 0, hide.ns=T) +
		facet_wrap(~group) +
		ylab("mcauprc") +
		theme(axis.text.x = element_text(size=9, angle = 30, vjust=0.75),
			axis.text.y = element_text(size=9),
			axis.title.y = element_text(size=9),
			axis.title.x = element_blank())
	ggsave(plot = p2, filename = paste0(plot_dir, "decoupler_bench_box_mcauprc_", ct, "_filtering.svg"), dpi=720, width=210, height=70, units="mm",device = "svg")
	
	
	
	p3 = ggviolin(data = data %>% filter(metric == "mcauroc") %>% mutate(filter = replace_na(filter, "no"),
        	                                                                group = factor(group, levels=c("S2Mb", "M2Kb", "S2Kb"))),
                	                x = "filter",
                        	        y = "score",
                                	add = "boxplot",
                                	add.params = list(fill = "white"))+
        	geom_pwc(method = "wilcox_test", label = "p.adj.format", ref.group="random",
                	        p.adjust.method = "fdr", vjust = 0.5, label.size = 2, bracket.nudge.y = 0.1, tip.length = 0, hide.ns=T) +
        	geom_pwc(method = "wilcox_test", label = "p.adj.format", ref.group="no",
                	        p.adjust.method = "fdr", vjust = 0.5, label.size = 2, bracket.nudge.y = 0.3, tip.length = 0, hide.ns=T) +
        	facet_wrap(~group) +
		ylab("mcauroc") +
        	theme(axis.text.x = element_text(size=9, angle = 30, vjust=0.75),
                	axis.text.y = element_text(size=9),
                	axis.title.y = element_text(size=9),
                	axis.title.x = element_blank())
	ggsave(plot = p3, filename = paste0(plot_dir, "decoupler_bench_box_mcauroc_", ct, "_filtering.svg"), dpi=720, width=210, height=70, units="mm",device = "svg")
} 



