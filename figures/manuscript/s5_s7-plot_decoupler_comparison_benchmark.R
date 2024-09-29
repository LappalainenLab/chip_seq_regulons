#------------------------------------------------------------------------------
# This script visualizes results of benchmarking of
# S2Mb, M2Kb, CollecTri, Dorothea, GTRD, ChIP-Atlas, RegNet, TRRUST regulons
#------------------------------------------------------------------------------

library(ggplot2); theme_set(theme_bw())
library(data.table)
library(dplyr)
library(ggpubr)
library(here)
library(see)
library(grid)
library(gridExtra)
library(tidyr)
library(cowplot)
library(foreach)

#---------------fixing-paths------------------------
setwd(here())

processed_data_dir = "data/2-plot_decoupler_comparison_benchmark_across_cells/"
plot_dir = "plots/s5_s7-plot_decoupler_comparison_benchmark/"
#---------------------------------------------------

cells = c("mcf7", "hepg2", "k562")

foreach(ct = cells) %do% {

	data = fread(paste0(processed_data_dir, ct, "_comparison_benchmark.tsv"))

	data %>% pivot_wider(names_from = metric, values_from = score, values_fn=mean) %>%
		filter(method == "consensus_estimate")-> to_plot

	to_plot$net = factor(to_plot$net, levels=c("S2Mb", "M100Kb", "S100Kb", "M2Kb", "CollecTri", "Dorothea", "ChIP-Atlas", "RegNet", "TRRUST"))
	p1 = ggplot(data = to_plot, aes(x = mcauroc, y = mcauprc, color = net)) +
		geom_point() +
		scale_color_okabeito(palette = "black_original", order=c(6, 2, 5, 4, 8, 1, 7, 3, 9)) +
		theme_pubr(legend = c(0.85, 0.25)) +
		theme(axis.text = element_text(size=9),
			legend.background = element_blank(),
                	axis.title = element_text(size=9),
			legend.text = element_text(size=9),
                	legend.title = element_text(size=9),
			legend.spacing.y = unit(0.1, 'pt')) +
		guides(colour=guide_legend(ncol=1,byrow=T,
			keywidth=3,
                        keyheight=3,
                        default.unit="mm"))
	ggsave(plot = p1, filename = paste0(plot_dir, "decoupler_bench_scatter_", ct, "_comparison.svg"), width=100, height=70, units="mm", dpi=720, device = "svg")
	
	data$net = factor(data$net, levels=c("S2Mb", "M100Kb", "S100Kb", "M2Kb", "CollecTri", "Dorothea", "GTRD", "ChIP-Atlas", "RegNet", "TRRUST"))
	
	p2 = ggviolin(data = data %>% filter(metric == "mcauprc"), 
			x = "net",
			y = "score",
			add = "boxplot", 
			add.params = list(fill = "white"),
			fill = "net") +
		scale_fill_oi(palette = "black_original", order=c(6, 2, 5, 4, 8, 1, 7, 3, 9)) +
		theme_pubr(legend = "none") +
		ylab("mcauprc") +
		theme(axis.text.x = element_text(size=9, angle = 30, vjust=0.75),
			axis.text.y = element_text(size=9),
			axis.title.y = element_text(size=9),
			axis.title.x = element_blank(),
			legend.key.size = unit(3, "pt"),
			legend.key.height = unit(2, "pt")) +
        	guides(colour=guide_legend(ncol=1,byrow=T,
			keywidth=3,
                        keyheight=2,
                        default.unit="mm"))

	ggsave(plot = p2, filename = paste0(plot_dir, "decoupler_bench_box_mcauprc_", ct, "_comparison.svg"), dpi=720, width=100, height=70, units="mm",device = "svg")
	p3 = ggviolin(data = data %>% filter(metric == "mcauroc"),
			x = "net",
                        y = "score",
                        add = "boxplot",
                        add.params = list(fill = "white"),
                        fill = "net") +
        	scale_fill_oi(palette = "black_original", order=c(6, 2, 5, 4, 8, 1, 7, 3, 9)) +
		theme_pubr(legend = "none") +
        	ylab("mcauroc") +
		theme(axis.text.x = element_text(size=9, angle = 30, vjust=0.75),
                	axis.text.y = element_text(size=9),
                	axis.title.y = element_text(size=9),
                	axis.title.x = element_blank(),
                	legend.key.size = unit(3, "pt"),
			legend.key.height = unit(2, "pt")) +
        	guides(colour=guide_legend(ncol=1,
                        keywidth=3,
                        keyheight=2,
                        default.unit="mm"))
	
	ggsave(plot = p3, filename = paste0(plot_dir, "decoupler_bench_box_mcauroc_", ct, "_comparison.svg"), dpi=720, width=100, height=70, units="mm",device = "svg") 

	p3 = ggviolin(data = data %>% filter(metric == "mcauroc"),
			x = "net",
                        y = "score",
                        add = "boxplot",
                        add.params = list(fill = "white"),
                        fill = "net") +
        	scale_fill_oi(palette = "black_original", order=c(6, 2, 5, 4, 8, 1, 7, 3, 9)) +
        	theme_pubr(legend = "right") +
        	ylab("mcauroc") +
        	theme(axis.text.x = element_text(size=9, angle = 30, vjust=0.75),
                	axis.text.y = element_text(size=9),
                	axis.title.y = element_text(size=9),
                	axis.title.x = element_blank(),
                	legend.key.size = unit(3, "pt"),
                	legend.key.height = unit(2, "pt")) +
        	guides(colour=guide_legend(ncol=1,
                        keywidth=3,
                        keyheight=2,
                        default.unit="mm"))

	leg = cowplot::get_legend(p3)

	svg(paste0(plot_dir, "decoupler_bench_legend.svg"), width=1.95, height=3.90)
	grid.newpage()
	grid.draw(leg)
	dev.off()
}

