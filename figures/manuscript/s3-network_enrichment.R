library(ggpubr)
library(cowplot)
library(data.table)
library(grid)
library(gridExtra)
library(here)
library(gtable)
library(stringr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(see)
library(ggbreak)
library(foreach)

setwd(here())

# Plot style
source("figures/fig_style.R")

#cell_line = c("HepG2", "MCF7", "GM12878", "K562")
cell_line = c("K562")
foreach(cl = cell_line) %do% {
	dodge <- position_dodge(width=0.9)

	input_file <- sprintf("data/s3-network_enrichment/enrich_scores_remap_all_tfs_%s.tsv", cl)

	enrich_scores = fread(input_file)

	enrich_scores$method = factor(enrich_scores$method, levels=c("method1",
                                                                "method5",
                                                                "method4",
                                                                "method2",
                                                                "method3",
                                                                "shuffled"),
                                                labels = c("S2Mb", "M100Kb", "S100Kb", "M2Kb", "S2Kb", "shuffled"))

	ppi_plot_file <- sprintf("plots/s3-network_enrichment/enrich_scores_remap_all_tfs_ppi_%s.svg", cl)
	coexp_plot_file <- sprintf("plots/s3-network_enrichment/enrich_scores_remap_all_tfs_coexpressed_%s.svg", cl)

	# ppi
	enrich_scores %>% filter(variable == "is_ppi") -> ppi
	p = ggboxplot(ppi, 
			x = "method", 
			y = "odds", 
			fill = "method", 
			outlier.shape=NA, 
			add.params = list(fill = "white"))+
		geom_pwc(method = "wilcox_test", label = "p.adj.signif",
                	        p.adjust.method = "fdr", #vjust = 1.0, 
				label.size = 2, bracket.nudge.y = -0.1, step.increase = 0.01, tip.length = 0, hide.ns=T) +
		ylab("Enrichment Score") +
		xlab("Method") +
		labs(fill='Method') +
		ylim(-1, 10) +
		scale_fill_oi(palette = "black_original", order=c(6, 5, 2, 4, 7)) +
		theme(axis.text.x = element_text(size=9),
			legend.position = "none",
                	axis.text.y = element_text(size=9),
	                axis.title.y = element_text(size=9),
	                axis.title.x = element_blank())

	ggsave(plot = p, filename = ppi_plot_file, device = "svg", dpi = 720, width=99, height=99, units="mm")

	# coexpressed

        input_file <- sprintf("data/s3-network_enrichment/enrich_scores_remap_all_tfs_coexp_%s.tsv", cl)

	coexpressed = fread(input_file)
	coexpressed$method = factor(coexpressed$method, levels=c("method1",
                                                                "method5",
                                                                "method4",
                                                                "method2",
                                                                "method3",
                                                                "shuffled"),
                                                labels = c("S2Mb", "M100Kb", "S100Kb", "M2Kb", "S2Kb", "shuffled"))

	#enrich_scores %>% filter(variable == "is_coexpressed") -> coexpressed
	p = ggboxplot(coexpressed,
        	        x = "method",
        	        y = "odds",
                	fill = "method",
	                outlier.shape = NA,
	                add.params = list(fill = "white"))+
		geom_pwc(method = "wilcox_test", label = "p.adj.signif",
	                        p.adjust.method = "fdr", vjust = 1.0,
	                        label.size = 2, bracket.nudge.y = -0.1, step.increase = 0.01, tip.length = 0, hide.ns=T) +
        	ylab("Enrichment Score") +
        	xlab("Method") +
		#ylim(c(-1, 10)) +
        	labs(fill='Method') +
        	scale_fill_oi(palette = "black_original", order=c(6, 5, 2, 4, 7)) + 
        	theme(axis.text.x = element_text(size=9), 
                	legend.position = "none",
               		axis.text.y = element_text(size=9),
	                axis.title.y = element_text(size=9),
	                axis.title.x = element_blank()) 

	ggsave(plot = p, filename = coexp_plot_file, device = "svg",  dpi = 720, width=99, height=99, units="mm")

}

p = ggviolin(ppi,   
       	        x = "method",
               	y = "odds",
                fill = "method",
       	        add = "boxplot",
               	add.params = list(fill = "white"))+
        geom_pwc(method = "wilcox_test", label = "p.adj.signif",
       	                p.adjust.method = "fdr", vjust = 1.0, label.size = 2, bracket.nudge.y = 0.1, step.increase = 0.05, tip.length = 0, hide.ns=T) +
        ylab("Enrichment Score") +
       	xlab("Method") +
        labs(fill='Method') +
       	scale_fill_oi(palette = "black_original", order=c(6, 5, 2, 4, 7)) +
        theme(axis.text.x = element_text(size=9),
       	        legend.position = "right",
               	axis.text.y = element_text(size=9),
                axis.title.y = element_text(size=9),
       	        axis.title.x = element_blank())



leg = cowplot::get_legend(p)

svg("plots/s3-network_enrichment/enrich_scores_remap_all_tfs_legend.svg", width=1.95, height=3.90)
grid.newpage()
grid.draw(leg)
dev.off()
