library(ggpubr)
library(cowplot)
library(data.table)
library(grid)
library(gridExtra)
library(here)
library(gtable)
library(jtools)
library(stringr)
library(biomaRt)
library(ggplot2)
library(fastglm)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(see)


#here::i_am("fig_style.R")

# Plot style
source("figures/fig_style.R")
enrich_scores = fread("data/enrich_analysis/enrich_scores_remap_all_tfs_K562_04_07_23.tsv")
enrich_scores_random = fread("data/enrich_analysis/enrich_scores_remap_all_tfs_K562_04_07_23_random.tsv")
dodge <- position_dodge(width=0.9)

bind_rows(enrich_scores, enrich_scores_random %>% filter(method == "method1_random")) -> enrich_scores

enrich_scores$method = factor(enrich_scores$method, levels=c("method1",
								"method2", 
								"method3",
								"method1_random"),
						labels = c("S2Mb", "M2Kb", "S2Kb", "Random"))

# ppi
enrich_scores %>% filter(variable == "is_ppi") -> ppi
p = ggviolin(ppi, 
		x = "method", 
		y = "odds", 
		fill = "method", 
		add = "boxplot", 
		add.params = list(fill = "white"))+
	geom_pwc(method = "wilcox_test", label = "p.adj.signif",
                        p.adjust.method = "fdr", vjust = 0.1, label.size = 2, bracket.nudge.y = 0.1, step.increase = 0.05, tip.length = 0, hide.ns=T) +
	ylab("Enrichment Score") +
	xlab("Method") +
	labs(fill='Method') +
	scale_fill_oi(palette = "black_original", order=c(6, 4, 7)) +
	theme(axis.text.x = element_text(size=9),
		legend.position = "none",
                axis.text.y = element_text(size=9),
                axis.title.y = element_text(size=9),
                axis.title.x = element_blank())

ggsave(plot = p, filename = "plots/enrich_scores_remap_all_tfs_ppi_K562_04_07_23.svg", device = "svg", dpi = 720, width=99, height=99, units="mm")

# coexpressed
enrich_scores %>% filter(variable == "is_coexpressed") -> coexpressed
p = ggviolin(coexpressed,
                x = "method",
                y = "odds",
                fill = "method",
                add = "boxplot",
                add.params = list(fill = "white"))+
        geom_pwc(method = "wilcox_test", label = "p.adj.signif",
			p.adjust.method = "fdr", vjust = 0.1, label.size = 2, bracket.nudge.y = 0.1, step.increase = 0.05, tip.length = 0, hide.ns=T) +
        ylab("Enrichment Score") +
        xlab("Method") +
        labs(fill='Method') +
        scale_fill_oi(palette = "black_original", order=c(6, 4, 7)) + 
        theme(axis.text.x = element_text(size=9), 
                legend.position = "none",
                axis.text.y = element_text(size=9),
                axis.title.y = element_text(size=9),
                axis.title.x = element_blank()) 

ggsave(plot = p, filename = "plots/enrich_scores_remap_all_tfs_coexpressed_K562_04_07_23.svg", device = "svg",  dpi = 720, width=99, height=99, units="mm")



p = ggviolin(ppi,   
                x = "method",
                y = "odds",
                fill = "method",
                add = "boxplot",
                add.params = list(fill = "white"))+
        geom_pwc(method = "wilcox_test", label = "p.adj.signif",
                        p.adjust.method = "fdr", vjust = 0.1, label.size = 2, bracket.nudge.y = 0.1, step.increase = 0.05, tip.length = 0, hide.ns=T) +
        ylab("Enrichment Score") +
        xlab("Method") +
        labs(fill='Method') +
        scale_fill_oi(palette = "black_original", order=c(6, 4, 7)) +
        theme(axis.text.x = element_text(size=9),
                legend.position = "right",
                axis.text.y = element_text(size=9),
                axis.title.y = element_text(size=9),
                axis.title.x = element_blank())



leg = cowplot::get_legend(p)

svg("plots/enrich_scores_remap_all_tfs_legend_04_07_23.svg", width=1.95, height=3.90)
grid.newpage()
grid.draw(leg)
dev.off()
