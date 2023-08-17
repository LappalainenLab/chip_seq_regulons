library(grid)
library(cowplot)
library(ggpubr)
library(see)
library(data.table)
library(here)
library(jtools)
library(stringr)
library(biomaRt)
library(ggplot2)
library(fastglm)
library(ggrepel)
library(dplyr)
library(RColorBrewer)


# Plot style
source("figures/fig_style.R")
enrich_scores = fread("data/enrich_analysis/enrich_scores_remap_all_tfs_K562_04_07_23.tsv")
#enrich_scores_random = fread("data/enrich_analysis/enrich_scores_remap_all_tfs_K562_04_07_23_random.tsv")
dodge <- position_dodge(width=0.9)

TFs = c("NFE2", "RUNX1", "GFI1B", "IKZF1") 
#bind_rows(enrich_scores, enrich_scores_random %>% filter(method == "method1_random")) -> enrich_scores

enrich_scores %>%  filter(((variable == "network") | 
				(variable == "is_ppi") |
				(variable == "is_coexpressed")) & 
				(TF %in% TFs)) -> enrich_scores


enrich_scores$method = factor(enrich_scores$method, levels=c("method1",
                                                                "method2",
                                                                "method3",
                                                                "method1_random"),
                                                labels = c("S2Mb", "M2Kb", "S2Kb", "Random"))

enrich_scores$pvalue = p.adjust(enrich_scores$pvalue, method = "fdr")

p = ggplot(enrich_scores,
               aes(y=as.numeric(unlist(odds)),
                   x=as.factor(unlist(access)),
                   group=interaction(as.factor(unlist(method)),
                                     as.factor(unlist(access))),
                   ymin=as.numeric(unlist(conf_int_l)),
                   ymax=as.numeric(unlist(conf_int_u)),
                   color=method,
		   alpha = pvalue < 0.01,
                   label=interaction(as.factor(unlist(method)), as.factor(unlist(access))))) +
        geom_pointrange(position = dodge, size=0.2) +
        gtex_v8_figure_theme() +
        geom_hline(yintercept = 0.0, color="black", linetype="dashed") +
        coord_flip() +
	#ylim(-2, 6) + 
        ylab("Enrichment Score") +
        xlab("Accession") +
        labs(color='Method') +
        guides(alpha = "none") + #,
		#shape = "none") +
	scale_color_oi(palette = "black_original", order=c(6, 4, 7, 1)) +
	theme_pubr() +
	facet_grid( unlist(TF) ~ variable, scales = "free", switch="y") +
        theme(panel.grid = element_line(color = "gray",
                                              size = 0.1,
                                              linetype = 2),
		legend.key = element_rect(colour = "transparent", fill = "white"),
		strip.text.y.left = element_text(size=9, angle=0),
                axis.title.y = element_blank(),
                strip.placement = "outside",
                strip.background.y = element_blank(),
		legend.position = "right", 
		strip.text.x = element_text(size=9),
		axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size=9),
		legend.text = element_text(size = 9),
		axis.text.y = element_blank())




leg = cowplot::get_legend(p)
p = p + theme(legend.position = 'none')
ggsave(plot = p, filename = "plots/enrich_scores_remap_all_tfs_K562_04_07_23.svg", device = "svg", dpi = 720, width=105, height=105, units = "mm")

 
svg("plots/enrich_scores_remap_all_tfs_K562_legend_04_07_23.svg", width=1.95, height=3.90)
grid.newpage()
grid.draw(leg)
dev.off()
