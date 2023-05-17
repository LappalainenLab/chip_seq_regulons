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

here::i_am("fig_style.R")


# Plot style
source("fig_style.R")

enrich_scores = fread("data/enrich_analysis/enrich_scores_remap_all_tfs_with_ppi.tsv")

dodge <- position_dodge(width=0.9)

enrich_scores %>% filter((variable == "is_in_msigdb") & (abs(odds) > 0.1 )) %>%
		select("TF") %>%
		distinct() -> TFs

samples = sample(TFs$TF, 20)
enrich_scores %>% filter(TF %in% samples) %>% filter((variable != "has_motif") & 
							(variable != "network") & 
							(variable != "is_ppi") & 
							(variable != "is_coexpressed") &
							(abs(conf_int_l) < 10) & 
							(abs(conf_int_u) < 10)) -> enrich_scores
enrich_scores
enrich_scores$pvalue = p.adjust(enrich_scores$pvalue, method = "fdr")

leg = ggplot(enrich_scores,
               aes(y=as.numeric(unlist(odds)),
                   x=as.factor(unlist(TF)),
                   ymin=as.numeric(unlist(conf_int_l)),
                   ymax=as.numeric(unlist(conf_int_u)),
                   shape=as.factor(unlist(method)))) +
        geom_pointrange(position = dodge, size=0.2) +
	scale_shape_discrete(labels = c("S2Mb", "M2Kb", "S2Kb"), name = "Method") +
        theme_bw() +
	theme(legend.text = element_text(size = 8),
		legend.title = element_text(size = 9))


plot = ggplot(enrich_scores,
               aes(y=as.numeric(unlist(odds)),
                   x=as.factor(unlist(TF)),
                   group=interaction(as.factor(unlist(method)),
                                     as.factor(unlist(TF))),
                   ymin=as.numeric(unlist(conf_int_l)),
                   ymax=as.numeric(unlist(conf_int_u)),
                   color=as.factor(unlist(TF)),
                   shape=as.factor(unlist(method)),
                   alpha=(as.numeric(unlist(pvalue)) < 0.05),
                   label=interaction(as.factor(unlist(method)), as.factor(unlist(TF))))) +
        geom_pointrange(position = dodge, size=0.2) +
        gtex_v8_figure_theme() +
        geom_hline(yintercept = 0.0, color="black", linetype="dashed") +
        coord_flip() +
        ylab("Enrichment Score") +
        xlab("TF") +
        labs(color='TF', shape="method") +
	scale_shape_discrete(labels = c("S2Mb", "M2Kb", "S2Kb")) +
        guides(alpha = "none") +
	facet_grid( unlist(TF) ~ variable, scales = "free") + 
        theme(panel.grid = element_line(color = "gray",
                                              size = 0.1,
                                              linetype = 2),
		strip.text.y = element_blank(),
		axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
		legend.position = "none")

leg = gtable_filter(ggplot_gtable(ggplot_build(p2)), "guide-box")

plotNew <- arrangeGrob(plot, leg,
          widths = unit.c(unit(1, "npc") - leg2$width, leg2$width), nrow = 1)


ggsave(plot = plotNew, filename = "plots/enrich_scores_remap_all_tfs_with_ppi.png", device = "png", dpi = 600, width=4, height=6)
