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

here::i_am("fig_style.R")


# Plot style
source("fig_style.R")

enrich_scores = fread("data/enrich_analysis/enrich_scores_remap_all_tfs_with_ppi.tsv")
dodge <- position_dodge(width=0.9)

TFs = c("NFE2", "RUNX1", "GFI1B", "IKZF1") 

enrich_scores %>%  filter(((variable == "network") | 
				(variable == "is_ppi") |
				(variable == "is_coexpressed")) & (TF %in% TFs)) -> enrich_scores


enrich_scores$method = factor(unlist(enrich_scores$method), levels = c("method1", "method2", "method3"), labels=c("S2Mb", "M2Kb", "S2Kb"))

enrich_scores$pvalue = p.adjust(enrich_scores$pvalue, method = "fdr")

ggplot(enrich_scores,
               aes(y=as.numeric(unlist(odds)),
                   x=as.factor(unlist(access)),
                   group=interaction(as.factor(unlist(method)),
                                     as.factor(unlist(access))),
                   ymin=as.numeric(unlist(conf_int_l)),
                   ymax=as.numeric(unlist(conf_int_u)),
                   color=method,
		   shape = method,
                   label=interaction(as.factor(unlist(method)), as.factor(unlist(access))))) +
        geom_pointrange(position = dodge, size=0.2) +
        gtex_v8_figure_theme() +
        geom_hline(yintercept = 0.0, color="black", linetype="dashed") +
        coord_flip() +
	#ylim(-2, 6) + 
        ylab("Enrichment Score") +
        xlab("Accession") +
        labs(color='Method') +
        guides(alpha = "none",
		shape = "none") +
        scale_colour_brewer(palette="Dark2", name = "Method") +
	facet_grid( unlist(TF) ~ variable, scales = "free") + 
        theme(panel.grid = element_line(color = "gray",
                                              size = 0.1,
                                              linetype = 2),
		axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
		legend.text = element_text(size = 8),
		axis.text.y = element_blank())


ggsave(filename = "plots/enrich_scores_remap_4_tfs_with_ppi.png", device = "png", dpi = 600, width=4, height=6)
