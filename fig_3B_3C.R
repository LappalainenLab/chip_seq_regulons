library(data.table)
library(ggplotify)
library(dplyr)
library(stringr)
library(here)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(UpSetR)
library(RColorBrewer)

here::i_am("fig_style.R")
source("fig_style.R")

to_plot = fread("data/enrich_analysis/homer_freqs_remap.tsv", header = T, sep="\t")


# Fig. 3C
to_plot %>% 
	drop_na() %>% 
	group_by(TF) %>% 
	summarise(maxtarg = max(`% target`), 
		maxtargmet = paste(method[which(`% target` ==  max(`% target`))])) -> to_plot_stats

position = position_dodge2()


ggplot(data.table(table(to_plot_stats$maxtargmet)), 
	aes(x=factor(V1, levels = c("method2", "method1", "method3"), 
		y=N, 
		fill=factor(V1, levels = c("method1", "method2", "method3"), 
		labels = c("S2Mb", "M2Kb", "S2Kb")))) +
        geom_col(position = position) +
        gtex_v8_figure_theme() +
	xlab("Method") +
        ylab("Number of TFs") +
        scale_fill_brewer(palette="Dark2", name = "Method")
ggsave(filename = "plots/homer_motif_in_tss_enrich_plot_remap_stats.png", device = "png", dpi = 600, width=4, height=4)


# Fig. 3B
pivot_wider(to_plot, names_from="method", values_from="% target") -> temp

temp %>% mutate(S2Mb = complete.cases(method1),
		S2Kb = complete.cases(method2),
		M2Kb = complete.cases(method3)) %>%
		select(TF, S2Mb, M2Kb, S2Kb) -> freq_table

freq_list = list(S2Mb = unique(freq_table[freq_table$S2Mb,]$TF),
	M2Kb = unique(freq_table[freq_table$M2Kb,]$TF),
	S2Kb = unique(freq_table[freq_table$S2Kb,]$TF),
	none = unique(freq_table[(!freq_table$S2Mb) & (!freq_table$M2Kb) & (!freq_table$S2Kb),]$TF))

png("plots/homer_motif_enrich_summary.png", res=600, width = 5, height = 5, units= "in")
upset(fromList(freq_list), order.by="freq", mainbar.y.label = "Experiments")
dev.off()
