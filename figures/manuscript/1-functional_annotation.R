library(SuperExactTest)
library(data.table)	
library(grid)
library(gtable)
library(foreach)
library(VennDiagram)
library(here)
library(jtools)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(UpSetR)
library(ggimage)

source("figures/fig_style.R")

here::i_am("extract_tss.R")

library(RColorBrewer)

tf_target_mapping = fread("data/remap2022/data/TF_target_mapping_filtered_merged_HepG2_with_motifs_with_ppi_with_dnase_with_atac_with_dist_score.tsv", nThread=10)

tf_target_mapping$is_motif_method1 = as.logical(tf_target_mapping$n_motifs_method1)
tf_target_mapping$is_motif_method2 = as.logical(tf_target_mapping$n_motifs_method2)
tf_target_mapping$is_motif_method3 = as.logical(tf_target_mapping$n_motifs_method3)

tf_target_mapping %>%
        distinct() %>%
        mutate(pair = paste0(tf, "_", gene_symbol)) -> to_plot

test_list = list("S2Mb" = to_plot%>% 
			  filter(is_method_1) %>% 
			  dplyr::select(pair) %>%
			  distinct() %>%
			  unlist(),
		"M2Kb" = to_plot %>%
                          filter(is_method_2) %>%
                          dplyr::select(pair) %>%
                          distinct() %>%
                          unlist(),
		"S2Kb" = to_plot %>%
                          filter(is_method_3) %>%
                          dplyr::select(pair) %>%
                          distinct() %>%
                          unlist(),
		"DNAse" = to_plot %>%
                          filter(is_dnase) %>%
                          dplyr::select(pair) %>%
                          distinct() %>%
                          unlist(),
		"ATAC" = to_plot %>%
                          filter(is_atac) %>%
                          dplyr::select(pair) %>%
                          distinct() %>%
                          unlist(),
		"S2Mb_motif" = to_plot %>%
                          filter(is_motif_method1) %>%
                          dplyr::select(pair) %>%
                          distinct() %>%
                          unlist(),
		"M2Kb_motif" = to_plot %>%
                          filter(is_motif_method2) %>%
                          dplyr::select(pair) %>%
                          distinct() %>%
                          unlist(),
		"S2Kb_motif" = to_plot %>%
                          filter(is_motif_method3) %>%
                          dplyr::select(pair) %>%
                          distinct() %>%
                          unlist()
)

svg(paste0("plots/overlap_annotations_HepG2.svg"))

p = upset(fromList(test_list), sets = c("DNAse", "ATAC", "S2Mb", "M2Kb", 
		"S2Kb", "S2Mb_motif", "M2Kb_motif", "S2Kb_motif"), 
	mb.ratio = c(0.7, 0.3),
        order.by="freq",
        mainbar.y.label = "Pairs intersection",
        sets.x.label = "TF-target pairs",
        show.numbers = F,
        text.scale = 1.5,
        scale.sets = "log10")
dev.off()

n = nrow(p$New_data)
n
# 5035767 unique interaction

p$New_data %>% filter((S2Mb == 1) & (S2Kb == 1) & (M2Kb == 1) & (ATAC == 1) & (DNAse ==1)) %>% nrow / n
# 0.3164567 percent of interactions confirmed by 3 approaches, DNAse and ATAC-Seq

p$New_data %>% filter((S2Mb == 1) & (S2Kb == 0) & (M2Kb == 0) & (ATAC == 0) & (DNAse == 0)) %>% nrow / n
# 0.2585541 percent of interactions confirmed only by S2Mb regulon

p$New_data %>% filter((S2Mb == 1) & (S2Kb == 0) & (M2Kb == 0) & (ATAC == 1) & (DNAse == 1)) %>% nrow / n
# 0.2100381 percent of interactions confirmed by S2Mb, DNAse and ATAC-Seq


p$New_data %>% filter(((S2Mb == 1) & (S2Mb_motif ==1)) | ((M2Kb == 1) & (M2Kb_motif ==1)) | ((S2Kb == 1) & (S2Kb_motif ==1))) %>% nrow  / n
# 0.1077216 percent of interactions enriched in TFBS
