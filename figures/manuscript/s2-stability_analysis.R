library(data.table)	
library(grid)
library(gtable)
library(foreach)
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

library(RColorBrewer)

v100 = fread("data/regulons/TF_target_mapping_filtered_merged_K562_enc_v100_with_ppi_with_dnase_with_atac.tsv", nThread=10)
v109 = fread("data/regulons/TF_target_mapping_filtered_merged_K562_with_ppi_with_dnase_with_atac.tsv", nThread=10)
v110 = fread("data/regulons/TF_target_mapping_filtered_merged_K562_enc_v110_with_ppi_with_dnase_with_atac.tsv", nThread=10)


v100 <- v100 %>% dplyr::mutate(pair = paste(tf, gene_symbol, sep = "_")) %>% filter(is_S2Mb | is_S2Kb | is_M2Kb)
v109 <- v109 %>% dplyr::mutate(pair = paste(tf, gene_symbol, sep = "_")) %>% filter(is_S2Mb | is_S2Kb | is_M2Kb)
v110 <- v110 %>% dplyr::mutate(pair = paste(tf, gene_symbol, sep = "_")) %>% filter(is_S2Mb | is_S2Kb | is_M2Kb)

	test_list = list("S2Mb_v100" = v100%>% 
            				  filter(is_S2Mb) %>% 
            				  dplyr::select(pair) %>%
            				  distinct() %>%
            				  unlist(),
			"S2Mb_v109" = v109%>%
                                          filter(is_S2Mb) %>%
                                          dplyr::select(pair) %>%
                                          distinct() %>%
                                          unlist(),
			"S2Mb_v110" = v110%>%
                                          filter(is_S2Mb) %>%
                                          dplyr::select(pair) %>%
                                          distinct() %>%
                                          unlist(),
			"M2Kb_v100" = v100%>%
                                          filter(is_M2Kb) %>%
                                          dplyr::select(pair) %>%
                                          distinct() %>%
                                          unlist(),
                        "M2Kb_v109" = v109%>%
                                          filter(is_M2Kb) %>%
                                          dplyr::select(pair) %>%
                                          distinct() %>%
                                          unlist(),
                        "M2Kb_v110" = v110%>%
                                          filter(is_M2Kb) %>%
                                          dplyr::select(pair) %>%
                                          distinct() %>%
                                          unlist(),
			"S2Kb_v100" = v100%>%
                                          filter(is_S2Kb) %>%
                                          dplyr::select(pair) %>%
                                          distinct() %>%
                                          unlist(),
                        "S2Kb_v109" = v109%>%
                                          filter(is_S2Kb) %>%
                                          dplyr::select(pair) %>%
                                          distinct() %>%
                                          unlist(),
                        "S2Kb_v110" = v110%>%
                                          filter(is_S2Kb) %>%
                                          dplyr::select(pair) %>%
                                          distinct() %>%
                                          unlist()
		)

svg("./plots/1-functional_annotation/s2-stability_analysis.svg")

p = upset(fromList(test_list), sets = c("S2Mb_v100", "S2Mb_v109", "S2Mb_v110",
					"M2Kb_v100", "M2Kb_v109", "M2Kb_v110",
					"S2Kb_v100", "S2Kb_v109", "S2Kb_v110"), 
	mb.ratio = c(0.7, 0.3),
       	order.by="freq",
        mainbar.y.label = "Number of TF-target gene pairs",
	sets.x.label = "TF-target\ninteractiosn",
        show.numbers = F,
       	text.scale = 2,
        scale.sets = "identity",
	nintersects = 15,
	intersection)
print(p)
dev.off()

