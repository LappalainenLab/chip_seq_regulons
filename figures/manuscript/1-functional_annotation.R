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
cell_line = c("K562", "HepG2", "MCF7", "GM12878")
foreach(cl = cell_line) %do% {
    #input_file <- sprintf("data/regulons/%s_all_regulons.tsv", cl)
	input_file <- sprintf("data/regulons/TF_target_mapping_filtered_merged_%s_with_ppi_with_dnase_with_atac_with_motifs_with_ccres.tsv", cl)

	to_plot = fread(input_file)
	
	to_plot <- to_plot %>% dplyr::mutate(pair = paste(tf, gene_symbol, sep = "_")) %>% filter(is_S2Mb | is_S2Kb | is_M2Kb)

	test_list = list("S2Mb" = to_plot%>% 
                    				  filter(is_S2Mb) %>% 
                    				  dplyr::select(pair) %>%
                    				  distinct() %>%
                    				  unlist(),
        			"M2Kb" = to_plot %>%
                        	          filter(is_M2Kb) %>%
                                	  dplyr::select(pair) %>%
        	                          distinct() %>%
                	                  unlist(),
        			"S2Kb" = to_plot %>%
                                	  filter(is_S2Kb) %>%
        	                          dplyr::select(pair) %>%
                	                  distinct() %>%
                        	          unlist(),
        			"is_motif" = to_plot %>%
                	                  filter(!is.na(n_motifs) | ! is.na(n_motifs_hoco)) %>%
                        	          dplyr::select(pair) %>%
                                	  distinct() %>%
        	                          unlist(),
	                        "is_PLS" = to_plot %>%
                                          filter(is_pls) %>%
                                          dplyr::select(pair) %>%
                                          distinct() %>%
                                          unlist(),
                                "is_pELS" = to_plot %>%
                                          filter(is_pels) %>%
                                          dplyr::select(pair) %>%
                                          distinct() %>%
                                          unlist(),
                                "is_dELS" = to_plot %>%
                                          filter(is_dels) %>%
                                          dplyr::select(pair) %>%
                                          distinct() %>%
                                          unlist()
			)

        plot_file <- sprintf("plots/1-functional_annotation/overlap_annotations_%s.svg", cl)
	svg(plot_file)

	p = upset(fromList(test_list), sets = c("S2Mb", "M2Kb", 
			"S2Kb", "is_motif", "is_PLS", "is_pELS", "is_dELS"), 
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
}



n = nrow(p$New_data)

n
# 5035767 unique interaction

p$New_data %>% filter((S2Mb == 1) & (S2Kb == 1) & (M2Kb == 1)) %>% nrow
# 0.3164567 percent of interactions confirmed by 3 approaches, DNAse and ATAC-Seq


