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

library(RColorBrewer)
cell_line = c("K562", "HepG2", "MCF7", "GM12878")
foreach(cl = cell_line) %do% {
        input_file <- sprintf("data/1-functional_annotation/%s_regulon_with_motifs.tsv", cl)

	to_plot = fread(input_file)

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

        plot_file <- sprintf("plots/1-functional_annotation/overlap_annotations_%s.svg", cl)
	svg(plot_file)

	p = upset(fromList(test_list), sets = c(#"DNAse", "ATAC", 
						"S2Mb", "M2Kb", 
			"S2Kb", "S2Mb_motif", "M2Kb_motif", "S2Kb_motif"), 
		mb.ratio = c(0.7, 0.3),
        	order.by="freq",
	        mainbar.y.label = "Number of TF-target gene pairs",
		sets.x.label = "TF-target\ninteractiosn",
	        show.numbers = F,
        	text.scale = 2,
	        scale.sets = "log10",
		nintersects = 12,
		intersection)
	print(p)
	dev.off()
}



n = nrow(p$New_data)

n
# 5035767 unique interaction

p$New_data %>% filter((S2Mb == 1) & (S2Kb == 1) & (M2Kb == 1)) %>% nrow
# 0.3164567 percent of interactions confirmed by 3 approaches, DNAse and ATAC-Seq

#p$New_data %>% filter((S2Mb == 1) & (S2Kb == 1) & (M2Kb == 1) & (ATAC == 1) & (DNAse ==1)) %>% nrow    
# 0.3164567 percent of interactions confirmed by 3 approaches, DNAse and ATAC-Seq

#p$New_data %>% filter((S2Mb == 1) & (S2Kb == 1) & (M2Kb == 1) & (ATAC == 1) & (DNAse ==1)) %>% nrow / n
# 0.3164567 percent of interactions confirmed by 3 approaches, DNAse and ATAC-Seq

#p$New_data %>% filter((S2Mb == 1) & (S2Kb == 0) & (M2Kb == 0) & (ATAC == 0) & (DNAse == 0)) %>% nrow / n
# 0.2585541 percent of interactions confirmed only by S2Mb regulon

#p$New_data %>% filter((S2Mb == 1) & (S2Kb == 0) & (M2Kb == 0) & (ATAC == 1) & (DNAse == 1)) %>% nrow / n
# 0.2100381 percent of interactions confirmed by S2Mb, DNAse and ATAC-Seq


p$New_data %>% filter(((S2Mb == 1) & (S2Mb_motif ==1)) | ((M2Kb == 1) & (M2Kb_motif ==1)) | ((S2Kb == 1) & (S2Kb_motif ==1))) %>% nrow  / n
# 0.1077216 percent of interactions enriched in TFBS

