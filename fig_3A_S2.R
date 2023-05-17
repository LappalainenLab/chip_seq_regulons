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
library(biomaRt)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(UpSetR)
library(ggimage)

source("fig_style.R")

here::i_am("extract_tss.R")

library(RColorBrewer)
myCol <- brewer.pal(3, "Dark2")




files = list.files("data/enrich_analysis/", pattern="overlap")

for (f in files){			
	sub_tf_target_mapping_acc = fread(paste0("data/enrich_analysis/", f))
	TF = sub_tf_target_mapping_acc$tf
	x = sub_tf_target_mapping_acc$access
	print(length(sub_tf_target_mapping_acc %>% filter(is_method_1) %>% dplyr::select(ensembl_gene_id) %>% distinct() %>% unlist()))
	png(paste0("plots/overlap_methods_", TF, "_", x, ".png"), res=600, width = 5, height = 5, units= "in")
	print(upset(fromList(list("S2Mb" = sub_tf_target_mapping_acc %>% 
				  filter(is_method_1) %>% 
				  dplyr::select(ensembl_gene_id) %>% 
				  distinct() %>% 
				  unlist(),
			"M2Kb" = sub_tf_target_mapping_acc %>% 
				     filter(is_method_2) %>% 
				     dplyr::select(ensembl_gene_id) %>% 
				     distinct() %>% 
				     unlist(),
			"S2Kb" = sub_tf_target_mapping_acc %>% 
				     filter(is_method_3) %>% 
				     dplyr::select(ensembl_gene_id) %>% 
				     distinct() %>%
				     unlist())), order.by="freq",  mainbar.y.label = "Number of target genes"))
	dev.off()
}	




