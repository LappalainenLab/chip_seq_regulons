library(data.table)
#library(SuperExactTest)
library(here)
library(jtools)
library(stringr)
library(dplyr)
library(doParallel)
library(foreach)
library(tidyr)
library(UpSetR)
library(ggbreak)
library(decoupleR)

registerDoParallel(cores=10)

here::i_am("README.md")


# Plot style
source("figures/fig_style.R")

cell_line = c("K562", "MCF7", "HepG2", "GM12878")

foreach(cl = cell_line) %do% {

	regulon_input_file = sprintf("data/regulons/%s_all_regulons.tsv", cl)
	chip_atlas_input_file = sprintf("data/regulons/ChIP-Atlas_target_genes_%s.tsv", cl)

	ga_normal = fread("data/regulons/Revised_Supplemental_Table_S3_Normal.csv", header = T)
	collectri = decoupleR::get_collectri(organism='human', split_complexes=FALSE)
	tf_target_mapping = fread(regulon_input_file, nThread=10)
	chip_atlas = fread(chip_atlas_input_file, header = T)

	tfs = intersect(intersect(intersect(tf_target_mapping$tf, ga_normal$TF), collectri$source), chip_atlas$tf)

	ga_normal %>% filter(TF %in% tfs) -> ga_normal
	collectri %>% filter(source %in% tfs) -> collectri
	chip_atlas %>% filter(tf %in% tfs) -> chip_atlas
	tf_target_mapping %>% filter(tf %in% tfs) -> tf_target_mapping

	tf_target_mapping %>%
		dplyr::distinct() %>%
		mutate(pair = paste0(tf, "_", gene_symbol)) -> to_plot

	test_list = list("S2Mb" = to_plot %>% 
					filter(is_S2Mb) %>% 
	                          	dplyr::select(pair) %>%
	                          	dplyr::distinct() %>%
					unlist(),
			"M2Kb" = to_plot %>%
	                                filter(is_M2Kb) %>%
	                                dplyr::select(pair) %>%
	                                dplyr::distinct() %>%
     		                        unlist(),
			"S2Kb" = to_plot %>%
                        	        filter(is_S2Kb) %>%
                                	dplyr::select(pair) %>%
                                	dplyr::distinct() %>%
                                	unlist(),
			"ChIP-Atlas" = chip_atlas %>%
        	                                dplyr::distinct() %>%
                	                        mutate(pair = paste0(tf, "_", Target_genes)) %>%
                        	                dplyr::select(pair) %>%
						dplyr::distinct() %>%	
                                        	unlist(),
			"Dorothea" = ga_normal %>%
						dplyr::select(TF, target) %>%
                	                        dplyr::distinct() %>%
                        	                mutate(pair = paste0(TF, "_", target)) %>%
                                	        dplyr::select(pair) %>%
						dplyr::distinct() %>% 
	                                        unlist(), 
			"CollecTri" =  collectri %>%
						dplyr::distinct() %>%
                        	                mutate(pair = paste0(source, "_", target)) %>%
                                	        dplyr::select(pair) %>% 
						dplyr::distinct() %>%
	                                        unlist())

	plot_file = sprintf("plots/2-regulon_overlap/overlap_with_db_%s.svg", cl)
	svg(plot_file, width = 8.268)

	p = upset(fromList(test_list), 
		sets = c("S2Mb", "M2Kb", "S2Kb", "ChIP-Atlas", "Dorothea", "CollecTri"), 
		mb.ratio = c(0.7, 0.3), 
		order.by="freq", 
		mainbar.y.label = "Pairs intersection", 
		sets.x.label = "TF-target pairs",
		show.numbers = T,
		text.scale = 1.5,
		scale.sets = "identity", 
		nintersects=12) 

	print(p)
	dev.off()
}

collectri %>% distinct()  %>% nrow
# 16767 average collectri regulon size

length(unique(collectri$source))
# 203 TFs 

#ChIP-Atlas

chip_atlas = fread("data/trans_networks/ChIP-Atlas_target_genes_K562.tsv", header = T)
chip_atlas %>% distinct()  %>% nrow
# 1529999 average ChIP-Atlas K-562 regulon size

length(unique(chip_atlas$tf))
# 400 TFs

chip_atlas = fread("data/trans_networks/ChIP-Atlas_target_genes_MCF7.tsv", header = T)
chip_atlas %>% distinct()  %>% nrow
# 2068087 average ChIP-Atlas MCF7 regulon size

length(unique(chip_atlas$tf))
# 366 TFs

chip_atlas = fread("data/trans_networks/ChIP-Atlas_target_genes_HepG2.tsv", header = T)
chip_atlas %>% distinct()  %>% nrow 
# 639027 average ChIP-Atlas HepG2 regulon size
 
length(unique(chip_atlas$tf))
# 182 TFs

chip_atlas = fread("data/trans_networks/ChIP-Atlas_target_genes_GM12878.tsv", header = T)
chip_atlas %>% distinct()  %>% nrow  
# 655040 average ChIP-Atlas GM-12878 regulon size
 
length(unique(chip_atlas$tf)) 
# 156 TFs


ga_normal %>% select(TF, target) %>% distinct()  %>% nrow
# 349005 average Dorothea regulon size

length(unique(ga_normal$TF))
# 203 TFs



p$New_data %>% dplyr::select(S2Mb, S2Kb, M2Kb) %>% filter((S2Mb ==1) & (S2Kb == 1) & (M2Kb == 1)) %>% nrow / nrow(p$New_data)
# 0.24414 of interactions confirmd by all 3 approaches in K-562
# 0.291794 of interactions confirmd by all 3 approaches in Hep-G2
# 0.3064423 of interactions confirmd by all 3 approaches in MCF-7
# 0.2678522 of interactions confirmd by all 3 approaches in GM-12878


p$New_data %>% filter((`ChIP-Atlas` == 1) & (S2Mb ==1) & (S2Kb == 1) & (M2Kb == 1)) %>% nrow / nrow(p$New_data)
# 0.1479902 of interactions confirmd by all 3 approaches and ChIP-Atlas in K-562
# 0.1968782 of interactions confirmd by all 3 approaches and ChIP-Atlas in Hep-G2
# 0.1763473 of interactions confirmd by all 3 approaches and ChIP-Atlas in MCF-7
# 0.1670479 of interactions confirmd by all 3 approaches and ChIP-Atlas in GM-12878

p$New_data %>% filter((Dorothea == 1) & ((S2Mb ==1) | (S2Kb == 1) | (M2Kb == 1) | (`ChIP-Atlas` == 1) | (CollecTri == 1))) %>% nrow / nrow(p$New_data)
# 0.0847815 of interactions confirmd by Dorothea and any other database in K-562
# 0.1024777 of interactions confirmd by Dorothea and any other database in Hep-G2
# 0.1312638 of interactions confirmd by Dorothea and any other database in MCF-7    
# 0.1131837 of interactions confirmd by Dorothea and any other database in GM-12878 

p$New_data %>% filter((Dorothea == 1) & (CollecTri == 1)) %>% nrow / nrow(p$New_data)
# 0.003277594 of interactions confirmd by Dorothea and CollecTri in K-562
# 0.003820505 of interactions confirmd by Dorothea and CollecTri in HepG2
# 0.00645772 of interactions confirmd by Dorothea and CollecTri in MCF-7
# 0.005587055 of interactions confirmd by Dorothea and CollecTri in K-562

p$New_data %>% rowSums() -> temp
sum(temp > 1) / nrow(p$New_data) 
# 0.4818014 interactions were present in at least two datasets (K-562)
# 0.6289997 interactions were present in at least two datasets (HepG2)
# 0.5600939 interactions were present in at least two datasets (MCF7)
# 0.5310304 interactions were present in at least two datasets (GM-12878)
