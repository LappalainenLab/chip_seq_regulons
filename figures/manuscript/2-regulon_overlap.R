library(data.table)
library(SuperExactTest)
library(here)
library(jtools)
library(stringr)
library(dplyr)
library(doParallel)
library(foreach)
library(tidyr)
library(UpSetR)

registerDoParallel(cores=10)

here::i_am("test_GA_enrichment.R")


# Plot style
source("figures/fig_style.R")

with_motif = F

tf_target_mapping = fread("data/remap2022/data/TF_target_mapping_filtered_merged_K562_with_motifs_with_ppi_with_coexpr_with_dnase_with_atac_with_dist_score.tsv", nThread=10)

ga_normal = fread("data/garcia_alonso/Revised_Supplemental_Table_S3_Normal.csv", header = T)
msigdb = fread("data/trans_networks/MSigDB_C3_TFT_GTRD_v2023.tsv", header = T)
chip_atlas = fread("data/trans_networks/ChIP-Atlas_target_genes_K562.tsv", header = T)

tfs = intersect(intersect(intersect(tf_target_mapping$tf, ga_normal$TF), msigdb$tf), chip_atlas$tf)

ga_normal %>% filter(TF %in% tfs) -> ga_normal
msigdb %>% filter(tf %in% tfs) -> msigdb
chip_atlas %>% filter(tf %in% tfs) -> chip_atlas
tf_target_mapping %>% filter(tf %in% tfs) -> tf_target_mapping

tf_target_mapping %>%
	dplyr::distinct() %>%
	mutate(pair = paste0(tf, "_", gene_symbol)) -> to_plot

test_list = list("S2Mb" = to_plot %>% 
				filter(is_method_1) %>% 
                          	dplyr::select(pair) %>%
                          	dplyr::distinct() %>%
				unlist(),
		"M2Kb" = to_plot %>%
                                filter(is_method_2) %>%
                                dplyr::select(pair) %>%
                                dplyr::distinct() %>%
                                unlist(),
		"S2Kb" = to_plot %>%
                                filter(is_method_3) %>%
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
		"GTRD" =  msigdb %>%
					dplyr::distinct() %>%
                                        mutate(pair = paste0(tf, "_", target)) %>%
                                        dplyr::select(pair) %>% 
					dplyr::distinct() %>%
                                        unlist())

svg('plots/overlap_with_db_K562.svg', width = 8.268)

p = upset(fromList(test_list), 
	sets = c("S2Mb", "M2Kb", "S2Kb", "ChIP-Atlas", "Dorothea", "GTRD"), 
	mb.ratio = c(0.7, 0.3), 
	order.by="freq", 
	mainbar.y.label = "Pairs intersection", 
	sets.x.label = "TF-target pairs",
	show.numbers = F,
	text.scale = 1.5,
	scale.sets = "log10", 
	nintersects=20)
p
dev.off()


msigdb %>% distinct()  %>% nrow
# 256331 average GTRD regulon size

length(unique(msigdb$tf))
# 502 TFs 

#ChIP-Atlas

chip_atlas = fread("data/trans_networks/ChIP-Atlas_target_genes_K562.tsv", header = T)
chip_atlas %>% distinct()  %>% nrow
# 1529999 average ChIP-Atlas K-562 regulon size

length(unique(chip_atlas$tf))
# 400 TFs

chip_atlas = fread("data/trans_networks/ChIP-Atlas_target_genes_HepG2.tsv", header = T)
chip_atlas %>% distinct()  %>% nrow
# 2068087 average ChIP-Atlas HepG2 regulon size

length(unique(chip_atlas$tf))
# 366 TFs

chip_atlas = fread("data/trans_networks/ChIP-Atlas_target_genes_MCF7.tsv", header = T)
chip_atlas %>% distinct()  %>% nrow 
# 639027 average ChIP-Atlas MCF-7 regulon size
 
length(unique(chip_atlas$tf))
# 182 TFs

chip_atlas = fread("data/trans_networks/ChIP-Atlas_target_genes_GM12878.tsv", header = T)
chip_atlas %>% distinct()  %>% nrow  
# 655040 average ChIP-Atlas GM-12878 regulon size
 
length(unique(chip_atlas$tf)) 
# 156 TFs


ga_normal %>% select(TF, target) %>% distinct()  %>% nrow
# 1076628 average GTRD regulon size

length(unique(ga_normal$TF))
# 1402 TFs



p$New_data %>% dplyr::select(S2Mb, S2Kb, M2Kb) %>% filter((S2Mb ==1) & (S2Kb == 1) & (M2Kb == 1)) %>% nrow / nrow(p$New_data)
# 0.1709467 of interactions confirmd by all 3 approaches in K-562
# 0.4185891 of interactions confirmd by all 3 approaches in Hep-G2
# 0.2641069 of interactions confirmd by all 3 approaches in MCF-7
# 0.1462951 of interactions confirmd by all 3 approaches in GM-12878


p$New_data %>% filter((`ChIP-Atlas` == 1) & (S2Mb ==1) & (S2Kb == 1) & (M2Kb == 1)) %>% nrow / nrow(p$New_data)
# 0.05322095 of interactions confirmd by all 3 approaches and ChIP-Atlas in K-562
# 0.234045 of interactions confirmd by all 3 approaches and ChIP-Atlas in Hep-G2
# 0.06009169 of interactions confirmd by all 3 approaches and ChIP-Atlas in MCF-7
# 0.06164239 of interactions confirmd by all 3 approaches and ChIP-Atlas in GM-12878

p$New_data %>% filter((Dorothea == 1) & ((S2Mb ==1) | (S2Kb == 1) | (M2Kb == 1) | (`ChIP-Atlas` == 1) | (GTRD == 1))) %>% nrow / nrow(p$New_data)
# 0.04156506 of interactions confirmd by Dorothea and any other database in K-562
# 0.02340525 of interactions confirmd by Dorothea and any other database in Hep-G2
# 0.02815413 of interactions confirmd by Dorothea and any other database in MCF-7    
# 0.05214271 of interactions confirmd by Dorothea and any other database in GM-12878 

p$New_data %>% filter((Dorothea == 1) & (GTRD == 1)) %>% nrow / nrow(p$New_data)
# 0.01067344 of interactions confirmd by Dorothea and GTRD in K-562

p$New_data %>% rowSums() -> temp
sum(temp > 1) / nrow(p$New_data) 
# 0.2728547 interactions were present in at least two datasets (K-562)
# 0.3571928 interactions were present in at least two datasets (MCF-7)
# 0.5837918 interactions were present in at least two datasets (HepG2)
# 0.2496658 interactions were present in at least two datasets (GM-12878)
