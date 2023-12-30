library(data.table)
library(here)
library(jtools)
library(stringr)
library(biomaRt)
library(fastglm)
library(dplyr)
library(doParallel)
library(foreach)
library(tidyr)

registerDoParallel(cores=10)

here::i_am("README.md")


# Plot style
source("figures/fig_style.R")

# Load regulons
tf_target_mapping = fread("data/regulons/K562_all_regulons.tsv", nThread=10)

shuffled_net = tf_target_mapping
shuffled_net[, c("gene_symbol", "tpm_total", "ensembl_gene_id")] = sample_n(shuffled_net[, c("gene_symbol", "tpm_total", "ensembl_gene_id")], nrow(shuffled_net), replace = TRUE)


networks = fread("data/trans_networks/STING_Seq_networks.tsv")
networks$is_network = T

tf_target_mapping %>% distinct(gene_symbol, tpm_total, ensembl_gene_id) -> tested_genes   

tfs = unique(tf_target_mapping$tf)


enrich_scores = foreach(TF = tfs) %dopar% {
                print(paste0("Started: ", TF))
                enrich_temp = list()

                tf_target_mapping %>% filter(tf == TF) -> sub_tf_target_mapping
                shuffled_net %>% filter(tf == TF) -> sub_shuffled_net

		networks %>% filter(tf == TF) -> sub_networks
		accessions = unique(sub_tf_target_mapping$accession)                

                enrich_temp= sapply(accessions, function(x){

                        tmp = c()
			
			#subset networks
                        sub_tf_target_mapping_acc = sub_tf_target_mapping[sub_tf_target_mapping$accession == x, ]
			
			sub_tf_target_mapping_acc %>% group_by(ensembl_gene_id) %>% 
				summarise(across(c(is_method_1, is_method_2, is_method_3,
							is_coexpressed, 
							is_ppi), function(x){as.logical(sum(x))})) -> sub_tf_target_mapping_acc
			left_join(tested_genes, sub_networks, by = c('gene_symbol')) -> sub_tested_genes
			left_join(sub_tested_genes, sub_tf_target_mapping_acc) -> sub_tested_genes                      	
			
			sub_tested_genes %>% replace(is.na(.), 0) -> sub_tested_genes
			
			# Subset shuffled net
                        sub_shuffled_net_acc = sub_shuffled_net[sub_shuffled_net$accession == x, ]

                        sub_shuffled_net_acc %>% group_by(ensembl_gene_id) %>%
                                summarise(across(c(is_method_1, is_method_2, is_method_3,
                                                        is_coexpressed,
                                                        is_ppi), function(x){as.logical(sum(x))})) -> sub_shuffled_net_acc
                        left_join(tested_genes, sub_networks, by = c('gene_symbol')) -> sub_tested_genes_shuffled
                        left_join(sub_tested_genes_shuffled, sub_shuffled_net_acc) -> sub_tested_genes_shuffled

                        sub_tested_genes_shuffled %>% replace(is.na(.), 0) -> sub_tested_genes_shuffled


			# S2Mb
                        log_model_target_m1 = glm("is_method_1 ~ is_coexpressed + tpm_total",
                                                data = sub_tested_genes, family="binomial")

                        odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1]))
                        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] - (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
                        conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] + (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
                        pvalue = as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,4])
                        log2(odds_ratio)
                        pvalue

                        tmp = cbind(tmp, c("is_coexpressed", "method1", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))


			log_model_target_m1 = glm("is_method_1 ~ is_ppi + tpm_total",
                                                data = sub_tested_genes, family="binomial")

                        odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1]))
                        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] - (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
                        conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] + (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
                        pvalue = as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,4])
                        log2(odds_ratio) 
                        pvalue 

                        tmp = cbind(tmp, c("is_ppi", "method1", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))


			log_model_target_m1 = glm(paste0("is_method_1 ~ is_network + tpm_total"),
        	                                      data = sub_tested_genes, family="binomial")
       	                odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1]))
               	        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] - (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
                       	conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] + (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
                       	pvalue = as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,4])
                        tmp = cbind(tmp, c("network", "method1", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))


			#---------------------------------------------------------------------------------------------------------------------------------------------------------------
			# M2Kb

                        log_model_target_m2 = glm("is_method_2 ~ is_coexpressed + tpm_total",
                                                data = sub_tested_genes, family="binomial") 
 
                        odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1])) 
                        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] - (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
                        conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] + (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
                        pvalue = as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,4]) 
                        log2(odds_ratio) 
                        pvalue 
 
                        tmp = cbind(tmp, c("is_coexpressed", "method2", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))


			log_model_target_m2 = glm("is_method_2 ~ is_ppi + tpm_total",
                                                data = sub_tested_genes, family="binomial")

                        odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1]))
                        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] - (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
                        conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] + (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
                        pvalue = as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,4])
                        log2(odds_ratio) 
                        pvalue 

                        tmp = cbind(tmp, c("is_ppi", "method2", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))
			
				
			log_model_target_m2 = glm(paste0("is_method_2 ~ is_network + tpm_total"),
                                                      data = sub_tested_genes, family="binomial")
                        odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1]))
                        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] - (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
                        conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] + (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
                        pvalue = as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,4])
                        tmp = cbind(tmp, c("network", "method2", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))

                        #---------------------------------------------------------------------------------------------------------------------------------------------------------------                        
			# S2Kb

                        log_model_target_m3 = glm("is_method_3 ~ is_coexpressed + tpm_total",
                                                data = sub_tested_genes, family="binomial") 
 
                        odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1])) 
                        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] - (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
                        conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] + (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
                        pvalue = as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,4]) 
                        log2(odds_ratio) 
                        pvalue 
 
                        tmp = cbind(tmp, c("is_coexpressed", "method3", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))


			log_model_target_m3 = glm("is_method_3 ~ is_ppi + tpm_total",
                                                data = sub_tested_genes, family="binomial")

                        odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1]))
                        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] - (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
                        conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] + (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
                        pvalue = as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,4])
                        log2(odds_ratio) 
                        pvalue 

                        tmp = cbind(tmp, c("is_ppi", "method3", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))


			log_model_target_m3 = glm(paste0("is_method_3 ~ is_network + tpm_total"),
                                                      data = sub_tested_genes, family="binomial")
                        odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1]))
                        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] - (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
                        conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] + (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
                        pvalue = as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,4])
                        tmp = cbind(tmp, c("network", "method3", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))

                        #---------------------------------------------------------------------------------------------------------------------------------------------------------------
                        # Shuffled network

                        log_model_target_shuf = glm("is_method_1 ~ is_coexpressed + tpm_total",
                                                data = sub_tested_genes_shuffled, family="binomial")

                        odds_ratio = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1]))
                        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1] - (data.table(summary(log_model_target_shuf)$coefficients)[2,2] * 1.96)))
                        conf_int_u = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1] + (data.table(summary(log_model_target_shuf)$coefficients)[2,2] * 1.96)))
                        pvalue = as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,4])
                        log2(odds_ratio)
                        pvalue

                        tmp = cbind(tmp, c("is_coexpressed", "Random", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))


                        log_model_target_shuf = glm("is_method_1 ~ is_ppi + tpm_total",
                                                data = sub_tested_genes_shuffled, family="binomial")

                        odds_ratio = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1]))
                        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1] - (data.table(summary(log_model_target_shuf)$coefficients)[2,2] * 1.96)))
                        conf_int_u = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1] + (data.table(summary(log_model_target_shuf)$coefficients)[2,2] * 1.96)))
                        pvalue = as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,4])
                        log2(odds_ratio)
                        pvalue

                        tmp = cbind(tmp, c("is_ppi", "Random", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))


                        log_model_target_shuf = glm(paste0("is_method_1 ~ is_network + tpm_total"),
                                                      data = sub_tested_genes_shuffled, family="binomial")
                        odds_ratio = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1]))
                        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1] - (data.table(summary(log_model_target_shuf)$coefficients)[2,2] * 1.96)))
                        conf_int_u = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1] + (data.table(summary(log_model_target_shuf)$coefficients)[2,2] * 1.96)))
                        pvalue = as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,4])
                        tmp = cbind(tmp, c("network", "Random", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))
                        #---------------------------------------------------------------------------------------------------------------------------------------------------------------
	})
	print(paste0("Done with ", TF))
        enrich_temp = matrix(enrich_temp, ncol=8, byrow=TRUE)
        as.data.frame(enrich_temp, stringsAsFactors=FALSE)
}

enrich_scores = rbindlist(enrich_scores)
colnames(enrich_scores) = c("variable", "method", "TF", "odds", "conf_int_l", "conf_int_u", "pvalue", "access")



fwrite(enrich_scores, "data/1-network_enrichment/enrich_scores_remap_all_tfs_K562.tsv", col.names = T, row.names = F, quote = F, sep="\t")
fwrite(enrich_scores, "data/s3-network_enrichment/enrich_scores_remap_all_tfs_K562.tsv", col.names = T, row.names = F, quote = F, sep="\t")


