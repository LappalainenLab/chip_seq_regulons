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

r2Log <- function(model) {

  summaryLog <- summary(model)
  1 - summaryLog$deviance / summaryLog$null.deviance

}

# Plot style
source("figures/fig_style.R")

# Load regulons
tf_target_mapping = fread("data/regulons/TF_target_mapping_filtered_merged_K562_with_ppi_with_dnase_with_atac_with_motifs_with_ccres_cleaned.tsv", nThread=10)
tf_target_mapping$is_motif = (!is.na(tf_target_mapping$n_motifs) | !is.na(tf_target_mapping$n_motifs_hoco))


networks = fread("data/trans_networks/STING_Seq_networks.tsv")
networks$is_network = T

tf_target_mapping %>% distinct(gene_symbol, tpm_total, ensembl_gene_id, is_coexpressed) %>% group_by(gene_symbol) %>% summarize(tpm_total = max(tpm_total)) -> tested_genes 

tfs = unique(tf_target_mapping$tf)


enrich_scores = foreach(TF = tfs) %dopar% {
            print(paste0("Started: ", TF))
            enrich_temp = list()
            tmp = c()
    
            tf_target_mapping %>% filter(tf == TF) -> sub_tf_target_mapping

            networks %>% filter(tf == TF) -> sub_networks  

            sub_tf_target_mapping %>% select(is_S2Mb, is_M2Kb, is_S2Kb, is_ppi, is_motif, gene_symbol, is_coexpressed) %>% distinct() -> sub_tf_target_mapping_acc
            sub_tf_target_mapping_acc %>% group_by(gene_symbol) %>% 
                                    summarise(across(c(is_S2Mb, 
                                                       is_M2Kb, 
                                                       is_S2Kb, 
                                                       is_ppi, 
                                                       is_motif,
						       is_coexpressed), function(x){as.logical(sum(x))})) -> sub_tf_target_mapping_acc
            left_join(tested_genes, sub_networks, by = c('gene_symbol')) -> sub_tested_genes
            left_join(sub_tested_genes, sub_tf_target_mapping_acc) -> sub_tested_genes                      	
            sub_tested_genes %>% select(-tf) -> sub_tested_genes
            sub_tested_genes %>% replace(is.na(.), 0) -> sub_tested_genes


		#S2Mb

            log_model_target_m1 = glm("is_S2Mb ~ is_coexpressed + tpm_total",
                                    data = sub_tested_genes, family="binomial")

            odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1]))
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] - (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] + (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,4])

            tmp = cbind(tmp, c("is_coexpressed", "method1", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))
            log2(odds_ratio)
            print("R^2")
            print(r2Log(log_model_target_m1))



            log_model_target_m1 = glm("is_S2Mb ~ is_motif + tpm_total",
                                    data = sub_tested_genes, family="binomial")

            odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1]))
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] - (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] + (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,4])
    
            tmp = cbind(tmp, c("is_motif", "method1", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))
            log2(odds_ratio)
            print("R^2")
            print(r2Log(log_model_target_m1))

			log_model_target_m1 = glm("is_S2Mb ~ is_ppi + tpm_total",
                                      data = sub_tested_genes, family="binomial")

            odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1]))
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] - (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] + (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,4])

            tmp = cbind(tmp, c("is_ppi", "method1", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))
            log2(odds_ratio)
            print("R^2")
            print(r2Log(log_model_target_m1))

			log_model_target_m1 = glm(paste0("is_S2Mb ~ is_network + tpm_total"),
        	                                      data = sub_tested_genes, family="binomial")
       	    odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1]))
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] - (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] + (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,4])
            tmp = cbind(tmp, c("network", "method1", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))
            log2(odds_ratio)
            print("R^2")
            print(r2Log(log_model_target_m1))

			#---------------------------------------------------------------------------------------------------------------------------------------------------------------
			# M2Kb
            log_model_target_m2 = glm("is_M2Kb ~ is_coexpressed + tpm_total",
                                    data = sub_tested_genes, family="binomial")

            odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1]))
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] - (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] + (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,4])

            tmp = cbind(tmp, c("is_coexpressed", "method2", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))
            log2(odds_ratio)
            print("R^2")
            print(r2Log(log_model_target_m2))



            log_model_target_m2 = glm("is_M2Kb ~ is_motif + tpm_total",
                                    data = sub_tested_genes, family="binomial") 

            odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1])) 
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] - (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] + (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,4]) 
            log2(odds_ratio) 
            pvalue 

            tmp = cbind(tmp, c("is_motif", "method2", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))

            print("R^2")
            print(r2Log(log_model_target_m2))

			log_model_target_m2 = glm("is_M2Kb ~ is_ppi + tpm_total",
                                                data = sub_tested_genes, family="binomial")

            odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1]))
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] - (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] + (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,4])
            log2(odds_ratio) 
            pvalue 

            tmp = cbind(tmp, c("is_ppi", "method2", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))

            print("R^2")
            print(r2Log(log_model_target_m2))
    
            log_model_target_m2 = glm(paste0("is_M2Kb ~ is_network + tpm_total"),
                                          data = sub_tested_genes, family="binomial")
            odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1]))
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] - (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] + (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,4])
            tmp = cbind(tmp, c("network", "method2", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))

            print("R^2")
            print(r2Log(log_model_target_m2))

                        #---------------------------------------------------------------------------------------------------------------------------------------------------------------                        
			# S2Kb
            log_model_target_m3 = glm("is_S2Kb ~ is_coexpressed + tpm_total",
                                    data = sub_tested_genes, family="binomial")

            odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1]))
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] - (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] + (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,4])

            tmp = cbind(tmp, c("is_coexpressed", "method3", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))
            log2(odds_ratio)
            print("R^2")
            print(r2Log(log_model_target_m3))




            log_model_target_m3 = glm("is_S2Kb ~ is_motif + tpm_total",
                                    data = sub_tested_genes, family="binomial") 

            odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1])) 
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] - (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] + (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,4]) 
            log2(odds_ratio) 
            pvalue 

            tmp = cbind(tmp, c("is_motif", "method3", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))

            print("R^2")
            print(r2Log(log_model_target_m3))

			log_model_target_m3 = glm("is_S2Kb ~ is_ppi + tpm_total",
                                                data = sub_tested_genes, family="binomial")

            odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1]))
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] - (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] + (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,4])
            log2(odds_ratio) 
            pvalue 

            tmp = cbind(tmp, c("is_ppi", "method3", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))

            print("R^2")
            print(r2Log(log_model_target_m3))

log_model_target_m3 = glm(paste0("is_S2Kb ~ is_network + tpm_total"),
                                          data = sub_tested_genes, family="binomial")
            odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1]))
            conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] - (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
            conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] + (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
            pvalue = as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,4])
            tmp = cbind(tmp, c("network", "method3", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))

            print("R^2")
            print(r2Log(log_model_target_m3))
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------
    	print(paste0("Done with ", TF))
        enrich_temp = matrix(tmp, ncol=7, byrow=TRUE)
        as.data.frame(enrich_temp, stringsAsFactors=FALSE)
}

enrich_scores = rbindlist(enrich_scores)
colnames(enrich_scores) = c("variable", "method", "TF", "odds", "conf_int_l", "conf_int_u", "pvalue")



fwrite(enrich_scores, "data/1-network_enrichment/enrich_scores_remap_all_tfs_K562.tsv", col.names = T, row.names = F, quote = F, sep="\t")
fwrite(enrich_scores, "data/s3-network_enrichment/enrich_scores_remap_all_tfs_K562.tsv", col.names = T, row.names = F, quote = F, sep="\t")


