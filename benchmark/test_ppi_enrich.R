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

cell_line = c("HepG2", "MCF7", "GM12878", "K562")

foreach(cl = cell_line) %do% {
	
	# Load regulons
	input_file = sprintf("data/regulons/TF_target_mapping_filtered_merged_%s_with_ppi_with_dnase_with_atac_with_motifs_with_ccres_cleaned.tsv", cl)
	tf_target_mapping = fread(input_file, nThread=10)
	
	tf_target_mapping %>% select(c(n_motifs, 
					n_motifs_hoco, 
					is_dnase, 
					is_S2Mb, 
					is_S2Kb, 
					is_M2Kb, 
					is_atac, 
					ensembl_gene_id, 
					tpm_total, 
					is_ppi, 
					tf, 
					gene_symbol, 
					accession)) %>% 
				distinct() -> tf_target_mapping	
	shuffled_net = tf_target_mapping
	shuffled_net[, c("gene_symbol", "tpm_total", "ensembl_gene_id")] = sample_n(shuffled_net[, c("gene_symbol", "tpm_total", "ensembl_gene_id")], nrow(shuffled_net), replace = TRUE)
	
	
	tf_target_mapping %>% distinct(gene_symbol, tpm_total, ensembl_gene_id) -> tested_genes   
	
	tfs = unique(tf_target_mapping$tf)
	
	
	enrich_scores = foreach(TF = tfs) %dopar% {
		                print(paste0("Started: ", TF))
		                enrich_temp = list()
		
	        	        tf_target_mapping %>% filter(tf == TF) -> sub_tf_target_mapping
	                	shuffled_net %>% filter(tf == TF) -> sub_shuffled_net
	
				accessions = unique(sub_tf_target_mapping$accession)                
		
		                enrich_temp= sapply(accessions, function(x){
		
	                        tmp = c()			
				#subset networks
                	        sub_tf_target_mapping_acc = sub_tf_target_mapping[sub_tf_target_mapping$accession == x, ]
				
				sub_tf_target_mapping_acc %>% group_by(ensembl_gene_id) %>% 
				summarise(across(c(is_S2Mb, is_M2Kb, is_S2Kb,
								is_ppi), function(x){as.logical(sum(x))})) -> sub_tf_target_mapping_acc
				left_join(tested_genes, sub_tf_target_mapping_acc) -> sub_tested_genes                      	
				
				sub_tested_genes %>% replace(is.na(.), 0) -> sub_tested_genes
				
				# Subset shuffled net
                	        sub_shuffled_net_acc = sub_shuffled_net[sub_shuffled_net$accession == x, ]
	
       		                sub_shuffled_net_acc %>% group_by(ensembl_gene_id) %>%
               		                summarise(across(c(is_S2Mb, is_M2Kb, is_S2Kb,
                               		                        is_ppi), function(x){as.logical(sum(x))})) -> sub_shuffled_net_acc
       		                left_join(tested_genes, sub_shuffled_net_acc) -> sub_tested_genes_shuffled
	
       		                sub_tested_genes_shuffled %>% replace(is.na(.), 0) -> sub_tested_genes_shuffled
		
	
				# S2Mb
	
				log_model_target_m1 = glm("is_S2Mb ~ is_ppi + tpm_total",
               		                                data = sub_tested_genes, family="binomial")
	
       		                odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1]))
               		        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] - (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
                       		conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,1] + (data.table(summary(log_model_target_m1)$coefficients)[2,2] * 1.96)))
	                        pvalue = as.numeric(data.table(summary(log_model_target_m1)$coefficients)[2,4])
       		                log2(odds_ratio) 
               		        pvalue 
	
       		                tmp = cbind(tmp, c("is_ppi", "S2Mb", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))
				
				print("R^2")
				print(r2Log(log_model_target_m1))
	
				#---------------------------------------------------------------------------------------------------------------------------------------------------------------
				# M2Kb
	
	
				log_model_target_m2 = glm("is_M2Kb ~ is_ppi + tpm_total",
               		                                data = sub_tested_genes, family="binomial")
	
       		                odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1]))
               		        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] - (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
                       		conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,1] + (data.table(summary(log_model_target_m2)$coefficients)[2,2] * 1.96)))
	                        pvalue = as.numeric(data.table(summary(log_model_target_m2)$coefficients)[2,4])
       		                log2(odds_ratio) 
               		        pvalue 
	
       		                tmp = cbind(tmp, c("is_ppi", "M2Kb", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))
				
                                print("R^2")
                                print(r2Log(log_model_target_m2))

        	                #---------------------------------------------------------------------------------------------------------------------------------------------------------------                        
				# S2Kb
	

				log_model_target_m3 = glm("is_S2Kb ~ is_ppi + tpm_total",
       	        	                                data = sub_tested_genes, family="binomial")
	
       		                odds_ratio = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1]))
               		        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] - (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
                       		conf_int_u = exp(as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,1] + (data.table(summary(log_model_target_m3)$coefficients)[2,2] * 1.96)))
	 	                pvalue = as.numeric(data.table(summary(log_model_target_m3)$coefficients)[2,4])
       	        	        log2(odds_ratio) 
               	        	pvalue 
	
       		                tmp = cbind(tmp, c("is_ppi", "S2Kb", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))

                                print("R^2")
                                print(r2Log(log_model_target_m3))

       		                #---------------------------------------------------------------------------------------------------------------------------------------------------------------
               		        # Shuffled network
	

        	                log_model_target_shuf = glm("is_S2Mb ~ is_ppi + tpm_total",
       	        	                                data = sub_tested_genes_shuffled, family="binomial")
	
       		                odds_ratio = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1]))
               		        conf_int_l = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1] - (data.table(summary(log_model_target_shuf)$coefficients)[2,2] * 1.96)))
                       		conf_int_u = exp(as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,1] + (data.table(summary(log_model_target_shuf)$coefficients)[2,2] * 1.96)))
	                        pvalue = as.numeric(data.table(summary(log_model_target_shuf)$coefficients)[2,4])
       		                log2(odds_ratio)
               		        pvalue
	
       		                tmp = cbind(tmp, c("is_ppi", "Random", TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue, x))

                                print("R^2")
                                print(r2Log(log_model_target_shuf))
				
				tmp		
		
	                        #---------------------------------------------------------------------------------------------------------------------------------------------------------------
		})
		print(paste0("Done with ", TF))
        	enrich_temp = matrix(enrich_temp, ncol=8, byrow=TRUE)
	        as.data.frame(enrich_temp, stringsAsFactors=FALSE)
	}

	enrich_scores = rbindlist(enrich_scores)
	colnames(enrich_scores) = c("variable", "method", "TF", "odds", "conf_int_l", "conf_int_u", "pvalue", "access")


	output_file = sprintf("data/s3-network_enrichment/enrich_scores_remap_all_tfs_%s.tsv", cl)
	fwrite(enrich_scores, output_file, col.names = T, row.names = F, quote = F, sep="\t")
}
