#!/usr/bin/env Rscript

# Load necessary libraries
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
library(predictmeans)

# Set up parallel processing with 10 cores
registerDoParallel(cores = 10)

# Set the working directory to the project root
setwd(here())

# Helper function to extract and format logistic model results
extract_model_results <- function(model, variable, method, TF) {
  summary_data <- data.table(summary(model)$coefficients)
  
  odds_ratio <- exp(summary_data[2, 1])
  conf_int_l <- exp(summary_data[2, 1] - summary_data[2, 2] * 1.96)
  conf_int_u <- exp(summary_data[2, 1] + summary_data[2, 2] * 1.96)
  pvalue <- as.numeric(summary_data[2, 4])
  
  return(c(variable, method, TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))
}

# Load the custom plot style
source("figures/fig_style.R")

# Load regulons and mark motif presence
tf_target_mapping <- fread("data/regulons/TF_target_mapping_filtered_merged_K562_with_ppi_with_dnase_with_atac_with_motifs_with_ccres_cleaned.tsv", nThread = 10)
tf_target_mapping$is_motif <- (!is.na(tf_target_mapping$n_motifs) | !is.na(tf_target_mapping$n_motifs_hoco))

# Load network data and mark them as networks
networks <- fread("data/trans_networks/STING_Seq_networks.tsv")
networks$is_network <- TRUE

# Filter for distinct tested genes based on gene_symbol
tested_genes <- tf_target_mapping %>%
  distinct(gene_symbol, tpm_total, ensembl_gene_id, is_coexpressed) %>%
  group_by(gene_symbol) %>%
  summarize(tpm_total = max(tpm_total))

# List of transcription factors (TFs) to analyze
tfs <- c("IKZF1", "GFI1B", "NFE2", "RUNX1")

# Calculate enrichment scores for each TF
enrich_scores = foreach(TF = tfs) %dopar% {
  print(paste0("Started: ", TF))
  tmp <- list()
  
  # Filter data for the specific TF
  sub_tf_target_mapping <- tf_target_mapping %>% filter(tf == TF)
  sub_networks <- networks %>% filter(tf == TF)
  
  # Select and clean relevant columns for model analysis
  sub_tf_target_mapping_acc <- sub_tf_target_mapping %>%
    select(is_S2Mb, is_M2Kb, is_S2Kb, is_S100Kb, is_M100Kb, is_ppi, is_motif, gene_symbol, is_coexpressed) %>%
    distinct() %>%
    group_by(gene_symbol) %>%
    summarize(across(c(is_S2Mb, is_M2Kb, is_S2Kb, is_S100Kb, is_M100Kb, is_ppi, is_motif, is_coexpressed),
                     ~ as.logical(sum(.))))
  
  # Join gene data with networks and tested genes
  sub_tested_genes <- tested_genes %>%
    left_join(sub_networks, by = "gene_symbol") %>%
    left_join(sub_tf_target_mapping_acc, by = "gene_symbol") %>%
    select(-tf) %>%
    replace(is.na(.), 0)
# -----------------------------------------------------------------------------------------------------------------------------
  # Logistic regression models for various annotations
  # S2Mb

  # Coexpression model for S2Mb
  s2mb_model <- glm("is_S2Mb ~ is_coexpressed + tpm_total", data = sub_tested_genes, family = "binomial")
  s2mb_results <- extract_model_results(s2mb_model, "is_coexpressed", "method1", TF)
  tmp <- cbind(tmp, s2mb_results)
  
  # Motif model for S2Mb
  s2mb_model <- glm("is_S2Mb ~ is_motif + tpm_total", data = sub_tested_genes, family = "binomial")
  motif_results <- extract_model_results(s2mb_model, "is_motif", "method1", TF)
  tmp <- cbind(tmp, motif_results)
  
  # PPI model for S2Mb
  s2mb_model <- glm("is_S2Mb ~ is_ppi + tpm_total", data = sub_tested_genes, family = "binomial")
  ppi_results <- extract_model_results(s2mb_model, "is_ppi", "method1", TF)
  tmp <- cbind(tmp, ppi_results)
  
  # Network model for S2Mb
  s2mb_model <- glm("is_S2Mb ~ is_network + tpm_total", data = sub_tested_genes, family = "binomial")
  network_results <- extract_model_results(s2mb_model, "is_network", "method1", TF)
  tmp <- cbind(tmp, network_results)

# -----------------------------------------------------------------------------------------------------------------------------

  # M2Kb

  # Coexpression model for M2Kb
  m2kb_model <- glm("is_M2Kb ~ is_coexpressed + tpm_total", data = sub_tested_genes, family = "binomial")
  m2kb_results <- extract_model_results(m2kb_model, "is_coexpressed", "method2", TF)
  tmp <- cbind(tmp, m2kb_results)
  
  # Motif model for M2Kb
  m2kb_model <- glm("is_M2Kb ~ is_motif + tpm_total", data = sub_tested_genes, family = "binomial")
  motif_results <- extract_model_results(m2kb_model, "is_motif", "method2", TF)
  tmp <- cbind(tmp, motif_results)
  
  # PPI model for M2Kb
  m2kb_model <- glm("is_M2Kb ~ is_ppi + tpm_total", data = sub_tested_genes, family = "binomial")
  ppi_results <- extract_model_results(m2kb_model, "is_ppi", "method2", TF)
  tmp <- cbind(tmp, ppi_results)
  
  # Network model for M2Kb
  m2kb_model <- glm("is_M2Kb ~ is_network + tpm_total", data = sub_tested_genes, family = "binomial")
  network_results <- extract_model_results(m2kb_model, "is_network", "method2", TF)
  tmp <- cbind(tmp, network_results)
# -----------------------------------------------------------------------------------------------------------------------------

  # S2Kb

  # Coexpression model for S2Kb
  s2kb_model <- glm("is_S2Kb ~ is_coexpressed + tpm_total", data = sub_tested_genes, family = "binomial")
  s2kb_results <- extract_model_results(s2kb_model, "is_coexpressed", "method3", TF)
  tmp <- cbind(tmp, s2kb_results)
  
  # Motif model for S2Kb
  s2kb_model <- glm("is_S2Kb ~ is_motif + tpm_total", data = sub_tested_genes, family = "binomial")
  motif_results <- extract_model_results(s2kb_model, "is_motif", "method3", TF)
  tmp <- cbind(tmp, motif_results)
  
  # PPI model for S2Kb
  s2kb_model <- glm("is_S2Kb ~ is_ppi + tpm_total", data = sub_tested_genes, family = "binomial")
  ppi_results <- extract_model_results(s2kb_model, "is_ppi", "method3", TF)
  tmp <- cbind(tmp, ppi_results)
  
  # Network model for S2Kb
  s2kb_model <- glm("is_S2Kb ~ is_network + tpm_total", data = sub_tested_genes, family = "binomial")
  network_results <- extract_model_results(s2kb_model, "is_network", "method3", TF)
  tmp <- cbind(tmp, network_results)

# -----------------------------------------------------------------------------------------------------------------------------

  # M100Kb

  # Coexpression model for M100Kb
  m100kb_model <- glm("is_M100Kb ~ is_coexpressed + tpm_total", data = sub_tested_genes, family = "binomial")
  m100kb_results <- extract_model_results(m100kb_model, "is_coexpressed", "method5", TF)
  tmp <- cbind(tmp, m100kb_results)
  
  # Motif model for M100Kb
  m100kb_model <- glm("is_M100Kb ~ is_motif + tpm_total", data = sub_tested_genes, family = "binomial")
  motif_results <- extract_model_results(m100kb_model, "is_motif", "method5", TF)
  tmp <- cbind(tmp, motif_results)
  
  # PPI model for M100Kb
  m100kb_model <- glm("is_M100Kb ~ is_ppi + tpm_total", data = sub_tested_genes, family = "binomial")
  ppi_results <- extract_model_results(m100kb_model, "is_ppi", "method5", TF)
  tmp <- cbind(tmp, ppi_results)
  
  # Network model for M100Kb
  m100kb_model <- glm("is_M100Kb ~ is_network + tpm_total", data = sub_tested_genes, family = "binomial")
  network_results <- extract_model_results(m100kb_model, "is_network", "method5", TF)
  tmp <- cbind(tmp, network_results)
# -----------------------------------------------------------------------------------------------------------------------------

  # S100Kb

  # Coexpression model for S100Kb
  s100kb_model <- glm("is_S100Kb ~ is_coexpressed + tpm_total", data = sub_tested_genes, family = "binomial")
  s100kb_results <- extract_model_results(s100kb_model, "is_coexpressed", "method4", TF)
  tmp <- cbind(tmp, s100kb_results)
  
  # Motif model for S100Kb
  s100kb_model <- glm("is_S100Kb ~ is_motif + tpm_total", data = sub_tested_genes, family = "binomial")
  motif_results <- extract_model_results(s100kb_model, "is_motif", "method4", TF)
  tmp <- cbind(tmp, motif_results)
  
  # PPI model for S100Kb
  s100kb_model <- glm("is_S100Kb ~ is_ppi + tpm_total", data = sub_tested_genes, family = "binomial")
  ppi_results <- extract_model_results(s100kb_model, "is_ppi", "method4", TF)
  tmp <- cbind(tmp, ppi_results)
  
  # Network model for S100Kb
  s100kb_model <- glm("is_S100Kb ~ is_network + tpm_total", data = sub_tested_genes, family = "binomial")
  network_results <- extract_model_results(s100kb_model, "is_network", "method4", TF)
  tmp <- cbind(tmp, network_results)

# -----------------------------------------------------------------------------------------------------------------------------

  print(paste0("Done with ", TF))
  enrich_temp = matrix(tmp, ncol=7, byrow=TRUE)
  return(as.data.frame(enrich_temp, stringsAsFactors=FALSE))
}

# Convert the result to a proper data.table with appropriate column names
enrich_scores = rbindlist(enrich_scores)
colnames(enrich_scores) = c("variable", "method", "TF", "odds", "conf_int_l", "conf_int_u", "pvalue")


# Write the results to a file
fwrite(enrich_scores, "data/1-network_enrichment/enrich_scores_remap_all_tfs_K562.tsv", col.names = T, row.names = F, quote = F, sep="\t")

cat("Analysis completed and results saved to data/1-network_enrichment/enrich_scores_remap_all_tfs_K562.tsv\n")

