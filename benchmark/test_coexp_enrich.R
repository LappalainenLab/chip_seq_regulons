#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(here)
library(dplyr)
library(foreach)
library(doParallel)
library(tidyr)
library(predictmeans)

# Register parallel backend with 3 cores
registerDoParallel(cores = 3)

# Set working directory to the project's root
setwd(here())

# Function to calculate pseudo-R2 for logistic regression
r2Log <- function(model) {
  summaryLog <- summary(model)
  1 - summaryLog$deviance / summaryLog$null.deviance
}

# Load plot style settings
source("figures/fig_style.R")

# Set seed for reproducibility
set.seed(123)

# Define cell lines to iterate over
cell_lines <- c("K562")

# Iterate over cell lines and perform enrichment analysis
foreach(cl = cell_lines) %do% {
  
  # Load regulon data for the current cell line
  input_file <- sprintf("data/regulons/TF_target_mapping_filtered_merged_%s_with_ppi_with_dnase_with_atac_with_motifs_with_ccres_cleaned.tsv", cl)
  tf_target_mapping <- fread(input_file, nThread = 10)

  # Filter and keep only relevant columns, removing duplicates
  tf_target_mapping <- tf_target_mapping %>%
    select(n_motifs, n_motifs_hoco, is_dnase, is_S2Mb, is_S2Kb, is_M2Kb, is_S100Kb, is_M100Kb,
           is_atac, ensembl_gene_id, tpm_total, is_coexpressed, tf, gene_symbol, accession) %>%
    distinct()

  # Create a list of tested genes (unique gene symbols, tpm_total, and ensembl_gene_id)
  tested_genes <- tf_target_mapping %>%
    distinct(gene_symbol, tpm_total, ensembl_gene_id)

  # Extract unique transcription factors (TFs)
  tfs <- unique(tf_target_mapping$tf)

  # Perform enrichment analysis for each TF in parallel
  enrich_scores <- foreach(TF = tfs) %dopar% {
    
    print(paste0("Started: ", TF))

    # Filter target mapping for the current TF
    sub_tf_target_mapping <- tf_target_mapping %>% filter(tf == TF)

    # Get unique accessions for the current TF
    accessions <- unique(sub_tf_target_mapping$accession)

    # Perform enrichment analysis for each accession
    enrich_temp <- sapply(accessions, function(accession) {
      
      # Subset by accession
      sub_tf_target_mapping_acc <- sub_tf_target_mapping %>%
        filter(accession == accession)

      # Summarize binary features for each gene
      sub_tf_target_mapping_acc <- sub_tf_target_mapping_acc %>%
        group_by(ensembl_gene_id) %>%
        summarise(across(c(is_S2Mb, is_M2Kb, is_S2Kb, is_M100Kb, is_S100Kb, is_coexpressed),
                         ~ as.logical(sum(.)))) %>%
        ungroup()

      # Merge summarized data with tested genes
      sub_tested_genes <- left_join(tested_genes, sub_tf_target_mapping_acc, by = "ensembl_gene_id") %>%
        replace(is.na(.), 0)

      tmp <- c()  # Placeholder for storing model results

      # List of models to evaluate (S2Mb, M2Kb, S2Kb, M100Kb, S100Kb)
      models <- list(
        "method1" = glm(is_S2Mb ~ is_coexpressed + tpm_total, data = sub_tested_genes, family = "binomial"),
        "method2" = glm(is_M2Kb ~ is_coexpressed + tpm_total, data = sub_tested_genes, family = "binomial"),
        "method3" = glm(is_S2Kb ~ is_coexpressed + tpm_total, data = sub_tested_genes, family = "binomial"),
        "method5" = glm(is_M100Kb ~ is_coexpressed + tpm_total, data = sub_tested_genes, family = "binomial"),
        "method4" = glm(is_S100Kb ~ is_coexpressed + tpm_total, data = sub_tested_genes, family = "binomial")
      )

      # Iterate over models and calculate odds ratios, p-values, and shuffled results for is_S2Mb
      for (model_name in names(models)) {
        model <- models[[model_name]]
        
        # Get model summary
        summary_model <- summary(model)
        coeff <- summary_model$coefficients[2, ]

        # Calculate odds ratio and confidence intervals
        odds_ratio <- exp(coeff[1])
        conf_int_l <- exp(coeff[1] - (coeff[2] * 1.96))
        conf_int_u <- exp(coeff[1] + (coeff[2] * 1.96))
        pvalue <- coeff[4]

        # Append results to tmp
        tmp <- rbind(tmp, c("is_coexpressed", model_name, TF, log2(odds_ratio), log2(conf_int_l), log2(conf_int_u), pvalue))

        # Special handling for is_S2Mb: run permutation model
        if (model_name == "method1") {
          permglm_s2mb <- permmodels(model, nperm = 1000, test.statistic = "Wald", ncore = 3)

          # Extract shuffled coefficients and calculate shuffled odds ratio
          coeffs <- sapply(permglm_s2mb$permlist, function(x) { x$coef[2] })
          coeffs_quantile <- quantile(coeffs, probs = c(0.05, 0.95))

          shuffled_odds_ratio <- exp(mean(coeffs))
          shuffled_conf_int_l <- exp(coeffs_quantile[1])
          shuffled_conf_int_u <- exp(coeffs_quantile[2])
          shuffled_pvalue <- as.numeric(permglm_s2mb$COEFFICENTS[2, 5])

          # Append shuffled results for is_S2Mb
          tmp <- rbind(tmp, c("is_coexpressed", "shuffled", TF, log2(shuffled_odds_ratio), 
                              log2(shuffled_conf_int_l), log2(shuffled_conf_int_u), shuffled_pvalue))
        }
      }

      # Return the results
      tmp
    }, simplify = FALSE)

    # Combine results into a data frame
    enrich_temp <- do.call(rbind, enrich_temp)
    as.data.frame(enrich_temp, stringsAsFactors = FALSE)
  }

  # Combine results for all TFs and save to file
  enrich_scores <- rbindlist(enrich_scores)
  colnames(enrich_scores) <- c("variable", "method", "TF", "odds", "conf_int_l", "conf_int_u", "pvalue")

  # Save results to a file
  output_file <- sprintf("data/s3-network_enrichment/enrich_scores_remap_all_tfs_%s_coex.tsv", cl)
  fwrite(enrich_scores, output_file, col.names = TRUE, row.names = FALSE, sep = "\t")
}
