#!/usr/bin/Rscript

library(limma)
library(dplyr)
library(tidyr)
library(data.table)
library(foreach)
library(WGCNA)
library(stringr)

here::i_am("README.md")

# Set working directory
setwd(here())

# Define a function to remove version from IDs
remove_version <- function(x) {
        unlist(str_split(x, pattern="[.]"))[1]
}

meta = fread("data/coexpression/metadata.tsv")
meta %>% filter((`Output type` == "gene quantifications") & (`File assembly` == "GRCh38") & (`Genome annotation` == "V29") & (`File analysis status` == "released")) -> temp

RNA = data.table()
foreach (acc = temp$`File accession`) %do% {
	RNA_temp = fread(paste0("data/coexpression/", acc, ".tsv"))
	RNA_temp$Sample = acc
	RNA = rbind(RNA, RNA_temp)
}

pivot_wider(RNA, id_cols = "gene_id", names_from = "Sample", values_from = "TPM") -> RNA_exp


RNA_exp %>% rowwise() %>% filter(mean(c_across(colnames(.)[2:ncol(.)]) != 0) >= 0.8) -> RNA_temp

gene_names = RNA_temp$gene_id

rownames(RNA_temp) = gene_names
RNA_temp %>% select(-c(gene_id)) -> RNA_temp

#transpose matrix to correlate genes in the following
WGCNA_matrix = t(RNA_temp)

s = abs(bicor(WGCNA_matrix))

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)

png("plots/power_plot.png")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')

dev.off()

# Create a connection to the Ensembl database
ensembl <- biomaRt::useEnsembl(
        biomart = "ensembl",
        dataset = "hsapiens_gene_ensembl",
        #host = "https://feb2023.archive.ensembl.org", #release 109
        #host = "https://apr2020.archive.ensembl.org", #release 100
        host = "https://jul2023.archive.ensembl.org", #release 110
        mirror = "uswest"
)


gene_names <- sapply(gene_names, remove_version)

gene_names = data.table(gene_names)
colnames(gene_names) = c("ensembl_gene_id")

mapping <- as.data.table(biomaRt::getBM(
               attributes = c("external_gene_name", "ensembl_gene_id"),
	       filter = c("ensembl_gene_id"),
               values = gene_names$ensembl_gene_id,
	       mart = ensembl
       ))

head(gene_names)

final_mapping <- gene_names %>%
  left_join(mapping, by = "ensembl_gene_id")

final_mapping %>% rowwise() %>% mutate(external_gene_name = ifelse(external_gene_name == "", ensembl_gene_id, external_gene_name)) -> final_mapping

nrow(s)
nrow(final_mapping)

colnames(s) = final_mapping$external_gene_name 
rownames(s) = final_mapping$external_gene_name

## Write corr tables
fwrite(s, "data/coexpression/Bicor_ENCODE_bulk_K562.tsv", sep="\t", col.names=T, row.names=F, quote=F)
