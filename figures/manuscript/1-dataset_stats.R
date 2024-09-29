library(see)
library(grid)
library(ggpubr)
library(data.table)
library(here)
library(jtools)
library(stringr)
library(biomaRt)
library(ggplot2)
library(fastglm)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(doParallel)

registerDoParallel(cores=10)

setwd(here())


# Plot style
source("figures/fig_style.R")



#GM12878



temp1 = fread("data/1-dataset_stats/per_gene_stats_S2Mb_GM12878.tsv")
temp2 = fread("data/1-dataset_stats/per_gene_stats_M2Kb_GM12878.tsv")
temp3 = fread("data/1-dataset_stats/per_gene_stats_S2Kb_GM12878.tsv")
temp4 = fread("data/1-dataset_stats/per_gene_stats_S100Kb_GM12878.tsv")
temp5 = fread("data/1-dataset_stats/per_gene_stats_M100Kb_GM12878.tsv")

temp1$method = "S2Mb"
temp2$method = "M2Kb"
temp3$method = "S2Kb"
temp4$method = "S100Kb"
temp5$method = "M100Kb"
bind_rows(temp1, temp2, temp3, temp4, temp5) -> temp

plot = ggplot(temp, aes(x = n_tfs, color = factor(method, levels = c("S2Mb", "M2Kb", "S2Kb", "S100Kb", "M100Kb")))) +
        geom_freqpoly() +
	labs(color = "Method") +
	scale_color_okabeito(palette = "black_original", order=c(6, 4, 7, 2, 5)) +
        geom_vline(xintercept = c(mean(temp1$n_tfs), mean(temp2$n_tfs), mean(temp3$n_tfs), mean(temp4$n_tfs), mean(temp5$n_tfs)),  
				linetype="dashed", color =  c("#0072B2", "#009E73", "#D55E00", "#E69F00", "#F5C710")) +
        theme_pubr(legend="none") + 
	theme(legend.position = "none", axis.title = element_text(size= 11), axis.text = element_text(size = 9)) + ylab("Genes") + xlab("Number of TFs per gene") +
        scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
	annotate("text", x = 90, y=1500, label=paste0("S2Kb: ", round(mean(temp3$n_tfs), 0))) +
	annotate("text", x = 90, y=2000, label=paste0("M2Kb: ", round(mean(temp2$n_tfs), 0))) +
        annotate("text", x = 90, y=2500, label=paste0("S2Mb: ", round(mean(temp1$n_tfs), 0))) +
	annotate("text", x = 90, y=1000, label=paste0("S100Kb: ", round(mean(temp4$n_tfs), 0))) +
	annotate("text", x = 90, y=750, label=paste0("M100Kb: ", round(mean(temp5$n_tfs), 0)))

ggsave(filename = "plots/1-dataset_stats/N_tfs_per_gene_GM12878.svg", plot = plot, device = "svg", dpi = 720, width=70, height=70, units="mm")

temp1 = fread("data/1-dataset_stats/per_tf_stats_S2Mb_GM12878.tsv")
temp2 = fread("data/1-dataset_stats/per_tf_stats_M2Kb_GM12878.tsv")
temp3 = fread("data/1-dataset_stats/per_tf_stats_S2Kb_GM12878.tsv")
temp4 = fread("data/1-dataset_stats/per_tf_stats_S100Kb_GM12878.tsv")
temp5 = fread("data/1-dataset_stats/per_tf_stats_M100Kb_GM12878.tsv")
temp1$method = "S2Mb"
temp2$method = "M2Kb"
temp3$method = "S2Kb"
temp4$method = "S100Kb"
temp5$method = "M100Kb"
bind_rows(temp1, temp2, temp3, temp4, temp5) -> temp

plot = ggplot(temp, aes(x = n_genes, color = factor(method, levels = c("S2Mb", "M2Kb", "S2Kb", "S100Kb", "M100Kb")))) +
        geom_freqpoly() +
	labs(color = "Method") +                    
        scale_color_okabeito(palette = "black_original", order=c(6, 4, 7, 2, 5)) +
	geom_vline(xintercept = c(mean(temp1$n_genes), mean(temp2$n_genes), mean(temp3$n_genes), mean(temp4$n_genes), mean(temp5$n_genes)),
                                linetype="dashed", color =  c("#0072B2", "#009E73", "#D55E00", "#E69F00", "#F5C710")) +
        theme_pubr(legend="none") + 
        scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
	theme(legend.position = "none", axis.title = element_text(size= 11), axis.text = element_text(size = 9)) + ylab("TFs") + xlab("Number of target genes per TF") +
        annotate("text", x = 24000, y=16, label=paste0("S2Kb: ", round(mean(temp3$n_genes), 0))) +
        annotate("text", x = 24000, y=19, label=paste0("M2Kb: ", round(mean(temp2$n_genes), 0))) +
        annotate("text", x = 24000, y=22, label=paste0("S2Mb: ", round(mean(temp1$n_genes), 0))) +
        annotate("text", x = 24000, y=13, label=paste0("S100Kb: ", round(mean(temp4$n_genes), 0))) +
        annotate("text", x = 24000, y=10, label=paste0("M100Kb: ", round(mean(temp5$n_genes), 0)))
ggsave(filename = "plots/1-dataset_stats/N_genes_per_tf_GM12878.svg", plot = plot, device = "svg", dpi = 720, width=70, height=70, units="mm")




#HepG2
 
temp1 = fread("data/1-dataset_stats/per_gene_stats_S2Mb_HepG2.tsv")
temp2 = fread("data/1-dataset_stats/per_gene_stats_M2Kb_HepG2.tsv")
temp3 = fread("data/1-dataset_stats/per_gene_stats_S2Kb_HepG2.tsv")
temp4 = fread("data/1-dataset_stats/per_gene_stats_S100Kb_HepG2.tsv")
temp5 = fread("data/1-dataset_stats/per_gene_stats_M100Kb_HepG2.tsv")

temp1$method = "S2Mb"
temp2$method = "M2Kb"
temp3$method = "S2Kb"
temp4$method = "S100Kb"
temp5$method = "M100Kb"
bind_rows(temp1, temp2, temp3, temp4, temp5) -> temp

plot = ggplot(temp, aes(x = n_tfs, color = factor(method, levels = c("S2Mb", "M2Kb", "S2Kb", "S100Kb", "M100Kb")))) +
        geom_freqpoly() +
        labs(color = "Method") +
        scale_color_okabeito(palette = "black_original", order=c(6, 4, 7, 2, 5)) +
        geom_vline(xintercept = c(mean(temp1$n_tfs), mean(temp2$n_tfs), mean(temp3$n_tfs), mean(temp4$n_tfs), mean(temp5$n_tfs)),
                                linetype="dashed", color =  c("#0072B2", "#009E73", "#D55E00", "#E69F00", "#F5C710")) +
        theme_pubr(legend="none") +
        theme(legend.position = "none", axis.title = element_text(size= 11), axis.text = element_text(size = 9)) + ylab("Genes") + xlab("Number of TFs per gene") +
        scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
        annotate("text", x = 180, y=1500, label=paste0("S2Kb: ", round(mean(temp3$n_tfs), 0))) +
        annotate("text", x = 180, y=2000, label=paste0("M2Kb: ", round(mean(temp2$n_tfs), 0))) +
        annotate("text", x = 180, y=2500, label=paste0("S2Mb: ", round(mean(temp1$n_tfs), 0))) +
        annotate("text", x = 180, y=1000, label=paste0("S100Kb: ", round(mean(temp4$n_tfs), 0))) +
        annotate("text", x = 180, y=750, label=paste0("M100Kb: ", round(mean(temp5$n_tfs), 0)))

ggsave(filename = "plots/1-dataset_stats/N_tfs_per_gene_HepG2.svg", plot = plot, device = "svg", dpi = 720, width=70, height=70, units="mm")

temp1 = fread("data/1-dataset_stats/per_tf_stats_S2Mb_HepG2.tsv")
temp2 = fread("data/1-dataset_stats/per_tf_stats_M2Kb_HepG2.tsv")
temp3 = fread("data/1-dataset_stats/per_tf_stats_S2Kb_HepG2.tsv")
temp4 = fread("data/1-dataset_stats/per_tf_stats_S100Kb_HepG2.tsv")
temp5 = fread("data/1-dataset_stats/per_tf_stats_M100Kb_HepG2.tsv")
temp1$method = "S2Mb"
temp2$method = "M2Kb"
temp3$method = "S2Kb"
temp4$method = "S100Kb"
temp5$method = "M100Kb"
bind_rows(temp1, temp2, temp3, temp4, temp5) -> temp

plot = ggplot(temp, aes(x = n_genes, color = factor(method, levels = c("S2Mb", "M2Kb", "S2Kb", "S100Kb", "M100Kb")))) +
        geom_freqpoly() +
        labs(color = "Method") +
        scale_color_okabeito(palette = "black_original", order=c(6, 4, 7, 2, 5)) +
        geom_vline(xintercept = c(mean(temp1$n_genes), mean(temp2$n_genes), mean(temp3$n_genes), mean(temp4$n_genes), mean(temp5$n_genes)),
                                linetype="dashed", color =  c("#0072B2", "#009E73", "#D55E00", "#E69F00", "#F5C710")) +
        theme_pubr(legend="none") +
        scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
        theme(legend.position = "none", axis.title = element_text(size= 11), axis.text = element_text(size = 9)) + ylab("TFs") + xlab("Number of target genes per TF") +
        annotate("text", x = 27000, y=29, label=paste0("S2Kb: ", round(mean(temp3$n_genes), 0))) +
        annotate("text", x = 27000, y=37, label=paste0("M2Kb: ", round(mean(temp2$n_genes), 0))) +
        annotate("text", x = 27000, y=45, label=paste0("S2Mb: ", round(mean(temp1$n_genes), 0))) +
        annotate("text", x = 27000, y=21, label=paste0("S100Kb: ", round(mean(temp4$n_genes), 0))) +
        annotate("text", x = 27000, y=13, label=paste0("M100Kb: ", round(mean(temp5$n_genes), 0)))
ggsave(filename = "plots/1-dataset_stats/N_genes_per_tf_HepG2.svg", plot = plot, device = "svg", dpi = 720, width=70, height=70, units="mm")



#K562

temp1 = fread("data/1-dataset_stats/per_gene_stats_S2Mb_K562.tsv")
temp2 = fread("data/1-dataset_stats/per_gene_stats_M2Kb_K562.tsv")
temp3 = fread("data/1-dataset_stats/per_gene_stats_S2Kb_K562.tsv")
temp4 = fread("data/1-dataset_stats/per_gene_stats_S100Kb_K562.tsv")
temp5 = fread("data/1-dataset_stats/per_gene_stats_M100Kb_K562.tsv")

temp1$method = "S2Mb"
temp2$method = "M2Kb"
temp3$method = "S2Kb"
temp4$method = "S100Kb"
temp5$method = "M100Kb"
bind_rows(temp1, temp2, temp3, temp4, temp5) -> temp

plot = ggplot(temp, aes(x = n_tfs, color = factor(method, levels = c("S2Mb", "M2Kb", "S2Kb", "S100Kb", "M100Kb")))) +
        geom_freqpoly() +
        labs(color = "Method") +
        scale_color_okabeito(palette = "black_original", order=c(6, 4, 7, 2, 5)) +
        geom_vline(xintercept = c(mean(temp1$n_tfs), mean(temp2$n_tfs), mean(temp3$n_tfs), mean(temp4$n_tfs), mean(temp5$n_tfs)),
                                linetype="dashed", color =  c("#0072B2", "#009E73", "#D55E00", "#E69F00", "#F5C710")) +
        theme_pubr(legend="none") +
        theme(legend.position = "none", axis.title = element_text(size= 11), axis.text = element_text(size = 9)) + ylab("Genes") + xlab("Number of TFs per gene") +
        scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
        annotate("text", x = 250, y=1500, label=paste0("S2Kb: ", round(mean(temp3$n_tfs), 0))) +
        annotate("text", x = 250, y=2000, label=paste0("M2Kb: ", round(mean(temp2$n_tfs), 0))) +
        annotate("text", x = 250, y=2500, label=paste0("S2Mb: ", round(mean(temp1$n_tfs), 0))) +
        annotate("text", x = 250, y=1000, label=paste0("S100Kb: ", round(mean(temp4$n_tfs), 0))) +
        annotate("text", x = 250, y=750, label=paste0("M100Kb: ", round(mean(temp5$n_tfs), 0)))

ggsave(filename = "plots/1-dataset_stats/N_tfs_per_gene_K562.svg", plot = plot, device = "svg", dpi = 720, width=70, height=70, units="mm")

temp1 = fread("data/1-dataset_stats/per_tf_stats_S2Mb_K562.tsv")
temp2 = fread("data/1-dataset_stats/per_tf_stats_M2Kb_K562.tsv")
temp3 = fread("data/1-dataset_stats/per_tf_stats_S2Kb_K562.tsv")
temp4 = fread("data/1-dataset_stats/per_tf_stats_S100Kb_K562.tsv")
temp5 = fread("data/1-dataset_stats/per_tf_stats_M100Kb_K562.tsv")
temp1$method = "S2Mb"
temp2$method = "M2Kb"
temp3$method = "S2Kb"
temp4$method = "S100Kb"
temp5$method = "M100Kb"
bind_rows(temp1, temp2, temp3, temp4, temp5) -> temp

plot = ggplot(temp, aes(x = n_genes, color = factor(method, levels = c("S2Mb", "M2Kb", "S2Kb", "S100Kb", "M100Kb")))) +
        geom_freqpoly() +
        labs(color = "Method") +
        scale_color_okabeito(palette = "black_original", order=c(6, 4, 7, 2, 5)) +
        geom_vline(xintercept = c(mean(temp1$n_genes), mean(temp2$n_genes), mean(temp3$n_genes), mean(temp4$n_genes), mean(temp5$n_genes)),
                                linetype="dashed", color =  c("#0072B2", "#009E73", "#D55E00", "#E69F00", "#F5C710")) +
        theme_pubr(legend="none") +
        scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
        theme(legend.position = "none", axis.title = element_text(size= 11), axis.text = element_text(size = 9)) + ylab("TFs") + xlab("Number of target genes per TF") +
        annotate("text", x = 35000, y=30, label=paste0("S2Kb: ", round(mean(temp3$n_genes), 0))) +
        annotate("text", x = 35000, y=40, label=paste0("M2Kb: ", round(mean(temp2$n_genes), 0))) +
        annotate("text", x = 35000, y=50, label=paste0("S2Mb: ", round(mean(temp1$n_genes), 0))) +
        annotate("text", x = 35000, y=20, label=paste0("S100Kb: ", round(mean(temp4$n_genes), 0))) +
        annotate("text", x = 35000, y=10, label=paste0("M100Kb: ", round(mean(temp5$n_genes), 0)))
ggsave(filename = "plots/1-dataset_stats/N_genes_per_tf_K562.svg", plot = plot, device = "svg", dpi = 720, width=70, height=70, units="mm")



#MCF7

temp1 = fread("data/1-dataset_stats/per_gene_stats_S2Mb_MCF7.tsv")
temp2 = fread("data/1-dataset_stats/per_gene_stats_M2Kb_MCF7.tsv")
temp3 = fread("data/1-dataset_stats/per_gene_stats_S2Kb_MCF7.tsv")
temp4 = fread("data/1-dataset_stats/per_gene_stats_S100Kb_MCF7.tsv")
temp5 = fread("data/1-dataset_stats/per_gene_stats_M100Kb_MCF7.tsv")

temp1$method = "S2Mb"
temp2$method = "M2Kb"
temp3$method = "S2Kb"
temp4$method = "S100Kb"
temp5$method = "M100Kb"
bind_rows(temp1, temp2, temp3, temp4, temp5) -> temp

plot = ggplot(temp, aes(x = n_tfs, color = factor(method, levels = c("S2Mb", "M2Kb", "S2Kb", "S100Kb", "M100Kb")))) +
        geom_freqpoly() +
        labs(color = "Method") +
        scale_color_okabeito(palette = "black_original", order=c(6, 4, 7, 2, 5)) +
        geom_vline(xintercept = c(mean(temp1$n_tfs), mean(temp2$n_tfs), mean(temp3$n_tfs), mean(temp4$n_tfs), mean(temp5$n_tfs)),
                                linetype="dashed", color =  c("#0072B2", "#009E73", "#D55E00", "#E69F00", "#F5C710")) +
        theme_pubr(legend="none") +
        theme(legend.position = "none", axis.title = element_text(size= 11), axis.text = element_text(size = 9)) + ylab("Genes") + xlab("Number of TFs per gene") +
        scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
        annotate("text", x = 120, y=1500, label=paste0("S2Kb: ", round(mean(temp3$n_tfs), 0))) +
        annotate("text", x = 120, y=2000, label=paste0("M2Kb: ", round(mean(temp2$n_tfs), 0))) +  
        annotate("text", x = 120, y=2500, label=paste0("S2Mb: ", round(mean(temp1$n_tfs), 0))) +
        annotate("text", x = 120, y=1000, label=paste0("S100Kb: ", round(mean(temp4$n_tfs), 0))) +
        annotate("text", x = 120, y=750, label=paste0("M100Kb: ", round(mean(temp5$n_tfs), 0)))

ggsave(filename = "plots/1-dataset_stats/N_tfs_per_gene_MCF7.svg", plot = plot, device = "svg", dpi = 720, width=70, height=70, units="mm")


temp1 = fread("data/1-dataset_stats/per_tf_stats_S2Mb_MCF7.tsv")
temp2 = fread("data/1-dataset_stats/per_tf_stats_M2Kb_MCF7.tsv")
temp3 = fread("data/1-dataset_stats/per_tf_stats_S2Kb_MCF7.tsv")
temp4 = fread("data/1-dataset_stats/per_tf_stats_S100Kb_MCF7.tsv")
temp5 = fread("data/1-dataset_stats/per_tf_stats_M100Kb_MCF7.tsv")
temp1$method = "S2Mb"
temp2$method = "M2Kb"
temp3$method = "S2Kb"
temp4$method = "S100Kb"
temp5$method = "M100Kb"
bind_rows(temp1, temp2, temp3, temp4, temp5) -> temp

plot = ggplot(temp, aes(x = n_genes, color = factor(method, levels = c("S2Mb", "M2Kb", "S2Kb", "S100Kb", "M100Kb")))) +
        geom_freqpoly() +
        labs(color = "Method") +
        scale_color_okabeito(palette = "black_original", order=c(6, 4, 7, 2, 5)) +
        geom_vline(xintercept = c(mean(temp1$n_genes), mean(temp2$n_genes), mean(temp3$n_genes), mean(temp4$n_genes), mean(temp5$n_genes)),
                                linetype="dashed", color =  c("#0072B2", "#009E73", "#D55E00", "#E69F00", "#F5C710")) +
        theme_pubr(legend="none") +
        scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
        theme(legend.position = "none", axis.title = element_text(size= 11), axis.text = element_text(size = 9)) + ylab("TFs") + xlab("Number of target genes per TF") +
        annotate("text", x = 35000, y=17, label=paste0("S2Kb: ", round(mean(temp3$n_genes), 0))) +
        annotate("text", x = 35000, y=20, label=paste0("M2Kb: ", round(mean(temp2$n_genes), 0))) +
        annotate("text", x = 35000, y=23, label=paste0("S2Mb: ", round(mean(temp1$n_genes), 0))) +
        annotate("text", x = 35000, y=14, label=paste0("S100Kb: ", round(mean(temp4$n_genes), 0))) +
        annotate("text", x = 35000, y=11, label=paste0("M100Kb: ", round(mean(temp5$n_genes), 0)))
ggsave(filename = "plots/1-dataset_stats/N_genes_per_tf_MCF7.svg", plot = plot, device = "svg", dpi = 720, width=70, height=70, units="mm")


plot = ggplot(temp, aes(x = n_genes, color = factor(method, levels = c("S2Mb", "M2Kb", "S2Kb", "S100Kb", "M100Kb")))) +
        geom_freqpoly() +
        labs(color = "Method") +
	theme_pubr() +
        scale_color_okabeito(palette = "black_original", order=c(6, 4, 7, 2, 5))


leg = cowplot::get_legend(plot)
svg("plots/1-dataset_stats/N_genes_tfs_legend.svg", width=1.95, height=3.90)
grid.newpage()
grid.draw(leg)
dev.off()

