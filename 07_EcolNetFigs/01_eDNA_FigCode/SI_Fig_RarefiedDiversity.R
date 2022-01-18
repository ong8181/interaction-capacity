####
#### CERrice2017 All data analysis
#### Fig. Network properties
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
load("../../02_EcolComAnalysis/12_TaxaSpecificISOut/12_TaxaSpecificISOut.RData")

# Ceate output director
# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.29
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.29
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.11
library(scales); packageVersion("scales") # 1.0.0, 2019.11.12
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.11
library(GGally); packageVersion("GGally") # 1.4.0, 2019.11.12
theme_set(theme_cowplot())

# Create output directory
fig_output <- "../00_RawFigs/01_Fig_eDNAts"
dir.create(fig_output)

#-------------------- Set threshold for DNA copy number --------------------#
# No.1 Global threshold
quantile(edna_all[,edna_var_coln][edna_all[,edna_var_coln]>0], na.rm = T)
total_div_all <- rowSums(edna_all[,edna_var_coln] > 0, na.rm = T)
total_div_001 <- rowSums(edna_all[,edna_var_coln] > 1, na.rm = T)
total_div_002 <- rowSums(edna_all[,edna_var_coln] > 10, na.rm = T)
total_div_003 <- rowSums(edna_all[,edna_var_coln] > 100, na.rm = T)
total_div_004 <- rowSums(edna_all[,edna_var_coln] > 1000, na.rm = T)
div_tab <- data.frame(total_div = total_div_all,
                      div_over_0001 = total_div_001,
                      div_over_0010 = total_div_002,
                      div_over_0100 = total_div_003,
                      div_over_1000 = total_div_004)

div_pairs1 <- ggpairs(div_tab, aes_string(alpha = 0.5)) + panel_border() +
  ggtitle("Global threshold for DNA copy number")

# No.2 DNA-region-specific thresholds
# No.2.1. Threshold based on the lowest concentration of standard DNA
pro_std_thr <- c(100000,50000,25000,10000,5000)[5]
euk_std_thr <- c(10000,5000,2500,1250,250)[5]
fun_std_thr <- c(400,200,100,50,10)[5]
inv_std_thr <- c(200,100,50,25,5)[5]
pro_var_colnames <- rownames(edna_tax3[edna_tax3$miseq_run == "RMR-076",])
euk_var_colnames <- rownames(edna_tax3[edna_tax3$miseq_run == "CMR-002",])
fun_var_colnames <- rownames(edna_tax3[edna_tax3$miseq_run == "RMR-078",])
inv_var_colnames <- rownames(edna_tax3[edna_tax3$miseq_run == "RMR-099",])

adjust_factor <- 100
specific_div_001 <- rowSums(edna_all[,pro_var_colnames] > pro_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,euk_var_colnames] > euk_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,fun_var_colnames] > fun_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,inv_var_colnames] > inv_std_thr/adjust_factor, na.rm = T)

adjust_factor <- 10
specific_div_002 <- rowSums(edna_all[,pro_var_colnames] > pro_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,euk_var_colnames] > euk_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,fun_var_colnames] > fun_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,inv_var_colnames] > inv_std_thr/adjust_factor, na.rm = T)

adjust_factor <- 1
specific_div_003 <- rowSums(edna_all[,pro_var_colnames] > pro_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,euk_var_colnames] > euk_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,fun_var_colnames] > fun_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,inv_var_colnames] > inv_std_thr/adjust_factor, na.rm = T)

adjust_factor <- 0.1
specific_div_004 <- rowSums(edna_all[,pro_var_colnames] > pro_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,euk_var_colnames] > euk_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,fun_var_colnames] > fun_std_thr/adjust_factor, na.rm = T) +
  rowSums(edna_all[,inv_var_colnames] > inv_std_thr/adjust_factor, na.rm = T)

div_spec_tab <- data.frame(total_div = total_div_all,
                      div_over_0001 = specific_div_001,
                      div_over_0010 = specific_div_002,
                      div_over_0100 = specific_div_003,
                      div_over_1000 = specific_div_004)

div_pairs2 <- ggpairs(div_spec_tab, aes_string(alpha = 0.5)) + panel_border() +
  ggtitle("MiSeq run-specific threshold for DNA copy number")


div_all <- plot_grid(ggmatrix_gtable(div_pairs1),
                     ggmatrix_gtable(div_pairs2),
                     ncol = 1, align = "hv", labels = "auto")

pdf(sprintf("%s/Fig_RarefiedDiversity_All.pdf", fig_output), width = 10, height = 16)
div_all; dev.off()
pdf(sprintf("%s/Fig_RarefiedDiversity.pdf", fig_output), width = 10, height = 8)
div_pairs1; dev.off()

saveRDS(div_pairs1, sprintf("%s/Fig_RarefiedDiversity.obj", fig_output))
