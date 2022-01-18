####
#### CERrice2017 additional data analysis
####

# ------------------------------------------ #
# Load workspace and objects
# ------------------------------------------ #
# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the data file
load("../../02_EcolComAnalysis/12_TaxaSpecificISOut/12_TaxaSpecificISOut.RData")


# ------------------------------------------ #
# Load libraries
# ------------------------------------------ #
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.7.27
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.7.27
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.7.27
library(ggsci); packageVersion("ggsci") # 2.9, 2021.7.27
theme_set(theme_cowplot())
# For eigenvalue calculations
library(Rcpp); packageVersion("Rcpp") # 1.0.7, 2020.7.29
library(RcppArmadillo); packageVersion("RcppArmadillo") # 0.10.5.0.0, 2020.7.29
# For dynamic stability
source("functions/DynamicStability.R")
sourceCpp("functions/RcppArmadillo_spEigen.cpp")
sourceCpp("functions/RcppArmadillo_Eigen.cpp")


# ------------------------------------------ #
# Generate output folder
# ------------------------------------------ #
(fn <- basename(rstudioapi::getSourceEditorContext()$path))
output_folder <- paste0(stringr::str_sub(fn, end = -3), "Out"); rm(fn)
dir.create(output_folder)


# ------------------------------------------ #
# Identify dominant community
# ------------------------------------------ #
## Identify abundant taxa
taxa_names_sorted <- edna_all[,edna_var_coln] %>%
  colSums(na.rm = TRUE) %>% sort(decreasing = TRUE) %>% names()

## Calculate the dynamic stability of subset community
## (Just extract the subset members and interactions)
stab_top10 <- c(); t10_o <- taxa_names_sorted[1:10]
stab_top50 <- c(); t50_o <- taxa_names_sorted[1:50]
stab_top100 <- c(); t100_o <- taxa_names_sorted[1:100]
stab_top200 <- c(); t200_o <- taxa_names_sorted[1:200]
stab_top400 <- c(); t400_o <- taxa_names_sorted[1:400]
stab_top500 <- c(); t500_o <- taxa_names_sorted[1:500]
stab_top800 <- c(); t800_o <- taxa_names_sorted[1:800]
stab_all <- c()

for(i in 1:length(compiled_smap_all)){
  # Set time
  start_time <- proc.time()[3]
  
  if(length(compiled_smap_all[[i]]) == 3) {
    ## Extract interaction matrix at time t
    J1_matrix_all_tmp <- compiled_smap_all[[i]]$J1_matrix_all
    J2_matrix_pre_all_tmp <- compiled_smap_all[[i]]$J2_matrix_pre_all
    
    ## Extract taxa names of sub-networks
    t10 <- t10_o[t10_o %in% colnames(J1_matrix_all_tmp[[1]])]
    t50 <- t50_o[t50_o %in% colnames(J1_matrix_all_tmp[[1]])]
    t100 <- t100_o[t100_o %in% colnames(J1_matrix_all_tmp[[1]])]
    t200 <- t200_o[t200_o %in% colnames(J1_matrix_all_tmp[[1]])]
    t400 <- t400_o[t400_o %in% colnames(J1_matrix_all_tmp[[1]])]
    t500 <- t500_o[t500_o %in% colnames(J1_matrix_all_tmp[[1]])]
    t800 <- t800_o[t800_o %in% colnames(J1_matrix_all_tmp[[1]])]

    if(min(length(t10), length(t50), length(t100),
           length(t200), length(t400), length(t500), length(t800)) >= 2){
      ## Extract sub-networks
      J1_top10 <- list(J1_matrix_all_tmp[[1]][t10,t10])
      J1_top50 <- list(J1_matrix_all_tmp[[1]][t50,t50])
      J1_top100 <- list(J1_matrix_all_tmp[[1]][t100,t100])
      J1_top200 <- list(J1_matrix_all_tmp[[1]][t200,t200])
      J1_top400 <- list(J1_matrix_all_tmp[[1]][t400,t400])
      J1_top500 <- list(J1_matrix_all_tmp[[1]][t500,t500])
      J1_top800 <- list(J1_matrix_all_tmp[[1]][t800,t800])
      J2_top10 <- list(J2_matrix_pre_all_tmp[[1]][t10,])
      J2_top50 <- list(J2_matrix_pre_all_tmp[[1]][t50,])
      J2_top100 <- list(J2_matrix_pre_all_tmp[[1]][t100,])
      J2_top200 <- list(J2_matrix_pre_all_tmp[[1]][t200,])
      J2_top400 <- list(J2_matrix_pre_all_tmp[[1]][t400,])
      J2_top500 <- list(J2_matrix_pre_all_tmp[[1]][t500,])
      J2_top800 <- list(J2_matrix_pre_all_tmp[[1]][t800,])
      
      ## Calculate dynamic stability
      dstab_top10 <- dynamic_stability(J1_top10, J2_top10, n_eigen = 5)
      dstab_top50 <- dynamic_stability(J1_top50, J2_top50, n_eigen = 5)
      dstab_top100 <- dynamic_stability(J1_top100, J2_top100, n_eigen = 5)
      dstab_top200 <- dynamic_stability(J1_top200, J2_top200, n_eigen = 5)
      dstab_top400 <- dynamic_stability(J1_top400, J2_top400, n_eigen = 5)
      dstab_top500 <- dynamic_stability(J1_top500, J2_top500, n_eigen = 5)
      dstab_top800 <- dynamic_stability(J1_top800, J2_top800, n_eigen = 5)
      dstab_all <- dynamic_stability(J1_matrix_all_tmp, J2_matrix_pre_all_tmp, n_eigen = 5)
      
      # Save the absolute values of the leading eigenvalue
      stab_top10 <- c(stab_top10, as.numeric(abs(dstab_top10[2])))
      stab_top50 <- c(stab_top50, as.numeric(abs(dstab_top50[2])))
      stab_top100 <- c(stab_top100, as.numeric(abs(dstab_top100[2])))
      stab_top200 <- c(stab_top200, as.numeric(abs(dstab_top200[2])))
      stab_top400 <- c(stab_top400, as.numeric(abs(dstab_top400[2])))
      stab_top500 <- c(stab_top500, as.numeric(abs(dstab_top500[2])))
      stab_top800 <- c(stab_top800, as.numeric(abs(dstab_top800[2])))
      stab_all <- c(stab_all, as.numeric(abs(dstab_all[2])))
    } else {
      stab_top10 <- c(stab_top10, NA)
      stab_top50 <- c(stab_top50, NA)
      stab_top100 <- c(stab_top100, NA)
      stab_top200 <- c(stab_top200, NA)
      stab_top400 <- c(stab_top400, NA)
      stab_top500 <- c(stab_top500, NA)
      stab_top800 <- c(stab_top800, NA)
      stab_all <- c(stab_all, NA)
    }
  } else {
    stab_top10 <- c(stab_top10, NA)
    stab_top50 <- c(stab_top50, NA)
    stab_top100 <- c(stab_top100, NA)
    stab_top200 <- c(stab_top200, NA)
    stab_top400 <- c(stab_top400, NA)
    stab_top500 <- c(stab_top500, NA)
    stab_top800 <- c(stab_top800, NA)
    stab_all <- c(stab_all, NA)
  }
  
  # Show messeges
  time_used <- round(proc.time()[3] - start_time, digits = 2)
  cat("[Process] Time ID", i, "/", length(compiled_smap_all), "finished;", time_used, "sec elapsed\n")
}


## Add calculated variables to "edna_prop"
edna_prop$dynamic_stab_t10 <- stab_top10
edna_prop$dynamic_stab_t50 <- stab_top50
edna_prop$dynamic_stab_t100 <- stab_top100
edna_prop$dynamic_stab_t200 <- stab_top200
edna_prop$dynamic_stab_t400 <- stab_top400
edna_prop$dynamic_stab_t500 <- stab_top500
edna_prop$dynamic_stab_t800 <- stab_top800
edna_prop$dynamic_stab_all <- stab_all


# ------------------------------------------ #
# Visualize patterns
# ------------------------------------------ #
edna_prop <- readRDS(sprintf("%s/edna_prop_all.obj", output_folder))
(g1 <- ggplot(edna_prop, aes(x = dynamic_stab, y = dynamic_stab_all)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_vline(xintercept = 1, linetype = 2) +
    xlab("Dynamic stability (Full community)") +
    ylab("Dynamic stability (Full community)") +
    xlim(0,8) + ylim(0,8) +
    NULL)
(g2 <- ggplot(edna_prop, aes(x = dynamic_stab, y = dynamic_stab_t800)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_vline(xintercept = 1, linetype = 2) +
    xlab("Dynamic stability (Full community)") +
    ylab("Dynamic stability (Top 800 taxa)") +
    xlim(0,8) + ylim(0,8) +
    NULL)
(g3 <- ggplot(edna_prop, aes(x = dynamic_stab, y = dynamic_stab_t500)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_vline(xintercept = 1, linetype = 2) +
    xlab("Dynamic stability (Full community)") +
    ylab("Dynamic stability (Top 500 taxa)") +
    xlim(0,8) + ylim(0,8) +
    NULL)
(g4 <- ggplot(edna_prop, aes(x = dynamic_stab, y = dynamic_stab_t400)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_vline(xintercept = 1, linetype = 2) +
    xlab("Dynamic stability (Full community)") +
    ylab("Dynamic stability (Top 400 taxa)") +
    xlim(0,8) + ylim(0,8) +
    NULL)
(g5 <- ggplot(edna_prop, aes(x = dynamic_stab, y = dynamic_stab_t200)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_vline(xintercept = 1, linetype = 2) +
    xlab("Dynamic stability (Full community)") +
    ylab("Dynamic stability (Top 200 taxa)") +
    xlim(0,8) + ylim(0,8) +
    NULL)
(g6 <- ggplot(edna_prop, aes(x = dynamic_stab, y = dynamic_stab_t100)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_vline(xintercept = 1, linetype = 2) +
    xlab("Dynamic stability (Full community)") +
    ylab("Dynamic stability (Top 100 taxa)") +
    xlim(0,8) + ylim(0,8) +
    NULL)
(g7 <- ggplot(edna_prop, aes(x = dynamic_stab, y = dynamic_stab_t50)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_vline(xintercept = 1, linetype = 2) +
    xlab("Dynamic stability (Full community)") +
    ylab("Dynamic stability (Top 50 taxa)") +
    xlim(0,8) + ylim(0,8) +
    NULL)
(g8 <- ggplot(edna_prop, aes(x = dynamic_stab, y = dynamic_stab_t10)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_vline(xintercept = 1, linetype = 2) +
    xlab("Dynamic stability (Full community)") +
    ylab("Dynamic stability (Top 10 taxa)") +
    xlim(0,8) + ylim(0,8) +
    NULL)

g_all <- plot_grid(g2, g3, g4, g5, g6, g7, g8,
                   ncol = 3, byrow = T, align = "hv",
                   labels = "AUTO")
g_all2 <- plot_grid(g2, g4, g5, g6,
                   ncol = 2, byrow = T, align = "hv",
                   labels = "AUTO")
g_list <- list(g2, g4, g5, g6)

# ------------------------------------------ #
# Save output
# ------------------------------------------ #
ggsave(file = sprintf("%s/Subset_Stability.pdf", output_folder),
       plot = g_all, width = 14, height = 14)
ggsave(file = sprintf("%s/Subset_Stability2.pdf", output_folder),
       plot = g_all2, width = 10, height = 10)

saveRDS(edna_prop, sprintf("%s/edna_prop_all.obj", output_folder))
saveRDS(g_list, sprintf("%s/Subset_Stability2.obj", output_folder))
