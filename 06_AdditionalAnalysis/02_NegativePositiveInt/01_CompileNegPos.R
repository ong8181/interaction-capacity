####
#### CERrice2017 additional data analysis
####

# ------------------------------------------ #
# Load workspace and objects
# ------------------------------------------ #
# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file
load("../../02_EcolComAnalysis/12_TaxaSpecificISOut/12_TaxaSpecificISOut.RData")


# ------------------------------------------ #
# Load libraries
# ------------------------------------------ #
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.7.27
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.7.27
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.7.27
library(ggsci); packageVersion("ggsci") # 2.9, 2021.7.27
theme_set(theme_cowplot())


# ------------------------------------------ #
# Generate output folder
# ------------------------------------------ #
(fn <- basename(rstudioapi::getSourceEditorContext()$path))
output_folder <- paste0(stringr::str_sub(fn, end = -3), "Out"); rm(fn)
dir.create(output_folder)


# ------------------------------------------ #
# Extract positive and negative interactions
# ------------------------------------------ #
mean_pos_int <- c()
mean_neg_int <- c()
n_pos_int <- c()
n_neg_int <- c()
n_int_per_sp <- matrix(NaN, ncol = nrow(edna_all), nrow = length(edna_var_coln))
rownames(n_int_per_sp) <- colnames(edna_all[edna_var_coln])

for(i in 1:length(compiled_smap_all)){
  if(length(compiled_smap_all[[i]]) == 3) {
    ## Extract interaction matrix at time t
    int_mat_tmp <- compiled_smap_all[[i]]$J1_matrix_all[[1]]
    ## Replace the diagonal elements and 0 values with NA
    diag(int_mat_tmp) <- NA
    int_mat_tmp[int_mat_tmp == 0] <- NA
    ## Identify the positive and negative elements
    pos_id <- int_mat_tmp > 0 & !is.na(int_mat_tmp)
    neg_id <- int_mat_tmp < 0 & !is.na(int_mat_tmp)
    int_id <- pos_id | neg_id
    n_int_per_sp[rownames(int_mat_tmp),i] <- rowSums(int_id, na.rm = T)

    ## Calculate the mean positive and negative interactions
    mean_pos_int <- c(mean_pos_int, mean(int_mat_tmp[pos_id]))
    mean_neg_int <- c(mean_neg_int, mean(int_mat_tmp[neg_id]))
    n_pos_int <- c(n_pos_int, sum(pos_id))
    n_neg_int <- c(n_neg_int, sum(neg_id))
  } else {
    ## Assign NA if there is no interaction matrix
    mean_pos_int <- c(mean_pos_int, NA)
    mean_neg_int <- c(mean_neg_int, NA)
    n_pos_int <- c(n_pos_int, NA)
    n_neg_int <- c(n_neg_int, NA)
  }
}

## Add calculated variables to "edna_prop"
edna_prop$mean_pos_int <- mean_pos_int
edna_prop$mean_neg_int <- mean_neg_int
edna_prop$n_pos_int <- n_pos_int
edna_prop$n_neg_int <- n_neg_int


# ------------------------------------------ #
# Check the number of interactions
# ------------------------------------------ #
# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the data file
best_E <- readRDS("../../02_EcolComAnalysis/00_CompileAllDataOut/BestE_eDNA_ts.obj")
## Calculating the difference between best_E - n_int_per_sp
## If this is negative, the number of interactions exceed the best E
## (= over-embedding and the estimated IS could be inaccurate)
E_int <- best_E - n_int_per_sp
sum(E_int < 0, na.rm = T)/sum(!is.na(E_int))
sum(n_int_per_sp < 14, na.rm = T)/sum(!is.na(n_int_per_sp))


# ------------------------------------------ #
# Visualize patterns
# ------------------------------------------ #
(g1 <- ggplot(edna_prop, aes(x = total_div, y = mean_pos_int)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_smooth(method = "gam", color = "red3") +
    xlab("ASV diversity") +
    ylab("Mean positive interaction strength") +
    NULL)
(g2 <- ggplot(edna_prop, aes(x = total_div, y = mean_neg_int)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_smooth(method = "gam", color = "red3") +
    xlab("ASV diversity") +
    ylab("Mean negative interaction strength") +
    NULL)
(g3 <- ggplot(edna_prop, aes(x = n_pos_int, y = mean_pos_int)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_smooth(method = "gam", color = "red3") +
    xlim(0, 1250) +
    xlab("The total number of positive interactions") +
    ylab("Mean positive interaction strength") +
    NULL)
(g4 <- ggplot(edna_prop, aes(x = n_neg_int, y = mean_neg_int)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_smooth(method = "gam", color = "red3") +
    xlim(0, 1250) +
    xlab("The total number of negative interactions") +
    ylab("Mean negative interaction strength") +
    NULL)

(g5 <- ggplot(edna_prop, aes(x = total_div, y = n_pos_int * mean_pos_int / total_div)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_smooth(method = "gam", color = "red3") +
    xlab("ASV diversity") +
    ylab("Positive interaction capacity") +
    NULL)
(g6 <- ggplot(edna_prop, aes(x = total_div, y = n_neg_int * mean_neg_int / total_div)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_smooth(method = "gam", color = "red3") +
    xlab("ASV diversity") +
    ylab("Negative interaction capacity") +
    NULL)

g_all <- plot_grid(g1, g2, g3, g4, g5, g6,
                   ncol = 3, byrow = F, align = "hv",
                   labels = "AUTO")

# ------------------------------------------ #
# Save output
# ------------------------------------------ #
ggsave(file = sprintf("%s/PositiveNegative.pdf", output_folder),
       plot = g_all, width = 14, height = 8)

