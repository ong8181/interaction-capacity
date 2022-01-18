####
#### CERrice2017 All data analysis
#### Fig. Network properties
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace of original data
load("../../02_EcolComAnalysis/09_NetworkPropCCMOut/09_NetworkPropCCMOut.RData")

# Load saved object of SI (surrogate data)
edna_prop_surr <- readRDS("../../02_EcolComAnalysis/SI/SI_01_RandomizedTS/SI_08_NetworkPropCompile_SurrogateOut/edna_prop.obj")

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.30
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.11
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.30
library(scales); packageVersion("scales") # 1.1.0, 2020.1.30
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.11
library(mgcv); packageVersion("mgcv") # 1.8.31, 2020.1.30
theme_set(theme_cowplot())

library(knitr); packageVersion("knitr") # 1.26, 2020.1.30
library(grid); packageVersion("grid") # 3.6.1, 2020.1.30

# Create output directory
fig_output <- "../00_RawFigs/02_Fig_EDMnet"
dir.create(fig_output)

#-------------------- Data exploration (original data) --------------------#
gt01 <- summary(gam(int_mean ~ s(total_div), data = edna_prop))
gt02 <- summary(gam(int_n ~ s(total_div), data = edna_prop))
gt03 <- summary(gam(int_per_sp ~ s(total_div), data = edna_prop))
gt04 <- summary(gam(temp_mean ~ s(total_div), data = edna_prop[edna_prop$total_div > 0,]))
gt05 <- summary(gam(connectance ~ s(total_div), data = edna_prop))
gt07 <- summary(gam(int_per_sp ~ s(temp_mean), data = edna_prop))
gt08 <- summary(gam(int_per_sp ~ s(total_dna), data = edna_prop))
gt09 <- summary(gam(connectance ~ s(temp_mean), data = edna_prop))
gt10 <- summary(gam(connectance ~ s(total_dna), data = edna_prop))
gt11 <- summary(gam(total_dna ~ s(temp_mean), data = edna_prop[edna_prop$total_dna > 0,]))
gt12 <- summary(gam(int_per_sp ~ s(connectance), data = edna_prop))
gt13 <- summary(gam(total_div ~ s(total_dna), data = edna_prop[edna_prop$total_dna > 0,]))
gt14 <- summary(gam(mean_cv ~ s(total_div), data = edna_prop))
gt15 <- summary(gam(int_max ~ s(total_div), data = edna_prop))
gt16 <- summary(gam(int_max ~ s(temp_mean), data = edna_prop)) # Non significant regression
gt17 <- summary(gam(dynamic_stab ~ s(mean_cv), data = edna_prop))
gt18 <- summary(gam(dynamic_stab ~ s(total_div), data = edna_prop)) # Non significant regression


#-------------------- Data exploration (surrogate data) --------------------#
st01 <- summary(gam(int_mean ~ s(total_div), data = edna_prop_surr))
st02 <- summary(gam(int_n ~ s(total_div), data = edna_prop_surr))
st03 <- summary(gam(int_per_sp ~ s(total_div), data = edna_prop_surr))
st04 <- summary(gam(total_div ~ s(temp_mean), data = edna_prop_surr[edna_prop_surr$total_div > 0,]))
st05 <- summary(gam(connectance ~ s(total_div), data = edna_prop_surr))
st07 <- summary(gam(int_per_sp ~ s(temp_mean), data = edna_prop_surr))
st08 <- summary(gam(int_per_sp ~ s(total_dna), data = edna_prop_surr))
st09 <- summary(gam(connectance ~ s(temp_mean), data = edna_prop_surr))
st10 <- summary(gam(connectance ~ s(total_dna), data = edna_prop_surr))
st11 <- summary(gam(total_dna ~ s(temp_mean), data = edna_prop_surr[edna_prop_surr$total_div > 0,]))
st12 <- summary(gam(int_per_sp ~ s(connectance), data = edna_prop_surr))
st13 <- summary(gam(total_div ~ s(total_dna), data = edna_prop_surr))
st14 <- summary(gam(mean_cv ~ s(total_div), data = edna_prop_surr))
st15 <- summary(gam(int_max ~ s(total_div), data = edna_prop_surr))
st16 <- summary(gam(int_max ~ s(temp_mean), data = edna_prop_surr))
st17 <- summary(gam(dynamic_stab ~ s(total_div), data = edna_prop_surr))


#-------------------- Export to tables --------------------#
add_summary_table_rows <- function(x, data_table = stat_summary_table){
  return(rbind(data_table,
               data.frame(Explained_variable = as.character(x$formula)[2],
                          Smooth_term = rownames(x$s.table),
                          Effective_df = round(x$s.table[,"edf"],1),
                          F_value = round(x$s.table[,"F"],1),
                          P_value = x$s.table[,"p-value"],
                          Adjusted_R2 = round(x$r.sq, 3))))
}

stat_summary_table <- data.frame(NULL)
stat_summary_table <- data.frame(Explained_variable = as.character(gt01$formula)[2],
                                 Smooth_term = rownames(gt01$s.table),
                                 Effective_df = round(gt01$s.table[,"edf"],1),
                                 F_value = round(gt01$s.table[,"F"],1),
                                 P_value = gt01$s.table[,"p-value"],
                                 Adjusted_R2 = round(gt01$r.sq, 3))

stat_summary_table <- add_summary_table_rows(gt02)
stat_summary_table <- add_summary_table_rows(gt03)
stat_summary_table <- add_summary_table_rows(gt04)
stat_summary_table <- add_summary_table_rows(gt05)
stat_summary_table <- add_summary_table_rows(gt07)
stat_summary_table <- add_summary_table_rows(gt08)
stat_summary_table <- add_summary_table_rows(gt09)
stat_summary_table <- add_summary_table_rows(gt10)
stat_summary_table <- add_summary_table_rows(gt11)
stat_summary_table <- add_summary_table_rows(gt12)
stat_summary_table <- add_summary_table_rows(gt13)
stat_summary_table <- add_summary_table_rows(gt14)
stat_summary_table <- add_summary_table_rows(gt15)
stat_summary_table <- add_summary_table_rows(gt16)
stat_summary_table <- add_summary_table_rows(gt17)
stat_summary_table <- add_summary_table_rows(gt18)

# Output table
Table_S1 <- as.matrix(stat_summary_table)
colnames(Table_S1) <- c("Explained variable", "Smooth term", "Effective d.f.",
                        "F value", "P value", "Adjusted R squared")
kable(Table_S1, format = "markdown")


#-------------------- Export surrogate summary table --------------------#
stat_surrogate_table <- data.frame(NULL)
stat_surrogate_table <- data.frame(Explained_variable = as.character(st01$formula)[2],
                                 Smooth_term = rownames(st01$s.table),
                                 Effective_df = round(st01$s.table[,"edf"],1),
                                 F_value = round(st01$s.table[,"F"],1),
                                 P_value = st01$s.table[,"p-value"],
                                 Adjusted_R2 = round(st01$r.sq, 3))

stat_surrogate_table <- add_summary_table_rows(st02, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st03, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st04, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st05, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st07, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st08, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st09, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st10, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st11, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st12, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st13, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st14, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st15, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st16, data_table = stat_surrogate_table)
stat_surrogate_table <- add_summary_table_rows(st17, data_table = stat_surrogate_table)

st_tables <- gridExtra::tableGrob(stat_surrogate_table)
grid.draw(st_tables)

