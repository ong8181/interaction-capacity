####
#### CERrice2017 All data analysis
#### No. 4 Compile surrogate results
####

# Load workspace
load("03_CCMSurrogateOut/03_CCMSurrogateOut_2_Clim.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder04 <- "04_CompileSurResOut"
dir.create(output_folder04)

# Load CCM surrogate results
ednacom_sur_res1_1 <- readRDS(sprintf("%s/ednacom_sur_res1_1.obj", output_folder03))
ednacom_sur_res1_2 <- readRDS(sprintf("%s/ednacom_sur_res1_2.obj", output_folder03))
ednacom_sur_res1_3 <- readRDS(sprintf("%s/ednacom_sur_res1_3.obj", output_folder03))

# Combine surrogate results of eDNA-eDNA
ednacom_sur_res1 <- rbind(ednacom_sur_res1_1, ednacom_sur_res1_2, ednacom_sur_res1_3)
dim(causal_dnaxdna)
dim(ednacom_sur_res1)
causal_dnaxdna <- cbind(causal_dnaxdna, ednacom_sur_res1)

# Delete redundant objects
rm(ednacom_sur_res1)
rm(ednacom_sur_res1_1)
rm(ednacom_sur_res1_2)
rm(ednacom_sur_res1_3)
rm(ednacom_sur_res2)
rm(ednacom_sur_res3)
rm(ednacom_sur_res)

# Delete temporal objects
rm(causal_data)
rm(causal_all)
rm(effect_all)
rm(causal_lib)
rm(y_surrogate)
rm(y_surrogate0)
rm(y_sur_res)
rm(y_sur_stats)
rm(y_sur_p_vals)

# Delete unused objects
rm(edna_plot1)
rm(edna_plot2)
rm(edna_plot3)
rm(edna_plot4)
rm(edna_plot5)

# Extract "significant" causality
# Criteria = "p_val_log_d_rmse" & "rmse_naive" (but naive is already checked)
# Check the number of causal relationships
sum(causal_dnaxdna[,"pval_log_d_rmse"] < 0.005) # expected # = 7164, but 19537 detected
sum(causal_dnaclim[,"pval_log_d_rmse"] < 0.005) # expected # = 454, but 1133 detected
sum(causal_climdna[,"pval_log_d_rmse"] < 0.005) # expected # = 454, detected only 93

# Exclude self loop (for causal_dnaxdna) and extract significant causality
causal_dnaxdna2 <- causal_dnaxdna[causal_dnaxdna$effect_var != causal_dnaxdna$cause_var & causal_dnaxdna[,"pval_log_d_rmse"] < 0.005,]
causal_dnaclim2 <- causal_dnaclim[causal_dnaclim[,"pval_log_d_rmse"] < 0.005,]
causal_climdna2 <- causal_climdna[causal_climdna[,"pval_log_d_rmse"] < 0.005,]

# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder04, output_folder04))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder04, substr(Sys.time(), 1, 10)))
