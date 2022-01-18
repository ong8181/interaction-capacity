####
#### CERrice2017 All data analysis
#### No.12 Testing how diversity is determined
####

# Load workspace and objects
load("10_NetworkPropSmapOut/10_NetworkPropSmapOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder11 <- "11_DiversityNetworkOut"
dir.create(output_folder11)

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.29
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.12
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.29
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.12
library(rEDM); packageVersion("rEDM") # 0.7.5, 2020.1.29
library(pforeach); packageVersion("pforeach") # 1.3, 2019.11.12
library(glmnet); packageVersion("glmnet") # 3.0.1, 2020.1.29

# Perform regularized S-map for diversity
theta_test <- c(0, 1e-3, 1e-2, 0.1, 0.5, 1, 2, 4)
lambda_test <- c(0, 1e-3, 1e-2, 0.1, 0.5, 1, 5, 10, 50, 100)

#----- Two-variable embedding -----#
#----- Diversity can be calculated from int_per_sp, int_mean and connectance (by definition) -----#
#----- Then there will be three questions: (1) How is interaction capacity determined? -----#
#-----                                     (2) How is mean IS determined? -----#
#-----                                     (3) How is connectance determined? -----#

#-------------------- (1) How is IS per species determined? ---------------------#
# Identify potential causal variable of the diversity
#head(ednaprop_causal)
int_cause <- ednaprop_causal2[ednaprop_causal2$effect_var == "int_per_sp",]
int_cause[,c("cause_var", "E_effect", "rho", "log_d_rmse", "ccm_tp", "pval_log_d_rmse")]

# Select potential variable(s)
int_cause_var <- c("total_dna", "temp_mean")
prop_smap_id1 <- match(sprintf("int_per_sp_<=_%s", int_cause_var), names(prop_smap_list))
hist(prop_smap_list[prop_smap_id1[1]][[1]][,int_cause_var[1]]) # +/-, total DNA
hist(prop_smap_list[prop_smap_id1[2]][[1]][,int_cause_var[2]]) # +, temperature
#-----------------------------------------------------------------------#

#-------------------- (2) How is mean IS determined? ---------------------#
# Identify potential causal variable of the diversity
mean_cause <- ednaprop_causal2[ednaprop_causal2$effect_var == "int_mean",]
mean_cause[,c("cause_var", "E_effect", "rho", "log_d_rmse", "ccm_tp", "pval_log_d_rmse")]

# Select potential variable(s)
mean_cause_var <- c("total_dna", "temp_mean")
prop_smap_id2 <- match(sprintf("int_mean_<=_%s", mean_cause_var), names(prop_smap_list))
hist(prop_smap_list[prop_smap_id2[1]][[1]][,mean_cause_var[1]]) # +, total DNA
hist(prop_smap_list[prop_smap_id2[2]][[1]][,mean_cause_var[2]]) # -, temperature
#-----------------------------------------------------------------------#

#-------------------- (3) How is connectance determined? ---------------------#
# Identify potential causal variable of the diversity
con_cause <- ednaprop_causal2[ednaprop_causal2$effect_var == "connectance",]
con_cause[,c("cause_var", "E_effect", "rho", "log_d_rmse", "ccm_tp", "pval_log_d_rmse")]

# Select potential variable(s)
con_cause_var <- c("total_dna", "temp_mean")
prop_smap_id3 <- match(sprintf("connectance_<=_%s", con_cause_var), names(prop_smap_list))
hist(prop_smap_list[prop_smap_id3[1]][[1]][,con_cause_var[1]]) # +/-, total DNA
hist(prop_smap_list[prop_smap_id3[2]][[1]][,con_cause_var[2]]) # +/-, temperature
#-----------------------------------------------------------------------#

#-------------------- Lastly, how is diversity determined? ---------------------#
# Quantify the influences of temperature and total DNA on the diversity (1)
div_cause <- ednaprop_causal2[ednaprop_causal2$effect_var == "total_div",]
div_cause[,c("cause_var", "E_effect", "rho", "log_d_rmse", "ccm_tp", "pval_log_d_rmse")]

# Select potential variable(s)
div_cause_var <- c("total_dna", "temp_mean", "connectance", "int_per_sp", "int_mean")
prop_smap_id4 <- match(sprintf("total_div_<=_%s", div_cause_var), names(prop_smap_list))
hist(prop_smap_list[prop_smap_id4[1]][[1]][,div_cause_var[1]]) # total DNA -
hist(prop_smap_list[prop_smap_id4[2]][[1]][,div_cause_var[2]]) # temperature +
hist(prop_smap_list[prop_smap_id4[3]][[1]][,div_cause_var[3]]) # connectance, -
hist(prop_smap_list[prop_smap_id4[4]][[1]][,div_cause_var[4]]) # int_per_sp, -
hist(prop_smap_list[prop_smap_id4[5]][[1]][,div_cause_var[5]]) # int_mean, -
#-----------------------------------------------------------------------#


#----- Multi-variable embedding -----#
# Check details of causal factors
tmp_cause1 <- ednaprop_causal[ednaprop_causal$effect_var == "int_per_sp",]
tmp_cause1[,c("cause_var", "E_effect", "rho", "log_d_rmse", "ccm_tp", "pval_log_d_rmse")]
tmp_cause2 <- ednaprop_causal2[ednaprop_causal2$effect_var == "int_per_sp",]
tmp_cause2[,c("cause_var", "E_effect", "rho", "log_d_rmse", "ccm_tp", "pval_log_d_rmse")]

# Make new block for diversity
# Estimating best E for the multivariate & regularized S-map
block_ic <- data.frame(int_per_sp = as.numeric(scale(edna_prop$int_per_sp)),
                       total_dna = c(rep(NaN, 6), as.numeric(scale(edna_prop$total_dna))[1:(nrow(edna_prop)-6)]),
                       temp_mean = c(as.numeric(scale(edna_prop$temp_mean))[2:nrow(edna_prop)], rep(NaN, 1)))
Eic <- bestE_bidirect_block_lnlp(block_ic, E_range = E_RANGE, lib = edna_lib, show_fig = T)
block_ic2 <- cbind(block_ic, make_block(block_ic[,1], max_lag = Eic)[,3:(Eic+1)])

best_par_block <- pforeach(k = 1:length(theta_test), .c=rbind)({
  # Prepare output
  best_par_check <- data.frame(theta = NA, lambda = NA, rmse = NA)
  theta_i <- theta_test[k]
  
  for(l in 1:length(lambda_test)){
    # Perform the regularized S-map
    lambda_i <- lambda_test[l]
    best_check_tmp <- try(extended_lnlp(block_ic2, lib = edna_lib, pred = edna_lib,
                                        tp = 1, target_column = 1,
                                        theta = theta_i, method = "s-map",
                                        regularized = TRUE,
                                        lambda = lambda_i,
                                        alpha = 0, # ridge S-map
                                        glmnet_parallel = FALSE)$stats$rmse,
                          silent = TRUE)
    
    if(class(best_check_tmp) == "try-error"){
      best_par_check[l,] <- data.frame(theta = theta_i, lambda = lambda_i, rmse = NA)
    }else{
      best_par_check[l,] <- data.frame(theta = theta_i, lambda = lambda_i, rmse = best_check_tmp)
    }
  }
  # Output results
  best_par_check
})

best_par <- best_par_block[which.min(best_par_block$rmse), c("theta", "lambda")]

regl_smap_ic <- extended_lnlp(block_ic2, lib = edna_lib, pred = edna_lib,
                              tp = 1, target_column = 1,
                              theta = best_par$theta, method = "s-map",
                              regularized = TRUE,
                              lambda = best_par$lambda, alpha = 0, # ridge S-map
                              glmnet_parallel = FALSE, save_smap_coefficients = TRUE)
colnames(regl_smap_ic$smap_coefficients) <- c("time", colnames(block_ic2), "c_0")

rm(best_par_block); rm(best_par)

hist(regl_smap_ic$smap_coefficients$total_dna) # +/-
plot(block_ic2$total_dna, regl_smap_ic$smap_coefficients$total_dna)
hist(regl_smap_ic$smap_coefficients$temp_mean) # ++
plot(block_ic2$temp_mean, regl_smap_ic$smap_coefficients$temp_mean)
#-----------------------------------------------------------------------#


# Quantify the influences of temperature and total DNA on connectance
# Check causal factors and time-lag
tmp_cause3 <- ednaprop_causal[ednaprop_causal$effect_var == "connectance",]
tmp_cause3[,c("cause_var", "E_effect", "rho", "log_d_rmse", "ccm_tp", "pval_log_d_rmse")]
tmp_cause4 <- ednaprop_causal2[ednaprop_causal2$effect_var == "connectance",]
tmp_cause4[,c("cause_var", "E_effect", "rho", "log_d_rmse", "ccm_tp", "pval_log_d_rmse")]

# Make new block for diversity
# Estimating best E for the multivariate & regularized S-map
block_con <- data.frame(connectance = as.numeric(scale(edna_prop$connectance)),
                        total_dna = c(rep(NaN, 6), as.numeric(scale(edna_prop$total_dna))[1:(nrow(edna_prop)-6)]),
                        temp_mean = c(as.numeric(scale(edna_prop$temp_mean))[2:nrow(edna_prop)], rep(NaN, 1)))
Econ <- bestE_bidirect_block_lnlp(block_con, E_range = E_RANGE, lib = edna_lib, show_fig = T)
block_con2 <- cbind(block_con, make_block(block_con[,1], max_lag = Econ)[,3:(Econ+1)])

best_par_block <- pforeach(k = 1:length(theta_test), .c=rbind)({
  # Prepare output
  best_par_check <- data.frame(theta = NA, lambda = NA, rmse = NA)
  theta_i <- theta_test[k]
  
  for(l in 1:length(lambda_test)){
    # Perform the regularized S-map
    lambda_i <- lambda_test[l]
    best_check_tmp <- try(extended_lnlp(block_con2, lib = edna_lib, pred = edna_lib,
                                        tp = 1, target_column = 1,
                                        theta = theta_i, method = "s-map",
                                        regularized = TRUE,
                                        lambda = lambda_i,
                                        alpha = 0, # ridge S-map
                                        glmnet_parallel = FALSE)$stats$rmse,
                          silent = TRUE)
    
    if(class(best_check_tmp) == "try-error"){
      best_par_check[l,] <- data.frame(theta = theta_i, lambda = lambda_i, rmse = NA)
    }else{
      best_par_check[l,] <- data.frame(theta = theta_i, lambda = lambda_i, rmse = best_check_tmp)
    }
  }
  # Output results
  best_par_check
})

best_par <- best_par_block[which.min(best_par_block$rmse), c("theta", "lambda")]

regl_smap_con <- extended_lnlp(block_con2, lib = edna_lib, pred = edna_lib,
                               tp = 1, target_column = 1,
                               theta = best_par$theta, method = "s-map",
                               regularized = TRUE,
                               lambda = best_par$lambda, alpha = 0, # ridge S-map
                               glmnet_parallel = FALSE, save_smap_coefficients = TRUE)
colnames(regl_smap_con$smap_coefficients) <- c("time", colnames(block_con2), "c_0")

rm(best_par_block); rm(best_par)

hist(regl_smap_con$smap_coefficients$total_dna) # - to ++
plot(block_con2$total_dna, regl_smap_con$smap_coefficients$total_dna)
hist(regl_smap_con$smap_coefficients$temp_mean) # +/-
plot(block_con2$temp_mean, regl_smap_con$smap_coefficients$temp_mean)
#-----------------------------------------------------------------------#


# Quantify the influences of temperature and total DNA on diversity
# Check causal factors and time-lag
tmp_cause5 <- ednaprop_causal[ednaprop_causal$effect_var == "total_div",]
tmp_cause5[,c("cause_var", "E_effect", "rho", "log_d_rmse", "ccm_tp", "pval_log_d_rmse")]
tmp_cause6 <- ednaprop_causal2[ednaprop_causal2$effect_var == "total_div",]
tmp_cause6[,c("cause_var", "E_effect", "rho", "log_d_rmse", "ccm_tp", "pval_log_d_rmse")]

# Make new block for diversity
# Estimating best E for the multivariate & regularized S-map
block_div <- data.frame(total_div = as.numeric(scale(edna_prop$total_div)),
                        total_dna = c(rep(NaN, 8), as.numeric(scale(edna_prop$total_dna))[1:(nrow(edna_prop)-8)]),
                        temp_mean = c(as.numeric(scale(edna_prop$temp_mean))[1:nrow(edna_prop)]))
Ediv <- bestE_bidirect_block_lnlp(block_div, E_range = E_RANGE, lib = edna_lib, show_fig = T)
block_div2 <- cbind(block_div, make_block(block_div[,1], max_lag = Ediv)[,3:(Ediv+1)])

best_par_block <- pforeach(k = 1:length(theta_test), .c=rbind)({
  # Prepare output
  best_par_check <- data.frame(theta = NA, lambda = NA, rmse = NA)
  theta_i <- theta_test[k]
  
  for(l in 1:length(lambda_test)){
    # Perform the regularized S-map
    lambda_i <- lambda_test[l]
    best_check_tmp <- try(extended_lnlp(block_div2, lib = edna_lib, pred = edna_lib,
                                        tp = 1, target_column = 1,
                                        theta = theta_i, method = "s-map",
                                        regularized = TRUE,
                                        lambda = lambda_i,
                                        alpha = 0, # ridge S-map
                                        glmnet_parallel = FALSE)$stats$rmse,
                          silent = TRUE)
    
    if(class(best_check_tmp) == "try-error"){
      best_par_check[l,] <- data.frame(theta = theta_i, lambda = lambda_i, rmse = NA)
    }else{
      best_par_check[l,] <- data.frame(theta = theta_i, lambda = lambda_i, rmse = best_check_tmp)
    }
  }
  # Output results
  best_par_check
})

best_par <- best_par_block[which.min(best_par_block$rmse), c("theta", "lambda")]

regl_smap_div <- extended_lnlp(block_div2, lib = edna_lib, pred = edna_lib,
                               tp = 1, target_column = 1,
                               theta = best_par$theta, method = "s-map",
                               regularized = TRUE,
                               lambda = best_par$lambda, alpha = 0, # ridge S-map
                               glmnet_parallel = FALSE, save_smap_coefficients = TRUE)
colnames(regl_smap_div$smap_coefficients) <- c("time", colnames(block_div2), "c_0")

rm(best_par_block); rm(best_par)

hist(regl_smap_div$smap_coefficients$total_dna) # -
plot(block_div2$total_dna, regl_smap_div$smap_coefficients$total_dna)
hist(regl_smap_div$smap_coefficients$temp_mean) # +
plot(block_div2$temp_mean, regl_smap_div$smap_coefficients$temp_mean)
#-----------------------------------------------------------------------#

#----- Predicting total diversity using connectance, interaction capacity and mean interaction strength -----#
# Mechanistic embedding
# Make new block for diversity
block_div3 <- as.data.frame(apply(edna_prop[,c("total_div", "connectance", "int_per_sp", "int_mean")], 2, function(x) as.numeric(scale(x))))

best_par_block <- pforeach(k = 1:length(theta_test), .c=rbind)({
  # Prepare output
  best_par_check <- data.frame(theta = NA, lambda = NA, rmse = NA)
  theta_i <- theta_test[k]
  
  for(l in 1:length(lambda_test)){
    # Perform the regularized S-map
    lambda_i <- lambda_test[l]
    best_check_tmp <- try(extended_lnlp(block_div3, lib = edna_lib, pred = edna_lib,
                                        tp = 0, target_column = 1, lib_column = 2:4,
                                        theta = theta_i, method = "s-map",
                                        regularized = TRUE,
                                        lambda = lambda_i,
                                        alpha = 0, # ridge S-map
                                        glmnet_parallel = FALSE)$stats$rmse,
                          silent = TRUE)
    
    if(class(best_check_tmp) == "try-error"){
      best_par_check[l,] <- data.frame(theta = theta_i, lambda = lambda_i, rmse = NA)
    }else{
      best_par_check[l,] <- data.frame(theta = theta_i, lambda = lambda_i, rmse = best_check_tmp)
    }
  }
  # Output results
  best_par_check
})

best_par <- best_par_block[which.min(best_par_block$rmse), c("theta", "lambda")]

regl_smap_div_mechnistic <- extended_lnlp(block_div3, lib = edna_lib, pred = edna_lib,
                               tp = 0, target_column = 1, lib_column = 2:4,
                               theta = best_par$theta, method = "s-map",
                               regularized = TRUE,
                               lambda = best_par$lambda, alpha = 0, # ridge S-map
                               glmnet_parallel = FALSE, save_smap_coefficients = TRUE)
colnames(regl_smap_div_mechnistic$smap_coefficients) <- c("time", colnames(block_div3)[2:4], "c_0")
rm(best_par_block); rm(best_par)

# Check patterns
hist(regl_smap_div_mechnistic$smap_coefficients$connectance) # -
hist(regl_smap_div_mechnistic$smap_coefficients$int_per_sp) # +
hist(regl_smap_div_mechnistic$smap_coefficients$int_mean) # -
#-----------------------------------------------------------------------#

# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder11, output_folder11))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder11, substr(Sys.time(), 1, 10)))

