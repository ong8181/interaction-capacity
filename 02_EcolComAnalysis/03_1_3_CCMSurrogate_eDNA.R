####
#### CERrice2017 All data analysis
#### No. 3 Perform suggogate analysis for potentially causal variables
####

# Load workspace
load("02_CompileCCMresOut/02_CompileCCMresOut_Minimum.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder03 <- "03_CCMSurrogateOut"
dir.create(output_folder03)

# Load library and functions
library(rEDM); packageVersion("rEDM") # 0.7.5, 2020.1.7
library(pforeach); packageVersion("pforeach") # 1.3, 2019.10.28
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.7
source("functions/CcmExact_v5_rEDM_0.7.4.R")
source("functions/BidirectSimplex_rEDM_0.7.4.R")
source("functions/BidirectBlocklnlp_rEDM_0.7.4.R")
source("functions/F02_HelperFunctions.R")
source("functions/Parallel_Sur_95CI.R")

# Set parameters  
total_cycle <- floor(nrow(causal_dnaxdna)/3)
cycle_set <- (floor(nrow(causal_dnaxdna)/3)*2+1):nrow(causal_dnaxdna)
N_SURR <- 1000
E_RANGE <- 1:14
cycle_n <- 1

# Set output objects
ednacom_sur_res1 <- data.frame()
#ednacom_sur_res2 <- data.frame()
#ednacom_sur_res3 <- data.frame()

# Main loop
for(data_set_id in 1){
  # Set temporal objects
  if(data_set_id == 1){
    causal_data <- causal_dnaxdna
    effect_all <- edna_all; causal_all <- edna_all; causal_lib <- edna_lib
  }else if(data_set_id == 2){
    causal_data <- causal_climdna
    effect_all <- clim_all; causal_all <- edna_all; causal_lib <- clim_lib
  }else if(data_set_id == 3){
    causal_data <- causal_dnaclim
    effect_all <- edna_all; causal_all <- clim_all; causal_lib <- edna_lib
  }
  
  # Set temporal output objects
  ednacom_sur_res <- data.frame()
  
  # Start calculation of p-values for each dataset
  for(i in cycle_set){
    start_time <- proc.time()[3]
    
    # Extract time series
    xs <- as.numeric(scale(effect_all[,causal_data$effect_var[i]]))
    ys <- as.numeric(scale(causal_all[,causal_data$cause_var[i]]))

    ## Generate surrogate
    lib_trim <- NULL
    for(j in 1:5) lib_trim <- c(lib_trim, ys[causal_lib[j,1]:causal_lib[j,2]])
    y_surrogate0 <- make_surrogate_data(lib_trim, method = "seasonal", num_surr = N_SURR, T_period = 122)
    y_surrogate0 <- as.data.frame(y_surrogate0)
    y_surrogate <- as.data.frame(matrix(rep(NaN, length(ys)*N_SURR), ncol = N_SURR))
    for(k in 1:5) y_surrogate[causal_lib[k,1]:causal_lib[k,2],] <- y_surrogate0[(length(lib_trim)/5 * (k-1) + 1):(length(lib_trim)/5 * k),]
    
    # Calculate CCM for all surrogates
    y_sur_res <- ccm_exact_surrogate(xs, ys, y_surrogate, surrogate = "cause",
                                     E_fix = causal_data$E_xmap[i],
                                     tp = causal_data$ccm_tp[i],
                                     simplex_criteria = "rmse",
                                     lib = causal_lib)
    
    # Calcuate additional statistics
    nonzero_mae_id <- (y_sur_res$mae > 0) & ((y_sur_res$mae - y_sur_res$d_mae) > 0)
    nonzero_rmse_id <- (y_sur_res$rmse > 0) & ((y_sur_res$rmse - y_sur_res$d_rmse) > 0)
    if(all(nonzero_rmse_id) & all(nonzero_mae_id)){
      y_sur_res$log_d_mae <- log(y_sur_res$mae/(y_sur_res$mae - y_sur_res$d_mae))
      y_sur_res$log_d_rmse <- log(y_sur_res$rmse/(y_sur_res$rmse - y_sur_res$d_rmse))
    }else{
      y_sur_res$log_d_mae <- y_sur_res$log_d_rmse <- NaN
      y_sur_res$log_d_mae[nonzero_mae_id] <- log(y_sur_res$maee[nonzero_mae_id]/(y_sur_res$maee[nonzero_mae_id] - y_sur_res$d_maee[nonzero_mae_id]))
      y_sur_res$log_d_rmse[nonzero_rmse_id] <- log(y_sur_res$rmse[nonzero_rmse_id]/(y_sur_res$rmse[nonzero_rmse_id] - y_sur_res$d_rmse[nonzero_rmse_id]))
    }
    
    # Extracting statistical indices
    y_sur_stats <- list(summarize_sur(y_sur_res))
    y_sur_p_vals <- data.frame(matrix(c(sum(y_sur_res[,"rho"] > causal_data[i,"rho"]),
                                        sum(y_sur_res[,"mae"] < causal_data[i,"mae"]),
                                        sum(y_sur_res[,"rmse"] < causal_data[i,"rmse"]),
                                        sum(y_sur_res[,"r2"] > causal_data[i,"r2"]),
                                        sum(y_sur_res[,"d_rho"] > causal_data[i,"d_rho"]),
                                        sum(y_sur_res[,"d_mae"] < causal_data[i,"d_mae"]),
                                        sum(y_sur_res[,"d_rmse"] < causal_data[i,"d_rmse"]),
                                        sum(y_sur_res[,"d_r2"] > causal_data[i,"d_r2"]),
                                        sum(y_sur_res[,"log_d_mae"] < causal_data[i,"log_d_mae"]),
                                        sum(y_sur_res[,"log_d_rmse"] < causal_data[i,"log_d_rmse"]))/N_SURR, nrow = 1))
    
    # Add colnames to y_sur_pvals
    colnames(y_sur_p_vals) <- c("pval_rho", "pval_mae", "pval_rmse", "pval_r2",
                                "pval_d_rho", "pval_d_mae", "pval_d_rmse", "pval_d_r2",
                                "pval_log_d_mae", "pval_log_d_rmse")
    y_sur_p_vals$xmap_name <- sprintf("%s_XMAP_%s_TP%s",
                                      causal_data$effect_var[i],
                                      causal_data$cause_var[i],
                                      causal_data$ccm_tp[i])
    
    # Bind results
    ednacom_sur_res <- bind_rows(ednacom_sur_res, y_sur_p_vals)
    
    # Message
    time_used <- round(proc.time()[3] - start_time, digits = 2)
    cat("Dateset", data_set_id,
        ": Cycle", cycle_n, "/", total_cycle, "finished;", time_used, "sec elapsed\n")
    cycle_n <- cycle_n + 1
    
    # Save temporal object at either of cycle seq(50000,total_cycle,50000)
    if(any(cycle_n == seq(30000,total_cycle,30000))) saveRDS(ednacom_sur_res, sprintf("%s/tmp%s.obj", output_folder03, i))
  }
  
  # Bind results
  if(data_set_id == 1){
    ednacom_sur_res1 <- ednacom_sur_res
    #saveRDS(ednacom_sur_res1, sprintf("%s/ednacom_sur_res1_1.obj", output_folder03))
    saveRDS(ednacom_sur_res1, sprintf("%s/ednacom_sur_res1_3.obj", output_folder03))
  }else if(data_set_id == 2){
    ednacom_sur_res2 <- ednacom_sur_res
    saveRDS(ednacom_sur_res2, sprintf("%s/ednacom_sur_res2.obj", output_folder03))
  }else if(data_set_id == 3){
    ednacom_sur_res3 <- ednacom_sur_res
    saveRDS(ednacom_sur_res3, sprintf("%s/ednacom_sur_res3.obj", output_folder03))
  }
  
}


# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s_1_3_eDNA.RData", output_folder03, output_folder03))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_1_3_eDNA_SessionInfo_%s.txt", output_folder03, substr(Sys.time(), 1, 10)))
