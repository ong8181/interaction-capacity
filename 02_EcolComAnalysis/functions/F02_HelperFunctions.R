####
#### Helper functions for CERrice2017 time series
#### Functions for empirical dynamic modeling
####

## Check best embedding dimension & best nonlinear parameter
# Embedding dimension
bestE <- function(ts, E = 1:30, tau = 1, tp = 1, lib = lib,
                  num_neibours = "e+1",
                  criteria = "mae",
                  show_fig = F,
                  save_raw_data = F){
  # Do simplex projection
  simp_res <- simplex(ts, E = E, tau = tau, tp = tp, lib = lib, num_neighbors = num_neibours, silent = T)
  
  if(criteria == "rho"){
    best_index <- which.max(simp_res[,criteria])
  }else{
    best_index <- which.min(simp_res[,criteria])
  }

  if(show_fig){
    plot(simp_res$E, simp_res[,criteria], type = "b",
         xlab = "E", ylab = "Forecasting skill")
    points(simp_res[best_index, "E"], simp_res[best_index, criteria], pch = 16)
  }
  
  if(!save_raw_data){
    return(simp_res[best_index, "E"])
  }else{
    return(list(best_E = simp_res[best_index, "E"],
                all_res = simp_res))
  }
}


# Nonlinear parameter
bestT <- function(ts, E = 1:20, tau = 1, tp = 1, lib = lib,
                  theta = c(0, 1e-04, 3e-04, 0.001,
                            0.003, 0.01, 0.03, 0.1,
                            0.3, 0.5, 0.75, 1, 1.5,
                            2, 3, 4, 6, 8),
                  criteria = "mae",
                  show_fig = F,
                  save_raw_data = F){
  # Do simplex projection
  smap_res <- s_map(ts, E = E, theta = theta, tau = tau, tp = tp, lib = lib, silent = T)
  
  if(criteria == "rho"){
    best_index <- which.max(smap_res[,criteria])
  }else{
    best_index <- which.min(smap_res[,criteria])
  }
  
  if(show_fig){
    plot(smap_res$theta+1e-05, smap_res[,criteria], type = "b",
         xlab = "Theta", ylab = "Forecasting skill", log = "x")
    points(smap_res[best_index, "theta"]+1e-05, smap_res[best_index, criteria], pch = 16)
  }
  
  if(!save_raw_data){
    return(smap_res[best_index, "theta"])
  }else{
    return(list(bestT = smap_res[best_index, "theta"],
                all_res = smap_res))
  }
}

## Extract strong causal variables
extract_causal_vars <- function(ccm_res,
                                causal_index = "log_d_rmse",
                                ccm_tp = -1,
                                select_tp_mode = "specific_tp", # or "lower_tp" or "higher_tp"
                                select_strongest = T, # for "specific_tp" mode, this option can be turned off
                                causal_threshold = 0,
                                check_naive = TRUE){
  require(tidyverse)
  
  # Extract causal variable names
  if(select_tp_mode == "specific_tp"){
    ccm_res_ntp <- ccm_res[ccm_res$ccm_tp == ccm_tp, ]
  }else if(select_tp_mode == "lower_tp"){
    ccm_res_ntp <- ccm_res[ccm_res$ccm_tp <= ccm_tp, ]
  }else if(select_tp_mode == "higher_tp"){
    ccm_res_ntp <- ccm_res[ccm_res$ccm_tp >= ccm_tp, ]
  }else{
    stop("select_tp_mode should be either of specific_tp, lower_tp, or higher_tp")
  }
  
  if(check_naive){
    if(!is.na(match(causal_index, c("d_rmse", "log_d_rmse")))){
      ccm_res_ntp <- ccm_res_ntp[ccm_res_ntp[,causal_index] < 0,]
      ccm_res_ntp <- ccm_res_ntp[ccm_res_ntp$rmse_naive > ccm_res_ntp$rmse,]
    }else if(!is.na(match(causal_index, c("d_mae", "log_d_mae")))){
      ccm_res_ntp <- ccm_res_ntp[ccm_res_ntp[,causal_index] < 0,]
      ccm_res_ntp <- ccm_res_ntp[ccm_res_ntp$mae_naive > ccm_res_ntp$mae,]
    }else{
      stop("Either of log_d_rmse, log_d_mae, d_rmse, or d_mae is recommended as causal_index!")
    }
  }else{
    if(!is.na(match(causal_index, c("d_rmse", "log_d_rmse", "d_mae", "log_d_mae")))){
      ccm_res_ntp <- ccm_res_ntp[ccm_res_ntp[,causal_index] < 0,]
    }else{
      stop("Either of log_d_rmse, log_d_mae, d_rmse, or d_mae is recommended as causal_index!")
    }
  }  
  
  ccm_res_causal_pairs_all <- ccm_res_ntp[ccm_res_ntp[,causal_index] <= causal_threshold, c("cause_var", "effect_var")]
  ccm_res_causal_pairs <- unique(sprintf("%s --> %s", ccm_res_causal_pairs_all[,1],  ccm_res_causal_pairs_all[,2]))
  causal_pairs <- strsplit(ccm_res_causal_pairs, " --> ")
  
  # Extract causal variables
  causal_vars <- data.frame(NULL)
  
  if(select_strongest){
    for(cause_i in causal_pairs){
      cond1 <- ccm_res_ntp[,causal_index] <= causal_threshold
      cond2 <- ccm_res_ntp$cause_var == cause_i[1]
      cond3 <- ccm_res_ntp$effect_var == cause_i[2]
      ccm_res_sub <- ccm_res_ntp[cond1 & cond2 & cond3,]
      ccm_max_cause <- ccm_res_sub[which.min(ccm_res_sub[,causal_index]),]
      causal_vars <- bind_rows(causal_vars, ccm_max_cause)
    }
  }else{
    causal_vars <- ccm_res_ntp
  }
  
  return(causal_vars)
}

