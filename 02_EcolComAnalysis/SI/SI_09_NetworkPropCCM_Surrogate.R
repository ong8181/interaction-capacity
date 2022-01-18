####
#### CERrice2017 All data analysis
#### No.9 Detecting causality between network properties
####

# Load workspace and objects
load("SI_08_NetworkPropCompile_SurrogateOut/SI_08_NetworkPropCompile_SurrogateOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder09 <- "SI_09_NetworkPropCCM_SurrogateOut"
dir.create(output_folder09)

# Load library
library(Rcpp); packageVersion("Rcpp") # 1.0.3, 2020.1.29
library(rEDM); packageVersion("rEDM") # 0.7.5, 2020.1.29
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.29
library(pforeach); packageVersion("pforeach") # 1.3, 2019.11.13

# Load functions
source("../../functions/BidirectSimplex_rEDM_0.7.4.R")
source("../../functions/BidirectBlocklnlp_rEDM_0.7.4.R")
source("../../functions/Extended_SSR.R")
source("../../functions/Extended_Smap.R")
source("../../functions/RegularizedSmapAll.R")

# Pre-determine best E for all edna_prop variables
Eprop <- c()
for(i in 3:ncol(edna_prop)){
  xs <- as.numeric(scale(edna_prop[,i]))
  Eprop[i-2] <- bestE_bidirect(xs, E_range = E_RANGE, lib = edna_lib)
  cat(i-2, "\n")
}

#---------- Determine the best embedding dimension for all eDNA properties v.s. eDNA properties ----------#
Eprop_prop <- matrix(NaN, nrow = ncol(edna_prop)-2, ncol = ncol(edna_prop)-2)
for(i in 3:ncol(edna_prop)){ # x = effect variable
  
  # Record the starting time
  start_time <- proc.time()[3]
  
  # Parallel computing the optimal embedding dimension
  Eprop_prop_tmp <- pforeach(j = 3:ncol(edna_prop), .c = c)({
    if(i == j){
      Eprop_prop_row_tmp <- Eprop[i-2]
    }else{
      prop_std_x <- as.numeric(scale(edna_prop[,i]))
      prop_std_y <- as.numeric(scale(edna_prop[,j]))
      Eprop_prop_row_tmp <- bestE_bidirect_block_lnlp(cbind(prop_std_x, prop_std_y),
                                                      E_range = E_RANGE,
                                                      lib = edna_lib,
                                                      criteria = "rmse",
                                                      show_fig = F)
      
    }
    names(Eprop_prop_row_tmp) <- colnames(edna_prop)[j]
    Eprop_prop_row_tmp
  })
  
  # Save the optimal E to Eedna_edna
  Eprop_prop[i-2,] <- Eprop_prop_tmp
  
  # Output message
  time_elapsed <- proc.time()[3] - start_time
  cat(sprintf("Best E_prop_prop cbind(%s,all) determined: %.02f sec\n", i, time_elapsed))
}

# CCM between stability and interaction network properties
CCM_TP <- c(-14:2)
CRITERIA <- "rmse"
N_SURR <- 1000
total_cycle <- length(3:ncol(edna_prop))^2
cycle_n <- 1

ednaprop_ccm_res <- data.frame()

#-------------------- Main loop --------------------#
for(i in 3:ncol(edna_prop)){
  for(j in 3:ncol(edna_prop)){
    # Set time
    start_time <- proc.time()[3]
    
    xs <- as.numeric(scale(edna_prop[,i]))
    ys <- as.numeric(scale(edna_prop[,j]))
    Exy <- Eprop_prop[i-2,j-2]
    Ex <- Eprop[i-2]
    Ey <- Eprop[j-2]
    
    # Refresh output object
    ednaprop_ccm_res0 <- data.frame()
    
    for(ccm_tp_i in CCM_TP){
      ccm_exact_res <- ccm_exact(cbind(xs, ys), E = Exy + 1, tp = ccm_tp_i, lib = edna_lib)
      
      # Forecasting skill and its improvement
      fore_skill <- ccm_exact_res[,c("rho", "mae", "rmse", "r2")][3,]
      fore_naive <- ccm_exact_res[,c("rho", "mae", "rmse", "r2")][1,]
      colnames(fore_naive) <- c("rho_naive", "mae_naive", "rmse_naive", "r2_naive")
      fore_improvement <- ccm_exact_res[,c("rho", "mae", "rmse", "r2")][3,] - ccm_exact_res[,c("rho", "mae", "rmse", "r2")][2,]
      colnames(fore_improvement) <- c("d_rho", "d_mae", "d_rmse", "d_r2")
      fore_skill_all <- data.frame(cbind(fore_naive, fore_skill, fore_improvement))
      fore_improvement_log <- log(ccm_exact_res[,c("mae", "rmse")][3,]) - log(ccm_exact_res[,c("mae", "rmse")][2,])
      fore_skill_all$log_d_mae <- as.numeric(fore_improvement_log[1])
      fore_skill_all$log_d_rmse <- as.numeric(fore_improvement_log[2])
      
      if(fore_skill_all$log_d_rmse < 0){
        ## Generate surrogate
        lib_trim <- NULL
        for(k in 1:5) lib_trim <- c(lib_trim, ys[edna_lib[k,1]:edna_lib[k,2]])
        na_ids <- is.na(lib_trim)
        lib_trim[is.na(lib_trim)] <- mean(lib_trim, na.rm = T) # tentative
        y_surrogate0 <- make_surrogate_data(lib_trim, method = "seasonal", num_surr = N_SURR, T_period = 122)
        y_surrogate0 <- as.data.frame(y_surrogate0)
        y_surrogate0[na_ids,] <- NaN
        y_surrogate <- as.data.frame(matrix(rep(NaN, length(ys)*N_SURR), ncol = N_SURR))
        for(k in 1:5) y_surrogate[edna_lib[k,1]:edna_lib[k,2],] <- y_surrogate0[(length(lib_trim)/5 * (k-1) + 1):(length(lib_trim)/5 * k),]
        
        # Calculate CCM for all surrogates
        y_sur_res <- ccm_exact_surrogate(xs, ys, y_surrogate, surrogate = "cause",
                                         E_fix = Exy + 1, tp = ccm_tp_i,
                                         simplex_criteria = CRITERIA,
                                         lib = edna_lib)
        
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
        y_sur_p_vals <- data.frame(matrix(c(sum(y_sur_res[,"rho"] > fore_skill_all$rho),
                                            sum(y_sur_res[,"mae"] < fore_skill_all$mae),
                                            sum(y_sur_res[,"rmse"] < fore_skill_all$rmse),
                                            sum(y_sur_res[,"r2"] > fore_skill_all$r2),
                                            sum(y_sur_res[,"d_rho"] > fore_skill_all$d_rho),
                                            sum(y_sur_res[,"d_mae"] < fore_skill_all$d_mae),
                                            sum(y_sur_res[,"d_rmse"] < fore_skill_all$d_rmse),
                                            sum(y_sur_res[,"d_r2"] > fore_skill_all$d_r2),
                                            sum(y_sur_res[,"log_d_mae"] < fore_skill_all$log_d_mae),
                                            sum(y_sur_res[,"log_d_rmse"] < fore_skill_all$log_d_rmse))/N_SURR, nrow = 1))
      }else{
        y_sur_p_vals <- data.frame(matrix(NaN, ncol = 10))
      }

      # Add colnames to y_sur_pvals
      colnames(y_sur_p_vals) <- c("pval_rho", "pval_mae", "pval_rmse", "pval_r2",
                                  "pval_d_rho", "pval_d_mae", "pval_d_rmse", "pval_d_r2",
                                  "pval_log_d_mae", "pval_log_d_rmse")
      y_sur_p_vals$xmap_name <- sprintf("%s_XMAP_%s_TP%s",
                                        colnames(edna_prop)[i],
                                        colnames(edna_prop)[j],
                                        ccm_tp_i)
      
      ednaprop_ccm_res_tmp <- data.frame(effect_var = colnames(edna_prop)[i],
                                         cause_var = colnames(edna_prop)[j],
                                         criteria = CRITERIA,
                                         E_xmap = Exy + 1,
                                         E_effect = Ex,
                                         E_cause = Ey,
                                         ccm_tp = ccm_tp_i,
                                         n_surr = N_SURR) %>% cbind(., fore_skill_all, y_sur_p_vals)
      
      # Combine CCM results
      ednaprop_ccm_res0 <- rbind(ednaprop_ccm_res0, ednaprop_ccm_res_tmp)
    }
    
    # Combine CCM results
    ednaprop_ccm_res <- rbind(ednaprop_ccm_res, ednaprop_ccm_res0)
    
    # Message
    time_used <- round(proc.time()[3] - start_time, digits = 2)
    cat("Cycle", cycle_n, "/", total_cycle, "finished;", time_used, "sec elapsed\n")
    cycle_n <- cycle_n + 1
  }
}
#-------------------- Main loop finished --------------------#

# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder09, output_folder09))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("SI_00_SessionInfo_Surrogate/%s_SessionInfo_%s.txt", output_folder09, substr(Sys.time(), 1, 10)))

