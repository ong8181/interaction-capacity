####
#### CERrice2017 All data analysis
#### No.10 Quantifying interactions between network properties
####

# Load workspace and objects
load("09_NetworkPropCCMOut/09_NetworkPropCCMOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder10 <- "10_NetworkPropSmapOut"
dir.create(output_folder10)

# Load library
library(Rcpp); packageVersion("Rcpp") # 1.0.3, 2020.1.28
library(rEDM); packageVersion("rEDM") # 0.7.5, 2020.1.28
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.28
library(glmnet); packageVersion("glmnet") # 3.0.1, 2020.1.28
library(pforeach); packageVersion("pforeach") # 1.3, 2019.11.11

# Load functions
source("functions/BidirectSimplex_rEDM_0.7.4.R")
source("functions/BidirectBlocklnlp_rEDM_0.7.4.R")
source("functions/Extended_SSR.R")
source("functions/Extended_Smap.R")

# Select causal relationships between properties
ednaprop_causal <- na.omit(ednaprop_ccm_res) %>% .[.$pval_log_d_rmse < 0.05 & .$ccm_tp < 1,]
ednaprop_causal <- ednaprop_causal[ednaprop_causal$rmse_naive > ednaprop_causal$rmse,]
ednaprop_causal[,c("effect_var", "cause_var")]
ednaprop_causal[ednaprop_causal$effect_var == "dynamic_stab",]
ednaprop_causal[ednaprop_causal$effect_var == "dynamic_stab", c("log_d_rmse", "cause_var")]
ednaprop_causal[ednaprop_causal$effect_var == "mean_cv", c("log_d_rmse", "cause_var")]
ednaprop_causal[ednaprop_causal$effect_var == "total_div", c("log_d_rmse", "cause_var")]

# Select the strongest time step
effect_var_list <- unique(ednaprop_causal$effect_var)
ednaprop_causal2 <- data.frame(NULL)

for(effect_var_i in effect_var_list){
  cause_var_list <- unique(ednaprop_causal[ednaprop_causal$effect_var == effect_var_i,"cause_var"])
  cause_var_list <- cause_var_list[cause_var_list != effect_var_i]
  
  for(cause_var_i in cause_var_list){
    ednaprop_selected <- ednaprop_causal[ednaprop_causal$effect_var == effect_var_i & ednaprop_causal$cause_var == cause_var_i,]
    ednaprop_strnog <- ednaprop_selected %>% .$log_d_rmse %>% which.min(.) %>% ednaprop_selected[.,]
    # Combine selected results
    ednaprop_causal2 <- rbind(ednaprop_causal2, ednaprop_strnog)
  }
}

ednaprop_causal2[ednaprop_causal2$effect_var == "mean_cv", c("rho", "log_d_rmse", "ccm_tp", "cause_var")]


# Perform S-map to quantify "net" influences between the network properties
theta_test <- c(0, 1e-3, 1e-2, 0.1, 0.5, 1, 2, 4)
lambda_test <- c(0, 1e-3, 1e-2, 0.1, 0.5, 1, 5, 10, 50, 100)
prop_smap_list <- list(NULL)
total_cycle <- nrow(ednaprop_causal2)
list_i <- 1

for(effect_var_i in as.character(effect_var_list)){
  # Extract subset data
  ednaprop_sub <- ednaprop_causal2[ednaprop_causal2$effect_var == effect_var_i,]
  cause_vars <- as.character(ednaprop_sub$cause_var)
  n_cause <- nrow(ednaprop_sub)
  lags_prop <- abs(ednaprop_sub$ccm_tp)
  
  # Use naive embedding (same with the interaction network reconstruction)
  for(cause_n in 1:length(cause_vars)){
    # Set time
    start_time <- proc.time()[3]
    E_prop <- unique(ednaprop_sub[ednaprop_sub$cause_var == cause_vars[cause_n],]$E_xmap) - 1
    
    # Identify lags
    block_effect <- matrix(c(NaN, as.numeric(scale(edna_prop[,effect_var_i])), rep(NaN, 14)), ncol = 1)
    
    # Combine all causal variables with lags
    block_effect <- cbind(block_effect, c(rep(NaN, lags_prop[cause_n]),
                                          as.numeric(scale(edna_prop[,cause_vars[cause_n]])),
                                          rep(NaN, 15-lags_prop[cause_n])))
    # Add colnames
    colnames(block_effect) <- c(effect_var_i, cause_vars[cause_n])
    
    # Make time-lags of the effect variable
    if(E_prop > 1){
      block_lags <- make_block(block_effect[,1], max_lag = E_prop)[,3:(E_prop+1)]
      if(!is.data.frame(block_lags)) block_lags <- as.data.frame(block_lags)
      colnames(block_lags) <- sprintf("%s_lag%s", effect_var_i, 1:(E_prop-1))
      block_smap <- cbind(block_effect, block_lags)
    }else{
      block_smap <- block_effect
    }
  
    # Perform multivariate S-map
    best_par_block <- pforeach(k = 1:length(theta_test), .c=rbind)({
      # Prepare output
      best_par_check <- data.frame(theta = NA, lambda = NA, rmse = NA)
      theta_i <- theta_test[k]
      
      for(l in 1:length(lambda_test)){
        # Perform the regularized S-map
        lambda_i <- lambda_test[l]
        best_check_tmp <- try(extended_lnlp(block_smap, lib = edna_lib + 1, pred = edna_lib + 1,
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
  
    if(all(is.na(best_par_block[,3]))){
      best_par <- data.frame(theta = NA, lambda = NA)
    }else{
      best_par <- best_par_block[which.min(best_par_block$rmse), c("theta", "lambda")]
    }
  
    regl_smap_res <- extended_lnlp(block_smap, lib = edna_lib + 1, pred = edna_lib + 1,
                                   tp = 1, target_column = 1,
                                   theta = best_par$theta, method = "s-map",
                                   regularized = TRUE,
                                   lambda = best_par$lambda,
                                   alpha = 0, # ridge S-map
                                   glmnet_parallel = FALSE, save_smap_coefficients = TRUE)
  
    # Add column names
    colnames(regl_smap_res$smap_coefficients) <- c("time", colnames(block_smap), "c_0")
  
    # Collect S-map coefficients
    prop_smap_list[list_i] <- list(regl_smap_res$smap_coefficients)
    names(prop_smap_list)[list_i] <- sprintf("%s_<=_%s", effect_var_i, cause_vars[cause_n])
  
    # Save progress messeges
    time_used <- round(proc.time()[3] - start_time, digits = 2)
    cat("Regularized S-map", list_i, "/", total_cycle, "cycle finished:", time_used, "sec elapsed\n",
        "(from", cause_vars[cause_n], "to", effect_var_i, ")\n")
    list_i <- list_i + 1
  }
}

# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder10, output_folder10))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder10, substr(Sys.time(), 1, 10)))

