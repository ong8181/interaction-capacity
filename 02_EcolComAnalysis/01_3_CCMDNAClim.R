####
#### CERrice2017 All data analysis
#### No. 1.3 CCM for eDNA-climate data
####

# Load workspace
load("00_CompileAllDataOut/00_CompileAllDataOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder01 <- "01_CCMEcolComClimOut"
output_rawccmres03 <- sprintf("%s/CCMDNAxClim_RawRes", output_folder01)
dir.create(output_rawccmres03)

# Load library and functions
library(rEDM); packageVersion("rEDM") # 0.7.5, 2020.1.7
library(pforeach); packageVersion("pforeach") # 1.3, 2019.10.25
source("functions/CcmExact_v5_rEDM_0.7.4.R")
source("functions/BidirectSimplex_rEDM_0.7.4.R")
source("functions/BidirectBlocklnlp_rEDM_0.7.4.R")
source("functions/F02_HelperFunctions.R")
source("functions/Parallel_Sur_95CI.R")

# Main loop for eDNA CCM v.s. climate variables
E_RANGE <- 1:14 # check important parameters
CCM_TP <- c(-1) # consider only tp = -1 in the dynamic stability analysis
CRITERIA <- "rmse"
total_cycle <- length(edna_var_coln) #^2 #*length(CCM_TP)
cycle_n <- 1

# Main CCM loop
for(col_i in 1:length(edna_var_coln)){
  # Refresh output object
  ednaclim_ccm_res <- data.frame()
  
  # Prepare a effect variable
  xs <- as.numeric(scale(edna_all[,edna_var_coln[col_i]]))
  xs_name <- colnames(edna_all)[edna_var_coln[col_i]]
  Ex <- Eedna[col_i]
  
  # Set time
  start_time <- proc.time()[3]
  
  ednaclim_ccm_res <- pforeach(col_j = 1:length(clim_var_coln), .c=rbind)({
    # Refresh output object
    ednaclim_ccm_res0 <- data.frame()
    
    # Potential cause variable
    ys <- as.numeric(scale(clim_all[,clim_var_coln[col_j]]))
    ys_name <- colnames(clim_all)[clim_var_coln[col_j]]
    Ey <- Eclim[col_j]
    
    # Add optimal embedding dimension by cross-mapping
    Exy <- Eedna_clim[col_i, col_j]
    
    for(ccm_tp_i in CCM_TP){
      ccm_exact_res <- ccm_exact(cbind(xs, ys), E = Exy + 1, tp = ccm_tp_i, lib = edna_lib)
      
      # Forecasting skill and its improvement
      fore_skill <- ccm_exact_res[,c("rho", "mae", "rmse", "r2")][3,]
      fore_naive <- ccm_exact_res[,c("rho", "mae", "rmse", "r2")][1,]
      colnames(fore_naive) <- c("rho_naive", "mae_naive", "rmse_naive", "r2_naive")
      fore_improvement <- ccm_exact_res[,c("rho", "mae", "rmse", "r2")][3,] - ccm_exact_res[,c("rho", "mae", "rmse", "r2")][2,]
      colnames(fore_improvement) <- c("d_rho", "d_mae", "d_rmse", "d_r2")
      fore_skill_all <- data.frame(cbind(fore_naive, fore_skill, fore_improvement))
      
      meta_info <- data.frame(effect_var = xs_name,
                              cause_var = ys_name,
                              criteria = CRITERIA,
                              E_xmap = Exy + 1,
                              E_effect = Ex,
                              E_cause = Ey,
                              ccm_tp = ccm_tp_i,
                              n_surr = NA)
      
      # Combine CCM results
      ednaclim_ccm_res_tmp <- cbind(meta_info, fore_skill_all)
      
      # Combine CCM results
      ednaclim_ccm_res0 <- rbind(ednaclim_ccm_res0, ednaclim_ccm_res_tmp)
    }
    # Combine objects
    ednaclim_ccm_res0
  })
  
  # Replace rownames
  rownames(ednaclim_ccm_res) <- 1:nrow(ednaclim_ccm_res)
  
  # Calcuate additional statistics
  nonzero_mae_id <- (ednaclim_ccm_res$mae > 0) & ((ednaclim_ccm_res$mae - ednaclim_ccm_res$d_mae) > 0)
  nonzero_rmse_id <- (ednaclim_ccm_res$rmse > 0) & ((ednaclim_ccm_res$rmse - ednaclim_ccm_res$d_rmse) > 0)
  if(all(nonzero_rmse_id) & all(nonzero_mae_id)){
    ednaclim_ccm_res$log_d_mae <- log(ednaclim_ccm_res$mae/(ednaclim_ccm_res$mae - ednaclim_ccm_res$d_mae))
    ednaclim_ccm_res$log_d_rmse <- log(ednaclim_ccm_res$rmse/(ednaclim_ccm_res$rmse - ednaclim_ccm_res$d_rmse))
  }else{
    ednaclim_ccm_res$log_d_mae <- ednaclim_ccm_res$log_d_rmse <- NaN
    ednaclim_ccm_res$log_d_mae[nonzero_mae_id] <- log(ednaclim_ccm_res$maee[nonzero_mae_id]/(ednaclim_ccm_res$maee[nonzero_mae_id] - ednaclim_ccm_res$d_maee[nonzero_mae_id]))
    ednaclim_ccm_res$log_d_rmse[nonzero_rmse_id] <- log(ednaclim_ccm_res$rmse[nonzero_rmse_id]/(ednaclim_ccm_res$rmse[nonzero_rmse_id] - ednaclim_ccm_res$d_rmse[nonzero_rmse_id]))
  }
  
  # Temporal saving of result objects
  saveRDS(ednaclim_ccm_res, sprintf("%s/%04d_DNAxClim_ccmres_%s.obj", output_rawccmres03, col_i, xs_name))
  
  # Message
  time_used <- round(proc.time()[3] - start_time, digits = 2)
  cat("Cycle", cycle_n, "/", total_cycle, "finished;", time_used, "sec elapsed\n")
  cycle_n <- cycle_n + 1
}

# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s_3_ComClim.RData", output_folder01, output_folder01))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_3_ComClim_SessionInfo_%s.txt", output_folder01, substr(Sys.time(), 1, 10)))
