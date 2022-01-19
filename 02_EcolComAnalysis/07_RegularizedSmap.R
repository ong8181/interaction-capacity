####
#### CERrice2017 All data analysis
#### No.7 Reconstruct interaction matrix and calculate dynamic stability using regularized S-map
####

# Load workspace
load("06_TaxaAssignUpdateOut/06_TaxaAssignUpdateOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder07 <- "07_RegularizedSmapOut"
dir.create(output_folder07)

# Load library
library(Rcpp); packageVersion("Rcpp") # 1.0.3, 2020.1.22
library(RcppArmadillo); packageVersion("RcppArmadillo") # 0.9.800.3.0, 2020.1.22
library(rEDM); packageVersion("rEDM") # 0.7.5, 2020.1.22
library(Matrix); packageVersion("Matrix") # 1.2.18, 2020.1.22
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.22
library(glmnet); packageVersion("glmnet") # 3.0.1, 2020.1.22
library(pforeach); packageVersion("pforeach") # 1.3, 2019.11.5

# Load defined functions
source("functions/BidirectSimplex_rEDM_0.7.4.R")
source("functions/BidirectBlocklnlp_rEDM_0.7.4.R")
source("functions/F03_HelperFunctions.R")
source("functions/Extended_SSR.R")
source("functions/Extended_Smap.R")
source("functions/RegularizedSmapAll.R")
source("functions/DynamicStability.R")
sourceCpp("functions/RcppArmadillo_spEigen.cpp")
sourceCpp("functions/RcppArmadillo_Eigen.cpp")

# Update taxa information
edna_tax3 <- read.delim(sprintf("%s/SelectedASV_merge_classigntax", output_folder06))
rownames(edna_tax3) <- edna_tax3$query
all(rownames(edna_tax2) == rownames(edna_tax3))
edna_tax3$query <- edna_tax2$query
edna_tax3$seq <- edna_tax2$seq
edna_tax3$seqlen <- edna_tax2$seqlen
edna_tax3$entropy <- edna_tax2$entropy
edna_tax3$miseq_run <- edna_tax2$miseq_run
#saveRDS(edna_tax3, "~/Desktop/edna_tax3.obj")

# No.1: Select subset at each time point
#       Select non-zero species (species currently present in the community)
taxa_list0 <- rownames(edna_tax3)
total_cycle <- nrow(edna_all)
e_value_all <- data.frame(NULL)
compiled_smap_all <- list(NULL)

# Start main loop
for(time_i in 1:nrow(edna_all)){
  # Set time
  start_time <- proc.time()[3]
  sub_lib <- edna_lib[edna_all[time_i,"plot"],]
  
  # Select subset
  nonzero_taxa <- colnames(edna_all[,taxa_list0][which(edna_all[time_i,taxa_list0] > 0)])
  
  if(length(nonzero_taxa) > 0){
    taxa_list <- nonzero_taxa # Use all non-zero taxa to calculate S-map coefficients and connectance
    
    # No.2: Determine the best embedding dimension using Regularized S-map
    bestPar_all <- det_bestPar_reglSmap(taxa_list,
                                        theta_test = c(0, 1e-3, 1e-2, 0.1, 0.5, 1, 2, 4),
                                        lambda_test = c(0, 1e-3, 1e-2, 0.1, 0.5, 1, 5, 10, 50, 100),
                                        lib = edna_lib,
                                        #pred = c(time_i-1, time_i+1),
                                        pred = sub_lib,
                                        embedding = "naive",
                                        alpha = 0, # Ridge S-map
                                        progress_folder = sprintf("%s_temp", output_folder07))
    
    # No.3: Perform multivariate S-map
    m_smap_all <- multi_reglSmap_all(taxa_list, bestPar_all,
                                     lib = edna_lib,
                                     #pred = c(time_i-1, time_i+1),
                                     pred = sub_lib,
                                     embedding = "naive",
                                     alpha = 0, # Ridge S-map
                                     progress_folder = sprintf("%s_temp", output_folder07))
    
    # Check prediction length
    prediction_success <- unlist(sapply(m_smap_all, function(x) nrow(x[[1]])))
    if(length(prediction_success) == length(taxa_list)){
      
      # No.4: Compile S-map results
      compiled_smap <- compile_reglSmap_res(taxa_list, m_smap_all, calc_lib_id = time_i)
      
      # Save compiled S-map results
      compiled_smap_all[[time_i]] <- compiled_smap
      
      # No.5: Calculate dynamic stability
      dstab_all <- dynamic_stability(compiled_smap$J1_matrix_all,
                                     compiled_smap$J2_matrix_pre_all,
                                     n_eigen = 10)
    }else{
      # Save compiled S-map results
      compiled_smap_all[[time_i]] <- NA
      
      # Predictions are not complete and cannot calculate dynamic stability
      # Prepare pseudo-results
      n_eigen <- 10
      dstab_all <- data.frame(time_index = time_i)
      for(j in 1:n_eigen) dstab_all[,j+1] <- as.complex(NaN)
      colnames(dstab_all)[2:(n_eigen+1)] <- sprintf("ev_%s", 1:n_eigen)
      cat("Cycle", time_i, ": Predictions are not complete and cannot calculate dynamic stability\n")
    }
  }else{
    # Save compiled S-map results
    compiled_smap_all[[time_i]] <- NA
    
    # Prepare pseudo-results
    n_eigen <- 10
    dstab_all <- data.frame(time_index = time_i)
    for(j in 1:n_eigen) dstab_all[,j+1] <- as.complex(NaN)
    colnames(dstab_all)[2:(n_eigen+1)] <- sprintf("ev_%s", 1:n_eigen)
  }
  
  # Save results
  dstab_all$time_index <- time_i
  e_value_all <- bind_rows(e_value_all, dstab_all)
  
  # Remove temporal files
  temp_files <- list.files(sprintf("%s_temp", output_folder07))
  file.remove(sprintf("%s/%s", sprintf("%s_temp", output_folder07), temp_files))
  file.remove(sprintf("%s_temp", output_folder07))
  
  # Show message
  time_used <- round(proc.time()[3] - start_time, digits = 2)
  cat("Quantifying local dynamic stability: Time", time_i, "/", total_cycle, "finished;", time_used, "sec elapsed\n\n")
}

# Quick check
plot(abs(e_value_all$ev_1), type = "l")
abline(h=1, lty = 2)

# Save RDS objects
saveRDS(e_value_all, sprintf("%s/dynamic_stability_nonzeroridge.obj", output_folder07))
saveRDS(compiled_smap_all, sprintf("%s/interaction_matrix_nonzeroridge.obj", output_folder07))

# Save workspace and ojcects
### Below codes are simply to keep RData file size small and not need to execute
#rm(causal_dnaclim); rm(causal_climdna); rm(causal_dnaxdna)
#rm(edna_all0)
#rm(plot1_clim_df); rm(plot1_edna_df)
#rm(plot2_clim_df); rm(plot2_edna_df)
#rm(plot3_clim_df); rm(plot3_edna_df)
#rm(plot4_clim_df); rm(plot4_edna_df)
#rm(plot5_clim_df); rm(plot5_edna_df)
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder07, output_folder07))
#save.image(sprintf("%s/%s.RData", output_folder07, output_folder07))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder07, substr(Sys.time(), 1, 10)))
