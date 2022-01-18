####
#### CERrice2017 All data analysis
#### SI: Surrogate Community Data: No.7 Reconstruct interaction matrix and calculate dynamic stability using regularized S-map
#### Single variable test to check whether the observed pattern is a statistical artifact
####

# Load workspace
load("../../06_TaxaAssignUpdateOut/06_TaxaAssignUpdateOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder07 <- "SI_07_RegularizedSmap_SurrogateOut"
dir.create(output_folder07)

# Load library
library(Rcpp); packageVersion("Rcpp") # 1.0.3, 2020.1.23
library(RcppArmadillo); packageVersion("RcppArmadillo") # 0.9.400.3.0, 2020.1.23
library(rEDM); packageVersion("rEDM") # 0.7.5, 2020.1.23
library(Matrix); packageVersion("Matrix") # 1.2.18, 2020.1.23
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.23
library(glmnet); packageVersion("glmnet") # 3.0.1, 2020.1.23
library(pforeach); packageVersion("pforeach") # 1.3, 2019.11.8
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.8

# Load defined functions
source("../../functions/BidirectSimplex_rEDM_0.7.4.R")
source("../../functions/BidirectBlocklnlp_rEDM_0.7.4.R")
source("../../functions/Extended_SSR.R")
source("../../functions/Extended_Smap.R")
#source("../../functions/RegularizedSmapAll.R")
source("functions/RegularizedSmapAll_Surrogate.R")
source("../../functions/DynamicStability.R")
sourceCpp("../../functions/RcppArmadillo_spEigen.cpp")
sourceCpp("../../functions/RcppArmadillo_Eigen.cpp")

# Update taxa information
edna_tax3 <- read.delim(sprintf("../../%s/SelectedASV_merge_classigntax", output_folder06))
rownames(edna_tax3) <- edna_tax3$query
all(rownames(edna_tax2) == rownames(edna_tax3))
edna_tax3$query <- edna_tax2$query
edna_tax3$seq <- edna_tax2$seq
edna_tax3$seqlen <- edna_tax2$seqlen
edna_tax3$entropy <- edna_tax2$entropy
edna_tax3$miseq_run <- edna_tax2$miseq_run

# Make surrogate causal matrix
#link_n_total <- aggregate(causal_dnaxdna2$cause_var, by = list(causal_dnaxdna2$effect_var), length)$x
set.seed(ran_seed)
n_asv <- length(edna_var_coln)
surrogate_mat_element <- sample(c(rep(0, n_asv*(n_asv-1)-nrow(causal_dnaxdna2)),
                                  rep(1, nrow(causal_dnaxdna2))), n_asv*(n_asv-1))
surrogate_mat <- matrix(0, ncol = n_asv, nrow = n_asv)
surrogate_mat[row(surrogate_mat) != col(surrogate_mat)] <- surrogate_mat_element
colnames(surrogate_mat) <- rownames(surrogate_mat) <- colnames(edna_all[,edna_var_coln])
surrogate_mat_df <- as.data.frame(surrogate_mat)
surrogate_mat_df$effect_var <- rownames(surrogate_mat_df)
causal_dnasurr <- melt(surrogate_mat_df, id.vars = c("effect_var"))
causal_dnasurr <- causal_dnasurr[causal_dnasurr$value > 0,]
colnames(causal_dnasurr) <- c("effect_var", "causal_var", "cause_assign")
causal_dnasurr <- causal_dnasurr[order(causal_dnasurr$effect_var),]
# Median = 14 causal variables per species

# Make surrogate community data
edna_surr <- edna_all[,edna_var_coln]
edna_surr[] <- NaN

for(i in 1:ncol(edna_surr)){
  xs <- edna_all[,edna_var_coln[i]]
  lib_trim <- NULL
  for(j in 1:5) lib_trim <- c(lib_trim, xs[edna_lib[j,1]:edna_lib[j,2]])
  set.seed(ran_seed) # Random shuffle surrogate
  x_surrogate0 <- make_surrogate_data(lib_trim, method = "random_shuffle", num_surr = 1)
  x_surrogate <- rep(NaN, length(xs))
  for(k in 1:5) x_surrogate[edna_lib[k,1]:edna_lib[k,2]] <- x_surrogate0[(length(lib_trim)/5 * (k-1) + 1):(length(lib_trim)/5 * k)]
  edna_surr[,i] <- x_surrogate
}

# Check diversity pattern
plot(apply(edna_surr, 1, function(x) sum(x > 0, na.rm = T)), type = "l")

# Pre-determine the best embedding dimension
if(F){
  Eedna_surr <- readRDS(sprintf("%s/Eedna_surr.obj", output_folder07))
}else{
  Eedna_surr <- c(NULL)
  for(i in 1:length(edna_var_coln)){
    Eedna_surr[i] <- bestE_bidirect(as.numeric(scale(edna_surr[,i])),
                                    lib = edna_lib,
                                    E_range = E_RANGE,
                                    show_fig = TRUE)
    cat("Cycle", i, "finished\n")
  }
  names(Eedna_surr) <- colnames(edna_surr)
  saveRDS(Eedna_surr, file = sprintf("%s/Eedna_surr.obj", output_folder07))
}

# Generation of random climate influences
causal_dnaclimsurr <- causal_dnaclim3
causal_dnaclimsurr$effect_var <- sort(sample(colnames(edna_all)[edna_var_coln], nrow(causal_dnaclimsurr)))

# No.1: Select subset at each time point
#       Select non-zero species (species currently present in the community)
taxa_list0 <- rownames(edna_tax3)
total_cycle <- nrow(edna_surr)
e_value_all <- data.frame(NULL)
compiled_smap_all <- list(NULL)

# Start main loop
for(time_i in 1:nrow(edna_all)){
  # Set time
  start_time <- proc.time()[3]
  sub_lib <- edna_lib[edna_all[time_i,"plot"],]

  # Select subset
  nonzero_taxa <- colnames(edna_surr[,taxa_list0][which(edna_surr[time_i,taxa_list0] > 0)])
  
  if(length(nonzero_taxa) > 0){
    #taxa_list <- extract_connected_vertex(nonzero_taxa)$connected_vertex
    taxa_list <- nonzero_taxa
    
    # No.2: Determine the best embedding dimension using Regularized S-map
    bestPar_all <- det_bestPar_reglSmap_surrogate(taxa_list,
                                               theta_test = c(0, 1e-3, 1e-2, 0.1, 0.5, 1, 2, 4),
                                               lambda_test = c(0, 1e-3, 1e-2, 0.1, 0.5, 1, 5, 10, 50, 100),
                                               lib = edna_lib,
                                               pred = sub_lib,
                                               embedding = "naive",
                                               alpha = 0, # Ridge S-map
                                               progress_folder = sprintf("%s_temp", output_folder07))
    
    # No.3: Perform multivariate S-map
    m_smap_all <- multi_reglSmap_all_surrogate(taxa_list,
                                               bestPar_all,
                                               lib = edna_lib,
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
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder07, output_folder07))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("SI_00_SessionInfo_Surrogate/%s_SessionInfo_%s.txt", output_folder07, substr(Sys.time(), 1, 10)))
