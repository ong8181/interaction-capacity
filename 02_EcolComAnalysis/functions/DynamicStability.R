####
#### CERrice2017 All data analysis
#### Compile S-map results and calculate dynamic stability
####


#-------------------- Compile regularized S-map results --------------------#
compile_reglSmap_res <- function(taxa_list,
                                 smap_res,
                                 calc_lib_id = pred_lib,
                                 calc_length = length(calc_lib_id)){
  # Extract S-map results
  # Prepare J1 matrix
  J1_matrix_all <- list(NULL)
  J1_matrix <- matrix(0, ncol = length(taxa_list), nrow = length(taxa_list))
  colnames(J1_matrix) <- rownames(J1_matrix) <- taxa_list
  
  # Prepare J2 matrix (time-delayed effects)
  J2_matrix_pre_all <- list(NULL)
  J2_matrix_pre <- matrix(0, ncol = max(E_RANGE)-1, nrow = length(taxa_list))
  rownames(J2_matrix_pre) <- taxa_list
  
  # Prepare intercept matrix
  C0_matrix_all <- list(NULL)
  C0_matrix <- matrix(0, ncol = 1, nrow = length(taxa_list))
  rownames(C0_matrix) <- taxa_list
  
  # Generate matrix for all time points
  for(i in 1:calc_length){
    J1_matrix_all[[i]] <- J1_matrix
    J2_matrix_pre_all[[i]] <- J2_matrix_pre
    C0_matrix_all[[i]] <- C0_matrix
  }
  
  # Extract and assign S-map results
  total_cycle <- calc_length
  cycle_n <- 1
  
  # Identify colnames in the interaction matrix
  # Take approx. 50 sec
  for(time_i in 1:calc_length){
    # Set time
    start_time <- proc.time()[3]
    
    for(spp_i in 1:length(taxa_list)){
      # Identify time index
      smap_coef_time <- which(smap_res[[spp_i]]$smap_coefficients$time == calc_lib_id[time_i])
      
      # Identify causal taxa
      cause_edna <- intersect(colnames(smap_res[[spp_i]]$smap_coefficients), taxa_list)
      edna_id <- match(cause_edna, taxa_list) # Identify eDNA indices
      intercept_id <- as.numeric(
        na.omit(match(c("temp_mean", "actvp_mean", "light_mean", "c_0"),
                      colnames(smap_res[[spp_i]]$smap_coefficients)))
      ) # Sum climate and c_0 values as intercept
      lag_id <- which(!is.na(sapply(strsplit(colnames(smap_res[[spp_i]]$smap_coefficients), "_"), "[", 3)))
      
      # Assign values to J1 matrix
      # rows are "effect variables" and cols are "causal variables"
      J1_matrix_all[[time_i]][spp_i, edna_id] <- as.matrix(smap_res[[spp_i]]$smap_coefficients[smap_coef_time, cause_edna], nrow = 1)
      
      # Assing values to J2 matrix (time-delayed effects)
      if(length(lag_id) > 0){
        lag_colnames <- sapply(strsplit(colnames(smap_res[[spp_i]]$smap_coefficients)[lag_id], "_"), "[", 3)
        lag_id2 <- as.numeric(substr(lag_colnames, 4, nchar(lag_colnames)))
        J2_matrix_pre_all[[time_i]][spp_i, lag_id2] <- as.matrix(smap_res[[spp_i]]$smap_coefficients[smap_coef_time, lag_id], nrow = 1)
      }
      
      # Assign values to C0 matrix
      if(length(intercept_id) > 1){
        # multiply climate values by their coefficients
        climate_id <- head(intercept_id, length(intercept_id)-1)
        climate_id_name <- colnames(smap_res[[spp_i]]$smap_coefficients)[climate_id]
        climate_val <- smap_res[[spp_i]]$smap_coefficients[smap_coef_time, climate_id] * clim_all[smap_coef_time, climate_id_name]
        C0_matrix_all[[time_i]][spp_i,1] <- sum(climate_val, smap_res[[spp_i]]$smap_coefficients[smap_coef_time, tail(intercept_id, 1)])
      }else{
        C0_matrix_all[[time_i]][spp_i,1] <- sum(smap_res[[spp_i]]$smap_coefficients[smap_coef_time, intercept_id])
      }
    }
    
    # Show messeges
    time_used <- round(proc.time()[3] - start_time, digits = 2)
    cat("Compile S-map results::Cycle", cycle_n, "/", total_cycle, "finished;", time_used, "sec elapsed\n")
    cycle_n <- cycle_n + 1
  }
  
  # Return results
  return(list(J1_matrix_all = J1_matrix_all,
              J2_matrix_pre_all = J2_matrix_pre_all,
              C0_matrix_all = C0_matrix_all))
}



#-------------------- Calculate dynamic stability --------------------#
dynamic_stability <- function(J1_matrix_all,
                              J2_matrix_pre_all,
                              n_eigen = 10,
                              edna_lib_id = pred_lib,
                              edna_data = edna_all,
                              save_eigenvectors = FALSE,
                              exclude_minor = FALSE,
                              exclude_factor = Inf){
  
  # Identify delayed effects
  delay_col <- rowSums(sapply(J2_matrix_pre_all, function(x) colSums(x, na.rm = T)), na.rm = T) != 0
  delay_nonzero_col <- sum(delay_col)
  
  # Set loop parameters
  total_cycle <- length(J1_matrix_all)
  cycle_n <- 1
  
  # Prepare output list for eigenvectors
  if(save_eigenvectors) e_vector_all <- list(NULL)
  
  # Prepare output dataframe for eigenvalues
  e_value_all <- data.frame(time_index = 1:total_cycle)
  for(i in 1:n_eigen) e_value_all <- cbind(e_value_all, as.complex(NaN))
  colnames(e_value_all)[2:(n_eigen + 1)] <- sprintf("ev_%s", 1:n_eigen)
  
  for(time_i in 1:length(J1_matrix_all)){
    # Set time
    start_time <- proc.time()[3]
    
    # Identify taxa that are absent at time = edna_lib_id[time_i]
    if(exclude_minor){
      RMR076_id <- rownames(edna_tax2[edna_tax2$miseq_run == "RMR-076",]) # lower limit of STD = 5000
      RMR078_id <- rownames(edna_tax2[edna_tax2$miseq_run == "RMR-078",]) # lower limit of STD = 10
      RMR099_id <- rownames(edna_tax2[edna_tax2$miseq_run == "RMR-099",]) # lower limit of STD = 5
      CMR002_id <- rownames(edna_tax2[edna_tax2$miseq_run == "CMR-002",]) # lower limit of STD = 250
      
      RMR076_lowlim <- 5000/exclude_factor
      RMR078_lowlim <- 10/exclude_factor
      RMR099_lowlim <- 5/exclude_factor
      CMR002_lowlim <- 250/exclude_factor
      
      taxa_included <- c(RMR076_id[!is.na(edna_all[time_i,RMR076_id]) & edna_all[time_i,RMR076_id] > RMR076_lowlim],
                         RMR078_id[!is.na(edna_all[time_i,RMR078_id]) & edna_all[time_i,RMR078_id] > RMR078_lowlim],
                         RMR099_id[!is.na(edna_all[time_i,RMR099_id]) & edna_all[time_i,RMR099_id] > RMR099_lowlim],
                         CMR002_id[!is.na(edna_all[time_i,CMR002_id]) & edna_all[time_i,CMR002_id] > CMR002_lowlim])
    }else{
      taxa_included <- colnames(edna_data[,colnames(J1_matrix_all[[time_i]])])
    }
    
    # Extraction of interaction matrix (J1)
    J1 <- J1_matrix_all[[time_i]][taxa_included, taxa_included]
    if(class(J1) != "matrix") J1 <- as.matrix(J1)
  
    # Construction of time-delay matrices (J2)
    if(length(taxa_included) == 1){
      J2 <- as.matrix(J2_matrix_pre_all[[time_i]][,1][taxa_included])
      if(delay_nonzero_col > 1){
        for(delay_i in 2:delay_nonzero_col){ # if delay matrix is all zero, it will be removed to reduce computation time
          J2 <- cbind(J2, as.matrix(J2_matrix_pre_all[[time_i]][,delay_i][taxa_included]))
        }
      }
    }else{
      J2 <- diag(J2_matrix_pre_all[[time_i]][,1][taxa_included])
      if(delay_nonzero_col > 1){
        for(delay_i in 2:delay_nonzero_col){ # if delay matrix is all zero, it will be removed to reduce computation time
          J2 <- cbind(J2, diag(J2_matrix_pre_all[[time_i]][,delay_i][taxa_included]))
        }
      }
    }
    
    # Construction of time-delay matrices (I), taxa length x (delay col length - 1)
    Unity <- diag(nrow(J1)*delay_nonzero_col)
    
    # Construction of time-delay matrices (O)
    Zero <- matrix(0, ncol = nrow(J1), nrow = nrow(J1)*delay_nonzero_col)
    
    # Combine J1, J2, I and O matrices (and convert it as a sparse matrix)
    if(delay_nonzero_col > 0){
      A <- rbind(cbind(J1, J2), cbind(Unity, Zero))
    }else{
      A <- J1
    }
    
    if(all(!is.na(A)) & dim(A)[1] > 0){
      # Eigenvalue calculation
      spA <- Matrix(A) # convert to a sparse matrix
      
      # Check sparsity
      sparsity <- sum(A!=0)/(dim(A)[1]*dim(A)[2]) < 0.5 # Calculation of the sparsity
      # n_eigen + 1 must be less than the number of rows of A
      if(sparsity){
        eigen_all <- sp_eigen_cpp(spA, min(n_eigen, dim(A)[1]-2))
      }else{
        eigen_all <- eigen_cpp(A)
        eigen_all[[1]] <- sort(eigen_all[[1]], decreasing = T)
      }
      # Save eigenvalues
      e_value <- as.complex(eigen_all[[1]])
      if(length(e_value) < n_eigen) e_value <- c(e_value, rep(as.complex(NaN), n_eigen - length(e_value)))
      e_value_all[time_i,2:(n_eigen+1)] <- e_value
      if(save_eigenvectors) e_vector_all[[time_i]] <- as.complex(eigen_all[[2]])
    }else{
      e_value_all[time_i,2:(n_eigen+1)] <- as.complex(NaN)
      if(save_eigenvectors) e_vector_all[[time_i]] <- as.complex(NaN)
    }
    
    # Show messeges
    time_used <- round(proc.time()[3] - start_time, digits = 2)
    cat("Quantifying dynamic stability: Cycle", cycle_n, "/", total_cycle, "finished;", time_used, "sec elapsed\n")
    cycle_n <- cycle_n + 1
  }
  
  if(save_eigenvectors){
    # Return results
    return(list(eigen_value = e_value_all,
                eigen_vector = e_vector_all))
  }else{
    return(e_value_all)
  }
}


