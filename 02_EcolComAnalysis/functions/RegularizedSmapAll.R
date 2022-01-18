####
#### CERrice2017 All data analysis
#### Regularrized S-map to reconstruct interaction matrix and calculate dynamic stability
####

#-------------------- Determine best theta --------------------#
det_bestPar_reglSmap <- function(taxa_list,
                                 theta_test = c(0, 1e-3, 1e-2, 0.1, 0.5, 1, 2, 4),
                                 lambda_test = c(0, 1e-3, 1e-2, 0.1, 0.5, 1, 5, 10, 50, 100),
                                 lib = edna_lib, pred = lib,
                                 embedding = "naive", # "intermediate", "limited_by_E"
                                 regularized = TRUE,
                                 alpha = 0,
                                 progress_folder = "99_progress"){
  # Create progress folder
  dir.create(progress_folder)
  
  # Set time (global)
  start_time_g <- proc.time()[3]
  
  # Extract taxa only included in the subset
  valid_id_dna <- !is.na(match(causal_dnaxdna2$effect_var, taxa_list)) & !is.na(match(causal_dnaxdna2$cause_var, taxa_list))
  valid_id_clm <- !is.na(match(causal_dnaclim3$effect_var, taxa_list))
  causal_dnaxdna_sub <- causal_dnaxdna2[valid_id_dna,]
  causal_dnaclim_sub <- causal_dnaclim3[valid_id_clm,]
  total_cycle <- length(taxa_list)

  # Most of climate xmap eDNA detected no causality
  # Also, it is unlikely that members of the ecological community drive climate processes
  # Thus, in the subsequent analyses, eDNAs (or ecological community) are assumed to have no influences on climate
  
  #for(i in 1:length(taxa_list)){
  best_par_block <- pforeach(i = 1:length(taxa_list), .c=rbind)({
    # Set time
    start_time <- proc.time()[3]
    
    # Extract information for S-map
    E_tax <- Eedna[taxa_list[i]]
    effect_tax <- taxa_list[i]
    causal_tax <- causal_dnaxdna_sub[causal_dnaxdna_sub$effect_var == effect_tax,]$cause_var
    causal_clim <- causal_dnaclim_sub[causal_dnaclim_sub$effect_var == effect_tax,]$cause_var
    
    # Calculate the number of causal variables
    n_cause <- length(causal_tax) + length(causal_clim)
    
    # Make block for S-map
    # Three options to make a block
    if(embedding == "naive"){
      # Apply "naive embedding" (use all lags and causal factors to reconstruct state space)
      if(E_tax == 1 & n_cause == 0){
        block_smap <- matrix(as.numeric(scale(edna_all[,c(effect_tax)])), ncol = 1)
      }else{
        # Include time-lag components if the number of causal variables is smaller than E
        block_smap <- apply(cbind(edna_all[,c(effect_tax, causal_tax)],
                                  clim_all[,c(causal_clim)]), 2, function(x) as.numeric(scale(x)))
        if(E_tax != 1){
          block_lags <- make_block(as.numeric(scale(edna_all[,effect_tax])), max_lag = E_tax)
          block_smap <- cbind(block_smap, block_lags[,3:(E_tax + 1)])
        }
      }
    }else if(embedding == "intermediate"){
      # Apply "intermediate embedding" (use all causal factors and some lags to reconstruct state space)
      if(E_tax == 1 & n_cause == 0){
        block_smap <- matrix(as.numeric(scale(edna_all[,c(effect_tax)])), ncol = 1)
      }else if(n_cause+1 >= E_tax){
        # Use naive embedding to reconstruct state space
        block_smap <- apply(cbind(edna_all[,c(effect_tax, causal_tax)],
                                  clim_all[,c(causal_clim)]), 2, function(x) as.numeric(scale(x)))
      }else{
        # Include time-lag components if the number of causal variables is smaller than E
        block_lags <- make_block(as.numeric(scale(edna_all[,effect_tax])), max_lag = E_tax - n_cause)
        block_smap <- apply(cbind(edna_all[,c(effect_tax, causal_tax)],
                                  clim_all[,c(causal_clim)]), 2, function(x) as.numeric(scale(x)))
        block_smap <- cbind(block_smap, block_lags[,3:(E_tax + 1 - n_cause)])
      }
    }else if(embedding == "limited_by_E"){
      # Apply "embedding limited by E" (the upper limit of the number of variables are E)
      if(E_tax == 1){
        block_smap <- as.numeric(scale(edna_all[,effect_tax]))
      }else if(n_cause+1 >= E_tax){
        # Select the most predictive causal variables based on "log_d_rmse"
        causal_set <- rbind(causal_dnaxdna_sub[causal_dnaxdna_sub$effect_var == effect_tax,],
                            causal_dnaclim_sub[causal_dnaclim_sub$effect_var == effect_tax,])
        causal_order <- causal_set[order(causal_set$log_d_rmse, decreasing = F),]$cause_var[1:(E_tax-1)]
        block_smap1 <- cbind(edna_all, clim_all)[,c(effect_tax, causal_order)]
        block_smap <- apply(block_smap1, 2, function(x) as.numeric(scale(x)))
      }else{
        # Include time-lag components if the number of causal variables is smaller than E
        block_lags <- make_block(as.numeric(scale(edna_all[,effect_tax])), max_lag = E_tax-n_cause)
        block_smap1 <- apply(cbind(edna_all[,c(effect_tax, causal_tax)],
                                   clim_all[,c(causal_clim)]), 2, function(x) as.numeric(scale(x)))
        block_smap <- cbind(block_smap1, block_lags[,3:(E_tax-n_cause+1)])
      }
    }
    
    # Determine the best nonlinear parameter for the block (parallel version)
    #start_time_g2 <- proc.time()[3]
    best_par_check <- data.frame(theta = NA, lambda = NA, rmse = NA)
    
    for(k in 1:length(theta_test)){
      for(l in 1:length(lambda_test)){
        # Perform the regularized S-map
        theta_i <- theta_test[k]
        lambda_i <- lambda_test[l]
        best_check_tmp <- try(extended_lnlp(block_smap, lib = lib, pred = pred,
                                            tp = 1, target_column = 1,
                                            theta = theta_i, method = "s-map",
                                            regularized = regularized,
                                            lambda = lambda_i, alpha = alpha, # lasso S-map
                                            glmnet_parallel = FALSE)$stats$rmse,
                              silent = TRUE)
        
        if(class(best_check_tmp) == "try-error"){
          best_par_check[((k-1)*length(lambda_test)+l),]  <- data.frame(theta = theta_i, lambda = lambda_i, rmse = NA)
        }else{
          best_par_check[((k-1)*length(lambda_test)+l),]  <- data.frame(theta = theta_i, lambda = lambda_i, rmse = best_check_tmp)
        }
      }
    }
    
    #proc.time()[3] - start_time_g2
    #------------------------------- edited until here ------------------------------#
    # Save progress messeges
    time_used <- round(proc.time()[3] - start_time, digits = 2)
    progress_message <- sprintf("Cycle %s/%s finished; %s sec elapsed",
                       i, total_cycle, time_used)
    write.csv(progress_message, sprintf("%s/progress_%04d.txt", progress_folder, i))
    
    # Store results
    if(all(is.na(best_par_check[,3]))){
      best_par <- data.frame(theta = NA, lambda = NA)
    }else{
      best_par <- best_par_check[which.min(best_par_check$rmse), c("theta", "lambda")]
    }
    best_par
  })
  
  # Add name to bestT_block
  rownames(best_par_block) <- taxa_list
  
  # Remove temporal files
  temp_files <- list.files(progress_folder)
  file.remove(sprintf("%s/%s", progress_folder, temp_files))
  file.remove(progress_folder)
  
  # Save progress messeges
  time_used_g <- round(proc.time()[3] - start_time_g, digits = 2)
  cat("Best nonlinear parameter determined:", total_cycle, "cycles;", time_used_g, "sec elapsed\n")
  
  # Return results
  return(best_par_block)
}


#-------------------- Do multivariate regularized S-map --------------------#
# Best nonlinear parameters are stored in "bestT_block"
multi_reglSmap_all <- function(taxa_list,
                               bestPar_block,
                               lib = edna_lib,
                               pred = lib,
                               embedding = "naive", # "intermediate", "limited_by_E"
                               regularized = TRUE,
                               alpha = 1,
                               progress_folder = "99_progress"){

  # Create progress folder
  dir.create(progress_folder)
  
  # Set time (global)
  start_time_g <- proc.time()[3]
  
  # Extract taxa only included in the subset
  valid_id_dna <- !is.na(match(causal_dnaxdna2$effect_var, taxa_list)) & !is.na(match(causal_dnaxdna2$cause_var, taxa_list))
  valid_id_clm <- !is.na(match(causal_dnaclim3$effect_var, taxa_list))
  causal_dnaxdna_sub <- causal_dnaxdna2[valid_id_dna,]
  causal_dnaclim_sub <- causal_dnaclim3[valid_id_clm,]
  
  # Performing multivariate S-map using best theta
  total_cycle <- length(taxa_list)
  #smap_res <- list(NULL)
  
  smap_res <- pforeach(i = 1:length(taxa_list), .combine = c)({
    # Set time
    start_time <- proc.time()[3]
    
    # Extract information for S-map
    E_tax <- Eedna[taxa_list[i]]
    effect_tax <- taxa_list[i]
    causal_tax <- causal_dnaxdna_sub[causal_dnaxdna_sub$effect_var == effect_tax,]$cause_var
    causal_clim <- causal_dnaclim_sub[causal_dnaclim_sub$effect_var == effect_tax,]$cause_var
    
    # Calculate the number of causal variables
    n_cause <- length(causal_tax) + length(causal_clim)
    
    # Make block for S-map
    # Three options to make a block
    if(embedding == "naive"){
      # Apply "naive embedding" (use all lags and causal factors to reconstruct state space)
      if(E_tax == 1 & n_cause == 0){
        block_smap <- matrix(as.numeric(scale(edna_all[,c(effect_tax)])), ncol = 1)
        colnames(block_smap) <- effect_tax
      }else{
        # Include time-lag components if the number of causal variables is smaller than E
        block_smap <- apply(cbind(edna_all[,c(effect_tax, causal_tax)],
                                  clim_all[,c(causal_clim)]), 2, function(x) as.numeric(scale(x)))
        colnames(block_smap) <- c(effect_tax, causal_tax, causal_clim)
        if(E_tax != 1){
          block_lags <- make_block(as.numeric(scale(edna_all[,effect_tax])), max_lag = E_tax)
          colnames(block_lags)[2:(E_tax+1)] <- c(effect_tax, sprintf("%s_lag%s", effect_tax, 1:(E_tax-1)))
          block_smap <- cbind(block_smap, block_lags[,3:(E_tax + 1)])
          if(E_tax == 2) colnames(block_smap)[ncol(block_smap)] <- colnames(block_lags)[E_tax+1]
        }
      }
    }else if(embedding == "intermediate"){
      # Apply "intermediate embedding" (use all causal factors and some lags to reconstruct state space)
      if(E_tax == 1 & n_cause == 0){
        block_smap <- matrix(as.numeric(scale(edna_all[,c(effect_tax)])), ncol = 1)
        colnames(block_smap) <- effect_tax
      }else if(n_cause+1 >= E_tax){
        # Use naive embedding to reconstruct state space
        block_smap <- apply(cbind(edna_all[,c(effect_tax, causal_tax)],
                                  clim_all[,c(causal_clim)]), 2, function(x) as.numeric(scale(x)))
        colnames(block_smap) <- c(effect_tax, causal_tax, causal_clim)
      }else{
        # Include time-lag components if the number of causal variables is smaller than E
        block_lags <- make_block(as.numeric(scale(edna_all[,effect_tax])), max_lag = E_tax - n_cause)
        block_smap <- apply(cbind(edna_all[,c(effect_tax, causal_tax)],
                                  clim_all[,c(causal_clim)]), 2, function(x) as.numeric(scale(x)))
        colnames(block_lags)[2:(E_tax+1)] <- c(effect_tax, sprintf("%s_lag%s", effect_tax, 1:(E_tax-1)))
        colnames(block_smap) <- c(effect_tax, causal_tax, causal_clim)
        block_smap <- cbind(block_smap, block_lags[,3:(E_tax+1-n_cause)])
        if(E_tax+1-n_cause == 3) colnames(block_smap)[ncol(block_smap)] <- colnames(block_lags)[E_tax+1-n_cause]
      }
    }else if(embedding == "limited_by_E"){
      # Apply "embedding limited by E" (the upper limit of the number of variables are E)
      if(E_tax == 1){
        block_smap <- as.numeric(scale(edna_all[,effect_tax]))
        colnames(block_smap) <- effect_tax
      }else if(n_cause+1 >= E_tax){
        # Select the most predictive causal variables based on "log_d_rmse"
        causal_set <- rbind(causal_dnaxdna_sub[causal_dnaxdna_sub$effect_var == effect_tax,],
                            causal_dnaclim_sub[causal_dnaclim_sub$effect_var == effect_tax,])
        causal_order <- causal_set[order(causal_set$log_d_rmse, decreasing = F),]$cause_var[1:(E_tax-1)]
        block_smap1 <- cbind(edna_all, clim_all)[,c(effect_tax, causal_order)]
        colnames(block_smap1) <- c(effect_tax, causal_order)
        block_smap <- apply(block_smap1, 2, function(x) as.numeric(scale(x)))
      }else{
        # Include time-lag components if the number of causal variables is smaller than E
        block_lags <- make_block(as.numeric(scale(edna_all[,effect_tax])), max_lag = E_tax-n_cause)
        colnames(block_lags)[2:(E_tax-n_cause+1)] <- c(effect_tax, sprintf("%s_lag%s", effect_tax, 1:(E_tax-n_cause-1)))
        block_smap1 <- apply(cbind(edna_all[,c(effect_tax, causal_tax)],
                                   clim_all[,c(causal_clim)]), 2, function(x) as.numeric(scale(x)))
        colnames(block_smap1) <- c(effect_tax, causal_tax, causal_clim)
        block_smap <- cbind(block_smap1, block_lags[,3:(E_tax-n_cause+1)])
        if(E_tax+1-n_cause == 3) colnames(block_smap)[ncol(block_smap)] <- colnames(block_lags)[E_tax+1-n_cause]
      }
    }
    
    # Check block_smap colnames
    if(any(!is.na(match(colnames(block_smap), "")))) stop("Column names of the block inappropriate.")
    
    # Perform the multivariate S-map and save S-map coefficients
    if(!is.na(bestPar_block[i,1])){
      smap_res_tmp  <- try(extended_lnlp(block_smap, lib = lib, pred = pred,
                                     tp = 1, target_column = 1,
                                     theta = bestPar_block[i,"theta"], method = "s-map",
                                     regularized = regularized,
                                     lambda = bestPar_block[i,"lambda"], alpha = alpha, # lasso S-map
                                     glmnet_parallel = FALSE,
                                     save_smap_coefficients = TRUE),
                           silent = TRUE)
      if(class(smap_res_tmp) == "try-error"){
        smap_res_tmp <- list(NA, NA, NA)
      }else{
        # Add colname to the smap_coefficient results
        colnames(smap_res_tmp$smap_coefficients)[2:(ncol(block_smap)+1)] <- colnames(block_smap)
      }
    }else{
      smap_res_tmp <- list(NA, NA, NA)
    }

    # Save progress messeges
    time_used <- round(proc.time()[3] - start_time, digits = 2)
    progress_message <- sprintf("Cycle %s/%s finished; %s sec elapsed",
                                i, total_cycle, time_used)
    write.csv(progress_message, sprintf("%s/progress_%04d.txt", progress_folder, i))
    
    list(smap_res_tmp)
  })
  
  # Re-assign list name
  for(i in 1:length(smap_res)) names(smap_res[[i]]) <- c("model_output", "stats", "smap_coefficients")
  names(smap_res) <- taxa_list

  # Remove temporal files
  temp_files <- list.files(progress_folder)
  file.remove(sprintf("%s/%s", progress_folder, temp_files))
  file.remove(progress_folder)
  
  # Save progress messeges
  time_used_g <- round(proc.time()[3] - start_time_g, digits = 2)
  cat("Multivariate regularized S-map finished:", total_cycle, "cycles;", time_used_g, "sec elapsed\n")
  
  # Return results
  return(smap_res)
}


