####
#### Helper functions for CERrice2017 time series
#### Functions for empirical dynamic modeling
####


#-------------------- Subset extraction function --------------------#
# Define function
extract_subset <- function(tax_row_i, edna_tax_data){
  edna_tax_tmp1 <- edna_tax_data[edna_tax_data$superkingdom == tax_summary[tax_row_i, 1],]
  edna_tax_tmp1 <- edna_tax_tmp1[edna_tax_tmp1$kingdom == tax_summary[tax_row_i, 2],]
  edna_tax_tmp1 <- edna_tax_tmp1[edna_tax_tmp1$phylum == tax_summary[tax_row_i, 3],]
  return(rownames(edna_tax_tmp1))
}
#--------------------------------------------------------------------#

# Function to extract connected vertecies
extract_connected_vertex <- function(taxa_list){
  # Extract taxa only included in the subset
  valid_id_dna <- !is.na(match(causal_dnaxdna2$effect_var, taxa_list)) & !is.na(match(causal_dnaxdna2$cause_var, taxa_list))
  causal_dnaxdna_sub <- causal_dnaxdna2[valid_id_dna,]
  
  if(nrow(causal_dnaxdna_sub) > 0){
    # Create a weighted (time-averaged) adjucent matrix
    AJ_weighted <- matrix(0, ncol = length(taxa_list), nrow = length(taxa_list))
    colnames(AJ_weighted) <- rownames(AJ_weighted) <- taxa_list
    for(i in 1:nrow(causal_dnaxdna_sub)){
      AJ_weighted[causal_dnaxdna_sub$effect_var[i],
                  causal_dnaxdna_sub$cause_var[i]] <- abs(causal_dnaxdna_sub$log_d_rmse[i])
    }
    
    # Convert to a unweighted (time-averaged) adjucent matrix
    AJ_unweighted <- AJ_weighted
    AJ_unweighted[AJ_weighted > 0] <- 1
    nolink_id <- !(colSums(AJ_unweighted) > 0) & !(rowSums(AJ_unweighted) > 0)
    connected_vertex <- rownames(AJ_unweighted)[!nolink_id]
    
    return(list(connected_vertex = connected_vertex, AJ_weighted = AJ_weighted, AJ_unweighted = AJ_unweighted))
  }else{
    return(list(connected_vertex = NULL, AJ_weighted = NULL, AJ_unweighted = NULL))
  }
}



#-------------------- Multivariate S-map, Network reconstruction and dynamic stability --------------------#
regl_smap_to_stability <- function(taxa_list,
                                   progress_output = "progress"){
  total_cycle <- nrow(edna_all)
  e_value_all <- data.frame(NULL)
  connected_vertex_all <- list(NULL)
  compiled_smap_all <- list(NULL)
  
  # Start main loop
  for(time_i in 1:nrow(edna_all)){
    # Set time
    start_time <- proc.time()[3]
    
    # Select subset
    nonzero_taxa <- colnames(edna_all[,taxa_list][which(edna_all[time_i,taxa_list] > 0)])
    if(length(nonzero_taxa) > 0){
      connected_vertex <- extract_connected_vertex(nonzero_taxa)$connected_vertex
    }else{
      connected_vertex <- NULL
    }
    
    if(length(connected_vertex) > 0){
      connected_vertex_all[[time_i]] <- connected_vertex
      
      # No.2: Determine the best embedding dimension using Regularized S-map
      bestT_all <- det_bestT_reglSmap(connected_vertex, embedding = "naive",
                                      lib = edna_lib, pred = c(time_i-1, time_i+1),
                                      alpha = 0, # Ridge S-map
                                      progress_folder = sprintf("%s_temp",  progress_output))
      
      # No.3: Perform multivariate S-map
      m_smap_all <- multi_reglSmap_all(connected_vertex, bestT_all, embedding = "naive",
                                       lib = edna_lib, pred = c(time_i-1, time_i+1),
                                       alpha = 0, # Ridge S-map
                                       progress_folder = sprintf("%s_temp", progress_output))
      
      # Check prediction length
      prediction_success <- unlist(sapply(m_smap_all, function(x) nrow(x[[1]])))
      if(length(prediction_success) == length(connected_vertex)){
        
        # No.4: Compile S-map results
        compiled_smap <- compile_reglSmap_res(connected_vertex,  m_smap_all, calc_lib_id = time_i)
        
        # Save compiled S-map results
        compiled_smap_all[[time_i]] <- compiled_smap
        
        # No.5: Calculate dynamic stability
        dstab_all <- dynamic_stability(compiled_smap$J1_matrix_all,
                                       compiled_smap$J2_matrix_pre_all,
                                       exclude_zero = FALSE,
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

    # Show message
    time_used <- round(proc.time()[3] - start_time, digits = 2)
    cat("Quantifying local dynamic stability: Time", time_i, "/", total_cycle, "finished;", time_used, "sec elapsed\n\n")
  }
  
  return(list(eigen_value = e_value_all,
              interaction_matrix = compiled_smap_all,
              connected_vertex = connected_vertex_all))
}



#-------------------- Axis name list --------------------#
network_prop_names <- c("dynamic_stab" = "Dynamic stability",
  "total_div" = "ASV diversity (total)",
  "int_n" = "N of interactions",
  "int_mean" = "Mean IS",
  "int_med" = "Median IS",
  "int_min" = "Minimum IS",
  "int_weak" = "Weak IS index",
  "int_sd" = "S.D. of IS",
  "int_max" = "Maximum IS",
  "Prok_dna_copy_sum" = "Toal DNA conc. (Prokaryote)",
  "prok_div" = "ASV diversity (Prokaryote)",
  "prok_int_mean" = "Mean IS (Prokaryote)",
  "prok_int_med" = "Median IS (Prokaryote)",
  "prok_int_n" = "N of interactions (Prokaryote)",
  "prok_int_weak" = "Weak IS index (Prokaryote)",
  
  "Euk_dna_copy_sum" = "Total DNA conc. (Eukaryote)",
  "euk_div" = "ASV diversity (Eukaryote)",
  "euk_int_mean" = "Mean IS (Eukaryote)",
  "euk_int_med" = "Median IS (Eukaryote)",
  "euk_int_n" = "N of interactions (Eukaryote)",
  "euk_int_weak" = "Weak IS index (Eukaryote)",
  
  "Fungi_dna_copy_sum" = "Total DNA conc. (Fungi)",
  "fungi_div" = "ASV diversity (Fungi)",
  "fungi_int_med" = "Median IS (Fungi)",
  "fungi_int_mean" = "Mean IS (Fungi)",
  "fungi_int_n" = "N of interactions (Fungi)",
  "fungi_int_weak" = "Weak IS index (Fungi)",
  
  "virid_div" = "ASV diversity (Viridiplantae)",
  "virid_int_weak" = "Weak IS index (Viridiplantae)",
  "virid_int_med" = "Median IS (Viridiplantae)",
  "virid_int_mean" = "Mean IS (Viridiplantae)",
  "virid_int_n" = "N of interactions (Viridiplantae)",
  
  "Inv_dna_copy_sum" = "Total DNA conc. (COI)",
  
  "metaz_div" = "ASV diversity (Metazoa)",
  "metaz_int_weak" = "Weak IS index (Metazoa)",
  "metaz_int_med" = "Median IS (Metazoa)",
  "metaz_int_mean" = "Mean IS (Metazoa)",
  "metaz_int_n" = "N of interactions (Metazoa)",
  
  "undet_div" = "ASV diversity (Undetermined)",
  "temp_mean" = "Mean air temperature",
  "actvp_mean" = "Mean actual vapor pressure")

