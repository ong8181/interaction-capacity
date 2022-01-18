####
#### CERrice2017 All data analysis
#### No.12 Calculating taxa-specific IS per sp
####

# Load workspace and objects
load("11_DiversityNetworkOut/11_DiversityNetworkOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder12 <- "12_TaxaSpecificISOut"
dir.create(output_folder12)

# Compile interaction matrix properties
int_mat <- readRDS("07_RegularizedSmapOut/interaction_matrix_nonzeroridge.obj")
edna_intn <- edna_ic <- edna_intmax <- edna_all[,edna_var_coln]
edna_intn[] <- edna_ic[] <- edna_intmax[] <- NA # Set all elements as zero
edna_cor <- rep(NA, dim(edna_all)[1])
edna_cor_abs <- rep(NA, dim(edna_all)[1])

for(time_i in 1:length(int_mat)){
  if(length(int_mat[[time_i]]) > 1){
    if(all(!is.na(int_mat[[time_i]]$J1_matrix_all[[1]]))){
      int_mat_tmp <- int_mat[[time_i]]$J1_matrix_all[[1]]
      diag(int_mat_tmp) <- 0
      
      # Effect link + Causal link = 2 * Total link
      edna_intn[time_i, colnames(int_mat_tmp)] <- rowSums(int_mat_tmp != 0, na.rm = T) + colSums(int_mat_tmp != 0, na.rm = T)
      # Effect IS + Causal IS = 2 * Total IS
      edna_ic[time_i, colnames(int_mat_tmp)] <- rowSums(abs(int_mat_tmp), na.rm = T) + colSums(abs(int_mat_tmp), na.rm = T)
      # Effect IS + Causal IS = 2 * Total IS
      col_max <- apply(abs(int_mat_tmp), 2, function(x) max(x, na.rm = T))
      row_max <- apply(abs(int_mat_tmp), 1, function(x) max(x, na.rm = T))
      edna_intmax[time_i, colnames(int_mat_tmp)] <- apply(cbind(col_max, row_max), 1, max)
      
      # Correlation between effect and causal strengths
      uptri <- int_mat_tmp[upper.tri(int_mat_tmp)]
      lwtri <- int_mat_tmp[lower.tri(int_mat_tmp)]
      nonzero_int <- (uptri != 0 | lwtri != 0)
      edna_cor[time_i] <- cor(uptri[nonzero_int], lwtri[nonzero_int], method = "pearson", use = "complete.obs")
      edna_cor_abs[time_i] <- cor(abs(uptri[nonzero_int]), abs(lwtri[nonzero_int]), method = "pearson", use = "complete.obs")
    }
  }
}

#hist(edna_cor); hist(edna_cor_abs)
valid_row_id <- c()
for(i in 1:nrow(edna_lib)) valid_row_id <- c(valid_row_id, edna_lib[i,1]:edna_lib[i,2])
edna_ic <- edna_ic[valid_row_id,]
edna_intn <- edna_intn[valid_row_id,]
edna_intmax <- edna_intmax[valid_row_id,]

edna_prop$int_max_per_com <- edna_prop$int_per_link <-
  edna_prop$int_cor <- edna_prop$int_cor_abs <- NaN
edna_prop$int_max_per_com[valid_row_id] <- apply(edna_ic, 1, function(x){if(any(!is.na(x))){max(x, na.rm = T)}else{NA}})
edna_prop$int_per_link[valid_row_id] <- rowSums(edna_ic, na.rm = T)/rowSums(edna_intn, na.rm = T)
edna_prop$int_cor[valid_row_id] <- edna_cor[valid_row_id]
edna_prop$int_cor_abs[valid_row_id] <- edna_cor_abs[valid_row_id]

# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder12, output_folder12))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder12, substr(Sys.time(), 1, 10)))
