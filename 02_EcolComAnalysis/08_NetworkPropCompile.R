####
#### CERrice2017 All data analysis
#### No.8 Compile network properties
####

# Load workspace and objects
load("07_RegularizedSmapOut/07_RegularizedSmapOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder08 <- "08_NetworkPropCompileOut"
dir.create(output_folder08)

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.28
library(Rcpp); packageVersion("Rcpp") # 1.0.3, 2020.1.28
library(RcppArmadillo); packageVersion("RcppArmadillo") # 0.9.800.3.0, 2020.1.28
library(Matrix); packageVersion("Matrix") # 1.2.18, 2020.1.28

# Load interaction matrix data
dy_stab <- readRDS(sprintf("%s/dynamic_stability_nonzeroridge.obj", output_folder07))
dy_stab <- dy_stab[pred_lib,]
int_mat <- readRDS(sprintf("%s//interaction_matrix_nonzeroridge.obj", output_folder07))

# Compile data
edna_prop <- edna_all[,c("date", "plot")]
edna_prop$total_dna <- rowSums(edna_all[,edna_var_coln], na.rm = T)
edna_prop$total_div <- rowSums(edna_all[,edna_var_coln] > 0, na.rm = T)

# Add stability values
edna_prop$dynamic_stab <- NaN
edna_prop$dynamic_stab[pred_lib] <- abs(dy_stab$ev_1)

# Add system productivity values
edna_prop$total_np <- c(NA, diff(edna_prop$total_dna, lag = 1, differences = 1))

# Add C.V. of community dynamics (calculate CV by the sliding window approach)
edna_cv <- edna_all[,edna_var_coln]
mean_div <- edna_prop$total_div
edna_cv[] <- NaN
mean_div[] <- NaN
window_width <- 3 # One week (day - 3 to day + 3) time window

for(time_window in window_width:(nrow(edna_all) - window_width)){
  edna_cv[time_window,] <- apply(edna_all[(time_window - window_width):(time_window + window_width), edna_var_coln], 2,
                                 function(x) sd(x, na.rm = T)/mean(x, na.rm = T))
  mean_div[time_window] <- mean(edna_prop$total_div[(time_window - window_width):(time_window + window_width)], na.rm = T)
}

edna_prop$mean_cv <- rowMeans(edna_cv[,1:length(edna_var_coln)], na.rm = T)
edna_prop$mean_cv[edna_prop$total_div == 0] <- NaN
edna_prop$mean_div4cv <- mean_div
plot(edna_prop$mean_cv ~ edna_prop$total_div)
plot(edna_prop$mean_cv ~ edna_prop$mean_div4cv)

# Identify taxa ID belonging to a specific taxa
prok_id <- which(edna_tax3$superkingdom == "Bacteria" |  edna_tax3$superkingdom == "Archaea")
euk_id <- which(edna_tax3$superkingdom == "Eukaryota")
undet_id <- which(edna_tax3$superkingdom == "")
fungi_id <- which(edna_tax3$kingdom == "Fungi")
metaz_id <- which(edna_tax3$kingdom == "Metazoa")
virid_id <- which(edna_tax3$kingdom == "Viridiplantae")

# Check taxa-specific properties
taxa_div_check <- F
taxa_abn_check <- F
taxa_np_check <- F
taxa_int_prop_check <- F

# Add taxa diversity
if(taxa_div_check){
  edna_prop$prok_div <- rowSums(edna_all[,edna_var_coln[prok_id]] > 0, na.rm = T)
  edna_prop$euk_div <- rowSums(edna_all[,edna_var_coln[euk_id]] > 0, na.rm = T)
  edna_prop$undet_div <- rowSums(edna_all[,edna_var_coln[undet_id]] > 0, na.rm = T)
  edna_prop$fungi_div <- rowSums(edna_all[,edna_var_coln[fungi_id]] > 0, na.rm = T)
  edna_prop$metaz_div <- rowSums(edna_all[,edna_var_coln[metaz_id]] > 0, na.rm = T)
  edna_prop$virid_div <- rowSums(edna_all[,edna_var_coln[virid_id]] > 0, na.rm = T)
}

# Add taxa abundance
if(taxa_abn_check){
  edna_prop$prok_abn <- rowSums(edna_all[,edna_var_coln[prok_id]], na.rm = T)
  edna_prop$euk_abn <- rowSums(edna_all[,edna_var_coln[euk_id]], na.rm = T)
  edna_prop$undet_abn <- rowSums(edna_all[,edna_var_coln[undet_id]], na.rm = T)
  edna_prop$fungi_abn <- rowSums(edna_all[,edna_var_coln[fungi_id]], na.rm = T)
  edna_prop$metaz_abn <- rowSums(edna_all[,edna_var_coln[metaz_id]], na.rm = T)
  edna_prop$virid_abn <- rowSums(edna_all[,edna_var_coln[virid_id]], na.rm = T)
}

# Add taxa-specific productivity
if(taxa_np_check){
  edna_prop$prok_np <- c(NA, diff(edna_prop$prok_abn))
  edna_prop$euk_np <- c(NA, diff(edna_prop$euk_abn))
  edna_prop$undet_np <- c(NA, diff(edna_prop$undet_abn))
  edna_prop$fungi_np <- c(NA, diff(edna_prop$fungi_abn))
  edna_prop$metaz_np <- c(NA, diff(edna_prop$metaz_abn))
  edna_prop$virid_np <- c(NA, diff(edna_prop$virid_abn))
}

# Calculate network properties
edna_prop$int_n <- edna_prop$int_mean <- edna_prop$int_med <-
  edna_prop$int_max <- edna_prop$int_min <- edna_prop$int_sd <-edna_prop$int_weak <-
  edna_prop$int_per_sp <- edna_prop$connectance <- NaN

if(taxa_int_prop_check){
  edna_prop$prok_int_n <- edna_prop$prok_int_mean <- edna_prop$prok_int_med <- edna_prop$prok_int_weak <-
    edna_prop$prok_int_per_sp <- edna_prop$prok_connectance <- 
    edna_prop$euk_int_n <- edna_prop$euk_int_mean <- edna_prop$euk_int_med <- edna_prop$euk_int_weak <-
    edna_prop$euk_int_per_sp <- edna_prop$euk_connectance <- 
    edna_prop$fungi_int_n <- edna_prop$fungi_int_mean <- edna_prop$fungi_int_med <- edna_prop$fungi_int_weak <-
    edna_prop$fungi_int_per_sp <- edna_prop$fungi_connectance <- 
    edna_prop$metaz_int_n <- edna_prop$metaz_int_mean <- edna_prop$metaz_int_med <- edna_prop$metaz_int_weak <-
    edna_prop$metaz_int_per_sp <- edna_prop$metaz_connectance <- 
    edna_prop$virid_int_n <- edna_prop$virid_int_mean <- edna_prop$virid_int_med <- edna_prop$virid_int_weak <- 
    edna_prop$virid_int_per_sp <- edna_prop$virid_connectance <- NaN
}

# Calculate and assign properties of interaction networks
for(i in 1:length(int_mat))
{
  if(!is.na(int_mat[i]))
  {
    int_temp <- int_mat[[i]]$J1_matrix_all[[1]]
    
    # Calculate network properties
    diag(int_temp) <- 0 # Remove self-regulation terms
    int_vals <- abs(c(int_temp[int_temp != 0 & !is.na(int_temp)])) # Remove 0 and NAs, and take absolute values

    if(length(int_vals) > 0)
    {
      edna_prop$int_n[i] <- length(int_vals) # The number of interactions
      edna_prop$int_mean[i] <- mean(int_vals)
      edna_prop$int_med[i] <- median(int_vals)
      edna_prop$int_max[i] <- max(int_vals)
      edna_prop$int_min[i] <- min(int_vals)
      edna_prop$int_sd[i] <- sd(int_vals)
      edna_prop$int_weak[i] <- edna_prop$int_med[i]/edna_prop$int_max[i]
      edna_prop$int_per_sp[i] <- 2*sum(int_vals)/edna_prop$total_div[i]
      edna_prop$connectance[i] <- edna_prop$int_n[i]/(edna_prop$total_div[i]^2)
    }
    
    if(taxa_int_prop_check){
      # Identify taxa-specific IDs
      prok_int_temp <- as.numeric(na.omit(match(rownames(edna_tax3)[prok_id], colnames(int_temp)))) %>% int_temp[,.]
      euk_int_temp <- as.numeric(na.omit(match(rownames(edna_tax3)[euk_id], colnames(int_temp)))) %>% int_temp[,.]
      fungi_int_temp <- as.numeric(na.omit(match(rownames(edna_tax3)[fungi_id], colnames(int_temp)))) %>% int_temp[,.]
      metaz_int_temp <- as.numeric(na.omit(match(rownames(edna_tax3)[metaz_id], colnames(int_temp)))) %>% int_temp[,.]
      virid_int_temp <- as.numeric(na.omit(match(rownames(edna_tax3)[virid_id], colnames(int_temp)))) %>% int_temp[,.]
      
      # Collect values
      prok_vals <- abs(c(prok_int_temp[prok_int_temp != 0 & !is.na(prok_int_temp)])) # Remove 0 and NAs, and take absolute values
      euk_vals <- abs(c(euk_int_temp[euk_int_temp != 0 & !is.na(euk_int_temp)])) # Remove 0 and NAs, and take absolute values
      fungi_vals <- abs(c(fungi_int_temp[fungi_int_temp != 0 & !is.na(fungi_int_temp)])) # Remove 0 and NAs, and take absolute values
      metaz_vals <- abs(c(metaz_int_temp[metaz_int_temp != 0 & !is.na(metaz_int_temp)])) # Remove 0 and NAs, and take absolute values
      virid_vals <- abs(c(virid_int_temp[virid_int_temp != 0 & !is.na(virid_int_temp)])) # Remove 0 and NAs, and take absolute values

      # Assign taxa-specific interation values
      if(length(prok_vals) > 0)
      {
        edna_prop$prok_int_n[i] <- length(prok_vals)
        edna_prop$prok_int_mean[i] <- mean(prok_vals)
        edna_prop$prok_int_med[i] <- median(prok_vals)
        edna_prop$prok_int_weak[i] <- median(prok_vals)/max(prok_vals)
        edna_prop$prok_int_per_sp[i] <- 2*sum(prok_vals)/edna_prop$prok_div[i]
        edna_prop$prok_connectance[i] <- length(prok_vals)/(edna_prop$prok_div[i]^2)
      }
      
      if(length(euk_vals) > 0)
      {
        edna_prop$euk_int_n[i] <- length(euk_vals)
        edna_prop$euk_int_mean[i] <- mean(euk_vals)
        edna_prop$euk_int_med[i] <- median(euk_vals)
        edna_prop$euk_int_weak[i] <- median(euk_vals)/max(euk_vals)
        edna_prop$euk_int_per_sp[i] <- 2*sum(euk_vals)/edna_prop$euk_div[i]
        edna_prop$euk_connectance[i] <- length(euk_vals)/(edna_prop$euk_div[i]^2)
      }
      
      if(length(fungi_vals) > 0)
      {
        edna_prop$fungi_int_n[i] <- length(fungi_vals)
        edna_prop$fungi_int_mean[i] <- mean(fungi_vals)
        edna_prop$fungi_int_med[i] <- median(fungi_vals)
        edna_prop$fungi_int_weak[i] <- median(fungi_vals)/max(fungi_vals)
        edna_prop$fungi_int_per_sp[i] <- 2*sum(fungi_vals)/edna_prop$fungi_div[i]
        edna_prop$fungi_connectance[i] <- length(fungi_vals)/(edna_prop$fungi_div[i]^2)
      }
      
      if(length(metaz_vals) > 0)
      {
        edna_prop$metaz_int_n[i] <- length(metaz_vals)
        edna_prop$metaz_int_mean[i] <- mean(metaz_vals)
        edna_prop$metaz_int_med[i] <- median(metaz_vals)
        edna_prop$metaz_int_weak[i] <- median(metaz_vals)/max(metaz_vals)
        edna_prop$metaz_int_per_sp[i] <- 2*sum(metaz_vals)/edna_prop$metaz_div[i]
        edna_prop$metaz_connectance[i] <- length(metaz_vals)/(edna_prop$metaz_div[i]^2)
      }
      
      if(length(virid_vals) > 0)
      {
        edna_prop$virid_int_n[i] <- length(virid_vals)
        edna_prop$virid_int_mean[i] <- mean(virid_vals)
        edna_prop$virid_int_med[i] <- median(virid_vals)
        edna_prop$virid_int_weak[i] <- median(virid_vals)/max(virid_vals)
        edna_prop$virid_int_per_sp[i] <- 2*sum(virid_vals)/edna_prop$virid_div[i]
        edna_prop$virid_connectance[i] <- length(virid_vals)/(edna_prop$virid_div[i]^2)
      }
    }
  }
}

# Add climate variables
edna_prop$temp_mean <- clim_all$temp_mean
#edna_prop$actvp_mean <- clim_all$actvp_mean
#edna_prop$light_mean <- clim_all$light_mean

# Delete unnecssary objects
#rm(AJ_unweighted)
#rm(AJ_weighted)
rm(int_mat)
rm(list=c("i", "j", "k"))

# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder08, output_folder08))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder08, substr(Sys.time(), 1, 10)))
