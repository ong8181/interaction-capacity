####
#### CERrice2017 All data analysis
#### No. 0 Compile all data (Rice, Climate and eDNA data)
####

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder00 <- "00_CompileAllDataOut"
dir.create("00_SessionInfo")
dir.create(output_folder00)

# Load library and functions
library(phyloseq); packageVersion("phyloseq") # 1.28.0, 2020.1.6
library(lubridate); packageVersion("lubridate") # 1.7.4, 2019.10.23
library(rEDM); packageVersion("rEDM") # 0.7.5, 2020.1.6
library(pforeach); packageVersion("pforeach") # 1.3, 2019.10.23
source("functions/F01_HelperFunctions.R")
source("functions/BidirectSimplex_rEDM_0.7.4.R")
source("functions/BidirectBlocklnlp_rEDM_0.7.4.R")

# Load main data (Climate and eDNA data)
clim_all0 <- readRDS("data_climate/clm_d_mean.obj")
edna_all0 <- readRDS("../02_TimeSeriesCompile/11_TSfilter02Out/ps_comb_filt.obj")

# Load supplement data
clim_max <- readRDS("data_climate/clm_d_max.obj")
clim_min <- readRDS("data_climate/clm_d_min.obj")

# Adjust data format (change column order)
clim_all0 <- cbind(clim_all0[,1:6], clim_all0[,17:21], clim_all0[,27:31], clim_all0[,7:11], clim_all0[,12:16])
clim_max <- cbind(clim_max[,1:6], clim_max[,17:21], clim_max[,12:16])
clim_min <- cbind(clim_min[,1:6], clim_min[,17:21], clim_min[,12:16])

# Set climate dataframe as a default
# Climate dataframe for each plot
date_for_df <- ymd(clim_all0[,"date"])
p1_col <- substr(colnames(clim_all0), 1, 5) == "D0202"
p2_col <- substr(colnames(clim_all0), 1, 5) == "D0203"
p3_col <- substr(colnames(clim_all0), 1, 5) == "D0204"
p4_col <- substr(colnames(clim_all0), 1, 5) == "D0205"
p5_col <- substr(colnames(clim_all0), 1, 5) == "D0206"
p1_col2 <- substr(colnames(clim_max), 1, 5) == "D0202"
p2_col2 <- substr(colnames(clim_max), 1, 5) == "D0203"
p3_col2 <- substr(colnames(clim_max), 1, 5) == "D0204"
p4_col2 <- substr(colnames(clim_max), 1, 5) == "D0205"
p5_col2 <- substr(colnames(clim_max), 1, 5) == "D0206"

plot1_clim_df <- cbind(date_for_df, rep(1, nrow(clim_all0)), clim_all0[,p1_col], clim_max[,p1_col2], clim_min[,p1_col2])
plot2_clim_df <- cbind(date_for_df, rep(2, nrow(clim_all0)), clim_all0[,p2_col], clim_max[,p2_col2], clim_min[,p2_col2])
plot3_clim_df <- cbind(date_for_df, rep(3, nrow(clim_all0)), clim_all0[,p3_col], clim_max[,p3_col2], clim_min[,p3_col2])
plot4_clim_df <- cbind(date_for_df, rep(4, nrow(clim_all0)), clim_all0[,p4_col], clim_max[,p4_col2], clim_min[,p4_col2])
plot5_clim_df <- cbind(date_for_df, rep(5, nrow(clim_all0)), clim_all0[,p5_col], clim_max[,p5_col2], clim_min[,p5_col2])
clim_colnames <- c("date", "plot", "temp_mean", "relhm_mean", "satdf_mean", "actvp_mean", "light_mean", "temp_max", "relhm_max", "light_max", "temp_min", "relhm_min", "light_min")
colnames(plot1_clim_df) <- colnames(plot2_clim_df) <- colnames(plot3_clim_df) <- colnames(plot4_clim_df) <- colnames(plot5_clim_df) <- clim_colnames
plot1_clim_df$plot <- as.factor(plot1_clim_df$plot)
plot2_clim_df$plot <- as.factor(plot2_clim_df$plot)
plot3_clim_df$plot <- as.factor(plot3_clim_df$plot)
plot4_clim_df$plot <- as.factor(plot4_clim_df$plot)
plot5_clim_df$plot <- as.factor(plot5_clim_df$plot)

# Add cumulative climate variables
plot1_clim_df <- calculate_cum_vals(plot1_clim_df, "temp_mean", cum_range = 2:14, cum_func = "sum")
plot1_clim_df <- calculate_cum_vals(plot1_clim_df, "relhm_mean", cum_range = 2:14, cum_func = "mean")
plot1_clim_df <- calculate_cum_vals(plot1_clim_df, "satdf_mean", cum_range = 2:14, cum_func = "sum")
plot1_clim_df <- calculate_cum_vals(plot1_clim_df, "actvp_mean", cum_range = 2:14, cum_func = "sum")
plot1_clim_df <- calculate_cum_vals(plot1_clim_df, "light_mean", cum_range = 2:14, cum_func = "sum")

plot2_clim_df <- calculate_cum_vals(plot2_clim_df, "temp_mean", cum_range = 2:14, cum_func = "sum")
plot2_clim_df <- calculate_cum_vals(plot2_clim_df, "relhm_mean", cum_range = 2:14, cum_func = "mean")
plot2_clim_df <- calculate_cum_vals(plot2_clim_df, "satdf_mean", cum_range = 2:14, cum_func = "sum")
plot2_clim_df <- calculate_cum_vals(plot2_clim_df, "actvp_mean", cum_range = 2:14, cum_func = "sum")
plot2_clim_df <- calculate_cum_vals(plot2_clim_df, "light_mean", cum_range = 2:14, cum_func = "sum")

plot3_clim_df <- calculate_cum_vals(plot3_clim_df, "temp_mean", cum_range = 2:14, cum_func = "sum")
plot3_clim_df <- calculate_cum_vals(plot3_clim_df, "relhm_mean", cum_range = 2:14, cum_func = "mean")
plot3_clim_df <- calculate_cum_vals(plot3_clim_df, "satdf_mean", cum_range = 2:14, cum_func = "sum")
plot3_clim_df <- calculate_cum_vals(plot3_clim_df, "actvp_mean", cum_range = 2:14, cum_func = "sum")
plot3_clim_df <- calculate_cum_vals(plot3_clim_df, "light_mean", cum_range = 2:14, cum_func = "sum")

plot4_clim_df <- calculate_cum_vals(plot4_clim_df, "temp_mean", cum_range = 2:14, cum_func = "sum")
plot4_clim_df <- calculate_cum_vals(plot4_clim_df, "relhm_mean", cum_range = 2:14, cum_func = "mean")
plot4_clim_df <- calculate_cum_vals(plot4_clim_df, "satdf_mean", cum_range = 2:14, cum_func = "sum")
plot4_clim_df <- calculate_cum_vals(plot4_clim_df, "actvp_mean", cum_range = 2:14, cum_func = "sum")
plot4_clim_df <- calculate_cum_vals(plot4_clim_df, "light_mean", cum_range = 2:14, cum_func = "sum")

plot5_clim_df <- calculate_cum_vals(plot5_clim_df, "temp_mean", cum_range = 2:14, cum_func = "sum")
plot5_clim_df <- calculate_cum_vals(plot5_clim_df, "relhm_mean", cum_range = 2:14, cum_func = "mean")
plot5_clim_df <- calculate_cum_vals(plot5_clim_df, "satdf_mean", cum_range = 2:14, cum_func = "sum")
plot5_clim_df <- calculate_cum_vals(plot5_clim_df, "actvp_mean", cum_range = 2:14, cum_func = "sum")
plot5_clim_df <- calculate_cum_vals(plot5_clim_df, "light_mean", cum_range = 2:14, cum_func = "sum")

# eDNA dataframe for each plot
edna_plot1 <- subset_samples(edna_all0, plot == 1)
edna_plot2 <- subset_samples(edna_all0, plot == 2)
edna_plot3 <- subset_samples(edna_all0, plot == 3)
edna_plot4 <- subset_samples(edna_all0, plot == 4)
edna_plot5 <- subset_samples(edna_all0, plot == 5)

plot1_edna_df <- adjust_edna_df_format(data.frame(sample_data(edna_plot1)), data.frame(otu_table(edna_plot1)), plot_number = 1)
plot2_edna_df <- adjust_edna_df_format(data.frame(sample_data(edna_plot2)), data.frame(otu_table(edna_plot2)), plot_number = 2)
plot3_edna_df <- adjust_edna_df_format(data.frame(sample_data(edna_plot3)), data.frame(otu_table(edna_plot3)), plot_number = 3)
plot4_edna_df <- adjust_edna_df_format(data.frame(sample_data(edna_plot4)), data.frame(otu_table(edna_plot4)), plot_number = 4)
plot5_edna_df <- adjust_edna_df_format(data.frame(sample_data(edna_plot5)), data.frame(otu_table(edna_plot5)), plot_number = 5)

# Check dimensions of dataframes
all(plot1_edna_df$date == plot1_clim_df$date)
all(plot2_edna_df$date == plot2_clim_df$date)
all(plot3_edna_df$date == plot3_clim_df$date)
all(plot4_edna_df$date == plot4_clim_df$date)
all(plot5_edna_df$date == plot5_clim_df$date)

# Combine all to make composite time series
clim_all <- rbind(plot1_clim_df, plot2_clim_df, plot3_clim_df, plot4_clim_df, plot5_clim_df)
edna_all <- rbind(plot1_edna_df, plot2_edna_df, plot3_edna_df, plot4_edna_df, plot5_edna_df)
edna_all$plot <- as.factor(edna_all$plot)
edna_tax <- tax_table(edna_all0)@.Data
edna_tax2 <- as.data.frame(edna_tax)

# Find column and row numbers to separate plots and data
colnames(clim_all); (clim_meta_coln <- 1:2); (clim_var_coln <- 3:ncol(clim_all))
colnames(edna_all); (edna_meta_coln <- 1:40); (edna_var_coln <- 41:ncol(edna_all))
clim_lib <- cbind(matrix(which(clim_all$date == "2017-05-01")+23, ncol = 1), matrix(which(clim_all$date == "2017-09-30")-8, ncol = 1))
edna_lib <- cbind(matrix(which(edna_all$date == "2017-05-01")+23, ncol = 1), matrix(which(edna_all$date == "2017-09-30")-8, ncol = 1))
all(edna_lib == clim_lib)
E_RANGE <- 1:14

#---------- Determine the best embedding dimension of all edna time series ----------#
edna_std <- as.numeric(scale(edna_all[,edna_var_coln[1]]))
(Eedna <- bestE_bidirect(edna_std, E_range = E_RANGE, lib = edna_lib, criteria = "rmse", show_fig = T))
for(i in edna_var_coln[-1]){
  (Eedna[(i-edna_var_coln[1]+1)] <- bestE_bidirect(as.numeric(scale(edna_all[,i])), E_range = E_RANGE,
                                                   lib = edna_lib, criteria = "rmse", show_fig = T))
}
names(Eedna) <- colnames(edna_all)[edna_var_coln]

#---------- Determine the best embedding dimension of all climate time series ----------#
clim_std <- as.numeric(scale(clim_all[,clim_var_coln[1]]))
(Eclim <- bestE_bidirect(clim_std, E_range = E_RANGE, lib = clim_lib, criteria = "rmse", show_fig = T))
for(i in clim_var_coln[-1]){
  (Eclim[(i-clim_var_coln[1]+1)] <- bestE_bidirect(as.numeric(scale(clim_all[,i])), E_range = E_RANGE,
                                                   lib = edna_lib, criteria = "rmse", show_fig = T))
}
names(Eclim) <- colnames(clim_all)[clim_var_coln]

#---------- Determine the best embedding dimension for all eDNA v.s. eDNA pairs ----------#
Eedna_edna <- matrix(NaN, nrow = length(edna_var_coln), ncol = length(edna_var_coln))
for(i in 1:length(edna_var_coln)){ # x = effect variable
  
  # Record the starting time
  start_time <- proc.time()[3]
  
  # Parallel computing the optimal embedding dimension
  Eedna_edna_tmp <- pforeach(j = 1:length(edna_var_coln), .c = c)({
    if(i == j){
      Eedna_edna_row_tmp <- Eedna[i]
    }else{
      edna_std_x <- as.numeric(scale(edna_all[,edna_var_coln[i]]))
      edna_std_y <- as.numeric(scale(edna_all[,edna_var_coln[j]]))
      Eedna_edna_row_tmp <- bestE_bidirect_block_lnlp(cbind(edna_std_x, edna_std_y),
                                                  E_range = E_RANGE,
                                                  lib = edna_lib,
                                                  criteria = "rmse",
                                                  show_fig = F)
      
    }
    names(Eedna_edna_row_tmp) <- colnames(edna_all[,edna_var_coln])[j]
    Eedna_edna_row_tmp
  })
  
  # Save the optimal E to Eedna_edna
  Eedna_edna[i,] <- Eedna_edna_tmp
  
  # Output message
  time_elapsed <- proc.time()[3] - start_time
  cat(sprintf("Best E_edna_edna cbind(%s,all) determined: %.02f sec\n", i, time_elapsed))
}

#---------- Determine the best embedding dimension for all eDNA xmap climate pairs ----------#
Eedna_clim <- matrix(NaN, nrow = length(edna_var_coln), ncol = length(clim_var_coln))
for(i in 1:length(edna_var_coln)){ # x = effect variable
  # Record the starting time
  start_time <- proc.time()[3]
  
  # Parallel computing the optimal embedding dimension
  Eedna_clim_tmp <- pforeach(j = 1:length(clim_var_coln), .c = c)({
    edna_std_x <- as.numeric(scale(edna_all[,edna_var_coln[i]]))
    clim_std_y <- as.numeric(scale(clim_all[,clim_var_coln[j]]))
    Eedna_clim_row_tmp <- bestE_bidirect_block_lnlp(cbind(edna_std_x, clim_std_y),
                                                    E_range = E_RANGE,
                                                    lib = edna_lib,
                                                    criteria = "rmse",
                                                    show_fig = F)
    
    names(Eedna_clim_row_tmp) <- colnames(clim_all[,clim_var_coln])[j]
    Eedna_clim_row_tmp
  })
  
  # Save the optimal E to Eedna_clim
  Eedna_clim[i,] <- Eedna_clim_tmp
  
  # Output message
  time_elapsed <- proc.time()[3] - start_time
  cat(sprintf("Best E_edna_clim, cbind(%s,all) determined: %.02f sec\n", i, time_elapsed))
}


#---------- Determine the best embedding dimension for all climate v.s. eDNA pairs ----------#
Eclim_edna <- matrix(NaN, nrow = length(clim_var_coln), ncol = length(edna_var_coln))
for(i in 1:length(clim_var_coln)){ # x = effect variable
  # Record the starting time
  start_time <- proc.time()[3]
  
  # Parallel computing the optimal embedding dimension
  Eclim_edna_tmp <- pforeach(j = 1:length(edna_var_coln), .c = c)({
    clim_std_x <- as.numeric(scale(clim_all[,clim_var_coln[i]]))
    edna_std_y <- as.numeric(scale(edna_all[,edna_var_coln[j]]))
    Eclim_edna_row_tmp <- bestE_bidirect_block_lnlp(cbind(clim_std_x, edna_std_y),
                                                    E_range = E_RANGE,
                                                    lib = edna_lib,
                                                    criteria = "rmse",
                                                    show_fig = F)
    
    names(Eclim_edna_row_tmp) <- colnames(edna_all[,edna_var_coln])[j]
    Eclim_edna_row_tmp
  })
  
  # Save the optimal E to Eedna_clim
  Eclim_edna[i,] <- Eclim_edna_tmp
  
  # Output message
  time_elapsed <- proc.time()[3] - start_time
  cat(sprintf("Best E_clim_edna, cbind(%s,all) determined: %.02f sec\n", i, time_elapsed))
}


#---------- Determine the best embedding dimension for all climate v.s. climate pairs ----------#
Eclim_clim <- matrix(NaN, nrow = length(clim_var_coln), ncol = length(clim_var_coln))
for(i in 1:length(clim_var_coln)){ # x = effect variable
  # Record the starting time
  start_time <- proc.time()[3]
  
  # Parallel computing the optimal embedding dimension
  Eclim_clim_tmp <- pforeach(j = 1:length(clim_var_coln), .c = c)({
    #for(j in 1:length(edna_var_coln)){ # j = causal variable
    if(i == j){
      Eclim_clim_row_tmp <- NaN
    }else{
      clim_std_x <- as.numeric(scale(clim_all[,clim_var_coln[i]]))
      clim_std_y <- as.numeric(scale(clim_all[,clim_var_coln[j]]))
      Eclim_clim_row_tmp <- bestE_bidirect_block_lnlp(cbind(clim_std_x, clim_std_y),
                                                      E_range = E_RANGE,
                                                      lib = edna_lib,
                                                      criteria = "rmse",
                                                      show_fig = F)
      
    }
    names(Eclim_clim_row_tmp) <- colnames(clim_all[,clim_var_coln])[j]
    Eclim_clim_row_tmp
  })
  
  # Save the optimal E to Eedna_clim
  Eclim_clim[i,] <- Eclim_clim_tmp
  
  # Output message
  time_elapsed <- proc.time()[3] - start_time
  cat(sprintf("Best E_clim_clim, cbind(%s,all) determined: %.02f sec\n", i, time_elapsed))
}

# Calculate nonlinearity
Tedna <- c(NULL)
total_cycle <- length(edna_var_coln)
cycle_n <- 1
for(i in 1:length(edna_var_coln)){
  # Set time
  start_time <- proc.time()[3]
  
  ts_tmp <- as.numeric(scale(edna_all[,edna_var_coln[i]]))
  ts_block <- make_block(ts_tmp, max_lag = Eedna[i], tau = 1)[,2:(Eedna[i]+1)]
  Tedna[i] <- bestT_bidirect_block_lnlp(ts_block, lib = edna_lib, tau = 1, tp_forward = 1,
                                        complete_case_only = F, criteria = "rmse", show_fig = F, save_stats = F)
  names(Tedna)[i] <- colnames(edna_all[,edna_var_coln])[i]
  
  # Show messeges
  time_used <- round(proc.time()[3] - start_time, digits = 2)
  cat("Cycle", cycle_n, "/", total_cycle, "finished;", time_used, "sec elapsed\n")
  cycle_n <- cycle_n + 1
}

# Save best E and theta
saveRDS(edna_tax2, sprintf("%s/eDNAtaxa_table.obj", output_folder00))
saveRDS(Eedna, sprintf("%s/BestE_eDNA_ts.obj", output_folder00))
saveRDS(Eclim, sprintf("%s/BestE_Clim_ts.obj", output_folder00))
saveRDS(Eedna_edna, sprintf("%s/BestE_edna_edna.obj", output_folder00))
saveRDS(Eedna_clim, sprintf("%s/BestE_edna_clim.obj", output_folder00))
saveRDS(Eclim_edna, sprintf("%s/BestE_clim_edna.obj", output_folder00))
saveRDS(Eclim_clim, sprintf("%s/BestE_clim_clim.obj", output_folder00))
saveRDS(Tedna, sprintf("%s/Tedna.obj", output_folder00))

# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/00_CompileAllDataOut.RData", output_folder00))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/00_SessionInfo_CompileAllData_%s.txt", substr(Sys.time(), 1, 10)))
