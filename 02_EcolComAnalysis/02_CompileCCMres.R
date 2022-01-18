####
#### CERrice2017 All data analysis
#### No. 2 Compile CCM results
####

# Load workspace
load("00_CompileAllDataOut/00_CompileAllDataOut.RData")
output_folder01 <- "01_CCMEcolComClimOut"

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder02 <- "02_CompileCCMresOut"
dir.create(output_folder02)

# Load library and functions
library(rEDM); packageVersion("rEDM") # 0.7.5, 2020.1.7
library(pforeach); packageVersion("pforeach") # 1.3, 2019.10.26
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.7
source("functions/CcmExact_v5_rEDM_0.7.4.R")
source("functions/BidirectSimplex_rEDM_0.7.4.R")
source("functions/BidirectBlocklnlp_rEDM_0.7.4.R")
source("functions/F02_HelperFunctions.R")
source("functions/Parallel_Sur_95CI.R")

# Load and combine CCM results
## Specify file path
folder_DNAxDNA <- sprintf("%s/CCMDNAxDNA_RawRes", output_folder01)
folder_ClimDNA <- sprintf("%s/CCMClimxDNA_RawRes", output_folder01)
folder_DNAClim <- sprintf("%s/CCMDNAxClim_RawRes", output_folder01)

files_DNAxDNA <- list.files(folder_DNAxDNA, pattern=".obj", full.names = T)
files_ClimDNA <- list.files(folder_ClimDNA, pattern=".obj", full.names = T)
files_DNAClim <- list.files(folder_DNAClim, pattern=".obj", full.names = T)

## Combine CCM results
ccmres_DNAxDNA <- ccmres_ClimDNA <- ccmres_DNAClim <- data.frame(NULL)
for(i in 1:length(files_DNAxDNA)){
  ccmres_tmp1 <- readRDS(files_DNAxDNA[i])
  ccmres_tmp1$effect_var <- as.character(ccmres_tmp1$effect_var)
  ccmres_tmp1$cause_var <- as.character(ccmres_tmp1$cause_var)
  ccmres_tmp1$criteria <- as.character(ccmres_tmp1$criteria)
  ccmres_DNAxDNA <- bind_rows(ccmres_DNAxDNA, ccmres_tmp1)
}
for(i in 1:length(files_ClimDNA)){
  ccmres_tmp2 <- readRDS(files_ClimDNA[i])
  ccmres_tmp2$effect_var <- as.character(ccmres_tmp2$effect_var)
  ccmres_tmp2$cause_var <- as.character(ccmres_tmp2$cause_var)
  ccmres_tmp2$criteria <- as.character(ccmres_tmp2$criteria)
  ccmres_ClimDNA <- bind_rows(ccmres_ClimDNA, ccmres_tmp2)
}
for(i in 1:length(files_DNAClim)){
  ccmres_tmp3 <- readRDS(files_DNAClim[i])
  ccmres_tmp3$effect_var <- as.character(ccmres_tmp3$effect_var)
  ccmres_tmp3$cause_var <- as.character(ccmres_tmp3$cause_var)
  ccmres_tmp3$criteria <- as.character(ccmres_tmp3$criteria)
  ccmres_DNAClim <- bind_rows(ccmres_DNAClim, ccmres_tmp3)
}

# Select causal variables
causal_dnaxdna <- extract_causal_vars(ccmres_DNAxDNA, causal_index = "log_d_rmse", select_strongest = F,
                                      select_tp_mode = "specific_tp", ccm_tp = -1, causal_threshold = 0, check_naive = T)
causal_climdna0 <- extract_causal_vars(ccmres_ClimDNA, causal_index = "log_d_rmse", select_strongest = F,
                                       select_tp_mode = "specific_tp", ccm_tp = -1, causal_threshold = 0, check_naive = T)
causal_dnaclim0 <- extract_causal_vars(ccmres_DNAClim, causal_index = "log_d_rmse", select_strongest = F,
                                       select_tp_mode = "specific_tp", ccm_tp = -1, causal_threshold = 0, check_naive = T)

# Rename causal objects and remove cumulative climate variables
select_climdna_id <- select_dnaclim_id <- c(NULL)
for(select_clim_var_i in clim_colnames[3:13]){
  select_climdna_id <- c(select_climdna_id, which(causal_climdna0$effect_var == select_clim_var_i))
  select_dnaclim_id <- c(select_dnaclim_id, which(causal_dnaclim0$cause_var == select_clim_var_i))
}
causal_climdna <- causal_climdna0[select_climdna_id,]
causal_dnaclim <- causal_dnaclim0[select_dnaclim_id,]

# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s_Full.RData", output_folder02, output_folder02))

rm(folder_DNAxDNA)
rm(folder_ClimDNA)
rm(folder_DNAClim)
rm(ccmres_DNAxDNA)
rm(ccmres_ClimDNA)
rm(ccmres_DNAClim)
rm(causal_climdna0)
rm(causal_dnaclim0)

# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s_Minimum.RData", output_folder02, output_folder02))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder02, substr(Sys.time(), 1, 10)))
