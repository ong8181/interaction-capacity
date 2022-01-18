####
#### CERrice2017 All data analysis
#### No. 5 Preparation for network reconstruction and stability calculation
####

# Load workspace
load("04_CompileSurResOut/04_CompileSurResOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder05 <- "05_NetworkPrepOut"
dir.create(output_folder05)

# Load library
library(Rcpp); packageVersion("Rcpp") # 1.0.3, 2020.1.21
library(rEDM); packageVersion("rEDM") # 0.7.5, 2020.1.21
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.21
library(igraph); packageVersion("igraph") # 1.2.4.2, 2020.1.21

# Compile S-map lib length
pred_length <- sum(apply(edna_lib, 1, function(x){x[2]-x[1]+1}))
pred_lib <- c(apply(edna_lib, 1, function(x) x[1]:x[2]))

# Select only representative climate variables
causal_climdna3 <- causal_climdna2[!is.na(match(causal_climdna2$effect_var, c("temp_mean"))),] #, "light_mean", "actvp_mean"
causal_dnaclim3 <- causal_dnaclim2[!is.na(match(causal_dnaclim2$cause_var, c("temp_mean"))),] #, "light_mean", "actvp_mean"

# Determine the number of subgraphs (are there any "isolated" species?)
dim(causal_dnaxdna2); colnames(causal_dnaxdna2)
length(unique(causal_dnaxdna2$effect_var)) # 1196
length(unique(causal_dnaxdna2$cause_var)) # 960
dna_igraph <- graph_from_data_frame(causal_dnaxdna2[,c("cause_var", "effect_var")])
comp_length <- c(NULL)
for(i in 1:nrow(edna_tax2)) comp_length[i] <- length(subcomponent(dna_igraph, i))
all(comp_length == nrow(edna_tax2)) # All vertexes are connected!


# Save workspace and ojcects
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder05, output_folder05))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder05, substr(Sys.time(), 1, 10)))
