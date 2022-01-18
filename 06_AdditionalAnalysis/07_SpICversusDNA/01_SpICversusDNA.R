####
#### CERrice2017 All data analysis
#### Species-level interaction capacity and phylogeny
####

# Library
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.12.6
library(ggsci); packageVersion("ggsci") # 2.9, 2021.12.7
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.12.6
library(scales); packageVersion("scales") # 1.1.1, 2021.12.8
theme_set(theme_cowplot())

# Load workspace and objects
# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the data file
load("../../02_EcolComAnalysis/12_TaxaSpecificISOut/12_TaxaSpecificISOut.RData")

# Generate output folder
od <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(od, end = -3), "Out")); rm(od)
dir.create(output_folder)


# ----------------------------------------------------- # 
# Compile data
# ----------------------------------------------------- # 
# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the data file
int_mat <- readRDS("../../02_EcolComAnalysis/07_RegularizedSmapOut/interaction_matrix_nonzeroridge.obj")
edna_connect <- edna_intn <- edna_ic <- edna_intmax <- edna_all[,edna_var_coln]
edna_connect[] <- edna_intn[] <- edna_ic[] <- edna_intmax[] <- NA # Set all elements as zero
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
      
      # Calculate connectance
      edna_connect[time_i, colnames(int_mat_tmp)] <- edna_intn[time_i, colnames(int_mat_tmp)]/(dim(int_mat_tmp)[1]*2-2)
      
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
edna_connect <- edna_connect[valid_row_id,]

edna_prop$int_max_per_com <- edna_prop$int_per_link <-
  edna_prop$int_cor <- edna_prop$int_cor_abs <- NaN
edna_prop$int_max_per_com[valid_row_id] <- apply(edna_ic, 1, function(x){if(any(!is.na(x))){max(x, na.rm = T)}else{NA}})
edna_prop$int_per_link[valid_row_id] <- rowSums(edna_ic, na.rm = T)/rowSums(edna_intn, na.rm = T)
edna_prop$int_cor[valid_row_id] <- edna_cor[valid_row_id]
edna_prop$int_cor_abs[valid_row_id] <- edna_cor_abs[valid_row_id]


# ----------------------------------------------------- # 
# Visualize patterns
# ----------------------------------------------------- # 
# Prepare sample names
taxa_id <- colnames(edna_all[,edna_var_coln])
edna_ic$Sample_Name2 <- edna_all$Sample_Name2[valid_row_id]
edna_all2 <- edna_all[valid_row_id,edna_var_coln]
edna_all2$Sample_Name2 <- edna_all$Sample_Name2[valid_row_id]
all(edna_all2$Sample_Name2 == edna_ic$Sample_Name2)

# Prepare edna_prop2
edna_prop2 <- edna_prop[valid_row_id, c("date", "total_div")]
edna_prop2$Sample_Name2 <- edna_all$Sample_Name2[valid_row_id]

# Make the data.frame long
edna_ic_long <- pivot_longer(edna_ic, cols = -c(Sample_Name2), names_to = "taxa", values_to = "capacity")
edna_all_long <- pivot_longer(edna_all2, cols = -c(Sample_Name2), names_to = "taxa", values_to = "abundance")
edna_all_long$superkingdom <- as.character(edna_tax2[match(edna_all_long$taxa, rownames(edna_tax2)), "superkingdom"])
edna_all_long$superkingdom[edna_all_long$superkingdom == ""] <- "Undetermined"
edna_all_long$order <- as.character(edna_tax2[match(edna_all_long$taxa, rownames(edna_tax2)), "order"])
edna_all_long$family <- as.character(edna_tax2[match(edna_all_long$taxa, rownames(edna_tax2)), "family"])
edna_all_long$genus <- as.character(edna_tax2[match(edna_all_long$taxa, rownames(edna_tax2)), "genus"])

# Add abundance information
edna_taxa_abun <- colSums(edna_all[,edna_var_coln], na.rm = T)
edna_all_long$total_dna <- edna_taxa_abun[match(edna_all_long$taxa, names(edna_taxa_abun))]

# Add sample-specific ID
edna_ic_long$sample_taxa <- paste0(edna_ic_long$Sample_Name2, "-", edna_ic_long$taxa)
edna_all_long$sample_taxa <- paste0(edna_all_long$Sample_Name2, "-", edna_all_long$taxa)
all(edna_ic_long$sample_taxa == edna_all_long$sample_taxa)

edna_con_abun <- edna_all_long %>%
  mutate(capacity = edna_ic_long$capacity) %>%
  na.omit() %>%
  filter(total_dna > 0 & capacity > 0) %>%
  filter(superkingdom == "Bacteria" | superkingdom == "Eukaryota")
dim(edna_con_abun)


# Visualize with ggplot2
g1 <- ggplot(edna_con_abun, aes(x = 1/total_dna, y = capacity, color = superkingdom)) +
  geom_point(alpha = 0.1, color = "gray", size = 0) +
  stat_smooth(method = "lm", aes(color = superkingdom)) +
  scale_x_continuous(trans = "log10",
                     breaks = c(1e-07, 1e-05, 1e-03), limits = c(1e-06, 2e-02),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))}) +
  scale_y_log10() +
  panel_border() +
  facet_wrap(~ superkingdom, ncol = 2) + 
  xlab("An index of body size (1/total DNA copy numbers)") + ylab("Interaction capacity") +
  scale_color_manual(values = pal_startrek(palette = "uniform")(4)[c(2,3)]) +
  coord_cartesian(ylim = c(0.03, 0.2)) + 
  NULL
g1$layers[[1]] <- NULL

g2 <- ggplot(edna_con_abun, aes(x = 1/total_dna, y = capacity, color = superkingdom)) +
  geom_point(alpha = 0.1, color = "gray") +
  stat_smooth(method = "lm", aes(color = superkingdom)) +
  panel_border() +
  facet_wrap(~ superkingdom, ncol = 2) + 
  xlab("An index of body size (1/total DNA copy numbers)") + ylab("Interaction capacity") +
  scale_color_manual(values = pal_startrek(palette = "uniform")(4)[c(2,3)]) +
  scale_x_continuous(trans = "log10",
                     breaks = c(1e-07, 1e-05, 1e-03), limits = c(1e-06, 2e-02),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))}) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e-05, 1e-02, 1e+01), limits = c(0.8e-05, 2e+01),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))}) +
  NULL
g2

g3 <- ggplot(edna_con_abun, aes(x = capacity, fill = superkingdom)) +
  geom_histogram(position = "identity", alpha = 0.6) +
  theme(legend.position = c(0.8, 0.8)) +
  scale_x_log10() +
  xlab("Capacity") + ylab("Count") +
  scale_fill_startrek() +
  ggtitle("Taxa-specific interaction capacity") +
  NULL

g_all <- plot_grid(g2, g1, ncol = 1, labels = c("(d)", "(e)"), align = "hv", axis = "lrbt")

ggsave(file = sprintf("%s/SpLevel_Capacity.pdf", output_folder),
       plot = g_all,
       width = 12, height = 8)
ggsave(file = sprintf("%s/SpLevel_Capacity.jpg", output_folder),
       plot = g_all,
       width = 12, height = 8)

saveRDS(g_all, file = sprintf("%s/SpLevel_Capacity.obj", output_folder))
