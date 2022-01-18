####
#### CERrice2017 All data analysis
#### Visualize the relationship between network properties
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace and objects
load("../../02_EcolComAnalysis/10_NetworkPropSmapOut/10_NetworkPropSmapOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
fig_output <- "../00_RawFigs/02_Fig_EDMnet"
dir.create(fig_output)

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.30
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.12
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.12
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.30
library(scales); packageVersion("scales") # 1.1.0, 2020.1.30
theme_set(theme_cowplot())

source("../../02_EcolComAnalysis/functions/F03_HelperFunctions.R")

# Summarize S-map coefficients
# Make property matrix
prop_names <- colnames(edna_prop)[3:18]
prop_mat_mean <- matrix(0, ncol = length(prop_names), nrow = length(prop_names))
prop_mat_med <- matrix(0, ncol = length(prop_names), nrow = length(prop_names))

# Collect names
n_total_cause <- length(prop_smap_list)
cause_prop_names <- sapply(strsplit(names(prop_smap_list), "_<=_"), `[`, 1)
effect_prop_names <- sapply(strsplit(names(prop_smap_list), "_<=_"), `[`, 2)
colnames(prop_mat_mean) <- rownames(prop_mat_mean) <- prop_names
colnames(prop_mat_med) <- rownames(prop_mat_med) <- prop_names

# Assign S-map values
for(i in 1:n_total_cause){
  # Calculate mean S-map coefficient values
  mean_smap_coef_tmp <- mean(prop_smap_list[[i]][,effect_prop_names[i]], na.rm = T)
  median_smap_coef_tmp <- median(prop_smap_list[[i]][,effect_prop_names[i]], na.rm = T)
  # Assign value
  prop_mat_mean[effect_prop_names[i], cause_prop_names[i]] <- mean_smap_coef_tmp
  prop_mat_med[effect_prop_names[i], cause_prop_names[i]] <- median_smap_coef_tmp
}

# Digitalize results
prop_mat_mean[prop_mat_mean != 0] <- 1
prop_mat_med[prop_mat_med != 0] <- 1

# Select visualized properties
prop_mat_mean <- prop_mat_mean[c(1:3,5,7:16),c(1:3,5,7:16)]
prop_mat_med <- prop_mat_med[c(1:3,5,7:16),c(1:3,5,7:16)]

# Visualize the results
mean_melt <- melt(prop_mat_mean[,1:(dim(prop_mat_mean)[2]-1)])
med_melt <- melt(prop_mat_med[,1:(dim(prop_mat_mean)[2]-1)])
colnames(mean_melt)[1:2] <- colnames(med_melt)[1:2] <- c("Cause", "Effect")
names(prop_smap_list)[4]; head(mean_melt)

# Prepare new labels
new_prop_labels <- c("Total DNA", "Total diversity", "Dynamic stability", "Mean C.V.",
                     "Connectance", "Interaction capacity", "Weak interaction index",
                     "S.D. of IS", "Minimum IS", "Maximum IS", "Median IS", "Mean IS",
                     "No. of interactions", "Mean temperature")
new_prop_labels2 <- c("Total DNA", "Total diversity", "Dynamic stability", "Mean C.V.",
                      "Connectance", "Interaction capacity", "Weak interaction index",
                      "S.D. of IS", "Minimum IS", "Maximum IS", "Median IS", "Mean IS",
                      "No. of interactions")

# Visualize results
#ordered(mean_melt$Effect, levels = rev(levels(mean_melt$Effect)))
g1 <- ggplot(mean_melt, aes(x = Cause, y = ordered(Effect, levels = rev(levels(Effect))), fill = value))
g1 <- g1 + geom_tile(color = "gray20") + scale_fill_gradient(low = "white", high = "gray40")
g1 <- g1 + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g1 <- g1 + scale_x_discrete(labels = new_prop_labels) + scale_y_discrete(labels = rev(new_prop_labels2))
g1 <- g1 + xlab("Causal property") + ylab("Effect property")

g2 <- ggplot(med_melt, aes(x = Cause, y = ordered(Effect, levels = rev(levels(Effect))), fill = value))
g2 <- g2 + geom_tile(color = "gray20") + scale_fill_gradient(low = "white", high = "gray40")
g2 <- g2 + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g2 <- g2 + scale_x_discrete(labels = new_prop_labels) + scale_y_discrete(labels = rev(new_prop_labels2))
g2 <- g2 + xlab("Causal property") + ylab("Effect property")


# Save figures
pdf(sprintf("%s/NetworkPropCCMmatrix_mean.pdf", fig_output), width = 8, height = 7.5)
g1; dev.off()
pdf(sprintf("%s/NetworkPropCCMmatrix_median.pdf", fig_output), width = 8, height = 7.5)
g2; dev.off()

saveRDS(g1, sprintf("%s/NetworkPropCCMmatrix_mean.obj", fig_output))
saveRDS(g2, sprintf("%s/NetworkPropCCMmatrix_median.obj", fig_output))
