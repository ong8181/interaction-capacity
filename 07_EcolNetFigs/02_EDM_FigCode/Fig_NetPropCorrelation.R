####
#### CERrice2017 All data analysis
#### Fig. Network properties
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace of original data
load("../../02_EcolComAnalysis/09_NetworkPropCCMOut/09_NetworkPropCCMOut.RData")

# Load saved object of SI (surrogate data)
edna_prop_surr <- readRDS("../../02_EcolComAnalysis/SI/SI_01_RandomizedTS/SI_08_NetworkPropCompile_SurrogateOut/edna_prop.obj")

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.17
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.11
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.17
library(scales); packageVersion("scales") # 1.1.0, 2020.1.17
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.11
library(mgcv); packageVersion("mgcv") # 1.8.31, 2020.1.17
theme_set(theme_cowplot())

# Create output directory
fig_output <- "../00_RawFigs/02_Fig_EDMnet"
dir.create(fig_output)

# General parameters
smooth_se <- T
reg_col <- "red3"

#-------------------- Data exploration (original data) --------------------#
g01 <- ggplot(edna_prop, aes(x = total_div, y = int_mean)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  geom_hline(yintercept = 0.03, linetype = 2) +
  xlab("ASV diversity") + ylab("Mean IS per link")
summary(gam(int_mean ~ s(total_div), data = edna_prop))

g02 <- ggplot(edna_prop, aes(x = total_div, y = int_n)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab("ASV diversity") + ylab("The number of interactions")
summary(gam(int_n ~ s(total_div), data = edna_prop))

g03 <- ggplot(edna_prop, aes(x = total_div, y = int_per_sp)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab("ASV diversity") + ylab("Mean interaction capacity")
summary(gam(int_per_sp ~ s(total_div), data = edna_prop))

g04 <- ggplot(edna_prop[edna_prop$total_div > 0,], aes(y = total_div, x = temp_mean)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("ASV diversity")
summary(gam(temp_mean ~ s(total_div), data = edna_prop[edna_prop$total_div > 0,]))

g05 <- ggplot(edna_prop, aes(y = connectance, x = total_div)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab("ASV diversity") + ylab("Connectance")
summary(gam(connectance ~ s(total_div), data = edna_prop))

g06 <- ggplot(edna_prop, aes(y = total_div, x = int_per_sp/0.03/connectance/2)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  xlab("Predicted potential diversity") + ylab("Observed diversity")

g07 <- ggplot(edna_prop, aes(y = int_per_sp, x = temp_mean)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("Mean interaction capacity")
summary(gam(int_per_sp ~ s(temp_mean), data = edna_prop))

g08 <- ggplot(edna_prop, aes(y = int_per_sp, x = total_dna)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab("Total DNA conc. (copies/ml)") + ylab("Mean interaction capacity") +
  scale_x_continuous(breaks = seq(0, 1e+07, 4e+06),
    label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
  summary(gam(int_per_sp ~ s(total_dna), data = edna_prop))
  
g09 <- ggplot(edna_prop, aes(y = connectance, x = temp_mean)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("Connectance")
summary(gam(connectance ~ s(temp_mean), data = edna_prop))

g10 <- ggplot(edna_prop, aes(y = connectance, x = total_dna)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab("Total DNA conc. (copies/ml)") + ylab("Connectance") +
  scale_x_continuous(breaks = seq(0, 1e+07, 4e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
summary(gam(connectance ~ s(total_dna), data = edna_prop))

g11 <- ggplot(edna_prop[edna_prop$total_dna > 0,], aes(y = total_dna, x = temp_mean)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("Total DNA conc. (copies/ml)")
summary(gam(total_dna ~ s(temp_mean), data = edna_prop[edna_prop$total_dna > 0,]))

g12 <- ggplot(edna_prop, aes(y = int_per_sp, x = connectance)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab("Connectance") + ylab("Mean interaction capacity")
summary(gam(int_per_sp ~ s(connectance), data = edna_prop))

g13 <- ggplot(edna_prop[edna_prop$total_dna > 0,], aes(y = total_div, x = total_dna)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab("Total DNA conc. (copies/ml)") + ylab("ASV diversity") + ylim(0,380) +
  scale_x_continuous(breaks = seq(0, 1e+07, 4e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
summary(gam(total_div ~ s(total_dna), data = edna_prop[edna_prop$total_dna > 0,]))

g14 <- ggplot(edna_prop, aes(x = total_div, y = mean_cv)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  #geom_smooth(method = "lm", color = "red3", se = smooth_se, size = 0.8) +
  xlab("ASV diversity") + ylab("C.V.")
summary(gam(mean_cv ~ s(total_div), data = edna_prop))

g15 <- ggplot(edna_prop, aes(x = total_div, y = int_max)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab("ASV diversity") + ylab("Maximum IS in a community") + ylim(0,2.5)
summary(gam(int_max ~ s(total_div), data = edna_prop))

g16 <- ggplot(edna_prop, aes(x = temp_mean, y = int_max)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  #geom_smooth(method = "gam", formula = y ~ s(x), color = "red3", se = smooth_se, size = 0.8) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("Maximum IS in a community") + ylim(0,2.5)
summary(gam(int_max ~ s(temp_mean), data = edna_prop)) # Non significant regression

g17 <- ggplot(edna_prop, aes(x = mean_cv, y = dynamic_stab)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = reg_col, se = smooth_se, size = 0.8) +
  xlab("C.V.") + ylab("Dynamic stability")
summary(gam(dynamic_stab ~ s(mean_cv), data = edna_prop))

g18 <- ggplot(edna_prop, aes(x = total_div, y = dynamic_stab)) +
  geom_point(alpha = 0.5, size = 2) + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  #geom_smooth(method = "gam", formula = y ~ s(x), color = "red3", se = smooth_se, size = 0.8) +
  xlab("ASV diversity") + ylab("Dynamic stability")
summary(gam(dynamic_stab ~ s(total_div), data = edna_prop)) # Non significant regression


#-------------------- Data exploration (surrogate data) --------------------#
s01 <- ggplot(edna_prop_surr, aes(x = total_div, y = int_mean)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  xlab("ASV diversity") + ylab("Mean IS per link") + ylim(0, 0.28)
summary(gam(int_mean ~ s(total_div), data = edna_prop_surr))

s02 <- ggplot(edna_prop_surr, aes(x = total_div, y = int_n)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  xlab("ASV diversity") + ylab("The number of interactions")
summary(gam(int_n ~ s(total_div), data = edna_prop_surr))

s03 <- ggplot(edna_prop_surr, aes(x = total_div, y = int_per_sp)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  xlab("ASV diversity") + ylab("Mean interaction capacity")
summary(gam(int_per_sp ~ s(total_div), data = edna_prop_surr))

s04 <- ggplot(edna_prop_surr[edna_prop_surr$total_div > 0,], aes(y = total_div, x = temp_mean)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("ASV diversity")
summary(gam(total_div ~ s(temp_mean), data = edna_prop_surr[edna_prop_surr$total_div > 0,]))

s05 <- ggplot(edna_prop_surr, aes(y = connectance, x = total_div)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  xlab("ASV diversity") + ylab("Connectance")
summary(gam(connectance ~ s(total_div), data = edna_prop_surr))

s06 <- ggplot(edna_prop_surr, aes(y = total_div, x = int_per_sp/0.03/connectance/2)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  xlab("Predicted potential diversity") + ylab("Observed diversity")

s07 <- ggplot(edna_prop_surr, aes(y = int_per_sp, x = temp_mean)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("Mean interaction capacity")
summary(gam(int_per_sp ~ s(temp_mean), data = edna_prop_surr))

s08 <- ggplot(edna_prop_surr, aes(y = int_per_sp, x = total_dna)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  xlab("Total DNA conc. (copies/ml)") + ylab("Mean interaction capacity") +
  scale_x_continuous(breaks = seq(0, 1e+07, 4e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
summary(gam(int_per_sp ~ s(total_dna), data = edna_prop_surr))

s09 <- ggplot(edna_prop_surr, aes(y = connectance, x = temp_mean)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("Connectance")
summary(gam(connectance ~ s(temp_mean), data = edna_prop_surr))

s10 <- ggplot(edna_prop_surr, aes(y = connectance, x = total_dna)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  xlab("Total DNA conc. (copies/ml)") + ylab("Connectance") +
  scale_x_continuous(breaks = seq(0, 1e+07, 4e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
summary(gam(connectance ~ s(total_dna), data = edna_prop_surr))

s11 <- ggplot(edna_prop_surr[edna_prop_surr$total_dna > 0,], aes(y = total_dna, x = temp_mean)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("Total DNA conc. (copies/ml)")
summary(gam(total_dna ~ s(temp_mean), data = edna_prop_surr[edna_prop_surr$total_div > 0,]))

s12 <- ggplot(edna_prop_surr, aes(y = int_per_sp, x = connectance)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  xlab("Connectance") + ylab("Mean interaction capacity")
summary(gam(int_per_sp ~ s(connectance), data = edna_prop_surr))

s13 <- ggplot(edna_prop_surr, aes(y = total_div, x = total_dna)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  xlab("Total DNA conc. (copies/ml)") + ylab("ASV diversity") +
  scale_x_continuous(breaks = seq(0, 1e+07, 4e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
summary(gam(total_div ~ s(total_dna), data = edna_prop_surr))

s14 <- ggplot(edna_prop_surr, aes(x = total_div, y = mean_cv)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  xlab("ASV diversity") + ylab("C.V.")
summary(gam(mean_cv ~ s(total_div), data = edna_prop_surr))

s15 <- ggplot(edna_prop_surr, aes(x = total_div, y = int_max)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "gray30", se = smooth_se, size = 0.8) +
  xlab("ASV diversity") + ylab("Maximum IS in a community")
summary(gam(int_max ~ s(total_div), data = edna_prop_surr))

s16 <- ggplot(edna_prop_surr, aes(x = temp_mean, y = int_max)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("Maximum IS in a community")
summary(gam(int_max ~ s(temp_mean), data = edna_prop_surr))

s17 <- ggplot(edna_prop_surr, aes(x = total_div, y = dynamic_stab)) +
  geom_point(alpha = 0.7, size = 2, color = "gray10", shape = 1) +
  xlab("ASV diversity") + ylab("Dynamic stability")
summary(gam(dynamic_stab ~ s(total_div), data = edna_prop_surr))


#---------- Generate combined figure objects ----------#
g_all <- plot_grid(g01, g03, g05, g14,
                   g07, g08, g09, g10, 
                   labels = "auto",
                   ncol = 4, align = "hv", hjust = -0.02)
s_all <- plot_grid(s01, s03, s05, s14,
                   s07, s08, s09, s10, 
                   labels = "auto",
                   ncol = 4, align = "hv", hjust = -0.02)

# Supplementary figures
g_all_SI <- plot_grid(g04, g13, g11, g02, g18, g17, ncol = 3, align = "hv", labels = "auto", hjust = -0.02)
g_all_nolabel <- plot_grid(g04, g13, g11, g14,g01, g03, g05, g18,g07, g08, g09, g10, ncol = 4, align = "hv")
s_all_nolabel <- plot_grid(s04, s13, s11, s14, s01, s03, s05, s17,s07, s08, s09, s10, ncol = 4, align = "hv")
s_subset <- plot_grid(s01, s07, s08, s10, labels = "auto", ncol = 2, align = "hv", hjust = -0.02)
gs_max <- plot_grid(g15, g16, s15, s16, labels = "auto", ncol = 2, align = "hv", hjust = -0.02)


#-------------------- Save all figures --------------------#
# Save figures
pdf(sprintf("%s/NetworkProperties.pdf", fig_output), width = 16, height = 12)
g_all; dev.off()
pdf(sprintf("%s/NetworkProperties_RandomSurrogate.pdf", fig_output), width = 16, height = 12)
s_all; dev.off()

saveRDS(g_all, sprintf("%s/NetworkProperties.obj", fig_output))
saveRDS(g_all_SI, sprintf("%s/NetworkProperties_SI.obj", fig_output))
saveRDS(s_all, sprintf("%s/NetworkProperties_RandomSurrogate.obj", fig_output))
saveRDS(s_subset, sprintf("%s/NetworkProperties_RandomSurrogate2.obj", fig_output))
saveRDS(gs_max, sprintf("%s/NetworkProperties_IntMax.obj", fig_output))




#-------------------- Comparison between original v.s. randomized surrogate --------------------#
gs01 <- ggplot(edna_prop, aes(x = total_div, y = int_mean)) +
  geom_point(size = 2, color = "royalblue", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "royalblue", se = FALSE, size = 0.8) +
  geom_point(aes(x = edna_prop_surr$total_div, y = edna_prop_surr$int_mean), fill = "gray90", color = "gray10", size = 2, shape = 21) +
  geom_hline(yintercept = 0.03, linetype = 2) +
  xlab("ASV diversity") + ylab("Mean IS per link")

gs02 <- ggplot(edna_prop, aes(x = total_div, y = int_n)) +
  geom_point(size = 2, color = "royalblue", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "royalblue", se = FALSE, size = 0.8) +
  geom_point(aes(x = edna_prop_surr$total_div, y = edna_prop_surr$int_n), fill = "gray90", color = "gray10", size = 2, shape = 21) +
  xlab("ASV diversity") + ylab("The number of interactions")

gs03 <- ggplot(edna_prop, aes(x = total_div, y = int_per_sp)) +
  geom_point(size = 2, color = "royalblue", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "royalblue", se = FALSE, size = 0.8) +
  geom_point(aes(x = edna_prop_surr$total_div, y = edna_prop_surr$int_per_sp), fill = "gray90", color = "gray10", size = 2, shape = 21) +
  xlab("ASV diversity") + ylab("Mean interaction capacity")

gs05 <- ggplot(edna_prop, aes(y = connectance, x = total_div)) +
  geom_point(size = 2, color = "royalblue", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "royalblue", se = FALSE, size = 0.8) +
  geom_point(aes(x = edna_prop_surr$total_div, y = edna_prop_surr$connectance), fill = "gray90", color = "gray10", size = 2, shape = 21) +
  xlab("ASV diversity") + ylab("Connectance") + ylim(0, 0.05)

gs07 <- ggplot(edna_prop, aes(y = int_per_sp, x = temp_mean)) +
  geom_point(size = 2, color = "royalblue", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "royalblue", se = FALSE, size = 0.8) +
  geom_point(aes(x = edna_prop_surr$temp_mean, y = edna_prop_surr$int_per_sp), fill = "gray90", color = "gray10", size = 2, shape = 21) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("Mean interaction capacity")

gs08 <- ggplot(edna_prop, aes(y = int_per_sp, x = total_dna)) +
  geom_point(size = 2, color = "royalblue", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "royalblue", se = FALSE, size = 0.8) +
  geom_point(aes(x = edna_prop_surr$total_dna, y = edna_prop_surr$int_per_sp), fill = "gray90", color = "gray10", size = 2, shape = 21) +
  xlab("Total DNA conc. (copies/ml)") + ylab("Mean interaction capacity") +
  scale_x_continuous(breaks = seq(0, 1e+07, 4e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})

gs09 <- ggplot(edna_prop, aes(y = connectance, x = temp_mean)) +
  geom_point(size = 2, color = "royalblue", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "royalblue", se = FALSE, size = 0.8) +
  geom_point(aes(x = edna_prop_surr$temp_mean, y = edna_prop_surr$connectance), fill = "gray90", color = "gray10", size = 2, shape = 21) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) + ylab("Connectance")

gs10 <- ggplot(edna_prop, aes(y = connectance, x = total_dna)) +
  geom_point(size = 2, color = "royalblue", shape = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "royalblue", se = FALSE, size = 0.8) +
  geom_point(aes(x = edna_prop_surr$total_dna, y = edna_prop_surr$connectance), fill = "gray90", color = "gray10", size = 2, shape = 21) +
  xlab("Total DNA conc. (copies/ml)") + ylab("Connectance") +
  scale_x_continuous(breaks = seq(0, 1e+07, 4e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})

gs_selected <- plot_grid(gs01, gs07, gs08, gs03, gs09, gs10, ncol = 3, labels = "auto", align = "hv", hjust = -0.02)
gs_selected_nolabel <- plot_grid(gs01, gs07, gs08, gs03, gs09, gs10, ncol = 3, align = "hv", hjust = -0.02)
pdf(sprintf("%s/NetworkProperties_Compare.pdf", fig_output), width = 10, height = 4)
gs_selected; dev.off()

saveRDS(gs_selected, sprintf("%s/NetworkProperties_Compare.obj", fig_output))
