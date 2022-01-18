####
#### CERrice2017 All data analysis
#### Fig. Network properties analyzed with EDM
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
load("../../02_EcolComAnalysis/11_DiversityNetworkOut/11_DiversityNetworkOut.RData")

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.8
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.12
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.8
library(scales); packageVersion("scales") # 1.1.0, 2020.1.8
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.12
library(mgcv); packageVersion("mgcv") # 1.8.31, 2020.1.8
theme_set(theme_cowplot())

# Create output directory
fig_output <- "../00_RawFigs/02_Fig_EDMnet"
dir.create(fig_output)

# Extracting EDM results on connectance, interaction capacity and total diversity
smap_div2 <- cbind(regl_smap_div_mechnistic$model_output, regl_smap_div_mechnistic$smap_coefficients[,-1])
smap_div <- cbind(regl_smap_div$model_output, regl_smap_div$smap_coefficients[,-1])
smap_con <- cbind(regl_smap_con$model_output, regl_smap_con$smap_coefficients[,-1])
smap_itc <- cbind(regl_smap_ic$model_output, regl_smap_ic$smap_coefficients[,-1])

# Rename column names
colnames(smap_div)[which(colnames(smap_div) == "total_dna")] <- 
  colnames(smap_con)[which(colnames(smap_con) == "total_dna")] <-
  colnames(smap_itc)[which(colnames(smap_itc) == "total_dna")] <- "total_dna_effect"

colnames(smap_div)[which(colnames(smap_div) == "temp_mean")] <- 
  colnames(smap_con)[which(colnames(smap_con) == "temp_mean")] <-
  colnames(smap_itc)[which(colnames(smap_itc) == "temp_mean")] <- "temp_mean_effect"

smap_div$temp_mean <- smap_con$temp_mean <- smap_itc$temp_mean <- edna_prop$temp_mean
smap_div$total_dna <- smap_con$total_dna <- smap_itc$total_dna <- edna_prop$total_dna

# Predicted potential diversity and observed diversity
g01 <- ggplot(smap_div2, aes(y = obs, x = pred)) +
  geom_point(alpha = 0.5, size = 2, color = "deepskyblue3") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  xlab("Predicted potential diversity") + ylab("Observed diversity") +
  ggtitle("Mechanistic embedding by C, IC, and mean IS")

# Effect of temperature and total DNA on diversity
div_melt <- melt(smap_div[,c("time", "temp_mean_effect", "total_dna_effect")], id.vars = c("time"))
con_melt <- melt(smap_con[,c("time", "temp_mean_effect", "total_dna_effect")], id.vars = c("time"))
itc_melt <- melt(smap_itc[,c("time", "temp_mean_effect", "total_dna_effect")], id.vars = c("time"))

div_melt$effect_var <- "Diversity"
con_melt$effect_var <- "Connectance"
itc_melt$effect_var <- "Interaction capacity"

all_melt <- rbind(div_melt, con_melt, itc_melt)
all_melt$effect_var <- factor(all_melt$effect_var, levels = c("Interaction capacity",
                                                              "Connectance",
                                                              "Diversity"))

g02 <- ggplot(all_melt, aes(y = value, x = variable, facet = effect_var)) +
  geom_boxplot(width = 0.3, outlier.shape = NULL, outlier.size = 0, outlier.colour = "white", colour = "red3") +
  geom_jitter(size = 1.5, width = 0.1, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab(NULL) + ylab("Effect of variable") +
  facet_wrap(. ~ effect_var) + scale_color_d3() + panel_border() +
  scale_x_discrete(labels = c("Temperature", "Total DNA")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                   legend.position = "none")

# More detailed patterns
g03 <- ggplot(smap_con, aes(x = total_dna, y = total_dna_effect)) +
  geom_point(alpha = 0.5, size = 2, color = "darkblue") + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "red3", se = T, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Total DNA conc.\n(copies/ml)") + ylab("Effects on connectance") + ylim(-0.18, 0.2) +
  scale_x_continuous(breaks = seq(0, 1e+07, 4e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
summary(gam(smap_con$total_dna_effect ~ s(smap_con$total_dna)))

g04 <- ggplot(smap_itc, aes(x = total_dna, y = total_dna_effect)) +
  geom_point(alpha = 0.5, size = 2, color = "darkblue") + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "red3", se = T, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Total DNA conc.\n(copies/ml)") + ylab("Effects on interaction capacity") + ylim(-0.18, 0.2) +
  scale_x_continuous(breaks = seq(0, 1e+07, 4e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
summary(gam(smap_itc$total_dna_effect ~ s(smap_itc$total_dna)))



g05 <- ggplot(smap_con, aes(x = temp_mean, y = temp_mean_effect)) +
  geom_point(alpha = 0.5, size = 2, color = "darkblue") + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "red3", se = T, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) +
  ylab("Effects on connectance") + ylim(-0.18, 0.2) + xlim(15, 29)
summary(gam(smap_con$temp_mean_effect ~ s(smap_con$temp_mean)))

g06 <- ggplot(smap_itc, aes(x = temp_mean, y = temp_mean_effect)) +
  geom_point(alpha = 0.5, size = 2, color = "darkblue") + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "red3", se = T, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) +
  ylab("Effects on interaction capacity") + ylim(-0.18, 0.2) + xlim(15, 29)
summary(gam(smap_itc$temp_mean_effect ~ s(smap_itc$temp_mean)))


# Visualizing net effects of total DNA and mean air temperature on diversity
g07 <- ggplot(smap_div, aes(x = temp_mean, y = temp_mean_effect)) +
  geom_point(alpha = 0.5, size = 2, color = "darkblue") + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "red3", se = T, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab(expression(paste("Temperature ("*degree*C*")"))) +
  ylab("Net effects on diversity") + ylim(-0.18, 0.2) + xlim(15, 29)
summary(gam(smap_div$temp_mean_effect ~ s(smap_div$temp_mean)))

g08 <- ggplot(smap_div, aes(x = total_dna, y = total_dna_effect)) +
  geom_point(alpha = 0.5, size = 2, color = "darkblue") + theme(legend.position = "none") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", trans = "date", midpoint = as.numeric(mean(edna_prop$date))) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "red3", se = T, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Total DNA conc.\n(copies/ml)") + ylab("Net effects on diversity") + ylim(-0.18, 0.2) +
  scale_x_continuous(breaks = seq(0, 1e+07, 4e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
summary(gam(smap_div$total_dna_effect ~ s(smap_div$total_dna)))


g_all <- plot_grid(g05, g06, g07, g03, g04, g08,
                   align = "hv", ncol = 6,
                   labels = c("i","j","k","l","m","n"))


# Figure output
pdf(sprintf("%s/EDMFig_EffectofDNATemp.pdf", fig_output), width = 16, height = 4)
g_all; dev.off()

saveRDS(g_all, sprintf("%s/EDMFig_EffectofDNATemp.obj", fig_output))
