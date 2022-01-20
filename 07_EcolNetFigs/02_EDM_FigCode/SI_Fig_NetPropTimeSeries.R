####
#### CERrice2017 All data analysis
#### Fig. Time series of network properties
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
#load("../../02_EcolComAnalysis/12_TaxaSpecificISOut/12_TaxaSpecificISOut.RData")
# Load workspace of original data
load("../../02_EcolComAnalysis/09_NetworkPropCCMOut/09_NetworkPropCCMOut.RData")

# Create output directory
fig_output <- "../00_RawFigs/02_Fig_EDMnet"
dir.create(fig_output)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.8
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.8
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.11
library(scales); packageVersion("scales") # 1.1.0, 2020.1.8
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.12
theme_set(theme_cowplot())

# Animation package
library(ggimage); packageVersion("ggimage") # 0.2.5, 2020.1.8
library(magick); packageVersion("magick") # 2.2, 2020.1.8

# Load examples of time-varying networks
network_output_folder2 <- "01_NetworkFigs_v2"
net01 <- image_scale(image_read(sprintf("%s/Network_63.jpg", network_output_folder2)), "400x400")
net02 <- image_scale(image_read(sprintf("%s/Network_95.jpg", network_output_folder2)), "400x400")
net03 <- image_scale(image_read(sprintf("%s/Network_142.jpg", network_output_folder2)), "400x400")

# Illustrations of time-varying network properties
edna_df <- edna_prop[!is.na(edna_prop$plot),]

# Dynamic stability
g01 <- ggplot(edna_df, aes(x = date, y = dynamic_stab, group = plot, shape = plot, fill = plot, color = plot))
g01 <- g01 + geom_point() + geom_line() + geom_hline(yintercept = 1, linetype = 2)
g01 <- g01 + scale_color_startrek() + scale_fill_startrek() + scale_shape_manual(values = c(21:25))
g01 <- g01 + xlim(as.Date("2017-05-22"), as.Date("2017-09-23")) + xlab(NULL) + ylab("Dynamic stability")
g01 <- g01 + ylim(0.7,8.5)

# Mean C.V.
g02 <- ggplot(edna_df, aes(x = date, y = mean_cv, group = plot, shape = plot, fill = plot, color = plot))
g02 <- g02 + geom_point() + geom_line()
g02 <- g02 + scale_color_startrek() + scale_fill_startrek() + scale_shape_manual(values = c(21:25))
g02 <- g02 + xlim(as.Date("2017-05-22"), as.Date("2017-09-23")) + ylab("Mean C.V.")
g02 <- g02 + ylim(1.2,2.1) + xlab(NULL)

# Mean interaction strength
g03 <- ggplot(edna_df, aes(x = date, y = int_mean, group = plot, shape = plot, fill = plot, color = plot))
g03 <- g03 + geom_point() + geom_line()
g03 <- g03 + scale_color_startrek() + scale_fill_startrek() + scale_shape_manual(values = c(21:25))
g03 <- g03 + xlim(as.Date("2017-05-22"), as.Date("2017-09-23")) + ylab("Mean interaction strength")
g03 <- g03 + xlab(NULL)


# Combine figures
net_figs <- plot_grid(as.ggplot(net01) + theme(plot.margin = unit(c(0,0.2,0,0.5),"cm")),
                      as.ggplot(net02) + theme(plot.margin = unit(c(0,0.2,0,0.5),"cm")),
                      as.ggplot(net03) + theme(plot.margin = unit(c(0,0.2,0,0.5),"cm")),
                      ncol = 3, labels = c("b", NULL, NULL), hjust = 0.02)
prop_figs <- plot_grid(g01, ncol = 1, labels = c("c"),
                       align = "hv", hjust = 0.02)
g_all <- plot_grid(net_figs,
                   prop_figs, rel_heights = c(1, 1), ncol = 1) +
  theme(plot.margin = unit(c(0,0,0,1),"cm"))

# Save figure
pdf(sprintf("%s/EDMFig_TimeVaryingNetProps.pdf", fig_output), width = 10, height = 8)
g_all; dev.off()

saveRDS(g_all, sprintf("%s/EDMFig_TimeVaryingNetProps.obj", fig_output))

