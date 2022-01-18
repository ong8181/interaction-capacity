####
#### CER eDNA study
#### FigCode: Re-formatting all figures
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the figure file(s)!

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.7.30
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.7.30
library(scales); packageVersion("scales") # 1.1.1, 2021.7.30
library(ggimage); packageVersion("ggimage") # 0.3.0, 2021.12.8
library(magick); packageVersion("magick") # 2.7.3, 2021.12.8
theme_set(theme_cowplot())

# Create output folder
dir.create("../00_ReformatFigs")

#-------------------- Load saved figure objects --------------------#
Image_ExpDesign <- "../00_RawFigs/00_ExpImages/ExpDesign_v3.jpg"
Fig_eDNAPattern <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_AllBar.obj")
Fig_AllDivTS <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_AllDivTS.obj")
Fig_Network <- readRDS("../00_RawFigs/02_Fig_EDMnet/EDMFig_NetworkTimeAvr_v2.obj")
Fig_NetworkLegend <- readRDS("../00_RawFigs/02_Fig_EDMnet/EDMFig_NetworkLegend.obj")
Fig_NetworkPropCor <- readRDS("../00_RawFigs/02_Fig_EDMnet/NetworkProperties.obj")
Fig_HowDivDet <- readRDS("../00_RawFigs/02_Fig_EDMnet/EDMFig_EffectofDNATemp.obj")
Image_SummaryFig <- "../00_RawFigs/04_SummaryImage/Fig4_SummaryImage.jpg"


#-------------------- Reformat figures --------------------#
# Figure 1
img <- Image_ExpDesign %>% image_read()
Fig_ExpDesign <- ggdraw() + draw_image(img)
Fig_eDNAPattern <- Fig_eDNAPattern + theme(legend.position = c(.1, .7))
Fig_01 <- plot_grid(plot_grid(Fig_ExpDesign, labels = "(a)", label_x = 0, label_y = .95),
                    plot_grid(Fig_eDNAPattern, Fig_AllDivTS,
                              nrow = 1, rel_widths = c(0.9, 1), labels = c("(b)", "(c)")),
                    ncol = 1, rel_heights = c(1, 0.6))

# Figure 2
ggsave(filename = "../00_ReformatFigs/Fig_02_network.jpg",
       plot = Fig_Network, dpi = 100, width = 26.5, height = 25)
ggsave(filename = "../00_ReformatFigs/Fig_02_legend.jpg",
       plot = as.ggplot(Fig_NetworkLegend), dpi = 100, width = 9, height = 1)

Fig_Network_img <- image_read("../00_ReformatFigs/Fig_02_network.jpg")
Fig_NetworkLegend_img <- image_read("../00_ReformatFigs/Fig_02_legend.jpg")
Fig_Network2 <- ggdraw() + draw_image(Fig_Network_img)
Fig_NetworkLegend2 <- ggdraw() + draw_image(Fig_NetworkLegend_img)
Fig_02 <- plot_grid(Fig_Network2, Fig_NetworkLegend2,
                    ncol = 1, rel_heights = c(1, 0.11))

# Figure 3
Fig_03_title01 <- ggdraw() + draw_label("(a-h) Correlation analysis", x = 0.01, hjust = 0, size = 20)
Fig_03_title02 <- ggdraw() + draw_label("(i-n) Causal influences quantified by EDM", x = 0.01, hjust = 0, size = 20)
Fig_03 <- plot_grid(Fig_03_title01,
                    Fig_NetworkPropCor,
                    Fig_03_title02,
                    Fig_HowDivDet,
                    nrow = 4, #axis = "hv",
                    rel_heights = c(0.1, 1, 0.1, 0.6), labels = NULL)

# Figure 4
img2 <- Image_SummaryFig %>% image_read()
Fig_04 <- ggdraw() + draw_image(img2)

#-------------------- Save figures --------------------#
# Figure 1
pdf("../00_ReformatFigs/Fig_01.pdf", width = 12, height = 12)
Fig_01; dev.off()

# Figure 2
ggsave(filename = "../00_ReformatFigs/Fig_02.jpg",
       plot = Fig_02,
       dpi = 100, width = 26, height = 25)

 # Figure 3
pdf("../00_ReformatFigs/Fig_03.pdf", width = 14, height = 12)
Fig_03; dev.off()

# Figure 4
pdf("../00_ReformatFigs/Fig_04.pdf", width = 10, height = 12)
Fig_04; dev.off()

