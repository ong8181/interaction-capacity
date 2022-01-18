####
#### CERrice2017 All data analysis
#### Fig. Network properties analyzed with EDM
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
load("../../03_MetaAnalysis/01_MetaAnalysisOut/01_MetaAnalysisOut.RData")

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.2.13
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.7.8
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.2.13
library(scales); packageVersion("scales") # 1.1.0, 2020.2.13
library(ggsci); packageVersion("ggsci") # 2.9, 2019.7.8
theme_set(theme_cowplot())

# Create output directory
fig_output <- "../00_RawFigs/03_Fig_MetaAnalysis"
dir.create(fig_output)

# Reformat figures
t1 <- geom_obspred(tara_ggdf, "TaraOceans expedition\nMetagenome") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(tara_r_adj, digits = 3))))
t2 <- geom_obspred(bahr_ggdf, "Global soil\nMicrobes") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(bahr_r_adj, digits = 3))))
t3 <- geom_obspred(maiz_ggdf, "Coastal region\nFish") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(maiz_r_adj, digits = 3))))
t4 <- geom_obspred(jplk_ggdf, "Japanese lakes\nProkaryotes") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(jplk_r_adj, digits = 3))))
t5 <- geom_obspred(suwa_ggdf, "Lake Suwa\nZooplankton") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(suwa_r_adj, digits = 3))))
t6 <- geom_obspred(lago_ggdf, "Freshwater lagoons\nMacroinvertebrates") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(lago_r_adj, digits = 3))))

t_all <- plot_grid(t1, t2, t3, t4, t5, t6,
                   labels = "auto", align = "hv")

saveRDS(t_all, sprintf("%s/DiversityMetaAnal.obj", fig_output))
