####
#### Compare the resuls of phyloFlash and qMiSeq
#### Visualize all figures
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.7.23
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.7.23
library(ggsci); packageVersion("ggsci") # 2.9, 2021.7.23
theme_set(theme_cowplot())

# Generate output folder
(fn <- basename(rstudioapi::getSourceEditorContext()$path))
output_folder <- paste0(str_sub(fn, end = -3), "Out"); rm(fn)
dir.create(output_folder)


# ---------------------------------------------- #
# Load ggplot objects
# ---------------------------------------------- #
g1 <- readRDS(sprintf("%s/g1_s001.obj", output_folder))
g2 <- readRDS(sprintf("%s/g1_s002.obj", output_folder))
g3 <- readRDS(sprintf("%s/g1_s003.obj", output_folder))
g4 <- readRDS(sprintf("%s/g1_s004.obj", output_folder))

g1 <- g1 + labs(title = "Plot 1\n2017/7/14", subtitle = "No. of detected taxa\nqMiSeq = 136\nShotgun = 118")
g2 <- g2 + labs(title = "Plot 1\n2017/7/15", subtitle = "No. of detected taxa\nqMiSeq = 129\nShotgun = 99")
g3 <- g3 + labs(title = "Plot 5\n2017/7/14", subtitle = "No. of detected taxa\nqMiSeq = 82\nShotgun = 116")
g4 <- g4 + labs(title = "Plot 5\n2017/7/15", subtitle = "No. of detected taxa\nqMiSeq = 91\nShotgun = 115")

# ---------------------------------------------- #
# Compile all figures
# ---------------------------------------------- #
g_legend <- get_legend(g1)
g5 <- g1 + theme(legend.position = "none", plot.margin = margin(.8,.8,.8,.8,"cm"))
g6 <- g2 + theme(legend.position = "none")
g7 <- g3 + theme(legend.position = "none")
g8 <- g4 + theme(legend.position = "none")

g_all <- plot_grid(g5, g6, g7, g8, NULL, g_legend,
                   ncol = 6, rel_widths = c(1,1,1,1,0.2,1),
                   labels = c("a","b","c","d",NA,NA),
                   align = "hv", axis = "lrbt")
g_obj <- list(g1, g2, g3, g4)

# ---------------------------------------------- #
# Save output
# ---------------------------------------------- #
ggsave(file = sprintf("%s/qMiSeq_Shotgun_All.pdf", output_folder),
       plot = g_all, width = 16, height = 7)
saveRDS(g_obj, sprintf("%s/qMiSeq_Shotgun_All.obj", output_folder))

