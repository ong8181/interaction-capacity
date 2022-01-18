####
#### Edit figures of the additional analyses
####


# ------------------------------------------ #
# Load packages
# ------------------------------------------ #
library(tidyverse); packageVersion("tidyverse") # 1.3.1
library(cowplot); packageVersion("cowplot") # 1.1.1
theme_set(theme_cowplot())


# ------------------------------------------ #
# Generate output folder
# ------------------------------------------ #
(fn <- basename(rstudioapi::getSourceEditorContext()$path))
output_folder <- paste0(stringr::str_sub(fn, end = -3), "Out"); rm(fn)
dir.create(output_folder)


# ------------------------------------------ #
# Load data
# ------------------------------------------ #
## Ratzke et al. 2020 Nature Ecology & Evolution
## DOI: 10.1038/s41559-020-1099-4
# You need to contact ong8181@gmail.com, or the original authors to obtain the data file
d_ratzke <- read_csv("data/Ratzke_Fig3b_comm.csv")
d_ratzke$interaction_strength <- factor(d_ratzke$interaction_strength, levels = c("weak", "moderate", "strong"))
## Yamawo et al. 2021 Functional Ecology
## DOI: 10.1111/1365-2435.13859
# You need to contact ong8181@gmail.com, or the original authors to obtain the data file
d_yamawo <- read_csv("data/Yamawo_2021_FuncEcol_data2.csv")
d_yamawo1 <- d_yamawo[,9:ncol(d_yamawo)]
d_yamawo1[d_yamawo1 == 0] <- NA
d_yamawo$mean_is <- rowMeans(d_yamawo1, na.rm = TRUE)
d_yamawo$total_is <- rowSums(d_yamawo1, na.rm = TRUE)


# ------------------------------------------ #
# Visualize results
# ------------------------------------------ #
(g1 <- ggplot(d_ratzke, aes(x = interaction_strength, y = diversity, color = community, shape = community)) +
        geom_boxplot(aes(group = interaction_strength), outlier.color = NA, outlier.size = 0) +
        geom_jitter(height = 0, width = 0.2, size = 2) +
        scale_color_manual(values = c("red3", "royalblue", "gray20"),
                           guide = guide_legend(override.aes = list(linetype = "blank"))) +
        xlab("Interaction strength") + ylab("Diversity index") +
        ylim(0, 9) +
        labs(color = "Community type", shape = "Community type") +
        NULL)
(g2 <- ggplot(d_yamawo, aes(x = mean_is, y = all_ant_sp, fill = season)) +
        geom_point(shape = 21, size = 2) + scale_fill_manual(values = c("white", "black")) +
        ylim(0,20) + xlim(0,20) +
        ylab("Ant species richness") + xlab("Mean interaction strength") +
        NULL)


# ------------------------------------------ #
# Save output
# ------------------------------------------ #
saveRDS(list(g1, g2), sprintf("%s/Ratzke_Yamawo.obj", output_folder))
