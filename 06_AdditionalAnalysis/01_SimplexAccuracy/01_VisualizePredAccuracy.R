####
#### CERrice2017 additional data analysis
#### Prediction accuracy
####

# ------------------------------------------ #
# Load libraries
# ------------------------------------------ #
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.7.27
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.7.27
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.7.27
library(ggsci); packageVersion("ggsci") # 2.9, 2021.7.27
theme_set(theme_cowplot())


# ------------------------------------------ #
# Generate output folder
# ------------------------------------------ #
(fn <- basename(rstudioapi::getSourceEditorContext()$path))
output_folder <- paste0(stringr::str_sub(fn, end = -3), "Out"); rm(fn)
dir.create(output_folder)


# ------------------------------------------ #
# Load prediction accuracy results
# ------------------------------------------ #
d <- readRDS("00_Robj/BestE_eDNA_ts_w_stats.obj")
best_E <- sapply(d, '[')

d_best <- d %>% map_dbl(function(x) x$E)
d_rho <- d %>% map_dbl(function(x) max(x$stat$rho))
d_mae <- d %>% map_dbl(function(x) min(x$stat$mae))
d_rmse <- d %>% map_dbl(function(x) min(x$stat$rmse))

d_all <- tibble(taxa = names(d), E = d_best, max_rho = d_rho, min_mae = d_mae, min_rmse = d_rmse)
d_long <- pivot_longer(d_all, cols = -taxa)

## Generate figures
(g1 <- d_long %>% filter(name != "E") %>%
        ggplot(aes(x = name, y = value)) +
        geom_jitter(width = 0.2, height = 0, alpha = 0.2) +
        geom_boxplot(width = 0.4, outlier.color = NA, outlier.shape = 0, fill = NA) +
        #scale_x_discrete(labels = c("A", "B", "C")) +
        xlab(NULL) + ylab("Prediction accuracy") +
        NULL)

(g2 <- ggplot(d_all, aes(x = max_rho)) +
        geom_histogram(alpha = 0.8) +
        #scale_x_discrete(labels = c("A", "B", "C")) +
        xlab(expression(rho)) + ylab("Frequency") +
        scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(-0.1, 1.1)) +
        #scale_y_log10() +
        geom_hline(yintercept = 0, linetype = 2) +
        ggtitle("Correlation coefficients between predicted and observed values\nby simplex projection") +
        NULL)

(g3 <- d_long %>% filter(name == "max_rho") %>% 
        ggplot(aes(x = name, y = value)) +
        #geom_violin() +
        geom_jitter(width = 0.2, height = 0, alpha = 0.2) +
        geom_boxplot(width = 0.4, outlier.color = NA, outlier.shape = 0, fill = "white", alpha = 0.4) +
        #scale_x_discrete(labels = c("A", "B", "C")) +
        xlab("Results of simplex projection") + ylab("Prediction accuracy") + theme(axis.text.x = element_blank()) +
        NULL)

ggsave(filename = sprintf("%s/SimplexProjection_rho.pdf", output_folder), plot = g2,
       width = 5, height = 4)
ggsave(filename = sprintf("%s/SimplexProjection_boxplot.pdf", output_folder), plot = g3,
       width = 5, height = 4)
