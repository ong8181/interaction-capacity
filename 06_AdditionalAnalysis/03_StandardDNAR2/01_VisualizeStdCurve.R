####
#### CERrice2017 additional data analysis
#### Combine all standard curve results
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
# Load workspace of Eukaryote
# ------------------------------------------ #
# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the data file
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Prokaryote/04_STDCheck_ProkOut/04_STDCheck_ProkOut.RData")
## Extract necessary objects
r2_pro <- r2_summary
std_tab_pro <- new_std_table2
keep_obj <- c("r2_pro", "std_tab_pro")
rm(list = ls()[!ls() %in% keep_obj])


# ------------------------------------------ #
# Load workspace of Eukaryote
# ------------------------------------------ #
# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the data file
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Fungi/04_STDCheck_FungiOut/04_STDCheck_FungiOut.RData")
## Extract necessary objects
r2_fun <- r2_summary
std_tab_fun <- new_std_table2
keep_obj <- c("r2_fun", "std_tab_fun", "r2_pro", "std_tab_pro")
rm(list = ls()[!ls() %in% keep_obj])


# ------------------------------------------ #
# Load workspace of Eukaryote
# ------------------------------------------ #
# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the data file
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Invertebrate/04_STDCheck_InvOut/04_STDCheck_InvOut.RData")
## Extract necessary objects
r2_inv <- r2_summary
std_tab_inv <- new_std_table2
keep_obj <- c("r2_inv", "std_tab_inv", "r2_fun", "std_tab_fun", "r2_pro", "std_tab_pro")
rm(list = ls()[!ls() %in% keep_obj])


# ------------------------------------------ #
# Load workspace of Eukaryote
# ------------------------------------------ #
# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the data file
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Eukaryote/04_STDCheck_EukOut/04_STDCheck_EukOut.RData")
## Extract necessary objects
r2_euk <- r2_summary
std_tab_euk <- new_std_table2
keep_obj <- c("r2_euk", "std_tab_euk", "r2_inv", "std_tab_inv", "r2_fun", "std_tab_fun", "r2_pro", "std_tab_pro")
rm(list = ls()[!ls() %in% keep_obj])


# ------------------------------------------ #
# Generate output folder
# ------------------------------------------ #
(fn <- basename(rstudioapi::getSourceEditorContext()$path))
output_folder <- paste0(stringr::str_sub(fn, end = -3), "Out"); rm(fn)
dir.create(output_folder)


# ------------------------------------------ #
# Compile data
# ------------------------------------------ #
d_r2 <- tibble("16S" = r2_pro,
               "ITS" = r2_fun,
               "COI" = r2_inv,
               "18S" = r2_euk)
d_r2_long <- pivot_longer(d_r2, cols = 1:4)

## Extract example standard curves
d_std_pro <- tibble(std_tab_pro)
lm_coef_fun <- function(x) summary(lm(as.numeric(x) ~ std_copy_n + 0))$coefficients[1]
std_copy_n <- c(100000,50000,25000,10000,5000)
coef_summary <- apply(d_std_pro, 1, lm_coef_fun)
max_id <- which.max(coef_summary); coef_summary[max_id]
med_id <- which(order(coef_summary) == round(length(coef_summary))/2); coef_summary[med_id]
min_id <- which(coef_summary == min(coef_summary[coef_summary > 0])); coef_summary[min_id]

d_std_example <- tibble(std_copy_n = rep(std_copy_n, 3),
                        reads = c(unlist(d_std_pro[max_id,]),
                                  unlist(d_std_pro[med_id,]),
                                  unlist(d_std_pro[min_id,])),
                        cat = c(rep("max", 5), rep("med", 5), rep("min", 5)))

                    
                    
# ------------------------------------------ #
# Visualize patterns
# ------------------------------------------ #
(g1 <- ggplot(d_std_example, aes(x = std_copy_n, y = reads, color = cat)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("red3", "darkred", "gray20"), name = NULL) +
    stat_smooth(method = "lm", formula = 'y ~ x + 0', se = FALSE, size = 0.5) +
    xlim(0, 100000) +
    theme(legend.position = c(0.1, 0.7)) +
    xlab(expression(paste("Standard DNA (copies/", mu, "l)"))) +
    ylab("Sequence reads") +
    ggtitle("Examples of standard curves") +
    NULL)

(g2 <- ggplot(d_r2_long, aes(x = value, fill = name, facet = name)) +
    geom_histogram(position = "identity", alpha = 1) +
    scale_fill_startrek() +
    facet_wrap(. ~ name) +
    panel_border() +
    scale_x_continuous(breaks = seq(0.5, 1, 0.1), limits = c(0.5, 1)) +
    xlab(expression(paste(R^2, " of linear regression"))) + ylab("Frequency") +
    ggtitle("The quality of standard curves") +
    NULL)

(g_all <- plot_grid(g1, g2, ncol = 2, labels = c("A", "B"), align = "hv", axis = "lrbt"))
g_obj <- list(g1, g2)

# ------------------------------------------ #
# Save output
# ------------------------------------------ #
ggsave(file = sprintf("%s/StdDNA_QualityCheck.pdf", output_folder),
       plot = g_all, width = 12, height = 6)
saveRDS(g_obj, sprintf("%s/StdDNA_QualityCheck.obj", output_folder))



