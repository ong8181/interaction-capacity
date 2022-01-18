####
#### qPCR and QuantIT calculation
#### Visualize summary
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.7.16
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.7.16
library(ggsci); packageVersion("ggsci") # 2.9, 2021.7.16
theme_set(theme_cowplot())

# Generate output folder
(fn <- basename(rstudioapi::getSourceEditorContext()$path))
output_folder <- paste0(str_sub(fn, end = -3), "Out"); rm(fn)
dir.create(output_folder)


# ---------------------------------------------- #
# Load data
# ---------------------------------------------- #
d_all <- read_csv("data/data_dna_quant.csv")
d_all$plot <- as.factor(d_all$plot)
## Exclude NA samples
na_sample <- is.na(d_all$qpcr_16s_ml) | is.na(d_all$quantit_ng_ml) |
  is.na(d_all$qmiseq_total_ml) | is.na(d_all$qmiseq_16s_ml) |
  is.na(d_all$miseq_16s_reads) | d_all$miseq_16s_reads < 10
d <- d_all[!na_sample,]


# ---------------------------------------------- #
# Perform Standardized Maximum Axis (SMA)
# ---------------------------------------------- #
## Quant-IT
qnit_sma <- d %>%
  smatr::sma("quantit_ng_ml ~ qmiseq_16s_ml", data = .)
qtit_sma_slope <- qnit_sma$coef[[1]][2,1]
qtit_sma_intct <- qnit_sma$coef[[1]][1,1]
## Quant-IT (log)
qnit_log <- d %>%
  smatr::sma("quantit_ng_ml ~ qmiseq_16s_ml", data = ., log = "xy")
qtit_log_slope <- qnit_log$coef[[1]][2,1]
qtit_log_intct <- qnit_log$coef[[1]][1,1]
## qPCR
qpcr_sma <- d %>%
  smatr::sma("qpcr_16s_ml ~ qmiseq_16s_ml", data = .)
qpcr_sma_slope <- qpcr_sma$coef[[1]][2,1]
qpcr_sma_intct <- qpcr_sma$coef[[1]][1,1]
## Quant-IT (log)
qpcr_log <- d %>%
  smatr::sma("qpcr_16s_ml ~ qmiseq_16s_ml", data = ., log = "xy")
qpcr_log_slope <- qpcr_log$coef[[1]][2,1]
qpcr_log_intct <- qpcr_log$coef[[1]][1,1]


# ---------------------------------------------- #
# Visualization
# ---------------------------------------------- #
label_func <- function(x) {
  ifelse(x == 0, "0", 
         parse(text = gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x))))
         )
  }
p_alpha = 0.75

(g1 <- ggplot(d, aes(x = qmiseq_total_ml, y = quantit_ng_ml, color = plot)) + 
    geom_point(alpha = p_alpha) +
    geom_abline(slope = qtit_log_slope, intercept = qtit_log_intct, linetype = 2) +
    scale_x_continuous(trans = "log10", label = label_func) +
    scale_y_continuous(trans = "log10") +
    scale_color_startrek() +
    xlab("qMiSeq (copies/ml water)") +
    ylab("Quant-IT (ng DNA/ml water)") +
    NULL)
(g2 <- ggplot(d_all, aes(x = qmiseq_16s_ml, y = qpcr_16s_ml, color = plot)) + 
    geom_point(alpha = p_alpha) +
    geom_abline(slope = qpcr_log_slope, intercept = qpcr_log_intct, linetype = 2) +
    scale_x_continuous(trans = "log10", label = label_func) +
    scale_y_continuous(trans = "log10", label = label_func) +
    scale_color_startrek() +
    xlab("qMiSeq (16S copies/ml water)") +
    ylab("qPCR (16S copies/ml water)") +
    NULL)
(g3 <- ggplot(d_all, aes(x = miseq_16s_reads, y = qpcr_16s_ml, color = plot)) + 
    geom_point(alpha = p_alpha) +
    scale_x_continuous(trans = "log10", label = label_func, limits = c(100,20000)) +
    scale_y_continuous(trans = "log10", label = label_func) +
    scale_color_startrek() +
    xlab("MiSeq 16S sequence reads") +
    ylab("qPCR (16S copies/ml water)") +
    NULL)

## Combine figures
g_all <- plot_grid(g1, g2, g3, ncol = 3, align = "hv")
a <- 0.7
g_all2 <- plot_grid(g1 + theme(legend.position = "none", plot.margin = unit(c(a,a,a,a), "cm")),
                    g2 + theme(legend.position = "none", plot.margin = unit(c(a,a,a,a), "cm")),
                    get_legend(g2),
                    ncol = 3, rel_widths = c(1,1,0.2),
                    labels = c("a", "b", NULL),
                    align = "hv")
g_obj <- list(g1, g1, g3)

# ---------------------------------------------- #
# Save output
# ---------------------------------------------- #
ggsave(file = sprintf("%s/qMiSeq_Comparison_all.pdf", output_folder),
       plot = g_all, width = 14, height = 4)
ggsave(file = sprintf("%s/qMiSeq_Comparison.pdf", output_folder),
       plot = g_all2, width = 9, height = 4)
saveRDS(g_obj, sprintf("%s/qMiSeq_Comparison_all.obj", output_folder))


