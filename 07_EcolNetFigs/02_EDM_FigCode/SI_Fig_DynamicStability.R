####
#### CERrice2017 All data analysis
#### Fig. Dynamic stability
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
load("../../02_EcolComAnalysis/07_RegularizedSmapOut/07_RegularizedSmapOut.RData")

# Ceate output directory
fig_output <- "../00_RawFigs/02_Fig_EDMnet"
dir.create(fig_output)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.30
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.30
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.11
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.11
theme_set(theme_cowplot())

# Extract dynamic stability values
ev_df <- data.frame(e_value_all)
ev_df$plot <- edna_all$plot
ev_df$date_id <- edna_all$date

ev_df <- melt(ev_df, id.vars = c("date_id", "plot"), measure.vars = "ev_1")
ev_df <- ev_df[!is.na(ev_df$plot),]

f1 <- ggplot(ev_df, aes(x = date_id, y = abs(value), group = plot, shape = plot, fill = plot, color = plot)) +
  geom_point() + geom_line() + geom_hline(yintercept = 1, linetype = 2) +
  scale_color_startrek() + scale_fill_startrek() + scale_shape_manual(values = c(21:25)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  xlim(as.Date("2017-05-22"), as.Date("2017-09-23")) + xlab("Date") + ylab("Dynamic stability") +
  ylim(0.8,8.5) + xlab(NULL)

pdf(sprintf("%s/EDMFig_DynamicStability.pdf", fig_output), width = 7, height = 5)
f1; dev.off()

saveRDS(f1, sprintf("%s/EDMFig_DynamicStability.obj", fig_output))
