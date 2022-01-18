####
#### CER eDNA study
#### FigCode: Overall diversity 1, Exclude non-target ASVs, all rare ASVs included
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
load("../../01_DNAtsCERrice2017/02_TimeSeriesCompile/08_TSfilterPrepOut/08_TSfilterPrepOut.RData")

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.29
library(lubridate); packageVersion("lubridate") # 1.7.4, 2019.11.11
library(phyloseq); packageVersion("phyloseq") # 1.30.0, 2020.1.29
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.29
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.11
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.11
theme_set(theme_cowplot())

source("../00_Fig_functions/F01_FigHelperFunctions.R")

# Create output folder
fig_output <- "../00_RawFigs"
dir.create(fig_output)

# Remove "5/22/27" samples
ps_pro_sample2 <- subset_samples(ps_pro_sample1, sample_data(ps_pro_sample1)[,"date"] != "5/22/17") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(., ceiling)
ps_fun_sample2 <- subset_samples(ps_fun_sample1, sample_data(ps_fun_sample1)[,"date"] != "5/22/17") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(., ceiling)
ps_inv_sample2 <- subset_samples(ps_inv_sample1, sample_data(ps_inv_sample1)[,"date"] != "5/22/17") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(., ceiling)
ps_euk_sample2 <- subset_samples(ps_euk_sample1, sample_data(ps_euk_sample1)[,"date"] != "5/22/17") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(., ceiling)

#---------- RMR-076 Prokaryote ----------#
# Compile phyloseq object
# Exclude non-prokayote DNAs
pro_tax_cond <- tax_table(ps_pro_sample2)[,"superkingdom"] == "Bacteria" | tax_table(ps_pro_sample2)[,"superkingdom"] == "Archaea"
pro_tax_extract <- taxa_names(tax_table(ps_pro_sample2)[pro_tax_cond])
ps_pro_sample3 <- prune_taxa(pro_tax_extract, ps_pro_sample2)

# Figures
sample_data(ps_pro_sample3)$date <- mdy(sample_data(ps_pro_sample3)$date)
f1 <- plot_richness(ps_pro_sample3, x = "date", color = "plot", measures = c("Observed"))
f1$layers <- f1$layers[-1]
f1 <- f1 + geom_line() + geom_point(size = 1.5) + ylab("No. of ASV")
f1 <- f1 + scale_color_startrek() + scale_shape_manual(values = c(21:25))
f1 <- f1 + theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
f1 <- f1 + ggtitle("Prokaryotic diversity")
#----------------------------------------#

#---------- RMR-078 Fungi ----------#
# Compile phyloseq object
# Exclude non-fungal DNAs
fun_tax_cond <- tax_table(ps_fun_sample2)[,"kingdom"] == "Fungi"
fun_tax_extract <- taxa_names(tax_table(ps_fun_sample2)[fun_tax_cond])
ps_fun_sample3 <- prune_taxa(fun_tax_extract, ps_fun_sample2)

# Figures
sample_data(ps_fun_sample3)$date <- mdy(sample_data(ps_fun_sample3)$date)
f2 <- plot_richness(ps_fun_sample3, x = "date", color = "plot", measures = c("Observed"))
f2$layers <- f2$layers[-1]
f2 <- f2 + geom_line() + geom_point(size = 1.5) + ylab("No. of ASV")
f2 <- f2 + scale_color_startrek() + scale_shape_manual(values = c(21:25))
f2 <- f2 + theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
f2 <- f2 + ggtitle("Fungal diversity")
#----------------------------------------#


#---------- RMR-099 Invertebrate ----------#
# Compile phyloseq object
# Exclude non-fungal DNAs
inv_tax_cond <- tax_table(ps_inv_sample2)[,"kingdom"] == "Metazoa"
inv_tax_extract <- taxa_names(tax_table(ps_inv_sample2)[inv_tax_cond])
ps_inv_sample3 <- prune_taxa(inv_tax_extract, ps_inv_sample2)

# Figures
sample_data(ps_inv_sample3)$date <- mdy(sample_data(ps_inv_sample3)$date)
f3 <- plot_richness(ps_inv_sample3, x = "date", color = "plot", measures = c("Observed"))
f3$layers <- f3$layers[-1]
f3 <- f3 + geom_line() + geom_point(size = 1.5) + ylab("No. of ASV")
f3 <- f3 + scale_color_startrek() + scale_shape_manual(values = c(21:25))
f3 <- f3 + theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
f3 <- f3 + ggtitle("Metazoan diversity")
#----------------------------------------#

#---------- CMR-002 Eukaryote ----------#
# Compile phyloseq object
# Exclude non-fungal DNAs
euk_tax_cond1 <- tax_table(ps_euk_sample2)[,"superkingdom"] == "Eukaryota"
euk_tax_cond2 <- tax_table(ps_euk_sample2)[,"kingdom"] != "Fungi"
euk_tax_cond3 <- tax_table(ps_euk_sample2)[,"kingdom"] != "Metazoa"
euk_tax_extract <- taxa_names(tax_table(ps_euk_sample2)[euk_tax_cond1 & euk_tax_cond2 & euk_tax_cond3])
ps_euk_sample3 <- prune_taxa(euk_tax_extract, ps_euk_sample2)

# Figures
sample_data(ps_euk_sample3)$date <- mdy(sample_data(ps_euk_sample3)$date)
f4 <- plot_richness(ps_euk_sample3, x = "date", color = "plot", measures = c("Observed"))
f4$layers <- f4$layers[-1]
f4 <- f4 + geom_line() + geom_point(size = 1.5) + ylab("No. of ASV")
f4 <- f4 + scale_color_startrek() + scale_shape_manual(values = c(21:25))
f4 <- f4 + theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
f4 <- f4 + ggtitle("Eukaryotic (non-fungi & non-metazoa) diversity")
#----------------------------------------#

# Combine all figures
legend_plot <- get_legend(f1)

f_all <- plot_grid(f1 + theme(legend.position = "none"),
                   f2 + theme(legend.position = "none"), legend_plot,
                   f3 + theme(legend.position = "none"),
                   f4 + theme(legend.position = "none"), legend_plot, ncol = 3,
                   labels = NULL, rel_widths = c(1, 1, 0.2))

pdf(sprintf("%s/01_Fig_eDNAts/Fig_DivNotFiltered.pdf", fig_output), width = 14, height = 9)
f_all; dev.off()

saveRDS(f_all, sprintf("%s/01_Fig_eDNAts/Fig_DivNotFiltered.obj", fig_output))
