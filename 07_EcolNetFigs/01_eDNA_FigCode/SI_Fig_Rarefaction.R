####
#### CER eDNA study
#### FigCode: Examine rarefaction curves
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.29
library(lubridate); packageVersion("lubridate") # 1.7.4, 2019.11.3
library(phyloseq); packageVersion("phyloseq") # 1.30.0, 2020.1.29
library(gridExtra); packageVersion("gridExtra") # 2.3, 2019.11.3
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.29
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.3
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.3
library(ape); packageVersion("ape") # 5.3, 2019.11.3
theme_set(theme_cowplot())

source("../00_Fig_functions/F01_FigHelperFunctions.R")

# Create output folder
fig_output <- "../00_RawFigs"
dir.create(fig_output)

## Load helper functions
scripts <- c("graphical_methods.R", "tree_methods.R", "plot_merged_trees.R",
             "specificity_methods.R", "ternary_plot.R", "richness.R",
             "edgePCA.R", "copy_number_correction.R", "import_frogs.R",
             "prevalence.R", "compute_niche.R")
urls <- paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/", scripts)
for (url in urls) source(url)


# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

#---------- Non-Filtered samples ----------#
# Load workspace
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Prokaryote/06_NCcheck_ProkOut/06_NCcheck_ProkOut.RData")
ps0_pro <- subset_samples(ps0, sample_nc == "sample")
ps0_pro <- subset_samples(ps0_pro, sample_data(ps0_pro)[,"date"] != "5/22/17")
ps0_pro <- subset_samples(ps0_pro, sample_data(ps0_pro)[,"NonSTD_all"] > 19) # 2 samples removed
ps0_pro <- subset_taxa(ps0_pro, std_or_field != "Standard DNA")
ps0_pro <- prune_taxa(taxa_sums(ps0_pro) > 0, ps0_pro)

load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Eukaryote/06_NCcheck_EukOut/06_NCcheck_EukOut.RData")
ps0_euk <- subset_samples(ps0, sample_nc == "sample")
ps0_euk <- subset_samples(ps0_euk, sample_data(ps0_euk)[,"date"] != "5/22/17")
ps0_euk <- subset_samples(ps0_euk, sample_data(ps0_euk)[,"NonSTD_all"] > 19) # 0 samples removed
ps0_euk <- subset_taxa(ps0_euk, std_or_field != "Standard DNA")
ps0_euk <- prune_taxa(taxa_sums(ps0_euk) > 0, ps0_euk)

load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Fungi/06_NCcheck_FungiOut/06_NCcheck_FungiOut.RData")
ps0_fun <- subset_samples(ps0, sample_nc == "sample")
ps0_fun <- subset_samples(ps0_fun, sample_data(ps0_fun)[,"date"] != "5/22/17")
ps0_fun <- subset_samples(ps0_fun, sample_data(ps0_fun)[,"NonSTD_all"] > 19) # 4 samples removed
ps0_fun <- subset_taxa(ps0_fun, std_or_field != "Standard DNA")
ps0_fun <- prune_taxa(taxa_sums(ps0_fun) > 0, ps0_fun)

load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Invertebrate/06_NCcheck_InvOut/06_NCcheck_InvOut.RData")
ps0_inv <- subset_samples(ps0, sample_nc == "sample")
ps0_inv <- subset_samples(ps0_inv, sample_data(ps0_inv)[,"date"] != "5/22/17")
ps0_inv <- subset_samples(ps0_inv, sample_data(ps0_inv)[,"NonSTD_all"] > 19) # 0 samples removed
ps0_inv <- subset_taxa(ps0_inv, std_or_field != "Standard DNA")
ps0_inv <- prune_taxa(taxa_sums(ps0_inv) > 0, ps0_inv)

# Add month information
sample_data(ps0_pro)$month <- factor(month(mdy(sample_data(ps0_pro)$date)))
sample_data(ps0_euk)$month <- factor(month(mdy(sample_data(ps0_euk)$date)))
sample_data(ps0_fun)$month <- factor(month(mdy(sample_data(ps0_fun)$date)))
sample_data(ps0_inv)$month <- factor(month(mdy(sample_data(ps0_inv)$date)))

sample_data(ps0_pro)$plot_label <- sprintf("Plot %s", sample_data(ps0_pro)$plot)
sample_data(ps0_euk)$plot_label <- sprintf("Plot %s", sample_data(ps0_euk)$plot)
sample_data(ps0_fun)$plot_label <- sprintf("Plot %s", sample_data(ps0_fun)$plot)
sample_data(ps0_inv)$plot_label <- sprintf("Plot %s", sample_data(ps0_inv)$plot)

# Illustrating rarefaction curve for non filtered samples
rr_pro <- ggrare(ps0_pro, step = 500, color = "month", se = FALSE, plot = FALSE) +
  xlab("Sequence reads") + ylab("No. of ASVs") + scale_color_startrek() +
  facet_wrap(.~ plot_label, ncol = 2) + panel_border() + ggtitle("Prokaryote (16S)") +
  theme(legend.position = c(0.8,0.1))

rr_euk <- ggrare(ps0_euk, step = 500, color = "month", se = FALSE, plot = FALSE) +
  xlab("Sequence reads") + ylab("No. of ASVs") + scale_color_startrek() +
  facet_wrap(.~ plot_label, ncol = 2) + panel_border() + ggtitle("Eukaryote (18S)") +
  theme(legend.position = c(0.8,0.1))

rr_fun <- ggrare(ps0_fun, step = 500, color = "month", se = FALSE, plot = FALSE) +
  xlab("Sequence reads") + ylab("No. of ASVs") + scale_color_startrek() +
  facet_wrap(.~ plot_label, ncol = 2) + panel_border() + ggtitle("Fungi (ITS)") +
  theme(legend.position = c(0.8,0.1))

rr_inv <- ggrare(ps0_inv, step = 500, color = "month", se = FALSE, plot = FALSE) +
  xlab("Sequence reads") + ylab("No. of ASVs") + scale_color_startrek() +
  facet_wrap(.~ plot_label, ncol = 2) + panel_border() + ggtitle("Animal (COI)") +
  theme(legend.position = c(0.8,0.1))


# Illustrating rarefaction curve for non filtered samples
rr_pro2 <- ggrare(ps0_pro, step = 500, color = "plot_label", se = FALSE, plot = FALSE) +
  xlab("Sequence reads") + ylab("No. of ASVs") + scale_color_startrek(name = "Plot", alpha = 0.7) +
  panel_border() + ggtitle("Prokaryote (16S)")

rr_euk2 <- ggrare(ps0_euk, step = 500, color = "plot_label", se = FALSE, plot = FALSE) +
  xlab("Sequence reads") + ylab("No. of ASVs") + scale_color_startrek(name = "Plot", alpha = 0.7) +
  panel_border() + ggtitle("Eukaryote (18S)")

rr_fun2 <- ggrare(ps0_fun, step = 500, color = "plot_label", se = FALSE, plot = FALSE) +
  xlab("Sequence reads") + ylab("No. of ASVs") + scale_color_startrek(name = "Plot", alpha = 0.7) +
  panel_border() + ggtitle("Fungi (ITS)")

rr_inv2 <- ggrare(ps0_inv, step = 500, color = "plot_label", se = FALSE, plot = FALSE) +
  xlab("Sequence reads") + ylab("No. of ASVs") + scale_color_startrek(name = "Plot", alpha = 0.7) +
  panel_border() + ggtitle("Animal (COI)")


# Combine and save figures
rr_all <- plot_grid(rr_pro, rr_euk, rr_fun, rr_inv,
                    labels = "auto",
                    ncol = 2, align = "hv", axis = "btlr")

pdf(sprintf("%s/01_Fig_eDNAts/Fig_Rarefaction_v1.pdf", fig_output), width = 14, height = 16)
rr_all; dev.off()

saveRDS(rr_all, file = "../00_RawFigs/01_Fig_eDNAts/Fig_Rarefaction_v1.obj")
saveRDS(rr_pro2, file = "../00_RawFigs/01_Fig_eDNAts/Fig_Rarefaction_pro.obj")
saveRDS(rr_euk2, file = "../00_RawFigs/01_Fig_eDNAts/Fig_Rarefaction_euk.obj")
saveRDS(rr_fun2, file = "../00_RawFigs/01_Fig_eDNAts/Fig_Rarefaction_fun.obj")
saveRDS(rr_inv2, file = "../00_RawFigs/01_Fig_eDNAts/Fig_Rarefaction_inv.obj")
