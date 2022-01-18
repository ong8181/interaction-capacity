####
#### CER eDNA study
#### FigCode: Overall pattern 1, Include all detected ASVs
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

# Remove "5/22/17" samples
ps_pro_sample2 <- subset_samples(ps_pro_sample1, sample_data(ps_pro_sample1)[,"date"] != "5/22/17") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_fun_sample2 <- subset_samples(ps_fun_sample1, sample_data(ps_fun_sample1)[,"date"] != "5/22/17") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_inv_sample2 <- subset_samples(ps_inv_sample1, sample_data(ps_inv_sample1)[,"date"] != "5/22/17") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps_euk_sample2 <- subset_samples(ps_euk_sample1, sample_data(ps_euk_sample1)[,"date"] != "5/22/17") %>%
  prune_taxa(taxa_sums(.) > 0, .)

#---------- RMR-076 Prokaryote ----------#
ps_pro_sample3 <- taxa_name_summarize(ps_pro_sample2, "phylum", top_n = 12)
ps_pro_m1 <- psmelt(ps_pro_sample3); ps_pro_m1$date <- mdy(ps_pro_m1$date)
ps_pro_m2 <- stats::aggregate(ps_pro_m1$Abundance, by=list(ps_pro_m1$date, ps_pro_m1$phylum), "sum") # Summed up to make phylum sum
ps_pro_m3 <- stats::aggregate(ps_pro_m1$Abundance, by=list(ps_pro_m1$date, ps_pro_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_pro_m2) <- c("date", "phylum", "dna_conc")
colnames(ps_pro_m3) <- c("date", "rep_tax", "dna_conc")
ps_pro_m2$phylum <- factor(ps_pro_m2$phylum, levels = c(levels(ps_pro_m2$phylum), "Undetermined"))
ps_pro_m2$phylum[ps_pro_m2$phylum == ""] <- "Undetermined"
ps_pro_m3$rep_tax <- factor(ps_pro_m3$rep_tax, levels = c(levels(ps_pro_m3$rep_tax)[c(1:8, 10:11, 13, 9, 12)]))

# Figures
f1 <- ggplot(ps_pro_m2, aes(x = as.Date(date), y = dna_conc/5, group = phylum, fill = phylum)) # Divide by 5 to make average copy number / plot
f1 <- f1 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
f1 <- f1 + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
f1 <- f1 + scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)")
f1 <- f1 + ggtitle("All phyla by RMR-076")

f2 <- ggplot(ps_pro_m3, aes(x = as.Date(date), y = dna_conc/5, group = rep_tax, fill = rep_tax)) # Divide by 5 to make average copy number / plot
f2 <- f2 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
f2 <- f2 + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
f2 <- f2 + scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)")
f2 <- f2 + ggtitle("Major phyla by RMR-076")
#----------------------------------------#

#---------- RMR-078 Fungi ----------#
ps_fun_sample3 <- taxa_name_summarize(ps_fun_sample2, "phylum", top_n = 9)
ps_fun_m1 <- psmelt(ps_fun_sample3); ps_fun_m1$date <- mdy(ps_fun_m1$date)
ps_fun_m2 <- stats::aggregate(ps_fun_m1$Abundance, by=list(ps_fun_m1$date, ps_fun_m1$phylum), "sum") # Summed up to make phylum sum
ps_fun_m3 <- stats::aggregate(ps_fun_m1$Abundance, by=list(ps_fun_m1$date, ps_fun_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_fun_m2) <- c("date", "phylum", "dna_conc")
colnames(ps_fun_m3) <- c("date", "rep_tax", "dna_conc")
ps_fun_m2$phylum <- factor(ps_fun_m2$phylum, levels = c(levels(ps_fun_m2$phylum), "Undetermined"))
ps_fun_m2$phylum[ps_fun_m2$phylum == ""] <- "Undetermined"
ps_fun_m3$rep_tax <- factor(ps_fun_m3$rep_tax, levels = c(levels(ps_fun_m3$rep_tax)[c(1:4, 6:7, 9, 5, 8)]))

# Figures
f3 <- ggplot(ps_fun_m2, aes(x = as.Date(date), y = dna_conc/5, group = phylum, fill = phylum)) # Divide by 5 to make average copy number / plot
f3 <- f3 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
f3 <- f3 + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
f3 <- f3 + scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)")
f3 <- f3 + ggtitle("All phyla by RMR-078")

f4 <- ggplot(ps_fun_m3, aes(x = as.Date(date), y = dna_conc/5, group = rep_tax, fill = rep_tax)) # Divide by 5 to make average copy number / plot
f4 <- f4 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
f4 <- f4 + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
f4 <- f4 + scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)")
f4 <- f4 + ggtitle("Major phyla by RMR-078")
#----------------------------------------#


#---------- RMR-099 Invertebrate ----------#
ps_inv_sample3 <- taxa_name_summarize(ps_inv_sample2, "order", top_n = 12)
ps_inv_m1 <- psmelt(ps_inv_sample3); ps_inv_m1$date <- mdy(ps_inv_m1$date)
ps_inv_m2 <- stats::aggregate(ps_inv_m1$Abundance, by=list(ps_inv_m1$date, ps_inv_m1$phylum), "sum") # Summed up to make phylum sum
ps_inv_m3 <- stats::aggregate(ps_inv_m1$Abundance, by=list(ps_inv_m1$date, ps_inv_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_inv_m2) <- c("date", "order", "dna_conc")
colnames(ps_inv_m3) <- c("date", "rep_tax", "dna_conc")
ps_inv_m2$order <- factor(ps_inv_m2$order, levels = c(levels(ps_inv_m2$order), "Undetermined"))
ps_inv_m2$order[ps_inv_m2$order == ""] <- "Undetermined"
ps_inv_m3$rep_tax <- factor(ps_inv_m3$rep_tax, levels = c(levels(ps_inv_m3$rep_tax)[c(1:7, 9:11, 8, 12)]))

# Figures
f5 <- ggplot(ps_inv_m2, aes(x = as.Date(date), y = dna_conc/5, group = order, fill = order)) # Divide by 5 to make average copy number / plot
f5 <- f5 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
f5 <- f5 + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
f5 <- f5 + scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)")
f5 <- f5 + ggtitle("All orders by RMR-099")

f6 <- ggplot(ps_inv_m3, aes(x = as.Date(date), y = dna_conc/5, group = rep_tax, fill = rep_tax)) # Divide by 5 to make average copy number / plot
f6 <- f6 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
f6 <- f6 + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
f6 <- f6 + scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)")
f6 <- f6 + ggtitle("Major orders by RMR-099")
#----------------------------------------#

#---------- CMR-002 Eukaryote ----------#
ps_euk_sample3 <- taxa_name_summarize(ps_euk_sample2, "phylum", top_n = 12)
ps_euk_m1 <- psmelt(ps_euk_sample3); ps_euk_m1$date <- mdy(ps_euk_m1$date)
ps_euk_m2 <- stats::aggregate(ps_euk_m1$Abundance, by=list(ps_euk_m1$date, ps_euk_m1$phylum), "sum") # Summed up to make phylum sum
ps_euk_m3 <- stats::aggregate(ps_euk_m1$Abundance, by=list(ps_euk_m1$date, ps_euk_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_euk_m2) <- c("date", "phylum", "dna_conc")
colnames(ps_euk_m3) <- c("date", "rep_tax", "dna_conc")
ps_euk_m2$phylum <- factor(ps_euk_m2$phylum, levels = c(levels(ps_euk_m2$phylum), "Undetermined"))
ps_euk_m2$phylum[ps_euk_m2$phylum == ""] <- "Undetermined"
ps_euk_m3$rep_tax <- factor(ps_euk_m3$rep_tax, levels = c(levels(ps_euk_m3$rep_tax)[c(1:8, 10, 11, 9, 12)]))

# Figures
f7 <- ggplot(ps_euk_m2, aes(x = as.Date(date), y = dna_conc/5, group = phylum, fill = phylum)) # Divide by 5 to make average copy number / plot
f7 <- f7 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
f7 <- f7 + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
f7 <- f7 + scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)")
f7 <- f7 + ggtitle("All phyla by CMR-002")

f8 <- ggplot(ps_euk_m3, aes(x = as.Date(date), y = dna_conc/5, group = rep_tax, fill = rep_tax)) # Divide by 5 to make average copy number / plot
f8 <- f8 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
f8 <- f8 + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
f8 <- f8 + scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)")
f8 <- f8 + ggtitle("Major phyla by CMR-002")
#----------------------------------------#

# Combine all figures
legend_pro <- get_legend(f1)
legend_fun <- get_legend(f3)
legend_inv <- get_legend(f5)
legend_euk <- get_legend(f7)

f_all1 <- plot_grid(f1 + theme(legend.position = "none"),
                     f3 + theme(legend.position = "none"),
                     f5 + theme(legend.position = "none"),
                     f7 + theme(legend.position = "none"), ncol = 1,
                     align = "hv", labels = NULL)
f_all2 <- plot_grid(legend_pro, legend_fun, legend_inv, legend_euk, ncol = 1, labels = NULL)
f_all <- plot_grid(f_all1, f_all2, ncol = 2, rel_widths = c(1, 0.3), labels = NULL)

# Save figures
pdf(sprintf("%s/01_Fig_eDNAts/Fig_BarAllRuns.pdf", fig_output), width = 22, height = 24)
f_all; dev.off()

pdf(sprintf("%s/01_Fig_eDNAts/Fig_ProkPhylumBars.pdf", fig_output), width = 16, height = 10)
plot_grid(f1, f2, ncol = 1, align = "hv"); dev.off()

pdf(sprintf("%s/01_Fig_eDNAts/Fig_FungiPhylumBars.pdf", fig_output), width = 12, height = 8)
plot_grid(f3, f4, ncol = 1, align = "hv"); dev.off()

pdf(sprintf("%s/01_Fig_eDNAts/Fig_InvPhylumBars.pdf", fig_output), width = 12, height = 8)
plot_grid(f5, f6, ncol = 1, align = "hv"); dev.off()

pdf(sprintf("%s/01_Fig_eDNAts/Fig_EukPhylumBars.pdf", fig_output), width = 16, height = 10)
plot_grid(f7, f8, ncol = 1, align = "hv"); dev.off()

# Save all figures as R objects
saveRDS(f1, "../00_RawFigs/01_Fig_eDNAts/Fig_ProkBarAll.obj")
saveRDS(f2, "../00_RawFigs/01_Fig_eDNAts/Fig_ProkBarAll2.obj")
saveRDS(f3, "../00_RawFigs/01_Fig_eDNAts/Fig_FungiBarAll.obj")
saveRDS(f4, "../00_RawFigs/01_Fig_eDNAts/Fig_FungiBarAll2.obj")
saveRDS(f5, "../00_RawFigs/01_Fig_eDNAts/Fig_InvBarAll.obj")
saveRDS(f6, "../00_RawFigs/01_Fig_eDNAts/Fig_InvBarAll2.obj")
saveRDS(f7, "../00_RawFigs/01_Fig_eDNAts/Fig_EukBarAll.obj")
saveRDS(f8, "../00_RawFigs/01_Fig_eDNAts/Fig_EukBarAll2.obj")
