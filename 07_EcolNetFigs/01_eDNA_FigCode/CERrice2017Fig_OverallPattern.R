####
#### CER eDNA study
#### FigCode: Overall pattern 1, Include all detected ASVs
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
load("../../01_DNAtsCERrice2017/02_TimeSeriesCompile/08_TSfilterPrepOut/08_TSfilterPrepOut.RData")

# Load library and functions
library(phyloseq); packageVersion("phyloseq") # 1.38.0, 2021.12.8
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.7.29
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.7.29
library(reshape2); packageVersion("reshape2") # 1.4.4, 2021.7.29
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.7.29
library(scales); packageVersion("scales") # 1.1.1, 2021.7.29
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.3
library(lubridate); packageVersion("lubridate") # 1.8.0, 2021.12.8
theme_set(theme_cowplot())
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))
palette_pro <- get_palette(13); palette_pro <- palette_pro[c(1,7,2,8,3,9,4,10,5,11,6,12,13)]; palette_pro[12:13] <- c("gray20","gray60")
palette_euk <- get_palette(12); palette_euk <- palette_euk[c(1,7,2,8,3,9,4,10,5,11,6,12)]; palette_euk[11:12] <- c("gray20","gray60")
palette_its <- get_palette(9); palette_its <- palette_its[c(1,6,2,7,3,8,4,9,5)]; palette_its[8:9]   <- c("gray20","gray60")
palette_coi <- get_palette(12); palette_coi <- palette_coi[c(1,7,2,8,3,9,4,10,5,11,6,12)]; palette_coi[11:12] <- c("gray20","gray60")

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
ps_pro_m1 <- speedyseq::psmelt(ps_pro_sample3); ps_pro_m1$date <- mdy(ps_pro_m1$date)
ps_pro_m2 <- stats::aggregate(ps_pro_m1$Abundance, by=list(ps_pro_m1$date, ps_pro_m1$phylum), "sum") # Summed up to make phylum sum
ps_pro_m3 <- stats::aggregate(ps_pro_m1$Abundance, by=list(ps_pro_m1$date, ps_pro_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_pro_m2) <- c("date", "phylum", "dna_conc")
colnames(ps_pro_m3) <- c("date", "rep_tax", "dna_conc")
ps_pro_m2$phylum <- factor(ps_pro_m2$phylum, levels = c(unique(ps_pro_m2$phylum), "Undetermined"))
ps_pro_m2$phylum[ps_pro_m2$phylum == ""] <- "Undetermined"
ps_pro_m3$rep_tax <- factor(ps_pro_m3$rep_tax, levels = c(unique(ps_pro_m3$rep_tax)[c(1:8, 10:11, 13, 9, 12)]))

# Figures
f1 <- ggplot(ps_pro_m2, aes(x = as.Date(date), y = dna_conc/5, group = phylum, fill = phylum)) + # Divide by 5 to make average copy number / plot
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") +
  ggtitle("All phyla by RMR-076")

f2 <- ggplot(ps_pro_m3, aes(x = as.Date(date), y = dna_conc/5, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  xlab(NULL) + ylab("DNA (copies/ml water)") +
  scale_fill_manual(values = palette_pro) +
  ggtitle("Major phyla detected by 16S sequencing") +
  scale_y_continuous(breaks = seq(0, 1e+07, 2e+06), limits = c(0, 6.1e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
#----------------------------------------#

#---------- RMR-078 Fungi ----------#
ps_fun_sample3 <- taxa_name_summarize(ps_fun_sample2, "phylum", top_n = 9)
ps_fun_m1 <- speedyseq::psmelt(ps_fun_sample3); ps_fun_m1$date <- mdy(ps_fun_m1$date)
ps_fun_m2 <- stats::aggregate(ps_fun_m1$Abundance, by=list(ps_fun_m1$date, ps_fun_m1$phylum), "sum") # Summed up to make phylum sum
ps_fun_m3 <- stats::aggregate(ps_fun_m1$Abundance, by=list(ps_fun_m1$date, ps_fun_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_fun_m2) <- c("date", "phylum", "dna_conc")
colnames(ps_fun_m3) <- c("date", "rep_tax", "dna_conc")
ps_fun_m2$phylum <- factor(ps_fun_m2$phylum, levels = c(unique(ps_fun_m2$phylum), "Undetermined"))
ps_fun_m2$phylum[ps_fun_m2$phylum == ""] <- "Undetermined"
ps_fun_m3$rep_tax <- factor(ps_fun_m3$rep_tax, levels = c(unique(ps_fun_m3$rep_tax)[c(1:5, 7, 9, 6, 8)]))

# Figures
f3 <- ggplot(ps_fun_m2, aes(x = as.Date(date), y = dna_conc/5, group = phylum, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") +
  ggtitle("All phyla by RMR-078")

f4 <- ggplot(ps_fun_m3, aes(x = as.Date(date), y = dna_conc/5, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  xlab(NULL) + ylab("DNA (copies/ml water)") +
  scale_fill_manual(values = palette_its) +
  ggtitle("Major phyla detected by ITS sequencing") +
  scale_y_continuous(breaks = seq(0, 3e+04, 1e+04), limits = c(0, 3.1e+04),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
#----------------------------------------#


#---------- RMR-099 Invertebrate ----------#
ps_inv_sample3 <- taxa_name_summarize(ps_inv_sample2, "order", top_n = 12)
ps_inv_m1 <- speedyseq::psmelt(ps_inv_sample3); ps_inv_m1$date <- mdy(ps_inv_m1$date)
ps_inv_m2 <- stats::aggregate(ps_inv_m1$Abundance, by=list(ps_inv_m1$date, ps_inv_m1$phylum), "sum") # Summed up to make phylum sum
ps_inv_m3 <- stats::aggregate(ps_inv_m1$Abundance, by=list(ps_inv_m1$date, ps_inv_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_inv_m2) <- c("date", "order", "dna_conc")
colnames(ps_inv_m3) <- c("date", "rep_tax", "dna_conc")
ps_inv_m2$order <- factor(ps_inv_m2$order, levels = c(unique(ps_inv_m2$order), "Undetermined"))
ps_inv_m2$order[ps_inv_m2$order == ""] <- "Undetermined"
ps_inv_m3$rep_tax <- factor(ps_inv_m3$rep_tax, levels = c(unique(ps_inv_m3$rep_tax)[c(1:8, 10:11, 9, 12)]))

# Figures
f5 <- ggplot(ps_inv_m2, aes(x = as.Date(date), y = dna_conc/5, group = order, fill = order)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") +
  ggtitle("All orders by RMR-099")

f6 <- ggplot(ps_inv_m3, aes(x = as.Date(date), y = dna_conc/5, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  xlab(NULL) + ylab("DNA (copies/ml water)") +
  scale_fill_manual(values = palette_coi) +
  ggtitle("Major orders detected by COI sequencing") +
  scale_y_continuous(breaks = seq(0, 1e+05, 2e+04), limits = c(0, 1e+05),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
#----------------------------------------#

#---------- CMR-002 Eukaryote ----------#
ps_euk_sample3 <- taxa_name_summarize(ps_euk_sample2, "phylum", top_n = 12)
ps_euk_m1 <- speedyseq::psmelt(ps_euk_sample3); ps_euk_m1$date <- mdy(ps_euk_m1$date)
ps_euk_m2 <- stats::aggregate(ps_euk_m1$Abundance, by=list(ps_euk_m1$date, ps_euk_m1$phylum), "sum") # Summed up to make phylum sum
ps_euk_m3 <- stats::aggregate(ps_euk_m1$Abundance, by=list(ps_euk_m1$date, ps_euk_m1$rep_tax), "sum") # Summed up to make phylum sum
colnames(ps_euk_m2) <- c("date", "phylum", "dna_conc")
colnames(ps_euk_m3) <- c("date", "rep_tax", "dna_conc")
ps_euk_m2$phylum <- factor(ps_euk_m2$phylum, levels = c(unique(ps_euk_m2$phylum), "Undetermined"))
ps_euk_m2$phylum[ps_euk_m2$phylum == ""] <- "Undetermined"
ps_euk_m3$rep_tax <- factor(ps_euk_m3$rep_tax, levels = c(unique(ps_euk_m3$rep_tax)[c(1:7, 9:11, 8, 12)]))

# Figures
f7 <- ggplot(ps_euk_m2, aes(x = as.Date(date), y = dna_conc/5, group = phylum, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("DNA (copies/ml water)") +
  ggtitle("All phyla by CMR-002")

f8 <- ggplot(ps_euk_m3, aes(x = as.Date(date), y = dna_conc/5, group = rep_tax, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  xlab(NULL) + ylab("DNA (copies/ml water)") +
  scale_fill_manual(values = palette_euk) +
  ggtitle("Major phyla detected by 18S sequencing") +
  scale_y_continuous(breaks = seq(0, 2.5e+06, 5e+05), limits = c(0, 2.3e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
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

