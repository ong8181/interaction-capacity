####
#### CER eDNA study
#### FigCode: Pattern visualization after data filtering
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
load("../../02_EcolComAnalysis/07_RegularizedSmapOut/07_RegularizedSmapOut.RData")

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.29
library(lubridate); packageVersion("lubridate") # 1.7.4, 2019.11.11
library(phyloseq); packageVersion("phyloseq") # 1.30.0, 2020.1.29
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.29
library(scales); packageVersion("scales") # 1.1.0, 2020.1.14
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.11
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.11
theme_set(theme_cowplot())
source("../00_Fig_functions/F01_FigHelperFunctions.R")

# Create output folder
fig_output <- "../00_RawFigs"
dir.create(fig_output)

# Update taxa information
ps_filt <- readRDS("../../01_DNAtsCERrice2017/02_TimeSeriesCompile/11_TSfilter02Out/ps_comb_filt.obj")
all(rownames(tax_table(ps_filt)) == rownames(edna_tax3))
dim(tax_table(ps_filt)); dim(tax_table(as.matrix(edna_tax3)))
colnames(tax_table(ps_filt)); colnames(tax_table(as.matrix(edna_tax3)))
tax_table(ps_filt) <- tax_table(as.matrix(edna_tax3))

# Superkingdom visualization
ps_filt_m1 <- psmelt(ps_filt)
ps_filt_m1$date <- mdy(ps_filt_m1$date)
ps_filt_m2 <- stats::aggregate(ps_filt_m1$Abundance, by=list(ps_filt_m1$date, ps_filt_m1$superkingdom), "sum") # Summed up to make superkingdom sum
colnames(ps_filt_m2) <- c("date", "superkingdom", "dna_conc")
ps_filt_m2$superkingdom <- as.character(ps_filt_m2$superkingdom)
ps_filt_m2[ps_filt_m2$superkingdom == "", "superkingdom"] <- "Undetermined"

# Divide by 5 to make average copy number / plot
f1 <- ggplot(ps_filt_m2, aes(x = as.Date(date), y = dna_conc/5, group = superkingdom, fill = superkingdom)) +
  geom_bar(stat = "identity", colour = NA) +
  scale_fill_startrek() +
  xlab(NULL) + ylab("DNA (copies/ml water)") + theme(legend.position = c(.2, .7)) +
  theme(legend.title = element_blank()) +
  scale_y_continuous(breaks = seq(0, 1e+07, 2e+06), limits = c(0, 6.1e+06),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})

# Diversity visualization
ps_filt_div <- ps_filt %>% transform_sample_counts(ceiling)
sample_data(ps_filt_div)$date <- mdy(sample_data(ps_filt_div)$date)
f2 <- plot_richness(ps_filt_div, x = "date", color = "plot", shape = "plot", measures = c("Observed"))
f2$layers <- f2$layers[-1]
f2 <- f2 + geom_line(alpha = 0.8) + geom_point(aes(fill = plot), size = 1.5, alpha = 0.8) +
  scale_color_startrek() + scale_fill_startrek() + scale_shape_manual(values = c(21:25)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  ylab("No. of ASV") + xlab(NULL) + ggtitle("ASV diversity (1196 ASVs)") # After filtered, No. of ASV = 1002
f2 <- f2 + ylim(0,420)

# Visualize diversity of subset taxa
f3 <- subset_taxa(ps_filt_div, superkingdom == "") %>% 
  plot_richness(x = "date", color = "plot", shape = "plot", measures = c("Observed"))
f3$layers <- f3$layers[-1]
f3 <- f3 + geom_line(alpha = 0.8) + geom_point(aes(fill = plot), size = 1.5, alpha = 0.8) +
  scale_color_startrek() + scale_fill_startrek() + scale_shape_manual(values = c(21:25)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  ylab("No. of ASV") + xlab(NULL) + ggtitle("Observed diversity (Undetermined; 127 ASVs)") +
  ylim(0,50)

f4 <- subset_taxa(ps_filt_div, superkingdom == "Bacteria" | superkingdom == "Archaea") %>%
  plot_richness(x = "date", color = "plot", shape = "plot", measures = c("Observed"))
f4$layers <- f4$layers[-1]
f4 <- f4 + geom_line(alpha = 0.8) + geom_point(aes(fill = plot), size = 1.5, alpha = 0.8) +
  scale_color_startrek() + scale_fill_startrek() + scale_shape_manual(values = c(21:25)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  ylab("No. of ASV") + xlab(NULL) + ggtitle("Prokaryotic diversity (total 483 ASVs)")
f4 <- f4 + ylim(0,190)

f5 <- subset_taxa(ps_filt_div, kingdom == "Fungi") %>%
  plot_richness(x = "date", color = "plot", shape = "plot", measures = c("Observed"))
f5$layers <- f5$layers[-1]
f5 <- f5 + geom_line(alpha = 0.8) + geom_point(aes(fill = plot), size = 1.5, alpha = 0.8)
f5 <- f5 + scale_color_startrek() + scale_fill_startrek() + scale_shape_manual(values = c(21:25))
f5 <- f5 + theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
f5 <- f5 + ylab("No. of ASV") + xlab(NULL) + ggtitle("Fungal diversity (total 239 ASVs)")
f5 <- f5 + ylim(0,120)

f6 <- subset_taxa(ps_filt_div, kingdom == "Metazoa") %>%
  plot_richness(x = "date", color = "plot", shape = "plot", measures = c("Observed"))
f6$layers <- f6$layers[-1]
f6 <- f6 + geom_line(alpha = 0.8) + geom_point(aes(fill = plot), size = 1.5, alpha = 0.8)
f6 <- f6 + scale_color_startrek() + scale_fill_startrek() + scale_shape_manual(values = c(21:25))
f6 <- f6 + theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
f6 <- f6 + ylab("No. of ASV") + xlab(NULL) + ggtitle("Metazoan diversity (total 59 ASVs)")
f6 <- f6 + ylim(0,35)

f7 <- subset_taxa(ps_filt_div, kingdom == "Viridiplantae") %>%
  plot_richness(x = "date", color = "plot", shape = "plot", measures = c("Observed"))
f7$layers <- f7$layers[-1]
f7 <- f7 + geom_line(alpha = 0.8) + geom_point(aes(fill = plot), size = 1.5, alpha = 0.8)
f7 <- f7 + scale_color_startrek() + scale_fill_startrek() + scale_shape_manual(values = c(21:25))
f7 <- f7 + theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
f7 <- f7 + ylab("No. of ASV") + xlab(NULL) + ggtitle("Viridiplantae diversity (total 57 ASVs)")
f7 <- f7 + ylim(0,30)

# Combine all figures
legend_bar <- get_legend(f1)
legend_div <- get_legend(f2)

f_all1 <- plot_grid(f1, f2, ncol = 2, labels = c("a", "b"), rel_widths = c(1,0.9),
                    align = "h", axis = "bt")
f_all2 <- plot_grid(f4 + theme(legend.position = "none"),
                    f5 + theme(legend.position = "none"),
                    f6 + theme(legend.position = "none"),
                    f7 + theme(legend.position = "none"),
                    ncol = 2, labels = c("c", "d", "e", "f"), align = "hv")
f_all <- plot_grid(f_all1, f_all2, ncol = 1, labels = NULL, rel_heights = c(1,1.5))

pdf(sprintf("%s/01_Fig_eDNAts/Fig_AbnDiv_summary.pdf", fig_output), width = 14, height = 5)
f_all1; dev.off()

pdf(sprintf("%s/01_Fig_eDNAts/Fig_AbnDiv_subset.pdf", fig_output), width = 14, height = 7)
f_all2; dev.off()

pdf(sprintf("%s/01_Fig_eDNAts/Fig_AbnDiv_all.pdf", fig_output), width = 14, height = 12)
f_all; dev.off()


# Save all figures as R objects
saveRDS(f_all1, "../00_RawFigs/01_Fig_eDNAts/Fig_AbnDiv_summary.obj")
saveRDS(f_all2, "../00_RawFigs/01_Fig_eDNAts/Fig_AbnDiv_subset.obj")
saveRDS(f_all, "../00_RawFigs/01_Fig_eDNAts/Fig_AbnDiv_all.obj")
saveRDS(f1, "../00_RawFigs/01_Fig_eDNAts/Fig_AllBar.obj")
saveRDS(f2, "../00_RawFigs/01_Fig_eDNAts/Fig_AllDivTS.obj")
saveRDS(f3, "../00_RawFigs/01_Fig_eDNAts/Fig_UndetDivTS.obj")
saveRDS(f4, "../00_RawFigs/01_Fig_eDNAts/Fig_ProkDivTS.obj")
saveRDS(f5, "../00_RawFigs/01_Fig_eDNAts/Fig_FungiDivTS.obj")
saveRDS(f6, "../00_RawFigs/01_Fig_eDNAts/Fig_MetazDivTS.obj")
saveRDS(f7, "../00_RawFigs/01_Fig_eDNAts/Fig_ViridiDivTS.obj")
