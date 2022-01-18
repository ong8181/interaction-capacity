####
#### CERrice2017 All data analysis
#### Sequence reads of negative and positive controls
####

# Ceate output director
# Load library and functions
library(phyloseq); packageVersion("phyloseq") # 1.30.1, 2020.1.29
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2020.1.29
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.3
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.29
library(scales); packageVersion("scales") # 1.1.0, 2020.1.29
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.3
library(GGally); packageVersion("GGally") # 1.4.0, 2019.11.3
library(lubridate); packageVersion("lubridate") # 1.7.4, 2019.11.25
theme_set(theme_cowplot())

# Create output directory
fig_output <- "../00_RawFigs/01_Fig_eDNAts"
dir.create(fig_output)


# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
#-------------------- Generate figures for Prokaryote 16S rRNA --------------------#
# Negative controls
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Prokaryote/06_NCcheck_ProkOut/06_NCcheck_ProkOut.RData")
prok_n1 <- ggplot(subset(ps_m3, sample_nc != "pc"), aes(x = sample_nc, y = value, group = variable:sample_nc, color = variable))
prok_n1 <- prok_n1 + geom_boxplot(colour = "black", outlier.colour = "white", outlier.size = 0) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL, labels = c("Standard DNA", "Non standard DNA"))
prok_n1 <- prok_n1 + geom_jitter(alpha = 0.3, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.8))
prok_n1 <- prok_n1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
prok_n1 <- prok_n1 + scale_x_discrete(limits = c("sample", "std_nc", "field_nc", "pcr_nc"),
                                      labels = c("field_nc" = "Field NC", "pcr_nc" = "PCR NC",
                                                 "sample" = "Sample", "std_nc" = "Standard NC"))
prok_n1 <- prok_n1 + ylab("Sequence reads/sample") + xlab(NULL) + ggtitle("MiSeq run: Prokaryote 16S rRNA")

prok_n2 <- ggplot(subset(ps_m10, Description != "PC-5_23-NC"), aes(x = Sample, y = Abundance, colour = phylum, fill = phylum))
prok_n2 <- prok_n2 + geom_bar(stat = "identity", colour = NA) + scale_fill_igv() + xlab(NULL) + ylab("Sequence reads")
prok_n2 <- prok_n2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ggtitle("Non standard DNA: Prokaryote 16S rRNA")

prok_n3 <- ggplot(subset(ps_m11, Description != "PC-5_23-NC"), aes(x = Sample, y = Abundance, colour = species, fill = species))
prok_n3 <- prok_n3 + geom_bar(stat = "identity", colour = NA) + scale_fill_igv() + xlab(NULL) + ylab("Sequence reads")
prok_n3 <- prok_n3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ggtitle("Standard DNA: Prokaryote 16S rRNA")

# Positive controls
dim(sample_sheet); dim(taxa_wo_std); dim(seqtab_conv); all(rownames(sample_sheet) == rownames(seqtab_conv)) # sample name check
taxa_wo_std$seq <- colnames(seqtab_conv) # save sequence info
colnames(seqtab_conv) <- rownames(taxa_wo_std) <- sprintf("Prok_%s", rownames(taxa_wo_std))
# Import data to phyloseq and extract PC samples
ps_pro_pc <- phyloseq(otu_table(seqtab_conv, taxa_are_rows=FALSE),
                       sample_data(sample_sheet),
                       tax_table(as.matrix(taxa_wo_std))) %>%
  subset_samples(sample_nc == "pc" & Description != "PC-5_23-NC") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_taxa(names(sort(taxa_sums(.), TRUE))[1:10], .)
# Add DNA extraction date of PC samples
sample_data(ps_pro_pc)$extraction_date <- ymd("2017-05-23", "2017-06-27", "2017-07-12", "2017-08-09", "2017-08-22",
                                              "2017-09-11", "2017-09-29", "2017-10-11", "2017-10-25", "2017-11-08")
sample_data(ps_pro_pc)$extraction_date2 <- factor("After three months", levels = c("Before three months", "After three months"))
sample_data(ps_pro_pc)$extraction_date2[sample_data(ps_pro_pc)$extraction_date < ymd("2017-09-01")] <- "Before three months"
pro_pc <- psmelt(ps_pro_pc) %>% .[.$Abundance > 0,]
pro_name_order <- aggregate(pro_pc$Abundance, by = list(pro_pc$OTU), median) %>% .[order(.$x, decreasing = T),1]
pro_pc$OTU <- factor(pro_pc$OTU, levels = pro_name_order)
prok_p1 <- ggplot(pro_pc, aes(x = OTU, y = Abundance + 0.5, colour = extraction_date2)) +
  geom_boxplot(colour = "gray20", outlier.colour = "white", outlier.size = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_colour_manual(values =c("red3", "royalblue"), name = NULL) +
  ylab("DNA copy number + 0.5 (copies/ml water)") + xlab(NULL)

saveRDS(prok_n1, "../00_RawFigs/01_Fig_eDNAts/SI_Fig_ProkNC01.obj")
saveRDS(prok_p1, "../00_RawFigs/01_Fig_eDNAts/SI_Fig_ProkPC01.obj")
rm(list = ls())


#-------------------- Generate figures for Eukaryote 18S rRNA --------------------#
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Eukaryote/06_NCcheck_EukOut/06_NCcheck_EukOut.RData")
euk_n1 <- ggplot(subset(ps_m3, sample_nc != "pc"), aes(x = sample_nc, y = value, group = variable:sample_nc, color = variable))
euk_n1 <- euk_n1 + geom_boxplot(colour = "black", outlier.colour = "white", outlier.size = 0) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL, labels = c("Standard DNA", "Non standard DNA"))
euk_n1 <- euk_n1 + geom_jitter(alpha = 0.3, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.8))
euk_n1 <- euk_n1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
euk_n1 <- euk_n1 + scale_x_discrete(limits = c("sample", "std_nc", "field_nc", "pcr_nc"),
                                      labels = c("field_nc" = "Field NC", "pcr_nc" = "PCR NC",
                                                 "sample" = "Sample", "std_nc" = "Standard NC"))
euk_n1 <- euk_n1 + ylab("Sequence reads/sample") + xlab(NULL) + ggtitle("MiSeq run: Eukaryote 18S rRNA")

euk_n2 <- ggplot(subset(ps_m10, Description != "PC-5_23-NC"), aes(x = Sample, y = Abundance, colour = phylum, fill = phylum))
euk_n2 <- euk_n2 + geom_bar(stat = "identity", colour = NA) + scale_fill_igv() + xlab(NULL) + ylab("Sequence reads")
euk_n2 <- euk_n2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ggtitle("Non standard DNA: Eukaryote 18S rRNA")

euk_n3 <- ggplot(subset(ps_m11, Description != "PC-5_23-NC"), aes(x = Sample, y = Abundance, colour = species, fill = species))
euk_n3 <- euk_n3 + geom_bar(stat = "identity", colour = NA) + scale_fill_igv() + xlab(NULL) + ylab("Sequence reads")
euk_n3 <- euk_n3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ggtitle("Standard DNA: Eukaryote 18S rRNA")

# Positive controls
dim(sample_sheet); dim(taxa_wo_std); dim(seqtab_conv); all(rownames(sample_sheet) == rownames(seqtab_conv)) # sample name check
taxa_wo_std$seq <- colnames(seqtab_conv) # save sequence info
colnames(seqtab_conv) <- rownames(taxa_wo_std) <- sprintf("Euk_%s", rownames(taxa_wo_std))
# Import data to phyloseq and extract PC samples
ps_euk_pc <- phyloseq(otu_table(seqtab_conv, taxa_are_rows=FALSE),
                      sample_data(sample_sheet),
                      tax_table(as.matrix(taxa_wo_std))) %>%
  subset_samples(sample_nc == "pc" & Description != "PC-5_23-NC") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_taxa(names(sort(taxa_sums(.), TRUE))[1:10], .)
# Add DNA extraction date of PC samples
sample_data(ps_euk_pc)$extraction_date <- ymd("2017-05-23", "2017-06-27", "2017-07-12", "2017-08-09", "2017-08-22",
                                              "2017-09-11", "2017-09-29", "2017-10-11", "2017-10-25", "2017-11-08")
sample_data(ps_euk_pc)$extraction_date2 <- factor("After three months", levels = c("Before three months", "After three months"))
sample_data(ps_euk_pc)$extraction_date2[sample_data(ps_euk_pc)$extraction_date < ymd("2017-09-01")] <- "Before three months"
euk_pc <- psmelt(ps_euk_pc) %>% .[.$Abundance > 0,]
euk_name_order <- aggregate(euk_pc$Abundance, by = list(euk_pc$OTU), median) %>% .[order(.$x, decreasing = T),1]
euk_pc$OTU <- factor(euk_pc$OTU, levels = euk_name_order)
euk_p1 <- ggplot(euk_pc, aes(x = OTU, y = Abundance + 0.5, colour = extraction_date2)) +
  geom_boxplot(colour = "gray20", outlier.colour = "white", outlier.size = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_colour_manual(values =c("red3", "royalblue"), name = NULL) +
  ylab("DNA copy number + 0.5 (copies/ml water)") + xlab(NULL)

saveRDS(euk_n1, "../00_RawFigs/01_Fig_eDNAts/SI_Fig_EukNC01.obj")
saveRDS(euk_p1, "../00_RawFigs/01_Fig_eDNAts/SI_Fig_EukPC01.obj")
rm(list = ls())


#-------------------- Generate figures for Fungal ITS --------------------#
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Fungi/06_NCcheck_FungiOut/06_NCcheck_FungiOut.RData")
fungi_n1 <- ggplot(subset(ps_m3, sample_nc != "pc"), aes(x = sample_nc, y = value, group = variable:sample_nc, color = variable))
fungi_n1 <- fungi_n1 + geom_boxplot(colour = "black", outlier.colour = "white", outlier.size = 0) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL, labels = c("Standard DNA", "Non standard DNA"))
fungi_n1 <- fungi_n1 + geom_jitter(alpha = 0.3, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.8))
fungi_n1 <- fungi_n1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fungi_n1 <- fungi_n1 + scale_x_discrete(limits = c("sample", "std_nc", "field_nc", "pcr_nc"),
                                    labels = c("field_nc" = "Field NC", "pcr_nc" = "PCR NC",
                                               "sample" = "Sample", "std_nc" = "Standard NC"))
fungi_n1 <- fungi_n1 + ylab("Sequence reads/sample") + xlab(NULL) + ggtitle("MiSeq run: Fungal ITS")

fungi_n2 <- ggplot(subset(ps_m10, Description != "PC-5_23-NC"), aes(x = Sample, y = Abundance, colour = phylum, fill = phylum))
fungi_n2 <- fungi_n2 + geom_bar(stat = "identity", colour = NA) + scale_fill_igv() + xlab(NULL) + ylab("Sequence reads")
fungi_n2 <- fungi_n2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ggtitle("Non standard DNA: Fungal ITS")

fungi_n3 <- ggplot(subset(ps_m11, Description != "PC-5_23-NC"), aes(x = Sample, y = Abundance, colour = species, fill = species))
fungi_n3 <- fungi_n3 + geom_bar(stat = "identity", colour = NA) + scale_fill_igv() + xlab(NULL) + ylab("Sequence reads")
fungi_n3 <- fungi_n3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ggtitle("Standard DNA: Fungal ITS")

# Positive controls
dim(sample_sheet); dim(taxa_wo_std); dim(seqtab_conv); all(rownames(sample_sheet) == rownames(seqtab_conv)) # sample name check
taxa_wo_std$seq <- colnames(seqtab_conv) # save sequence info
colnames(seqtab_conv) <- rownames(taxa_wo_std) <- sprintf("Fungi_%s", rownames(taxa_wo_std))
# Import data to phyloseq and extract PC samples
ps_fungi_pc <- phyloseq(otu_table(seqtab_conv, taxa_are_rows=FALSE),
                      sample_data(sample_sheet),
                      tax_table(as.matrix(taxa_wo_std))) %>%
  subset_samples(sample_nc == "pc" & Description != "PC-5_23-NC") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_taxa(names(sort(taxa_sums(.), TRUE))[1:10], .)
# Add DNA extraction date of PC samples
sample_data(ps_fungi_pc)$extraction_date <- ymd("2017-05-23", "2017-06-27", "2017-07-12", "2017-08-09", "2017-08-22",
                                              "2017-09-11", "2017-09-29", "2017-10-11", "2017-10-25", "2017-11-08")
sample_data(ps_fungi_pc)$extraction_date2 <- factor("After three months", levels = c("Before three months", "After three months"))
sample_data(ps_fungi_pc)$extraction_date2[sample_data(ps_fungi_pc)$extraction_date < ymd("2017-09-01")] <- "Before three months"
fungi_pc <- psmelt(ps_fungi_pc) %>% .[.$Abundance > 0,]
fungi_name_order <- aggregate(fungi_pc$Abundance, by = list(fungi_pc$OTU), median) %>% .[order(.$x, decreasing = T),1]
fungi_pc$OTU <- factor(fungi_pc$OTU, levels = fungi_name_order)
fungi_p1 <- ggplot(fungi_pc, aes(x = OTU, y = Abundance + 0.5, colour = extraction_date2)) +
  geom_boxplot(colour = "gray20", outlier.colour = "white", outlier.size = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_colour_manual(values =c("red3", "royalblue"), name = NULL) +
  ylab("DNA copy number + 0.5 (copies/ml water)") + xlab(NULL)

saveRDS(fungi_n1, "../00_RawFigs/01_Fig_eDNAts/SI_Fig_FungiNC01.obj")
saveRDS(fungi_p1, "../00_RawFigs/01_Fig_eDNAts/SI_Fig_FungiPC01.obj")
rm(list = ls())


#-------------------- Generate figures for Animal COI --------------------#
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Invertebrate/06_NCcheck_InvOut/06_NCcheck_InvOut.RData")
inv_n1 <- ggplot(subset(ps_m3, sample_nc != "pc"), aes(x = sample_nc, y = value, group = variable:sample_nc, color = variable))
inv_n1 <- inv_n1 + geom_boxplot(colour = "black", outlier.colour = "white", outlier.size = 0) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL, labels = c("Standard DNA", "Non standard DNA"))
inv_n1 <- inv_n1 + geom_jitter(alpha = 0.3, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.8))
inv_n1 <- inv_n1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
inv_n1 <- inv_n1 + scale_x_discrete(limits = c("sample", "std_nc", "field_nc", "pcr_nc"),
                                        labels = c("field_nc" = "Field NC", "pcr_nc" = "PCR NC",
                                                   "sample" = "Sample", "std_nc" = "Standard NC"))
inv_n1 <- inv_n1 + ylab("Sequence reads/sample") + xlab(NULL) + ggtitle("MiSeq run: Animal COI")

inv_n2 <- ggplot(subset(ps_m10, Description != "PC-5_23-NC"), aes(x = Sample, y = Abundance, colour = phylum, fill = phylum))
inv_n2 <- inv_n2 + geom_bar(stat = "identity", colour = NA) + scale_fill_igv() + xlab(NULL) + ylab("Sequence reads")
inv_n2 <- inv_n2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ggtitle("Non standard DNA: Animal COI")

inv_n3 <- ggplot(subset(ps_m11, Description != "PC-5_23-NC"), aes(x = Sample, y = Abundance, colour = species, fill = species))
inv_n3 <- inv_n3 + geom_bar(stat = "identity", colour = NA) + scale_fill_igv() + xlab(NULL) + ylab("Sequence reads")
inv_n3 <- inv_n3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ggtitle("Standard DNA: Animal COI")

# Positive controls
dim(sample_sheet); dim(taxa_wo_std); dim(seqtab_conv); all(rownames(sample_sheet) == rownames(seqtab_conv)) # sample name check
taxa_wo_std$seq <- colnames(seqtab_conv) # save sequence info
colnames(seqtab_conv) <- rownames(taxa_wo_std) <- sprintf("Inv_%s", rownames(taxa_wo_std))
# Import data to phyloseq and extract PC samples
ps_inv_pc <- phyloseq(otu_table(seqtab_conv, taxa_are_rows=FALSE),
                        sample_data(sample_sheet),
                        tax_table(as.matrix(taxa_wo_std))) %>%
  subset_samples(sample_nc == "pc" & Description != "PC-5_23-NC") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_taxa(names(sort(taxa_sums(.), TRUE))[1:10], .)
sample_data(ps_inv_pc)$extraction_date <- ymd("2017-05-23", "2017-06-27", "2017-07-12", "2017-08-09", "2017-08-22",
                                                "2017-09-11", "2017-09-29", "2017-10-11", "2017-10-25", "2017-11-08")
sample_data(ps_inv_pc)$extraction_date2 <- factor("After three months", levels = c("Before three months", "After three months"))
sample_data(ps_inv_pc)$extraction_date2[sample_data(ps_inv_pc)$extraction_date < ymd("2017-09-01")] <- "Before three months"
inv_pc <- psmelt(ps_inv_pc) %>% .[.$Abundance > 0,]
inv_name_order <- aggregate(inv_pc$Abundance, by = list(inv_pc$OTU), median) %>% .[order(.$x, decreasing = T),1]
inv_pc$OTU <- factor(inv_pc$OTU, levels = inv_name_order)
inv_p1 <- ggplot(inv_pc, aes(x = OTU, y = Abundance + 0.5, colour = extraction_date2)) +
  geom_boxplot(colour = "gray20", outlier.colour = "white", outlier.size = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_colour_manual(values =c("red3", "royalblue"), name = NULL) +
  ylab("DNA copy number + 0.5 (copies/ml water)") + xlab(NULL)

saveRDS(inv_n1, "../00_RawFigs/01_Fig_eDNAts/SI_Fig_InvNC01.obj")
saveRDS(inv_p1, "../00_RawFigs/01_Fig_eDNAts/SI_Fig_InvPC01.obj")
