####
#### CER eDNA study
#### No.6 RMR-099 Invertebrate: Negative & Positive control check
####

# Load workspace and data
load("05_ReadNSummary_InvOut/05_ReadNSummary_InvOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder06 <- "06_NCcheck_InvOut"
dir.create(output_folder06)

# Load library and functions
library(ggplot2); packageVersion("ggplot2") #3.2.1, 2020.1.6
library(cowplot); packageVersion("cowplot") #1.0.0, 2020.1.6
library(reshape2); packageVersion("reshape2") #1.4.3, 2019.10.22
library(lubridate); packageVersion("lubridate") #1.7.4, 2019.10.22
library(phyloseq); packageVersion("phyloseq") #1.28.0, 2020.1.6
library(ggsci); packageVersion("ggsci") #2.9, 2019.10.22
theme_set(theme_cowplot())

# Preparetion to import to phyloseq
dim(sample_sheet); dim(tax_claident); dim(seqtab_nochim)
all(rownames(sample_sheet) == rownames(seqtab_nochim)) # sample name check

tax_claident2$seq <- colnames(seqtab_nochim) # save sequence info
tax_claident2$seqlen <- nchar(colnames(seqtab_nochim)) # calculate sequence length
colnames(seqtab_nochim) <- rownames(tax_claident2) # change col name
all(rownames(tax_claident2) == colnames(seqtab_nochim)) # taxa name check

tax_claident2$std_or_field <- "Field DNA"
tax_claident2[substr(tax_claident2$species, 1, nchar(std_seq_head)) == std_seq_head,]$std_or_field <- "Standard DNA"

# Import data to phyloseq
ps0 <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE),
                sample_data(sample_sheet),
                tax_table(as.matrix(tax_claident2)))

# Generating table for manual NC check
sample_sheet_na <- matrix(rep(NaN, ncol(tax_claident2)*ncol(sample_sheet)), ncol=ncol(sample_sheet))
colnames(sample_sheet_na) <- colnames(sample_sheet)
sample_sheet_comb <- rbind(sample_sheet_na, as.matrix(sample_sheet))
tax_table_comb <- rbind(t(as.matrix(tax_claident2)), seqtab_nochim)
table_for_nc_check <- cbind(sample_sheet_comb, tax_table_comb)
table_for_nc_check[,"std_validity"][1:ncol(tax_claident2)] <- colnames(tax_claident2)
write.csv(table_for_nc_check, sprintf("%s/SummarizedTable_for_NCcheck.csv", output_folder06))

# Visualize pattern
ps_m3 <- melt(sample_summary, measure.vars = c("STD_all", "NonSTD_all"),
              id.vars = c("sample_nc"))

p2 <- ggplot(ps_m3, aes(x = sample_nc, y = value, fill = variable, group = variable:sample_nc))
p2 <- p2 + geom_boxplot(colour = "black") + scale_fill_nejm() #+ scale_y_log10()
p2 <- p2 + ylab("Sequence reads/sample") + xlab(NULL)
p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2 <- p2 + scale_x_discrete(limits = c("sample", "pc", "std_nc", "field_nc", "pcr_nc"),
                            labels = c("field_nc" = "Field NC",
                                       "pc" = "PC",
                                       "pcr_nc" = "PCR NC",
                                       "sample" = "Sample",
                                       "std_nc" = "Standard NC"))
ggsave(sprintf("%s/SequenceReads_Barplot.pdf", output_folder06), plot = p2, width = 6, height = 6)

# Examination of ASVs that were detected in field NC, PCR NC and standard NC
# Extract taxa that are detected from NC samples
# Barplot of NC-detected taxa (Optional)
ps_nctax <- prune_taxa(taxa_sums(prune_samples(sample_data(ps0)$sample_nc != "sample" & sample_data(ps0)$sample_nc != "pc", ps0)) > 0, ps0)
p4 <- plot_bar(ps_nctax, x = "sample_nc", fill = "phylum")
p4 <- p4 + geom_bar(stat = "identity", colour = NA) + scale_fill_igv()

# Boxplot for NC-detected taxa
ps_m4 <- psmelt(ps_nctax)
ps_m5 <- ps_m4[ps_m4$sample_nc != "pc" & ps_m4$sample_nc != "sample",]
ps_m5 <- ps_m5[ps_m5$std_or_field == "Field DNA" & ps_m5$Abundance > 0,]
p5 <- ggplot(ps_m5, aes(x = sample_nc, y = Abundance, group = phylum:sample_nc, fill = phylum))
p5 <- p5 + geom_boxplot(width = 0.6, position = position_dodge(width=0.7)) + scale_y_log10() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p5 <- p5 + scale_fill_igv() + xlab(NULL) + ylab("Log10(Sequence reads)")
p5 <- p5 + scale_x_discrete(limits = c("field_nc", "std_nc", "pcr_nc"),
                            labels = c("field_nc" = "Field NC",
                                       "pcr_nc" = "PCR NC",
                                       "std_nc" = "Standard NC"))
ggsave(sprintf("%s/SequenceReads_NCtaxaBoxplot.pdf", output_folder06), plot = p5, width = 8, height = 6)

# Examination of positive controls and field negative controls
# (Time series plot)
# (Individual sample plot: Field negative controls)
ps_nctax2 <- prune_taxa(taxa_sums(prune_samples(sample_data(ps0)$sample_nc == "field_nc" | sample_data(ps0)$sample_nc == "pc", ps0)) > 0, ps0)
ps_m7 <- psmelt(ps_nctax2)
if(class(ps_m7$date) != "Date") ps_m7$date <- mdy(ps_m7$date)
ps_m8 <- ps_m7[ps_m7$sample_nc == "field_nc" & ps_m7$std_or_field == "Field DNA",]
p7 <- ggplot(ps_m8, aes(x = date, y = Abundance, fill = phylum))
p7 <- p7 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_igv()
p8 <- p7 + xlim(as.Date(c("2017-05-23", "2017-06-07"))) + ylim(0,15000)
p9 <- p7 + xlim(as.Date(c("2017-06-08", "2017-09-23"))) + ylim(0,15000)

ps_m9 <- ps_m7[ps_m7$sample_nc == "field_nc" & ps_m7$std_or_field == "Standard DNA",]
p10 <- ggplot(ps_m9, aes(x = date, y = Abundance, fill = species))
p10 <- p10 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_igv()
p11 <- p10 + xlim(as.Date(c("2017-05-23", "2017-06-07"))) + ylim(0,30000)
p12 <- p10 + xlim(as.Date(c("2017-06-08", "2017-09-23"))) + ylim(0,30000)

field_nc_reads <- plot_grid(p8, p9, p11, p12, ncol = 2, labels = "auto", align = "hv")
ggsave(sprintf("%s/SequenceReads_FieldNCreads.pdf", output_folder06), plot = field_nc_reads, width = 12, height = 8)


# (Individual sample plot: Positive control)
ps_m10 <- ps_m7[ps_m7$sample_nc == "pc" & ps_m7$std_or_field == "Field DNA",]
p13 <- ggplot(ps_m10, aes(x = Sample, y = Abundance, colour = phylum, fill = phylum))
p13 <- p13 + geom_bar(stat = "identity", colour = NA) + scale_fill_igv()
p13 <- p13 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylim(0,30000)

ps_m11 <- ps_m7[ps_m7$sample_nc == "pc" & ps_m7$std_or_field == "Standard DNA",]
p14 <- ggplot(ps_m11, aes(x = Sample, y = Abundance, colour = species, fill = species))
p14 <- p14 + geom_bar(stat = "identity", colour = NA) + scale_fill_igv()
p14 <- p14 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylim(0,12000)

pc_reads <- plot_grid(p13, p14, ncol = 1, align = "hv")
ggsave(sprintf("%s/SequenceReads_PCreads.pdf", output_folder06), plot = pc_reads, width = 6, height = 6)


# Examination of PCR negative controls and standard negative controls
# (Individual sample plot: PCR negative controls)
ps_nctax3 <- prune_taxa(taxa_sums(prune_samples(sample_data(ps0)$sample_nc == "std_nc" | sample_data(ps0)$sample_nc == "pcr_nc", ps0)) > 0, ps0)
ps_m12 <- psmelt(ps_nctax3)
if(class(ps_m12$date) != "Date") ps_m12$date <- mdy(ps_m12$date)
ps_m13 <- ps_m12[ps_m12$sample_nc == "pcr_nc" & ps_m12$std_or_field == "Field DNA",]
p15 <- ggplot(ps_m13, aes(x = Sample, y = Abundance, fill = phylum))
p15 <- p15 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_igv() + ylim(0,15000)
p15 <- p15 + ggtitle("PCR negative controls: Field DNA reads")

ps_m14 <- ps_m12[ps_m12$sample_nc == "pcr_nc" & ps_m12$std_or_field == "Standard DNA",]
p16 <- ggplot(ps_m14, aes(x = Sample, y = Abundance, fill = species))
p16 <- p16 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_igv() + ylim(0,15000)
p16 <- p16 + ggtitle("PCR negative controls: Standard DNA reads")

pcr_nc_reads <- plot_grid(p15, p16, ncol = 1, labels = "auto", align = "hv")
ggsave(sprintf("%s/SequenceReads_PCRNCreads.pdf", output_folder06), plot = pcr_nc_reads, width = 8, height = 8)


# (Individual sample plot: Standard negative controls)
ps_m15 <- ps_m12[ps_m12$sample_nc == "std_nc" & ps_m12$std_or_field == "Field DNA",]
p17 <- ggplot(ps_m15, aes(x = Sample, y = Abundance, fill = phylum))
p17 <- p17 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_igv() + ylim(0,15000)
p17 <- p17 + ggtitle("Standard negative controls: Field DNA reads")

ps_m16 <- ps_m12[ps_m12$sample_nc == "std_nc" & ps_m12$std_or_field == "Standard DNA",]
p18 <- ggplot(ps_m16, aes(x = Sample, y = Abundance, fill = species))
p18 <- p18 + geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90)) + scale_fill_igv() + ylim(0,30000)
p18 <- p18 + ggtitle("Standard negative controls: Standard DNA reads")

std_nc_reads <- plot_grid(p17, p18, ncol = 1, labels = "auto", align = "hv")
ggsave(sprintf("%s/SequenceReads_STDNCreads.pdf", output_folder06), plot = std_nc_reads, width = 8, height = 8)


# Calculating ratio
dna_summary <- aggregate(sample_summary$NonSTD_all, list(sample_summary$sample_nc), mean)
dna_summary$prop <- dna_summary$x/dna_summary[dna_summary$Group.1 == "sample","x"]
write.csv(dna_summary, sprintf("%s/DNAreads_ContaminationLevel.csv", output_folder06))

# Save and output results
#save(list = ls(all.names = TRUE),
#     file = sprintf("%s/06_NCcheck_InvOut.RData", output_folder06))
save.image(sprintf("%s/06_NCcheck_InvOut.RData", output_folder06))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/06_SessionInfo_NCcheck_%s.txt", substr(Sys.time(), 1, 10)))

