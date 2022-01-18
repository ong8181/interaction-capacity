####
#### CER eDNA study
#### No.8 Filtering time series preparation (Shannon entropy + DNA concentration)
####

# Load workspace
load("07_PhyloseqImportOut/07_PhyloseqImportOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder08 <- "08_TSfilterPrepOut"
dir.create(output_folder08)

# Load library and functions
library(ggplot2); packageVersion("ggplot2") #3.2.1, 2020.1.6
library(cowplot); packageVersion("cowplot") #1.0.0, 2020.1.6
library(reshape2); packageVersion("reshape2") #1.4.3, 2019.10.23
library(lubridate); packageVersion("lubridate") #1.7.4, 2019.10.23
library(phyloseq); packageVersion("phyloseq") #1.28.0, 2020.1.6
library(ggsci); packageVersion("ggsci") #2.9, 2019.10.23
library(infotheo); packageVersion("infotheo") #1.2.0, 2019.10.23
theme_set(theme_cowplot())

# Extracting data as time series
ps_pro_sample1 <- prune_taxa(taxa_sums(ps_pro_sample) > 0, ps_pro_sample)
ps_fun_sample1 <- prune_taxa(taxa_sums(ps_fun_sample) > 0, ps_fun_sample)
ps_inv_sample1 <- prune_taxa(taxa_sums(ps_inv_sample) > 0, ps_inv_sample)
ps_euk_sample1 <- prune_taxa(taxa_sums(ps_euk_sample) > 0, ps_euk_sample)

# Remove "5/22/17" samples because it is preliminary sampling
ps_pro_sample2 <- subset_samples(ps_pro_sample1, sample_data(ps_pro_sample1)[,"date"] != "5/22/17")
ps_fun_sample2 <- subset_samples(ps_fun_sample1, sample_data(ps_fun_sample1)[,"date"] != "5/22/17")
ps_inv_sample2 <- subset_samples(ps_inv_sample1, sample_data(ps_inv_sample1)[,"date"] != "5/22/17")
ps_euk_sample2 <- subset_samples(ps_euk_sample1, sample_data(ps_euk_sample1)[,"date"] != "5/22/17")

ps_pro_sample2 <- prune_taxa(taxa_sums(ps_pro_sample2) > 0, ps_pro_sample2)
ps_fun_sample2 <- prune_taxa(taxa_sums(ps_fun_sample2) > 0, ps_fun_sample2)
ps_inv_sample2 <- prune_taxa(taxa_sums(ps_inv_sample2) > 0, ps_inv_sample2)
ps_euk_sample2 <- prune_taxa(taxa_sums(ps_euk_sample2) > 0, ps_euk_sample2)

# Check the number of taxa detected
(tax_n_pro <- ncol(otu_table(ps_pro_sample2))) # 2620 (before DADA2 revise on 23 Oct 2019; 2979)
(tax_n_fun <- ncol(otu_table(ps_fun_sample2))) # 7412 (before DADA2 revise on 23 Oct 2019; 4210)
(tax_n_inv <- ncol(otu_table(ps_inv_sample2))) # 2508 (before DADA2 revise on 23 Oct 2019; 2517)
(tax_n_euk <- ncol(otu_table(ps_euk_sample2))) # 1405 (before DADA2 revise on 23 Oct 2019; 1344)
sum(tax_n_pro, tax_n_fun, tax_n_inv, tax_n_euk)
# 13945 taxa from 610 samples * 30-100 ml water!
# (before 23 Oct 2019; 11050 taxa from 610 samples * 30-100 ml water!)

# Define entropy estimation function
entropy_ps <- function(ts_object){
  n_bin_for_ts <- ceiling(sqrt(length(ts_object)))
  binned_ts <- discretize(c(ts_object), disc = "equalwidth", nbins = n_bin_for_ts)
  entropy_estimation <- infotheo::entropy(binned_ts)
  return(entropy_estimation)
}

pro_ent <- apply(data.frame(otu_table(ps_pro_sample2)), 2, entropy_ps)
fun_ent <- apply(data.frame(otu_table(ps_fun_sample2)), 2, entropy_ps)
inv_ent <- apply(data.frame(otu_table(ps_inv_sample2)), 2, entropy_ps)
euk_ent <- apply(data.frame(otu_table(ps_euk_sample2)), 2, entropy_ps)

# Set the lowest DNA copy numbers of standard DNAs to illustrate them
ll_pro <- 5000
ll_fun <- 10
ll_inv <- 5
ll_euk <- 250

# Calculate quantiles of entropy
quantile(pro_ent, probs = seq(0,1,0.05))
quantile(fun_ent, probs = seq(0,1,0.05))
quantile(inv_ent, probs = seq(0,1,0.05))
quantile(euk_ent, probs = seq(0,1,0.05))
quantile(c(pro_ent, fun_ent, inv_ent, euk_ent), probs = seq(0,1,0.05))
q90 <- quantile(c(pro_ent, fun_ent, inv_ent, euk_ent), probs = seq(0,1,0.05))['90%']
q95 <- quantile(c(pro_ent, fun_ent, inv_ent, euk_ent), probs = seq(0,1,0.05))['95%']

# Compile entropy information
g8_1 <- ggplot(NULL, aes(x = taxa_sums(ps_pro_sample2)+0.5, y = pro_ent))
g8_1 <- g8_1 + geom_point(size = 2, alpha = 0.3) + scale_x_log10(limits = c(0.03, 2.5e+7)) + ylim(0,1.7)
g8_1 <- g8_1 + geom_vline(xintercept = c(ll_pro/100)*nrow(sample_data(ps_pro_sample2)), linetype = 2)
g8_1 <- g8_1 + geom_hline(yintercept = q90, linetype = 2)
g8_1 <- g8_1 + xlab("Log10(Sum of DNA copy estimate + 0.5)") + ylab("Shannon entropy")
g8_1 <- g8_1 + ggtitle("Prokaryote (N = 2620):\nDNA copy and entropy")

g8_2 <- ggplot(NULL, aes(x = taxa_sums(ps_fun_sample2)+0.5, y = fun_ent))
g8_2 <- g8_2 + geom_point(size = 2, alpha = 0.3) + scale_x_log10(limits = c(0.03, 2.5e+7)) + ylim(0,1.7)
g8_2 <- g8_2 + geom_vline(xintercept = c(ll_fun/100)*nrow(sample_data(ps_fun_sample2)), linetype = 2)
g8_2 <- g8_2 + geom_hline(yintercept = q90, linetype = 2)
g8_2 <- g8_2 + xlab("Log10(Sum of DNA copy estimate + 0.5)") + ylab("Shannon entropy")
g8_2 <- g8_2 + ggtitle("Fungi (N=7412):\nDNA copy and entropy")

g8_3 <- ggplot(NULL, aes(x = taxa_sums(ps_inv_sample2)+0.5, y = inv_ent))
g8_3 <- g8_3 + geom_point(size = 2, alpha = 0.3) + scale_x_log10(limits = c(0.03, 2.5e+7)) + ylim(0,1.7)
g8_3 <- g8_3 + geom_vline(xintercept = c(ll_inv/100)*nrow(sample_data(ps_inv_sample2)), linetype = 2)
g8_3 <- g8_3 + geom_hline(yintercept = q90, linetype = 2)
g8_3 <- g8_3 + xlab("Log10(Sum of DNA copy estimate + 0.5)") + ylab("Shannon entropy")
g8_3 <- g8_3 + ggtitle("Invertebrate (N=2508):\nDNA copy and entropy")

g8_4 <- ggplot(NULL, aes(x = taxa_sums(ps_euk_sample2)+0.5, y = euk_ent))
g8_4 <- g8_4 + geom_point(size = 2, alpha = 0.3) + scale_x_log10(limits = c(0.03, 2.5e+7)) + ylim(0,1.7)
g8_4 <- g8_4 + geom_vline(xintercept = c(ll_euk/100)*nrow(sample_data(ps_euk_sample2)), linetype = 2)
g8_4 <- g8_4 + geom_hline(yintercept = q90, linetype = 2)
g8_4 <- g8_4 + xlab("Log10(Sum of DNA copy estimate + 0.5)") + ylab("Shannon entropy")
g8_4 <- g8_4 + ggtitle("Eukaryote (N=1405):\nDNA copy and entropy")

g8_all <- plot_grid(g8_1, g8_2, g8_3, g8_4, ncol = 2, labels = "auto", align = "hv")

ggsave(sprintf("%s/DNAcopyEntropy.pdf", output_folder08), plot = g8_all, width = 10, height = 10)

# Save and output results
#save(list = ls(all.names = TRUE),
#     file = sprintf("%s/08_TSfilterPrepOut.RData", output_folder08))
save.image(sprintf("%s/08_TSfilterPrepOut.RData", output_folder08))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/08_SessionInfo_TSfilterPrep_%s.txt", substr(Sys.time(), 1, 10)))
