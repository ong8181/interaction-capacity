####
#### CER eDNA study
#### No.10 Combine all phyloseq objects
####

# Load workspace
load("09_TSfilter01Out/09_TSfilter01Out.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder10 <- "10_TSfilterCombineOut"
dir.create(output_folder10)

# Load library
library(phyloseq); packageVersion("phyloseq") #1.28.0, 2020.1.6

# Combine all phyloseq object
# Prepare sample sheet combined
sample_meta_all <- as.data.frame(sample_data(ps_pro_sample3)[,c(3,9:15)])
sample_meta_pro <- as.data.frame(sample_data(ps_pro_sample3)[,16:23])
sample_meta_fun <- as.data.frame(sample_data(ps_fun_sample3)[,16:23])
sample_meta_inv <- as.data.frame(sample_data(ps_inv_sample3)[,16:23])
sample_meta_euk <- as.data.frame(sample_data(ps_euk_sample3)[,16:23])

meta_colnames <- colnames(sample_sheet_prok)[16:23]
colnames(sample_meta_pro) <- sprintf("Prok_%s", meta_colnames)
colnames(sample_meta_fun) <- sprintf("Fungi_%s", meta_colnames)
colnames(sample_meta_inv) <- sprintf("Inv_%s", meta_colnames)
colnames(sample_meta_euk) <- sprintf("Euk_%s", meta_colnames)

sample_combined <- cbind(sample_meta_all, sample_meta_pro, sample_meta_fun, sample_meta_inv, sample_meta_euk)

# Pre-combined
ps_combined0 <- merge_phyloseq(ps_pro_sample3,
                               ps_fun_sample3,
                               ps_inv_sample3,
                               ps_euk_sample3)

# Re-order taxtable information
potential_taxcol_names <- c("query",
                            "superkingdom", "kingdom", "subkingdom",
                            "phylum",
                            "class", "subclass", "infraclass",
                            "superorder", "order", "suborder", "infraorder", "parvorder",
                            "superfamily", "family","subfamily",
                            "tribe", "subtribe",
                            "genus", "subgenus",
                            "species", "subspecies", "forma",
                            "seq", "seqlen", "entropy", "miseq_run")
tax_table(ps_combined0) <- tax_table(ps_combined0)[,potential_taxcol_names]

# Re-merge phyloseq object (this process is done to keep sample information)
ps_combined <- phyloseq(otu_table(ps_combined0),
                        sample_data(sample_combined), # Keep all sample information
                        tax_table(ps_combined0))

# Save and output results
#save(list = ls(all.names = TRUE),
#     file = sprintf("%s/10_TSfilterCombineOut.RData", output_folder10))
save.image(sprintf("%s/10_TSfilterCombineOut.RData", output_folder10))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/10_SessionInfo_TSfilterCombine_%s.txt", substr(Sys.time(), 1, 10)))
