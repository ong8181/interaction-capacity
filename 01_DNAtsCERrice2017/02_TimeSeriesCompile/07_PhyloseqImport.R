####
#### CER Rice eDNA study
#### No.7 Importing all MiSeq run data to phyloseq objects
####

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder07 <- "07_PhyloseqImportOut"
dir.create(output_folder07)
dir.create("00_SessionInfo")

# Load library and functions
library(ggplot2); packageVersion("ggplot2") #3.2.1, 2020.1.6
library(cowplot); packageVersion("cowplot") #1.0.0, 2020.1.6
library(reshape2); packageVersion("reshape2") #1.4.3, 2019.10.23
library(lubridate); packageVersion("lubridate") #1.7.4, 2019.10.23
library(phyloseq); packageVersion("phyloseq") #1.28.0, 2020.1.6
library(ggsci); packageVersion("ggsci") #2.9, 2019.10.23
theme_set(theme_cowplot())

# Specify data folders
# RMR-076 Prokaryote
prok_data_folder <- "../01_SequenceProcessing/CERrice2017_Prokaryote/05_ReadNSummary_ProkOut/"
# RMR-078 Fungi
fungi_data_folder <- "../01_SequenceProcessing/CERrice2017_Fungi/05_ReadNSummary_FungiOut/"
# RMR-099 Invertebrate
inv_data_folder <- "../01_SequenceProcessing/CERrice2017_Invertebrate/05_ReadNSummary_InvOut/"
# CMR-002 Eukaryote
euk_data_folder <- "../01_SequenceProcessing/CERrice2017_Eukaryote/05_ReadNSummary_EukOut/"

#------ RMR-076 Prokaryote ------#
# Load objects
seqtab_conv_prok <- readRDS(sprintf("%s/seqtab_conv_Prok.obj", prok_data_folder))
sample_sheet_prok <- readRDS(sprintf("%s/sample_sheet_Prok.obj", prok_data_folder))
taxa_wo_std_prok <- readRDS(sprintf("%s/taxa_wo_std_Prok.obj", prok_data_folder))
# Preparetion to import to phyloseq
dim(sample_sheet_prok); dim(taxa_wo_std_prok); dim(seqtab_conv_prok)
all(rownames(sample_sheet_prok) == rownames(seqtab_conv_prok)) # sample name check
taxa_wo_std_prok$seq <- colnames(seqtab_conv_prok) # save sequence info
taxa_wo_std_prok$seqlen <- nchar(colnames(seqtab_conv_prok)) # calculate sequence length
taxa_wo_std_prok$miseq_run <- "RMR-076"
# Change taxa name to group-specific name
taxa_name_prok <- sprintf("Prok_%s", rownames(taxa_wo_std_prok))
colnames(seqtab_conv_prok) <- rownames(taxa_wo_std_prok) <- taxa_name_prok # change col name
all(rownames(taxa_wo_std_prok) == colnames(seqtab_conv_prok)) # taxa name check
# Import data to phyloseq
ps_pro_all <- phyloseq(otu_table(seqtab_conv_prok, taxa_are_rows=FALSE),
                       sample_data(sample_sheet_prok),
                       tax_table(as.matrix(taxa_wo_std_prok)))

#------ RMR-078 Fungi ------#
# Load objects
seqtab_conv_fungi <- readRDS(sprintf("%s/seqtab_conv_Fungi.obj", fungi_data_folder))
sample_sheet_fungi <- readRDS(sprintf("%s/sample_sheet_Fungi.obj", fungi_data_folder))
taxa_wo_std_fungi <- readRDS(sprintf("%s/taxa_wo_std_Fungi.obj", fungi_data_folder))
# Preparetion to import to phyloseq
dim(sample_sheet_fungi); dim(taxa_wo_std_fungi); dim(seqtab_conv_fungi)
all(rownames(sample_sheet_fungi) == rownames(seqtab_conv_fungi)) # sample name check
taxa_wo_std_fungi$seq <- colnames(seqtab_conv_fungi) # save sequence info
taxa_wo_std_fungi$seqlen <- nchar(colnames(seqtab_conv_fungi)) # calculate sequence length
taxa_wo_std_fungi$miseq_run <- "RMR-078"
# Change taxa name to group-specific name
taxa_name_fungi <- sprintf("Fungi_%s", rownames(taxa_wo_std_fungi))
colnames(seqtab_conv_fungi) <- rownames(taxa_wo_std_fungi) <- taxa_name_fungi # change col name
all(rownames(taxa_wo_std_fungi) == colnames(seqtab_conv_fungi)) # taxa name check
# Import data to phyloseq
ps_fun_all <- phyloseq(otu_table(seqtab_conv_fungi, taxa_are_rows=FALSE),
                       sample_data(sample_sheet_fungi),
                       tax_table(as.matrix(taxa_wo_std_fungi)))

#------ RMR-099 Invertebrate ------#
# Load objects
seqtab_conv_inv <- readRDS(sprintf("%s/seqtab_conv_Inv.obj", inv_data_folder))
sample_sheet_inv <- readRDS(sprintf("%s/sample_sheet_Inv.obj", inv_data_folder))
taxa_wo_std_inv <- readRDS(sprintf("%s/taxa_wo_std_Inv.obj", inv_data_folder))
# Preparetion to import to phyloseq
dim(sample_sheet_inv); dim(taxa_wo_std_inv); dim(seqtab_conv_inv)
all(rownames(sample_sheet_inv) == rownames(seqtab_conv_inv)) # sample name check
taxa_wo_std_inv$seq <- colnames(seqtab_conv_inv) # save sequence info
taxa_wo_std_inv$seqlen <- nchar(colnames(seqtab_conv_inv)) # calculate sequence length
taxa_wo_std_inv$miseq_run <- "RMR-099"
# Change taxa name to group-specific name
taxa_name_inv <- sprintf("Inv_%s", rownames(taxa_wo_std_inv))
colnames(seqtab_conv_inv) <- rownames(taxa_wo_std_inv) <- taxa_name_inv # change col name
all(rownames(taxa_wo_std_inv) == colnames(seqtab_conv_inv)) # taxa name check
# Import data to phyloseq
ps_inv_all <- phyloseq(otu_table(seqtab_conv_inv, taxa_are_rows=FALSE),
                         sample_data(sample_sheet_inv),
                         tax_table(as.matrix(taxa_wo_std_inv)))

#------ CMR-002 Eukaryote ------#
# Load objects
seqtab_conv_euk <- readRDS(sprintf("%s/seqtab_conv_Euk.obj", euk_data_folder))
sample_sheet_euk <- readRDS(sprintf("%s/sample_sheet_Euk.obj", euk_data_folder))
taxa_wo_std_euk <- readRDS(sprintf("%s/taxa_wo_std_Euk.obj", euk_data_folder))
# Preparetion to import to phyloseq
dim(sample_sheet_euk); dim(taxa_wo_std_euk); dim(seqtab_conv_euk)
all(rownames(sample_sheet_euk) == rownames(seqtab_conv_euk)) # sample name check
taxa_wo_std_euk$seq <- colnames(seqtab_conv_euk) # save sequence info
taxa_wo_std_euk$seqlen <- nchar(colnames(seqtab_conv_euk)) # calculate sequence length
taxa_wo_std_euk$miseq_run <- "CMR-002"
# Change taxa name to group-specific name
taxa_name_euk <- sprintf("Euk_%s", rownames(taxa_wo_std_euk))
colnames(seqtab_conv_euk) <- rownames(taxa_wo_std_euk) <- taxa_name_euk # change col name
all(rownames(taxa_wo_std_euk) == colnames(seqtab_conv_euk)) # taxa name check
# Import data to phyloseq
ps_euk_all <- phyloseq(otu_table(seqtab_conv_euk, taxa_are_rows=FALSE),
                       sample_data(sample_sheet_euk),
                       tax_table(as.matrix(taxa_wo_std_euk)))


#------ Extracting "sample" only ------#
# Removing field negative/pcr negative/standard negative/positive control samples
ps_pro_sample <- prune_samples(sample_data(ps_pro_all)$sample_nc == "sample", ps_pro_all)
ps_fun_sample <- prune_samples(sample_data(ps_fun_all)$sample_nc == "sample", ps_fun_all)
ps_inv_sample <- prune_samples(sample_data(ps_inv_all)$sample_nc == "sample", ps_inv_all)
ps_euk_sample <- prune_samples(sample_data(ps_euk_all)$sample_nc == "sample", ps_euk_all)


# Save and output results
#save(list = ls(all.names = TRUE),
#     file = sprintf("%s/07_PhyloseqImportOut.RData", output_folder07))
save.image(sprintf("%s/07_PhyloseqImportOut.RData", output_folder07))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/07_SessionInfo_PhyloseqImport_%s.txt", substr(Sys.time(), 1, 10)))

