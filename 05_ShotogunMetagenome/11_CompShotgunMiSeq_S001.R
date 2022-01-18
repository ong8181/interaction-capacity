####
#### Compare the resuls of phyloFlash and qMiSeq
#### S003
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2021.7.23
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.7.23
library(ggsci); packageVersion("ggsci") # 2.9, 2021.7.23
library(phyloseq); packageVersion("phyloseq") # 1.36.0, 2021.7.23
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.7.23
theme_set(theme_cowplot())
get_palette <- colorRampPalette(brewer.pal(8, "Paired"))


# Generate output folder
(fn <- basename(rstudioapi::getSourceEditorContext()$path))
output_folder <- paste0(str_sub(fn, end = -3), "Out"); rm(fn)
dir.create(output_folder)


# ---------------------------------------------- #
# Set parameters
# ---------------------------------------------- #
Sxxx <- "S001"
Sample_Description <- "P1-7_14"
Sample_Date <- "2017/7/14"
Sample_Plot <- "1"
n_taxa_col <- 18
min_reads <- 3

# Generate palette
palette_manual <- get_palette(n_taxa_col)
palette_manual[length(palette_manual)] <- "#808080"

# ---------------------------------------------- #
# Load taxa-rename data
# ---------------------------------------------- #
rename_df <- read_csv("data/rename_taxa1.csv")


# ---------------------------------------------- #
# Load phyloFlash data
# ---------------------------------------------- #
data_folder <- sprintf("01_phyloFlash/phyloflash_%s_out/%s.phyloFlash/", Sxxx, Sxxx)
d_pf <- read_csv(sprintf("%s/%s.phyloFlash.extractedSSUclassifications.csv", data_folder, Sxxx))
d_ntu_full <- read_csv(sprintf("%s/%s.phyloFlash.NTUfull_abundance.csv", data_folder, Sxxx),
                       col_names = FALSE)
d_taxa1 <- as.data.frame(t(sapply(str_split(d_ntu_full$X1, pattern = ";"), "[")))
colnames(d_taxa1) <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
d_otu1 <- data.frame(Sxxx = d_ntu_full$X2)
rownames(d_otu1) <- rownames(d_taxa1) <- sprintf("taxa_%03d", 1:nrow(d_ntu_full))

## Import to phyloseq
ps_pf0 <- phyloseq(otu_table(d_otu1, taxa_are_rows = TRUE),
                   tax_table(as.matrix(d_taxa1)))
ps_pf_r <- subset_taxa(ps_pf0, taxa_sums(ps_pf0) >= min_reads) %>% # Remove spurious seqs
  subset_taxa(order != "Chloroplast") %>%
  subset_taxa(superkingdom == "Bacteria" | superkingdom == "Archaes") %>%
  transform_sample_counts(function(x) x/sum(x))
ps_pf <- subset_taxa(ps_pf0, taxa_sums(ps_pf0) >= min_reads) %>% # Remove spurious seqs
  subset_taxa(order != "Chloroplast") %>%
  subset_taxa(superkingdom == "Bacteria" | superkingdom == "Archaes")


# ---------------------------------------------- #
# Load qMiSeq data
# ---------------------------------------------- #
ps_qm0 <- readRDS("data/CER2017_RiceField_all.obj") %>%
  subset_samples(Description == Sample_Description)
ps_qm_r <- subset_taxa(ps_qm0, taxa_sums(ps_qm0) > 0) %>%
  subset_taxa(superkingdom == "Bacteria" | superkingdom == "Archaes") %>%
  transform_sample_counts(function(x) x/sum(x))
ps_qm <- subset_taxa(ps_qm0, taxa_sums(ps_qm0) > 0) %>%
  subset_taxa(superkingdom == "Bacteria" | superkingdom == "Archaes")


# ---------------------------------------------- #
# Merge phyloFlash and qMiSeq
# ---------------------------------------------- #
## phyloFlash data
ps_pf_merge <- ps_pf
ps_pf_merge <- tax_glom(ps_pf_merge, taxrank = "phylum", NArm = TRUE)
tax_table(ps_pf_merge) <- as.data.frame(tax_table(ps_pf_merge)) %>%
  select(superkingdom, phylum) %>% as.matrix()

## qMiSeq data
ps_qm_merge <- ps_qm
tax_table(ps_qm_merge)[,"query"] <- NA
tax_table(ps_qm_merge)[,"phylum"][tax_table(ps_qm_merge)[,"phylum"] == ""] <- "Others"
ps_qm_merge <- tax_glom(ps_qm_merge, taxrank = "phylum", NArm = TRUE)
tax_table(ps_qm_merge) <- as.data.frame(tax_table(ps_qm_merge)) %>%
  select(superkingdom, phylum) %>% as.matrix()

# Re-extract data and mege data
## Sample data
df_s1 <- data.frame(Sample_ID = sprintf("%s_Shotgun", Sxxx),
                    date = Sample_Date, plot = Sample_Plot)
df_s2 <- data.frame(Sample_ID = sprintf("%s_qMiSeq", Sxxx),
                    date = Sample_Date, plot = Sample_Plot)

## OTU tables
df_o1 <- as.data.frame(otu_table(ps_pf_merge))
df_o2 <- as.data.frame(t(otu_table(ps_qm_merge)))
pf_o1_name <- as.character(tax_table(ps_pf_merge)[,"phylum"])
pf_o2_name <- as.character(tax_table(ps_qm_merge)[,"phylum"])
pf_o1_name[which("(Bacteria)" == pf_o1_name)] <- "Others"

## Taxa tables
df_t1 <- as.data.frame(tax_table(ps_pf_merge))
df_t1$phylum[which("(Bacteria)" == df_t1$phylum)] <- "Others"
df_t2 <- as.data.frame(tax_table(ps_qm_merge))
rownames(df_t1) <- rownames(df_o1) <- pf_o1_name
rownames(df_t2) <- rownames(df_o2) <- pf_o2_name
colnames(df_o1) <- rownames(df_s1) <- df_s1$Sample_ID
colnames(df_o2) <- rownames(df_s2) <- df_s2$Sample_ID

## Sort
df_t1 <- df_t1[sort(rownames(df_t1)),]
df_t2 <- df_t2[sort(rownames(df_t2)),]
df_o1 <- as.data.frame(df_o1[sort(rownames(df_o1)),])
rownames(df_o1) <- sort(pf_o1_name)
colnames(df_o1) <- df_s1$Sample_ID
df_o2 <- as.data.frame(df_o2[sort(rownames(df_o2)),])
rownames(df_o2) <- sort(pf_o2_name)
colnames(df_o2) <- df_s2$Sample_ID
pf_o1_name <- sort(pf_o1_name)
pf_o2_name <- sort(pf_o2_name)

### Rename taxa ###
rename_id1 <- match(pf_o1_name, rename_df$pf)
pf_tmp_name <- rename_df$qm[rename_id1]
pf_tmp_name[is.na(pf_tmp_name)] <- pf_o1_name[which(is.na(rename_id1))]
pf_o1_name <- pf_tmp_name
rownames(df_t1) <- rownames(df_o1) <- df_t1$phylum <- pf_o1_name
### Rename taxa ###

### Add zero taxa ###
no_taxa_id1 <- which(is.na(match(rename_df$qm, pf_o1_name)))
df_t1_tmp <- data.frame(superkingdom = rep("Bacteria", length(no_taxa_id1)),
                        phylum = rename_df$qm[no_taxa_id1])
df_o1_tmp <- data.frame(S000 = rep(0, length(no_taxa_id1)))
rownames(df_t1_tmp) <- df_t1_tmp$phylum
df_t1 <- df_t1 %>% add_row(df_t1_tmp)
rownames(df_t1)[(nrow(df_t1)-nrow(df_t1_tmp)+1):nrow(df_t1)] <- rownames(df_t1_tmp)
colnames(df_o1_tmp) <- colnames(df_o1)
df_o1 <- df_o1 %>% add_row(df_o1_tmp)
rownames(df_o1)[(nrow(df_o1)-nrow(df_o1_tmp)+1):nrow(df_o1)] <- rownames(df_t1_tmp)

no_taxa_id2 <- which(is.na(match(rename_df$qm, pf_o2_name)))
df_t2_tmp <- data.frame(superkingdom = rep("Bacteria", length(no_taxa_id2)),
                        phylum = rename_df$qm[no_taxa_id2])
df_o2_tmp <- data.frame(S000 = rep(0, length(no_taxa_id2)))
rownames(df_t2_tmp) <- df_t2_tmp$phylum
df_t2 <- df_t2 %>% add_row(df_t2_tmp)
rownames(df_t2)[(nrow(df_t2)-nrow(df_t2_tmp)+1):nrow(df_t2)] <- rownames(df_t2_tmp)
colnames(df_o2_tmp) <- colnames(df_o2)
df_o2 <- df_o2 %>% add_row(df_o2_tmp)
rownames(df_o2)[(nrow(df_o2)-nrow(df_o2_tmp)+1):nrow(df_o2)] <- rownames(df_t2_tmp)


## Re-import to phyloseq
ps1 <- phyloseq(otu_table(df_o1, taxa_are_rows = TRUE),
                sample_data(df_s1),
                tax_table(as.matrix(df_t1)))
ps2 <- phyloseq(otu_table(df_o2, taxa_are_rows = TRUE),
                sample_data(df_s2),
                tax_table(as.matrix(df_t2)))
ps3 <- merge_phyloseq(ps1, ps2); ps3
ps3 <- ps3 %>% transform_sample_counts(function(x) x/sum(x))

## Reset phylum levels
ps3_phylum_u <- sort(unique((tax_table(ps3)[,"phylum"])))
ps3_phylum_levels <- c(ps3_phylum_u[ps3_phylum_u != "Others"], "Others")

## Visualize
g1 <- plot_bar(ps3, fill = "phylum") +
  scale_fill_manual(values = palette_manual, breaks = ps3_phylum_levels) +
    labs(title = sprintf("Plot %s, %s", Sample_Plot, Sample_Date),
         subtitle = sprintf("qMiSeq = %s, Shotgun = %s",
                            ntaxa(ps_qm), ntaxa(ps_pf))) +
  xlab(NULL) + ylab("Relative abundance") + 
  theme(axis.text.x = element_text(vjust = 0.5)) +
  NULL

g1$data$phylum <- factor(g1$data$phylum, levels = ps3_phylum_levels)
g1


# ---------------------------------------------- #
# Save output
# ---------------------------------------------- #
saveRDS(g1, file = sprintf("15_VisualizeAllOut/g1_%s.obj", Sxxx))
ggsave(file = sprintf("%s/g1_%s.pdf", output_folder, Sxxx),
       plot = g1, width = 6, height = 8)

