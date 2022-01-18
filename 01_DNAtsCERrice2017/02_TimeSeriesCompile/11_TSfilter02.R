####
#### CER eDNA study
#### No.11 Filter time series using correlation between variables
####

# Load workspace
load("10_TSfilterCombineOut/10_TSfilterCombineOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder11 <- "11_TSfilter02Out"
dir.create(output_folder11)

# Load library
library(phyloseq); packageVersion("phyloseq") #1.28.0, 2020.1.6

# Calculate correlation coefficients
ps_cor <- cor(as.data.frame(otu_table(ps_combined)@.Data))
diag(ps_cor) <- 0
ps_cor[upper.tri(ps_cor)] <- 0

# Find highly correlated ASV columns
cor_set <- data.frame(NULL)
for(i in 1:nrow(ps_cor)){
  cor_set_j <- data.frame(NULL)
  j0 <- which(ps_cor[,i] > 0.90)
  if(length(j0) > 0 & length(j0) < 2){
    cor_set_j <- data.frame(taxa1 = i, taxa2 = j0, cor = ps_cor[j0,i])
  }else if(length(j0) > 1){
    for(j1 in 1:length(j0)){
      cor_set_j <- rbind(cor_set_j, data.frame(taxa1 = i, taxa2 = j0[j1], cor = ps_cor[j0[j1],i]))
    }
  }
  cor_set <- rbind(cor_set, cor_set_j)
}

# Add taxa names
cor_set_taxa1 <- tax_table(ps_combined)[cor_set$taxa1, c("superkingdom", "phylum", "family", "genus", "miseq_run")]@.Data
cor_set_taxa2 <- tax_table(ps_combined)[cor_set$taxa2, c("superkingdom", "phylum", "family", "genus", "miseq_run")]@.Data
cor_set$taxa_name1 <- sprintf("%s_%s_%s_%s_%s", cor_set_taxa1[,5], cor_set_taxa1[,1], cor_set_taxa1[,2], cor_set_taxa1[,3], cor_set_taxa1[,4])
cor_set$taxa_name2 <- sprintf("%s_%s_%s_%s_%s", cor_set_taxa2[,5], cor_set_taxa2[,1], cor_set_taxa2[,2], cor_set_taxa2[,3], cor_set_taxa2[,4])
cor_set$taxa_id1 <- taxa_names(ps_combined)[cor_set$taxa1]
cor_set$taxa_id2 <- taxa_names(ps_combined)[cor_set$taxa2]
rownames(cor_set) <- 1:nrow(cor_set)

# Compare taxa information
cor_set$assign_same_taxa <- NaN
# Check NA
for(cor_tax_i in 1:nrow(cor_set)){
  comp_id <- intersect(which(cor_set_taxa1[cor_tax_i,1:4] != ""), which(cor_set_taxa2[cor_tax_i,1:4] != ""))
  cor_set$assign_same_taxa[cor_tax_i] <- all(cor_set_taxa1[cor_tax_i,comp_id] == cor_set_taxa2[cor_tax_i,comp_id])
}

# Output CSV file
write.csv(cor_set, sprintf("%s/cor_set.csv", output_folder11), row.names = F)

# Check correlation figures
pdf(sprintf("%s/cor_all.pdf", output_folder11), width = 50, height = 50)
op <- par(mfrow=c(20,20))
for(i in 1:nrow(cor_set)){
  plot(otu_table(ps_combined)[,cor_set[i,1]], otu_table(ps_combined)[,cor_set[i,2]],
       xlab = cor_set[i,4], ylab = cor_set[i,5], main = cor_set[i,3])
  abline(0,1)
}
par(op)
dev.off()

# Identify taxa names with more than r = 0.90
sum(cor_set$cor > 0.9)
cor_set_filter <- cor_set[cor_set$cor > 0.90 & cor_set$assign_same_taxa == 1,]
cor_set_filter$taxa1_sum <- taxa_sums(ps_combined)[cor_set_filter$taxa1]
cor_set_filter$taxa2_sum <- taxa_sums(ps_combined)[cor_set_filter$taxa2]
cor_set_filter$keep <- NaN
cor_set_filter$remove <- NaN
keep_taxa1 <- cor_set_filter$taxa1_sum - cor_set_filter$taxa2_sum >= 0
keep_taxa2 <- cor_set_filter$taxa1_sum - cor_set_filter$taxa2_sum < 0
cor_set_filter$keep[keep_taxa1] <- taxa_names(ps_combined)[cor_set_filter$taxa1[keep_taxa1]]
cor_set_filter$keep[keep_taxa2] <- taxa_names(ps_combined)[cor_set_filter$taxa2[keep_taxa2]]
cor_set_filter$remove[!keep_taxa1] <- taxa_names(ps_combined)[cor_set_filter$taxa1[!keep_taxa1]]
cor_set_filter$remove[!keep_taxa2] <- taxa_names(ps_combined)[cor_set_filter$taxa2[!keep_taxa2]]
remove_cor_taxa <- unique(cor_set_filter$remove)

filtered_taxa <- taxa_names(ps_combined)[!(taxa_names(ps_combined) %in% remove_cor_taxa)]
ps_filt <- prune_taxa(filtered_taxa, ps_combined)
ps_filt

# Save phyloseq objects
saveRDS(ps_filt, sprintf("%s/ps_comb_filt.obj", output_folder11))

# Save and output results
#save(list = ls(all.names = TRUE),
#     file = sprintf("%s/11_TSfilter02Out.RData", output_folder11))
save.image(sprintf("%s/11_TSfilter02Out.RData", output_folder11))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/11_SessionInfo_TSfilter02_%s.txt", substr(Sys.time(), 1, 10)))

# Check No. of ASV
subset_taxa(ps_filt, miseq_run == "RMR-076")
subset_taxa(ps_filt, miseq_run == "RMR-078")
subset_taxa(ps_filt, miseq_run == "RMR-099")
subset_taxa(ps_filt, miseq_run == "CMR-002")
