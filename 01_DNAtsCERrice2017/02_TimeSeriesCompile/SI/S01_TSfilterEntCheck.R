####
#### CER eDNA study
#### No.S01 Entropy estimation
####

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folderS01 <- "S01_TSfilterEntCheckOut"
dir.create(output_folderS01)

# Load library and functions
library(phyloseq); packageVersion("phyloseq") #1.20.0, 2018.1.16
library(infotheo); packageVersion("infotheo") #1.2.0

# Load workspace
load("../07_PhyloseqImportOut/07_PhyloseqImportOut.RData")

# Extracting data as time series
ps_prok_sample <- prune_taxa(taxa_sums(ps_prok_sample) > 0, ps_prok_sample)
ps_fungi_sample <- prune_taxa(taxa_sums(ps_fungi_sample) > 0, ps_fungi_sample)
ps_euk_sample <- prune_taxa(taxa_sums(ps_euk_sample) > 0, ps_euk_sample)
ps_inv_sample <- prune_taxa(taxa_sums(ps_inv_sample) > 0, ps_inv_sample)

# Function to test the dependency of entropy estimation on n_bin
ent_bin_test <- function(ps_object){
  ent10 <- c(NULL); ent20 <- c(NULL); ent30 <- c(NULL)
  ent50 <- c(NULL); ent100 <- c(NULL); ent200 <- c(NULL)
  for(i in 1:length(otu_table(ps_object)[,1])){
    xx <- c(otu_table(ps_object)[,i])
    i10 <- infotheo::entropy(discretize(xx, disc = "equalwidth", nbins = 10))
    i20 <- infotheo::entropy(discretize(xx, disc = "equalwidth", nbins = 20))
    i30 <- infotheo::entropy(discretize(xx, disc = "equalwidth", nbins = 30))
    i50 <- infotheo::entropy(discretize(xx, disc = "equalwidth", nbins = 50))
    i100 <- infotheo::entropy(discretize(xx, disc = "equalwidth", nbins = 100))
    i200 <- infotheo::entropy(discretize(xx, disc = "equalwidth", nbins = 200))
    
    ent10 <- c(ent10, i10)
    ent20 <- c(ent20, i20)
    ent30 <- c(ent30, i30)
    ent50 <- c(ent50, i50)
    ent100 <- c(ent100, i100)
    ent200 <- c(ent200, i200)
  }
  
  ent_df <- data.frame(ent_bin_010 = ent10,
                       ent_bin_020 = ent20,
                       ent_bin_030 = ent30,
                       ent_bin_050 = ent50,
                       ent_bin_100 = ent100,
                       ent_bin_200 = ent200)
  return(ent_df)
}

prok_ent_test <- ent_bin_test(ps_prok_sample)
fungi_ent_test <- ent_bin_test(ps_fungi_sample)
inv_ent_test <- ent_bin_test(ps_inv_sample)
euk_ent_test <- ent_bin_test(ps_euk_sample)

# Save prok entropy
pdf(file = sprintf("%s/Prok_1_Scatter.pdf", output_folderS01), width = 8, height = 8)
plot(prok_ent_test)
dev.off()

pdf(file = sprintf("%s/Prok_2_RankScatter.pdf", output_folderS01), width = 8, height = 8)
prok_rank <- as.data.frame(apply(prok_ent_test, 2, rank))
plot(prok_rank)
dev.off()

pdf(file = sprintf("%s/Prok_3_BinDependence.pdf", output_folderS01), width = 9, height = 6)
op <- par(mfrow = c(1,2))
n_bin <- c(10, 20, 30, 50, 100, 200)
plot(as.numeric((prok_ent_test/prok_ent_test$ent_bin_200)[1,]) ~ n_bin,
     type = "l", ylim=c(0,1.2),
     ylab = "Scaled shannon entropy",
     xlab = "No. of binning")
for(i in 2:nrow(prok_ent_test)) lines(as.numeric((prok_ent_test/prok_ent_test$ent_bin_200)[i,]) ~ n_bin, lwd=0.1)
plot(as.numeric((prok_ent_test)[1,]) ~ n_bin, type = "l", ylim=c(0,3),
     ylab = "Shannon entropy",
     xlab = "No. of binning")
for(i in 2:nrow(prok_ent_test)) lines(as.numeric(prok_ent_test[i,]) ~ n_bin, lwd=0.1)
par(op)
dev.off()


# Save fungi entropy
pdf(file = sprintf("%s/Fungi_1_Scatter.pdf", output_folderS01), width = 8, height = 8)
plot(fungi_ent_test)
dev.off()

pdf(file = sprintf("%s/Fungi_2_RankScatter.pdf", output_folderS01), width = 8, height = 8)
fungi_rank <- as.data.frame(apply(fungi_ent_test, 2, rank))
plot(fungi_rank)
dev.off()

pdf(file = sprintf("%s/Fungi_3_BinDependence.pdf", output_folderS01), width = 9, height = 6)
op <- par(mfrow = c(1,2))
n_bin <- c(10, 20, 30, 50, 100, 200)
plot(as.numeric((fungi_ent_test/fungi_ent_test$ent_bin_200)[1,]) ~ n_bin,
     type = "l", ylim=c(0,1.2),
     ylab = "Scaled shannon entropy",
     xlab = "No. of binning")
for(i in 2:nrow(fungi_ent_test)) lines(as.numeric((fungi_ent_test/fungi_ent_test$ent_bin_200)[i,]) ~ n_bin, lwd=0.1)
plot(as.numeric((fungi_ent_test)[1,]) ~ n_bin, type = "l", ylim=c(0,3),
     ylab = "Shannon entropy",
     xlab = "No. of binning")
for(i in 2:nrow(fungi_ent_test)) lines(as.numeric(fungi_ent_test[i,]) ~ n_bin, lwd=0.1)
par(op)
dev.off()


# Save invertebrate entropy
pdf(file = sprintf("%s/Inv_1_Scatter.pdf", output_folderS01), width = 8, height = 8)
plot(inv_ent_test)
dev.off()

pdf(file = sprintf("%s/Inv_2_RankScatter.pdf", output_folderS01), width = 8, height = 8)
inv_rank <- as.data.frame(apply(inv_ent_test, 2, rank))
plot(inv_rank)
dev.off()

pdf(file = sprintf("%s/Inv_3_BinDependence.pdf", output_folderS01), width = 9, height = 6)
op <- par(mfrow = c(1,2))
n_bin <- c(10, 20, 30, 50, 100, 200)
plot(as.numeric((inv_ent_test/inv_ent_test$ent_bin_200)[1,]) ~ n_bin,
     type = "l", ylim=c(0,1.2),
     ylab = "Scaled shannon entropy",
     xlab = "No. of binning")
for(i in 2:nrow(inv_ent_test)) lines(as.numeric((inv_ent_test/inv_ent_test$ent_bin_200)[i,]) ~ n_bin, lwd=0.1)
plot(as.numeric((inv_ent_test)[1,]) ~ n_bin, type = "l", ylim=c(0,3),
     ylab = "Shannon entropy",
     xlab = "No. of binning")
for(i in 2:nrow(inv_ent_test)) lines(as.numeric(inv_ent_test[i,]) ~ n_bin, lwd=0.1)
par(op)
dev.off()


# Save eukaryote entropy
pdf(file = sprintf("%s/Euk_1_Scatter.pdf", output_folderS01), width = 8, height = 8)
plot(euk_ent_test)
dev.off()

pdf(file = sprintf("%s/Euk_2_RankScatter.pdf", output_folderS01), width = 8, height = 8)
euk_rank <- as.data.frame(apply(euk_ent_test, 2, rank))
plot(euk_rank)
dev.off()

pdf(file = sprintf("%s/Euk_3_BinDependence.pdf", output_folderS01), width = 9, height = 6)
op <- par(mfrow = c(1,2))
n_bin <- c(10, 20, 30, 50, 100, 200)
plot(as.numeric((euk_ent_test/euk_ent_test$ent_bin_200)[1,]) ~ n_bin,
     type = "l", ylim=c(0,1.2),
     ylab = "Scaled shannon entropy",
     xlab = "No. of binning")
for(i in 2:nrow(euk_ent_test)) lines(as.numeric((euk_ent_test/euk_ent_test$ent_bin_200)[i,]) ~ n_bin, lwd=0.1)
plot(as.numeric((euk_ent_test)[1,]) ~ n_bin, type = "l", ylim=c(0,3),
     ylab = "Shannon entropy",
     xlab = "No. of binning")
for(i in 2:nrow(euk_ent_test)) lines(as.numeric(euk_ent_test[i,]) ~ n_bin, lwd=0.1)
par(op)
dev.off()


# Save and output results
save.image(sprintf("%s/S01_TSfilterEntCheckOut.RData", output_folderS01))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("%s/S01_SessionInfo_TSfilter_%s.txt", output_folderS01, substr(Sys.time(), 1, 10)))
