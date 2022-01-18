####
#### Rice field samples (from 2017.5.23-9.22, 5 plots)
#### No.1 CMR-002 Eukaryote sequence analysis by DADA2
####

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)
dir.create("00_SessionInfo")
dir.create("01_CERrice2017_EukSeqOut")

# Load library and functions
library(dada2); packageVersion("dada2") #1.11.5, 2019.10.17
library(ShortRead); packageVersion("ShortRead") #1.40.0, 2019.10.17
library(ggplot2); packageVersion("ggplot2") #3.1.1, 2019.10.17
source("../functions/F01_HelperFunctions.R")

# Load sequence reads
# CMR-002-Eukaryote
path <- "/XXXXXX/01_SequenceProcessing/seqdata_claident/CMR-002-Eukaryote"
fnFs <- sort(list.files(path, pattern=".forward.fastq", full.names = T)) # Forward read files
fnRs <- sort(list.files(path, pattern=".reverse.fastq", full.names = T)) # Reverse read files

# Identify primers
FWD <- "GTACACACCGCCCGTC" # Euk1391f
REV <- "TGATCCTTCTGCAGGTTCACCTAC" # EukBr
FWD_orients <- AllOrients(FWD)
REV_orients <- AllOrients(REV)

# Pre-filtering to remove Ns
fnFs_filtN <- file.path(path, "01_filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs_filtN <- file.path(path, "01_filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs_filtN, fnRs, fnRs_filtN, maxN = 0, multithread = TRUE)

# Identify primers
rbind(FWD.ForwardReads = sapply(FWD_orients, PrimerHits, fn = fnFs_filtN[[4]]),
      FWD.ReverseReads = sapply(FWD_orients, PrimerHits, fn = fnRs_filtN[[4]]), 
      REV.ForwardReads = sapply(REV_orients, PrimerHits, fn = fnFs_filtN[[4]]), 
      REV.ReverseReads = sapply(REV_orients, PrimerHits, fn = fnRs_filtN[[4]]))

# Tell the path to cutadapt command
cutadapt <- "XXXXXX/cutadapt"
system2(cutadapt, args = "--version")

# Remove primers
path_cut <- file.path(path, "02_cutadapt")
if(!dir.exists(path_cut)) dir.create(path_cut)
fnFs_cut <- file.path(path_cut, basename(fnFs))
fnRs_cut <- file.path(path_cut, basename(fnRs))

FWD_RC <- dada2:::rc(FWD)
REV_RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1_flags <- paste("-g", FWD, "-a", REV_RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2_flags <- paste("-G", REV, "-A", FWD_RC)

# Run Cutadapt
# Cutadapt commands may crash from R console. In that case, excute it from terminal.
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-j 72", # Multithred option
                             R1_flags, R2_flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs_cut[i], "-p", fnRs_cut[i], # output files
                             fnFs_filtN[i], fnRs_filtN[i])) # input files
}

# Sanity check
rbind(FWD.ForwardReads = sapply(FWD_orients, PrimerHits, fn = fnFs_cut[[4]]), 
      FWD.ReverseReads = sapply(FWD_orients, PrimerHits, fn = fnRs_cut[[4]]), 
      REV.ForwardReads = sapply(REV_orients, PrimerHits, fn = fnFs_cut[[4]]), 
      REV.ReverseReads = sapply(REV_orients, PrimerHits, fn = fnRs_cut[[4]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path_cut, pattern = ".forward.fastq", full.names = TRUE))
cutRs <- sort(list.files(path_cut, pattern = ".reverse.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
sample_names <- unname(sapply(cutFs, function(x) strsplit(basename(x), "_")[[1]][2]))
head(sample_names)
# Visualize quality
#plotQualityProfile(fnFs[101:103])
#plotQualityProfile(fnRs[101:103])

# Perform quality filtering
filtFs <- file.path(path, "03_filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "03_filtered", paste0(sample_names, "_R_filt.fastq.gz"))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     # Remove Ns and prokaryote primers; Euk1391f = 18 bp, EukBr = 24 bp
                     # Output sequences of Claident already trimmed Ns and primers
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=F, minLen = 20,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Exclude 0 seq samples, rename filtFs and filtRs
if(length(sample_names[out[,2]<1 | out[,1]<1]) > 0){
  filtFs <- file.path(path, "03_filtered", paste0(sample_names[out[,2]>0 & out[,1]>0], "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "03_filtered", paste0(sample_names[out[,2]>0 & out[,1]>0], "_R_filt.fastq.gz"))
}

# Learn the error rates
min_nbases <- 5e+09
errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)

# Visualize errors
#plotErrors(errF, nominalQ = T)
#plotErrors(errR, nominalQ = T)

# Dereplicatin
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample_names[out[,2] > 0 & out[,1] > 0]
names(derepRs) <- sample_names[out[,2] > 0 & out[,1] > 0]

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merging paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, trimOverhang = TRUE, minOverlap = 20, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Cutting unexpected length sequences
seqtab2 <- seqtab # No filtering at this point
table(nchar(getSequences(seqtab2)))

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab2, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab_nochim)
sum(seqtab_nochim)/sum(seqtab)

# Track reads thourhg the pipeline
out2 <- out[out[,2] > 0 & out[,1] > 0,]
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab2), rowSums(seqtab_nochim),  rowSums(seqtab_nochim)/out2[,1])
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "tabled2", "nonchim", "prop(last/first)")
rownames(track) <- sample_names[out[,2] > 0 & out[,1] > 0]
head(track)

# Taxa output for claident tax assginment
seqs <- colnames(seqtab_nochim)
seqs_out <- as.matrix(c(rbind(sprintf(">Taxa%05d", 1:length(seqs)), seqs)), ncol = 1)
write.table(seqs_out, "01_CERrice2017_EukSeqOut/EukASV_seqs.fa", col.names = FALSE, row.names = FALSE, quote = FALSE)

#save.image("01_CERrice2017_EukSeqOut/01_CERrice2017_EukSeq_AllOut.RData")

# Reduce size
rm(dadaFs)
rm(dadaRs)
rm(derepFs)
rm(derepRs)

#save(list = ls(all.names = TRUE),
#     file = "01_CERrice2017_EukSeqOut/01_CERrice2017_EukSeqOut.RData")
save.image("01_CERrice2017_EukSeqOut/01_CERrice2017_EukSeqOut.RData")

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/01_SessionInfo_EukSeq_%s.txt", substr(Sys.time(), 1, 10)))
