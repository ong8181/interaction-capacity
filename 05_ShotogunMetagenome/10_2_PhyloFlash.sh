####
#### Analyzing Shotgun Metagenome by PhyloFlash
#### (http://hrgv.github.io/phyloFlash/)
####

# ------------------------------------------ #
# Preparation
# ------------------------------------------ #
## Activate virtual environment
conda activate pf
## Generate alias
alias pf34='perl ~/phyloFlash-pf3.4/phyloFlash.pl'
## Check environment
pf34 -check_env


# ------------------------------------------ #
# Set working directory and analyze data
# Simple phyloFlash analysis
# ------------------------------------------ #
## Set working directory
OUT_DIR="PhyloseqOut"
cd ~/Desktop/${OUT_DIR}
mkdir phyloflash_S001_out
mkdir phyloflash_S002_out
mkdir phyloflash_S003_out
mkdir phyloflash_S004_out

## Sample S001
## Run SPAdes (skip EMIRGE) and produce all optional outputs (recommended)
cd phyloflash_S001_out
pf34 -lib S001 -almosteverything -read1 ../seqdata_flt/S001_S1_L001_R1_flt.fastq.gz -read2 ../seqdata_flt/S001_S1_L001_R2_flt.fastq.gz

## Sample S002
## Run SPAdes (skip EMIRGE) and produce all optional outputs (recommended)
cd ../phyloflash_S002_out
pf34 -lib S002 -almosteverything -read1 ../seqdata_flt/S002_S2_L001_R1_flt.fastq.gz -read2 ../seqdata_flt/S002_S2_L001_R2_flt.fastq.gz

## Sample S003
## Run SPAdes (skip EMIRGE) and produce all optional outputs (recommended)
cd ../phyloflash_S003_out
pf34 -lib S003 -almosteverything -read1 ../seqdata_flt/S003_S3_L001_R1_flt.fastq.gz -read2 ../seqdata_flt/S003_S3_L001_R2_flt.fastq.gz

## Sample S004
## Run SPAdes (skip EMIRGE) and produce all optional outputs (recommended)
cd ../phyloflash_S004_out
pf34 -lib S004 -almosteverything -read1 ../seqdata_flt/S004_S4_L001_R1_flt.fastq.gz -read2 ../seqdata_flt/S004_S4_L001_R2_flt.fastq.gz


