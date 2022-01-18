####
#### Analyzing Shotgun Metagenome
#### Quality control by fastp (https://github.com/OpenGene/fastp)
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
# Global quality filtering using fastp
# ------------------------------------------ #
## Preparation
SEQ_DATA="PhyloseqOut/seqdata" # "seqdata" should include the shotgun metagenome data
cd ~/Desktop/${SEQ_DATA}
OUTPUT_DIR="../seqdata_flt"
mkdir $OUTPUT_DIR

## Average quality > 30
## Length > 150 bp
## No Nextera adaptor
for file in *_R1_001.fastq.gz; do
fastp \
--average_qual=30 \
--adapter_sequence=CTGTCTCTTATACACATCT \
--adapter_sequence_r2=CTGTCTCTTATACACATCT \
--in1 ${file%_R1_001.fastq.gz}_R1_001.fastq.gz \
--in2 ${file%_R1_001.fastq.gz}_R2_001.fastq.gz \
--out1 ${OUTPUT_DIR}/${file%_R1_001.fastq.gz}_R1_flt.fastq.gz \
--out2 ${OUTPUT_DIR}/${file%_R1_001.fastq.gz}_R2_flt.fastq.gz \
-j ${OUTPUT_DIR}/report_fastp.json \
-h ${OUTPUT_DIR}/report_${file%_R1_001.fastq.gz}_fastp.html
done

# Delete fastp summary
#rm *fastp.html
#rm *fastp.json

