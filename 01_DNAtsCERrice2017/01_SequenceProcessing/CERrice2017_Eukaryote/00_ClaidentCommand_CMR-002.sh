####
#### Commands to demultiplex bcl to fastq
#### (For 18S rRNA sequences)
####

# Set parameters
BASECALL_FOLDER="181130_M00962_0005_000000000-C5HNN"
RUN_FOLDER="181130_M00962_0005_000000000-C5HNN"
FASTQ_OUT_FOLDER="fastq_out"
FASTQ_OUT_FOLDER2="CMR-002_v2_fastq_out"
RUNNAME="CMR-002"
F_PRIMER="../TagPrimerFiles/CMR-002_tagfile/CMR-002_F_primer.txt"
R_PRIMER="../TagPrimerFiles/CMR-002_tagfile/CMR-002_R_primer.txt"
I7_INDEX="../TagPrimerFiles/CMR-002_tagfile/CMR-002_i7_index.txt"
I5_INDEX="../TagPrimerFiles/CMR-002_tagfile/CMR-002_i5_index.txt"
DEMULTIPLEX_OUT="CMR-002_v2_ClaidentDemultiplexed"

# Convert Bcl to Fastq (bcl2fastq2 v2.18)
# !!! Rename "SampleSheet.csv" before this command (e.g., "SampleSheet_rename.csv") !!!
bcl2fastq --processing-threads 72 --use-bases-mask Y150n,I8,I8,Y150n --create-fastq-for-index-reads --runfolder-dir $RUN_FOLDER --output-dir $FASTQ_OUT_FOLDER2

# Demultiplexing
cd $FASTQ_OUT_FOLDER2
clsplitseq --runname=$RUNNAME --index1file=$I7_INDEX --index2file=$I5_INDEX --primerfile=$F_PRIMER --reverseprimerfile=$R_PRIMER --minqualtag=30 --numthreads=72 --truncateN=enable *_R1_001.fastq.gz *_I1_001.fastq.gz *_I2_001.fastq.gz *_R2_001.fastq.gz $DEMULTIPLEX_OUT

# Move undetermined fastq to another folder
cd $DEMULTIPLEX_OUT
mkdir Undetermined
mv *undetermined.*.fastq.gz Undetermined

# ---------- These processes should be done after picking ASV by DADA2 ---------- #
# Taxa assignment using claident
# Check overall_class
clmakecachedb --blastdb=overall_class --numthreads=72 EukASV_seqs.fa EukASV_overall_cache
clidentseq --blastdb=EukASV_overall_cache --numthreads=72 EukASV_seqs.fa EukASV_overall_clidentseq
classigntax --taxdb=overall_class --maxpopposer=0.05 --minsoratio=19 EukASV_overall_clidentseq EukASV_overall_classigntax

# Overall genus
clmakecachedb --blastdb=overall_genus --numthreads=72 EukASV_seqs.fa EukASV_overallg_cache
clidentseq --blastdb=EukASV_overallg_cache --numthreads=72 EukASV_seqs.fa EukASV_overallg_clidentseq
classigntax --taxdb=overall_genus --maxpopposer=0.05 --minsoratio=19 EukASV_overallg_clidentseq EukASV_overallg_classigntax

# Merge identification results (overall_class + overall_genus)
clmergeassign --priority=descend EukASV_overallg_classigntax EukASV_overall_classigntax EukASV_merge_classigntax


