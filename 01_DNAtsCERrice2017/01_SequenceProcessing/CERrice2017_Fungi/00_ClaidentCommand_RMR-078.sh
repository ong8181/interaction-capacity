####
#### Command to demultiplex bcl to fastq
#### (For ITS sequences)
####

# Set parameters
BASECALL_FOLDER="180302_M04195_0084_000000000-BP7V4"
RUN_FOLDER="180302_M04195_0084_000000000-BP7V4"
FASTQ_OUT_FOLDER="fastq_out"
FASTQ_OUT_FOLDER2="RMR-078_v2_fastq_out"
RUNNAME="RMR-078"
F_PRIMER="../TagPrimerFiles/RMR-078_tagfile/RMR-078_F_primer.txt"
R_PRIMER="../TagPrimerFiles/RMR-078_tagfile/RMR-078_R_primer.txt"
I7_INDEX="../TagPrimerFiles/RMR-078_tagfile/RMR-078_i7_index.txt"
I5_INDEX="../TagPrimerFiles/RMR-078_tagfile/RMR-078_i5_index.txt"
DEMULTIPLEX_OUT="RMR-078_v2_ClaidentDemultiplexed"

# Convert Bcl to Fastq (bcl2fastq2 v2.18)
bcl2fastq --processing-threads 72 --use-bases-mask Y300n,I8,I8,Y300n --create-fastq-for-index-reads --runfolder-dir $RUN_FOLDER --output-dir $FASTQ_OUT_FOLDER2

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
clmakecachedb --blastdb=overall_class --numthreads=72 FungiASV_seqs.fa FungiASV_overall_cache
clidentseq --blastdb=FungiASV_overall_cache --numthreads=72 FungiASV_seqs.fa FungiASV_overall_clidentseq
classigntax --taxdb=overall_class --maxpopposer=0.05 --minsoratio=19 FungiASV_overall_clidentseq FungiASV_overall_classigntax

# Overall genus
clmakecachedb --blastdb=overall_genus --numthreads=72 FungiASV_seqs.fa FungiASV_overallg_cache
clidentseq --blastdb=FungiASV_overallg_cache --numthreads=72 FungiASV_seqs.fa FungiASV_overallg_clidentseq
classigntax --taxdb=overall_genus --maxpopposer=0.05 --minsoratio=19 FungiASV_overallg_clidentseq FungiASV_overallg_classigntax

# Merge identification results (overall_class + overall_genus)
clmergeassign --priority=descend FungiASV_overallg_classigntax FungiASV_overall_classigntax FungiASV_merge_classigntax
