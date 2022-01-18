####
#### Standard STD identification
####
####

cd "XXXXXX/01_SequenceProcessing/CERrice2017_Eukaryote/"

DBPATH=/XXXXXX/DADA2_DB/STDseqs/EukaryoteSTD/Euk1391f_STD
QUERYPATH=01_CERrice2017_EukSeqOut/EukASV_seqs.fa
OUTPUT=02_ident_EukSTD_BLASTnOut/EukSTD_out.txt
EVALUE_SET=1e-50

mkdir 02_ident_EukSTD_BLASTnOut
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT}

