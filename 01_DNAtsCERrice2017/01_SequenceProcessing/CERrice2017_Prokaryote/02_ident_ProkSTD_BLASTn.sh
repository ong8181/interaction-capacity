####
#### Standard STD identification
####
####

cd "XXXXXX/01_SequenceProcessing/CERrice2017_Prokaryote/"

DBPATH=/XXXXXX/DADA2_DB/STDseqs/ProkaryoteSTD/ProkSTD_515F
QUERYPATH=01_CERrice2017_ProkSeqOut/ProkASV_seqs.fa
OUTPUT=02_ident_ProkSTD_BLASTnOut/ProkSTD_out.txt
EVALUE_SET=1e-100

mkdir 02_ident_ProkSTD_BLASTnOut
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT}

