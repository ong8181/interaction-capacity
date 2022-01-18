####
#### Standard STD identification
####
####

cd "/XXXXXX/01_SequenceProcessing/CERrice2017_Fungi/"

DBPATH=/XXXXXX/DADA2_DB/STDseqs/FungiSTD/FungiSTD_ITSKYO
QUERYPATH=01_CERrice2017_FungiSeqOut/FungiASV_seqs.fa
OUTPUT=02_ident_FungiSTD_BLASTnOut/FungiSTD_out.txt
EVALUE_SET=1e-118

mkdir 02_ident_FungiSTD_BLASTnOut
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT}


