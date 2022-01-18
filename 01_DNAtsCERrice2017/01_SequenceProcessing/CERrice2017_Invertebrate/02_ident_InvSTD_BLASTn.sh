####
#### Standard STD identification
####
####

cd "XXXXXX/01_SequenceProcessing/CERrice2017_Invertebrate/"

DBPATH=/XXXXXX/DADA2_DB/STDseqs/InvertebrateSTD/InvertebrateSTD_mlCOI
QUERYPATH=01_CERrice2017_InvSeqOut/InvASV_seqs.fa
OUTPUT=02_ident_InvSTD_BLASTnOut/InvSTD_out.txt
EVALUE_SET=1e-167

mkdir 02_ident_InvSTD_BLASTnOut
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT}

