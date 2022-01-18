####
#### CERrice2017 All data analysis
#### Generating sample meta-data tables
####

# Ceate output director
# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.29
library(lubridate); packageVersion("lubridate") # 1.7.4, 2019.11.3

# Create output directory
fig_output <- "../00_RawFigs/01_Fig_eDNAts"
dir.create(fig_output)

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
# Generate tables for Prokaryote 16S rRNA
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Prokaryote/06_NCcheck_ProkOut/06_NCcheck_ProkOut.RData")
all(rownames(sample_sheet) == rownames(track))
colnames(sample_sheet)
sample_sheet$date <- mdy(sample_sheet$date)
write_csv(sample_sheet, "../00_RawFigs/01_Fig_eDNAts/SI_Table_1_ProkMetaData.csv")
rm(list = ls())

# Generate tables for Eukaryote 18S rRNA
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Eukaryote/06_NCcheck_EukOut/06_NCcheck_EukOut.RData")
all(rownames(sample_sheet) == rownames(track))
colnames(sample_sheet)
sample_sheet$date <- mdy(sample_sheet$date)
write_csv(sample_sheet, "../00_RawFigs/01_Fig_eDNAts/SI_Table_2_EukMetaData.csv")
rm(list = ls())

# Generate tables for Fungal ITS
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Fungi/06_NCcheck_FungiOut/06_NCcheck_FungiOut.RData")
dim(sample_sheet); dim(track)
all(rownames(sample_sheet) == rownames(track))
colnames(sample_sheet)
sample_sheet$date <- mdy(sample_sheet$date)
write_csv(sample_sheet, "../00_RawFigs/01_Fig_eDNAts/SI_Table_3_FungiMetaData.csv")
rm(list = ls())

# Generate tables for Animal COI
load("../../01_DNAtsCERrice2017/01_SequenceProcessing/CERrice2017_Invertebrate/06_NCcheck_InvOut/06_NCcheck_InvOut.RData")
all(rownames(sample_sheet) == rownames(track))
colnames(sample_sheet)
sample_sheet$date <- mdy(sample_sheet$date)
write_csv(sample_sheet, "../00_RawFigs/01_Fig_eDNAts/SI_Table_4_InvMetaData.csv")
rm(list = ls())


## Generate summary table for all
prok_table <- read_csv("../00_RawFigs/01_Fig_eDNAts/SI_Table_1_ProkMetaData.csv")
euk_table <- read_csv("../00_RawFigs/01_Fig_eDNAts/SI_Table_2_EukMetaData.csv")
fungi_table <- read_csv("../00_RawFigs/01_Fig_eDNAts/SI_Table_3_FungiMetaData.csv")
inv_table <- read_csv("../00_RawFigs/01_Fig_eDNAts/SI_Table_4_InvMetaData.csv")

dim(prok_table); dim(euk_table); dim(fungi_table); dim(inv_table)

summary_table <- data.frame(Sample_Code = prok_table$Sample_Name2,
                            Description = prok_table$Description,
                            Date = prok_table$date,
                            Plot = prok_table$plot,
                            Sample_Category = prok_table$sample_nc,
                            Filter045_ml = prok_table$filt045_ml,
                            Filter022_ml = prok_table$filt022_ml,
                            DNA_TotalReads_STD = prok_table$STD_all +
                              euk_table$STD_all +
                              fungi_table$STD_all +
                              inv_table$STD_all,
                            DNA_TotalReads_NonSTD = prok_table$NonSTD_all +
                              euk_table$NonSTD_all +
                              fungi_table$NonSTD_all +
                              inv_table$NonSTD_all,
                            DNA_CopySum = prok_table$dna_copy_sum +
                              euk_table$dna_copy_sum +
                              fungi_table$dna_copy_sum +
                              inv_table$dna_copy_sum,
                            Standard_curve_validity_16S = prok_table$std_validity,
                            Standard_curve_validity_18S = euk_table$std_validity,
                            Standard_curve_validity_ITS = fungi_table$std_validity,
                            Standard_curve_validity_COI = inv_table$std_validity,
                            DADA2_Input_16S = prok_table$dada2_input,
                            DADA2_Input_18S = euk_table$dada2_input,
                            DADA2_Input_ITS = fungi_table$dada2_input,
                            DADA2_Input_COI = inv_table$dada2_input,
                            DADA2_OutputProp_16S = prok_table$dada2_prop,
                            DADA2_OutputProp_18S = euk_table$dada2_prop,
                            DADA2_OutputProp_ITS = fungi_table$dada2_prop,
                            DADA2_OutputProp_COI = inv_table$dada2_prop,
                            DNA_CopySum_16S = prok_table$dna_copy_sum,
                            DNA_CopySum_18S = euk_table$dna_copy_sum,
                            DNA_CopySum_ITS = fungi_table$dna_copy_sum,
                            DNA_CopySum_COI = inv_table$dna_copy_sum
                            )

write_csv(summary_table, "../00_RawFigs/01_Fig_eDNAts/SI_Table_MetaDataSummary.csv")

                            

