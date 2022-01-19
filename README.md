# R and shell scripts for Ushio (in press) _Proceedings of the Royal Society B: Biological Sciences_
R and shell codes to reproduce the results in Ushio M (in press) "Interaction capacity as a potential driver of community diversity" _Proceedings of the Royal Society B: Biological Sciences_.

(Preprint: "Interaction capacity underpins community diversity" _bioRxiv_ https://doi.org/10.1101/2020.04.08.032524)

**Note**: In this repository, only core result files are included due to the file size limitation. You may reproduce the main figures using the information available in this repository, but if you want to reproduce the full results, please re-run all analyses using DNA sequence data (deposited in DDBJ; see the manuscript) or contact the author (ong8181@gmail.com).

# License
see LICENSE

# Core softwares
- `Claident` (https://www.claident.org/)
- `DADA2` (https://benjjneb.github.io/dada2/index.html)
- `rEDM` (https://ha0ye.github.io/rEDM/index.html)
- Regularized S-map (custom scripts) (https://github.com/ong8181/random-scripts/tree/master/02_RegularizedSmap)

# Analysis workflow
## No.1: Sequence analysis (`01_DNAtsCERrice2017` folder)
R and shell scripts for analyzing DNA sequences  
1. Sequence data (fastq files) generated by MiSeq was analyzed using `Claident` and `dada2`
    - Amplicon Sequence Variants (ASVs) were picked up by DADA2
    - Taxa were assigned by `Claident`
2. Sequence reads were converted to DNA copy numbers
    - `02_TimeSeriesCompile` folder
 
## No.2: Empirical dynamic modeling (`02_EcolComAnalysis` folder)
R scripts for nonlinear time series analysis (Empirical dynamic modeling).

## No.3: Meta-analysis (`03_MetaAnalysis` folder)
R scripts for the meta-analysis. The data for the meta-analysis were collected by contacting authors of the original publications, and thus are not included in the repository.

# Codes to revise the paper
## No.4: Quantitative metabarcoding v.s. other quantification methods (`04_QunatifyDNA` folder)
R scripts to compare the quantitative capacity of quantitative metabarcoding v.s. other quantification methods.

## No.5: Quantitative metabarcoding v.s. Shotgun metagenome analysis (`05_ShotgunMetagenome` folder)
Compare the community compositions detected by quantitative metabarcoding and shotgun metagenome analysis.

## No.6: Additional analyses (`06_AdditionalAnalysis` folder)
Additional statistical analyses.

## No.7: Visualization (`07_EcoNetFigs` folder)
Codes to create figures.

## data (`data_compiled` folder)
Climate data, DNA-based community data, sequence summary, and RData file are included.




