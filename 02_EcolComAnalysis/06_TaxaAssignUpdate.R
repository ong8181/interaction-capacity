####
#### CERrice2017 All data analysis
#### No. 6 Update of the assigned taxa
####

# Load workspace
load("05_NetworkPrepOut/05_NetworkPrepOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder06 <- "06_TaxaAssignUpdateOut"
dir.create(output_folder06)

# Update taxa assignment?
taxa_update <- F

if(taxa_update){
  # Output selected sequences
  write.table(as.matrix(c(rbind(sprintf(">%s", rownames(edna_tax2)),
                                as.character(edna_tax2$seq))), ncol = 1),
              sprintf("%s/SelectedASV.fa", output_folder06),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  
  # Perform Claident commonds from R
  setwd("06_TaxaAssignUpdateOut")
  
  # Search neighbor sequences using overall_genus and overall_taxa
  system("clmakecachedb --blastdb=overall_genus --numthreads=72 SelectedASV.fa SelectedASV_overallg_cache")
  system("clidentseq --blastdb=SelectedASV_overallg_cache --numthreads=72 SelectedASV.fa SelectedASV_overallg_clidentseq")
  system("clmakecachedb --blastdb=overall_class --numthreads=72 SelectedASV.fa SelectedASV_overallc_cache")
  system("clidentseq --blastdb=SelectedASV_overallc_cache --numthreads=72 SelectedASV.fa SelectedASV_overallc_clidentseq")
  
  # Assgin taxa name
  system("classigntax --taxdb=overall_genus --maxpopposer=0.05 --minsoratio=19 SelectedASV_overallg_clidentseq SelectedASV_overallg_classigntax")
  system("classigntax --taxdb=overall_class --maxpopposer=0.05 --minsoratio=19 SelectedASV_overallc_clidentseq SelectedASV_overallc_classigntax")
  #system("clmergeassign --priority=equal --preferlower SelectedASV_overallg_classigntax SelectedASV_overallc_classigntax SelectedASV_merge_classigntax")
  system("clmergeassign --priority=descend SelectedASV_overallg_classigntax SelectedASV_overallc_classigntax SelectedASV_merge_classigntax")
  
  # Save workspace and ojcects
  setwd("..")
}

save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder06, output_folder06))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder06, substr(Sys.time(), 1, 10)))
