####
#### CER eDNA study
#### No.9 Filtering time series preparation (Shannon entropy + DNA concentration)
####

# Load workspace
load("08_TSfilterPrepOut/08_TSfilterPrepOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder09 <- "09_TSfilter01Out"
dir.create(output_folder09)

# Load library and functions
library(phyloseq); packageVersion("phyloseq") #1.28.0, 2020.1.6

# Add entropy values to phyloseq objects
pro_taxtab_tmp <- data.frame(tax_table(ps_pro_sample2)@.Data)
fun_taxtab_tmp <- data.frame(tax_table(ps_fun_sample2)@.Data)
inv_taxtab_tmp <- data.frame(tax_table(ps_inv_sample2)@.Data)
euk_taxtab_tmp <- data.frame(tax_table(ps_euk_sample2)@.Data)

pro_taxtab_tmp$entropy <- pro_ent
fun_taxtab_tmp$entropy <- fun_ent
inv_taxtab_tmp$entropy <- inv_ent
euk_taxtab_tmp$entropy <- euk_ent

tax_table(ps_pro_sample2) <- tax_table(as.matrix(pro_taxtab_tmp))
tax_table(ps_fun_sample2) <- tax_table(as.matrix(fun_taxtab_tmp))
tax_table(ps_inv_sample2) <- tax_table(as.matrix(inv_taxtab_tmp))
tax_table(ps_euk_sample2) <- tax_table(as.matrix(euk_taxtab_tmp))

# Filter time series based on DNA copy numbers and Shannon entropy
ll_pro; ll_fun; ll_inv; ll_euk
# Exclude taxa of which average DNA copy number is lower than the lower limit/100
pro_cond1 <- taxa_sums(ps_pro_sample2) > nrow(sample_data(ps_pro_sample2))*ll_pro/100
fun_cond1 <- taxa_sums(ps_fun_sample2) > nrow(sample_data(ps_fun_sample2))*ll_fun/100
inv_cond1 <- taxa_sums(ps_inv_sample2) > nrow(sample_data(ps_inv_sample2))*ll_inv/100
euk_cond1 <- taxa_sums(ps_euk_sample2) > nrow(sample_data(ps_euk_sample2))*ll_euk/100

pro_cond2 <- pro_ent > q90; fun_cond2 <- fun_ent > q90
inv_cond2 <- inv_ent > q90; euk_cond2 <- euk_ent > q90

sum(pro_cond1); sum(fun_cond1); sum(inv_cond1); sum(euk_cond1) # 772; 787; 413; 432
sum(pro_cond2); sum(fun_cond2); sum(inv_cond2); sum(euk_cond2) # 570; 404; 154; 267
sum(pro_cond1 & pro_cond2) # 516
sum(fun_cond1 & fun_cond2) # 400
sum(inv_cond1 & inv_cond2) # 147
sum(euk_cond1 & euk_cond2) # 241

ps_pro_sample3 <- prune_taxa(pro_cond1 & pro_cond2, ps_pro_sample2)
ps_fun_sample3 <- prune_taxa(fun_cond1 & fun_cond2, ps_fun_sample2)
ps_inv_sample3 <- prune_taxa(inv_cond1 & inv_cond2, ps_inv_sample2)
ps_euk_sample3 <- prune_taxa(euk_cond1 & euk_cond2, ps_euk_sample2)

# Save and output results
#save(list = ls(all.names = TRUE),
#     file = sprintf("%s/09_TSfilter01Out.RData", output_folder09))
save.image(sprintf("%s/09_TSfilter01Out.RData", output_folder09))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/09_SessionInfo_TSfilter01_%s.txt", substr(Sys.time(), 1, 10)))
