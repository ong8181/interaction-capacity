####
#### CERrice2017 All data analysis
#### No.1 Meta-analysis on the diversity-temperature-abundance relationship
####

# Load workspace and objects
load("../02_EcolComAnalysis/12_TaxaSpecificISOut/12_TaxaSpecificISOut.RData")

# Set random seeds (for reproduction)
ran_seed <- 8181
set.seed(ran_seed)
output_folder01 <- "01_MetaAnalysisOut"
dir.create(output_folder01)
dir.create("00_SessionInfo")

# Load library
library(tidyverse); packageVersion("tidyverse") # 1.2.1, 2019.11.21
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.21
library(cowplot); packageVersion("cowplot") # 0.9.4, 2019.11.21
library(scales); packageVersion("scales") # 1.0.0, 2019.11.21
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.21
library(untb); packageVersion("untb") # 1.7.4, 2019.11.21
library(mgcv); packageVersion("mgcv") # 1.8.28, 2019.11.21
source("functions/F04_HelperFunctions.R")


#------------------------- No.1 -------------------------#
#---------- Test prediction using TaraOcean data ----------#
# downloaded from http://ocean-microbiome.embl.de/companion.html on 2019.6.6
tara_d <- read.csv("data_meta_analysis/data_taraocean.csv")
tara_d <- tara_d[,c("silva_richness", "silva_chao", "MeanTemp",
                    "heterotroph", "autotroph", "bacteria", "picoeukaryote")]
tara_d <- na.omit(tara_d)
tara_d$total_biomass <- tara_d$heterotrop + tara_d$autotroph
tara_d$log_div <- log(tara_d$silva_richness)
tara_d$log_temp <- log(tara_d$MeanTemp + 273.15)
tara_d$log_abun <- log(tara_d$total_biomass)
tara_gam_t <- gam(log_div ~ s(log_temp), data = tara_d)
tara_gam_a <- gam(log_div ~ s(log_abun), data = tara_d)
tara_gam <- gam(log_div ~ s(log_temp) + s(log_abun), data = tara_d)
tara_ggdf <- data.frame(observed = tara_d$log_div, predicted = predict(tara_gam))


#------------------------- No.2 -------------------------#
#---------- Test prediction using Bahram et al. (2018) Nature ----------#
# https://doi.org/10.1038/s41586-018-0386-6
bahram_d <- read.csv("data_meta_analysis/data_bahram.csv")
bahram_d <- na.omit(bahram_d[,c("Bacterial_taxonomic_diversity","Fungal_taxonomic_diversity",
                                "Bacterial_biomass","Fungal_biomass", "MAT")])
bahram_d$log_div <- log(bahram_d$Bacterial_taxonomic_diversity + bahram_d$Fungal_taxonomic_diversity)
bahram_d$log_temp <- log(bahram_d$MAT + 273.15)
bahram_d$log_abun <- log(bahram_d$Bacterial_biomass + bahram_d$Fungal_biomass)
bahr_gam_t <- gam(log_div ~ s(log_temp), data = bahram_d)
bahr_gam_a <- gam(log_div ~ s(log_abun), data = bahram_d)
bahr_gam <- gam(log_div ~ s(log_temp) + s(log_abun), data = bahram_d)
bahr_ggdf <- data.frame(observed = bahram_d$log_div, predicted = predict(bahr_gam))


#------------------------- No.3 -------------------------#
#---------- Test prediction using Maizuru data ----------#
maizuru_d <- read.csv("data_meta_analysis/data_maizuru.csv")
fish_col <- 3:100
maizuru_d$maizuru_sp <- rowSums(maizuru_d[,fish_col] > 0)
maizuru_d$total_abn <- rowSums(maizuru_d[,fish_col])
maizuru_d$log_div <- log(maizuru_d$maizuru_sp)
maizuru_d$log_temp <- log(maizuru_d$surf.t + 273.15)
maizuru_d$log_abun <- log(maizuru_d$total_abn)
maiz_gam_t <- gam(log_div ~ s(log_temp), data = maizuru_d)
maiz_gam_a <- gam(log_div ~ s(log_abun), data = maizuru_d)
maiz_gam <- gam(log_div ~ s(log_temp) + s(log_abun), data = maizuru_d)
maiz_ggdf <- data.frame(observed = maizuru_d$log_div, predicted = predict(maiz_gam))


#------------------------- No.4 -------------------------#
#---------- Test prediction using Lake Biwa data ----------#
jplk_d <- read.csv("data_meta_analysis/data_jplakes.csv")
jplk_d <- na.omit(jplk_d)
jplk_d$log_div <- log(jplk_d$otu_richness)
jplk_d$log_abun <- log(jplk_d$dapi_abundance)
jplk_d$log_temp <- log(jplk_d$temp + 273.15)
jplk_gam_t <- gam(log_div ~ s(log_temp), data = jplk_d)
jplk_gam_a <- gam(log_div ~ s(log_abun), data = jplk_d)
jplk_gam <- gam(log_div ~ s(log_temp) + s(log_abun), data = jplk_d)
jplk_ggdf <- data.frame(observed = jplk_d$log_div, predicted = predict(jplk_gam))


#------------------------- No.5 -------------------------#
#---------- Test prediction using Sakamoto et al. (2017) Ecol Res (Data Paper) ----------#
# https://doi.org/10.1007/s11284-017-1528-2
suwa_d <- read.csv("data_meta_analysis/data_sakamoto.csv")
suwa_d$log_div <- log(apply(suwa_d[,3:(ncol(suwa_d)-1)], 1, function(x) sum(x > 0, na.rm = T)))
suwa_d$log_temp <- log(suwa_d$WT + 273.15)
suwa_d$log_abun <- log(rowSums(suwa_d[,3:(ncol(suwa_d)-1)]))
suwa_gam_t <- gam(log_div ~ s(log_temp), data = suwa_d)
suwa_gam_a <- gam(log_div ~ s(log_abun), data = suwa_d)
suwa_gam <- gam(log_div ~ s(log_temp) + s(log_abun), data = suwa_d)
suwa_ggdf <- data.frame(observed = suwa_d$log_div, predicted = predict(suwa_gam))


#------------------------- No.6 -------------------------#
#---------- Test prediction using Okano et al. 2018 Limnology ----------#
# https://doi.org/10.1007/s10201-017-0530-2
okano_d <- read.csv("data_meta_analysis/data_okano.csv")
okano_d$log_div <- log(okano_d$benthos_div)
okano_d$log_temp <- log(okano_d$Temp + 273.15)
okano_d$log_abun <- log(okano_d$total_abundance)
lago_gam_t <- gam(log_div ~ s(log_temp), data = okano_d)
lago_gam_a <- gam(log_div ~ s(log_abun), data = okano_d)
lago_gam <- gam(log_div ~ s(log_temp) + s(log_abun), data = okano_d)
lago_ggdf <- data.frame(observed = okano_d$log_div, predicted = predict(lago_gam))


#------------------------- Check AIC, BIC and R2 -------------------------#
AIC(tara_gam, tara_gam_a, tara_gam_t) # temp & abundance
AIC(bahr_gam, bahr_gam_a, bahr_gam_t) # temp & abundance
AIC(maiz_gam, maiz_gam_a, maiz_gam_t) # temp & abundance
AIC(jplk_gam, jplk_gam_a, jplk_gam_t) # temp & abundance
AIC(suwa_gam, suwa_gam_a, suwa_gam_t) # temp & abundance
AIC(lago_gam, lago_gam_a, lago_gam_t) # temp & abundance

BIC(tara_gam, tara_gam_a, tara_gam_t) # temp
BIC(bahr_gam, bahr_gam_a, bahr_gam_t) # temp & abundance
BIC(maiz_gam, maiz_gam_a, maiz_gam_t) # temp & abundance
BIC(jplk_gam, jplk_gam_a, jplk_gam_t) # temp & abundance
BIC(suwa_gam, suwa_gam_a, suwa_gam_t) # temp & abundance
BIC(lago_gam, lago_gam_a, lago_gam_t) # temp & abundance or abundance

(tara_r_adj <- summary(tara_gam)$r.sq)
(bahr_r_adj <- summary(bahr_gam)$r.sq)
(maiz_r_adj <- summary(maiz_gam)$r.sq)
(jplk_r_adj <- summary(jplk_gam)$r.sq)
(suwa_r_adj <- summary(suwa_gam)$r.sq)
(lago_r_adj <- summary(lago_gam)$r.sq)


#------------------------- Add figure title and R2 values -------------------------#
t1 <- geom_obspred(tara_ggdf, "TaraOceans expedition, Metagenome") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(tara_r_adj, digits = 3))))
t2 <- geom_obspred(bahr_ggdf, "Global soil, Microbes") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(bahr_r_adj, digits = 3))))
t3 <- geom_obspred(maiz_ggdf, "Coastal region, Fish") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(maiz_r_adj, digits = 3))))
t4 <- geom_obspred(jplk_ggdf, "Japanese lakes, Microbes") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(jplk_r_adj, digits = 3))))
t5 <- geom_obspred(suwa_ggdf, "Lake Suwa, Zooplanktons") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(suwa_r_adj, digits = 3))))
t6 <- geom_obspred(lago_ggdf, "Lagoons, Macroinvertebrates") +
  theme(plot.title = element_text(size = 12)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -2,
           label = bquote(R^2 == .(format(lago_r_adj, digits = 3))))

# Output figure
t_all <- plot_grid(t1, t2, t3, t4, t5, t6,
                   labels = "auto", align = "hv")
pdf(sprintf("%s/DiversityMetaAnal.pdf", output_folder01), width = 12, height = 8)
t_all; dev.off()

# Save workspace and ojcects
#save.image(sprintf("%s/%s.RData", output_folder01, output_folder01))
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", output_folder01, output_folder01))

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/%s_SessionInfo_%s.txt", output_folder01, substr(Sys.time(), 1, 10)))
