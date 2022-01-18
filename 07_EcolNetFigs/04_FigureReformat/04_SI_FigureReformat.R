####
#### CER eDNA study
#### FigCode: Re-formatting all supplementary figures
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the figure file(s)!

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2021.7.29
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2021.7.29
library(cowplot); packageVersion("cowplot") # 1.1.1, 2021.7.29
library(scales); packageVersion("scales") # 1.1.1, 2021.7.29
library(ggimage); packageVersion("ggimage") # 0.3.0, 2021.12.8
library(magick); packageVersion("magick") # 2.7.3, 2021.12.8
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.12
theme_set(theme_cowplot())
# Create output folder
dir.create("../00_ReformatFigs")


#-------------------------------------------------#
# Load saved figure objects
#-------------------------------------------------#
## Experimental design
# Plot images
Fig_PlotImage <- image_read("../00_RawFigs/00_ExpImages/PlotImages.png")

## eDNA data 
# Rarefaction curve
Fig_ProRarefaction <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_Rarefaction_pro.obj")
Fig_EukRarefaction <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_Rarefaction_euk.obj")
Fig_FunRarefaction <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_Rarefaction_fun.obj")
Fig_InvRarefaction <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_Rarefaction_inv.obj")
Legend_Rarefaction <- get_legend(Fig_ProRarefaction + theme(legend.position = "right", legend.justification = "center", legend.title = element_blank()))
# Negative and positive controls
Fig_ProkNC <- readRDS("../00_RawFigs/01_Fig_eDNAts/SI_Fig_ProkNC01.obj")
Fig_ProkPC <- readRDS("../00_RawFigs/01_Fig_eDNAts/SI_Fig_ProkPC01.obj")
Fig_EukNC <- readRDS("../00_RawFigs/01_Fig_eDNAts/SI_Fig_EukNC01.obj")
Fig_EukPC <- readRDS("../00_RawFigs/01_Fig_eDNAts/SI_Fig_EukPC01.obj")
Fig_FungiNC <- readRDS("../00_RawFigs/01_Fig_eDNAts/SI_Fig_FungiNC01.obj")
Fig_FungiPC <- readRDS("../00_RawFigs/01_Fig_eDNAts/SI_Fig_FungiPC01.obj")
Fig_InvNC <- readRDS("../00_RawFigs/01_Fig_eDNAts/SI_Fig_InvNC01.obj")
Fig_InvPC <- readRDS("../00_RawFigs/01_Fig_eDNAts/SI_Fig_InvPC01.obj")
# Diversity pattern
Fig_DivTS_Prok <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_ProkDivTS.obj")
Fig_DivTS_Fungi <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_FungiDivTS.obj")
Fig_DivTS_Metaz <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_MetazDivTS.obj")
Fig_DivTS_Virid <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_ViridiDivTS.obj")
# Major phyla
Fig_ProkBar <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_ProkBarAll2.obj") + guides(fill=guide_legend(title="Major phyla"))
Fig_FungiBar <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_FungiBarAll2.obj") + guides(fill=guide_legend(title="Major phyla"))
Fig_InvBar <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_InvBarAll2.obj") + guides(fill=guide_legend(title="Major phyla"))
Fig_EukBar <- readRDS("../00_RawFigs/01_Fig_eDNAts/Fig_EukBarAll2.obj") + guides(fill=guide_legend(title="Major phyla"))

## Time series analysis
# Time-varying network properties
Fig_TimeVaryingProp <- readRDS("../00_RawFigs/02_Fig_EDMnet/EDMFig_TimeVaryingNetProps2.obj")
# Random shuffle time series
Fig_RandomSurProp <- readRDS("../00_RawFigs/02_Fig_EDMnet/NetworkProperties_Compare.obj")
# Species-level IS patterns
Fig_SpIntCom <- as.ggplot(image_read("../00_RawFigs/02_Fig_EDMnet/EDMFig_ISpropertiesPerTaxa.jpg"))
# Species-level connectance
Fig_SpConnect <- as.ggplot(image_read("../../06_AdditionalAnalysis/06_SpLevelConnectance/01_SpLevelConnectanceOut/SpLevel_Connectance.jpg"))
# Species-level IC v.s. 1/abundance
Fig_SpICDNA <- as.ggplot(image_read("../../06_AdditionalAnalysis/07_SpICversusDNA/01_SpICversusDNAOut/SpLevel_Capacity.jpg"))
# Max interaction strength (community- and species-level)
Fig_MaxIntCom <- readRDS("../00_RawFigs/02_Fig_EDMnet/NetworkProperties_IntMax.obj")
Fig_NetPropSI <- readRDS("../00_RawFigs/02_Fig_EDMnet/NetworkProperties_SI.obj")
# EDM between network properties
Fig_NetPropCCM <- readRDS("../00_RawFigs/02_Fig_EDMnet/NetworkPropCCMmatrix_mean.obj")

## Meta analysis
Fig_MetaAnalysis <- readRDS("../00_RawFigs/03_Fig_MetaAnalysis/DiversityMetaAnal.obj")
Fig_MetaAnalysis2 <- readRDS("../../06_AdditionalAnalysis/05_OtherDataSets/01_EditFigsOut/Ratzke_Yamawo.obj")


#-------------------------------------------------#
# Reformat figures
#-------------------------------------------------#
## Positive and negative controls
NC_legend <- get_legend(Fig_ProkNC + theme(legend.position = "right", legend.justification = "center"))
Fig_AllNCRare <- plot_grid(Fig_ProRarefaction + theme(legend.position = "none"),
                           Fig_EukRarefaction + theme(legend.position = "none"),
                           Fig_FunRarefaction + theme(legend.position = "none"),
                           Fig_InvRarefaction + theme(legend.position = "none"),
                           Legend_Rarefaction,
                           Fig_ProkNC + theme(legend.position = "none"),
                           Fig_EukNC + theme(legend.position = "none"),
                           Fig_FungiNC + theme(legend.position = "none"),
                           Fig_InvNC + theme(legend.position = "none"),
                           NC_legend,
                           ncol = 5, rel_widths = c(1,1,1,1,0.4),
                           labels = c("(a)","(b)","(c)","(d)",NA,
                                      "(e)","(f)","(g)","(h)",NA),
                           align = "h", axis = "lr")

PC_legend <- get_legend(Fig_ProkPC + theme(legend.position = "right", legend.justification = "left"))
Fig_AllPC0 <- plot_grid(Fig_ProkPC + theme(legend.position = "none"),
                        Fig_EukPC + theme(legend.position = "none"),
                        Fig_FungiPC + theme(legend.position = "none"),
                        Fig_InvPC + theme(legend.position = "none"),
                        rel_widths = c(1,1,1,1),
                        ncol = 4, labels = c("(i)","(j)","(k)","(l)"),
                        align = "hv", axis = "btlr")
Fig_AllPC <- plot_grid(Fig_AllPC0, PC_legend, rel_widths = c(1,0.1))
Fig_AllNCPCRare <- plot_grid(Fig_AllNCRare, Fig_AllPC, ncol = 1, rel_heights = c(1,0.6), labels = NULL)

## Barplots
DivTS_legend <- get_legend(Fig_DivTS_Prok + theme(legend.position = "bottom", legend.justification = "center"))
Fig_DivTS_subset <- plot_grid(Fig_DivTS_Prok + ggtitle("Prokaryotes") + theme(legend.position = "none"),
                              Fig_DivTS_Virid + ggtitle("Viridiplantae") + theme(legend.position = "none"),
                              Fig_DivTS_Fungi + ggtitle("Fungi") + theme(legend.position = "none"),
                              Fig_DivTS_Metaz + ggtitle("Metazoan") + theme(legend.position = "none"),
                              DivTS_legend, labels = c("(a)","(c)","(e)","(g)",NA),
                              ncol = 1, rel_heights = c(1,1,1,1,0.2), hjust = -0.02)
Fig_BarTS_major <- plot_grid(Fig_ProkBar + ggtitle("Major phyla (16S sequencing)"),
                             Fig_EukBar + ggtitle("Major phyla (18S sequencing)"),
                             Fig_FungiBar + ggtitle("Major phyla (ITS sequencing)"),
                             Fig_InvBar + ggtitle("Major orders (COI sequencing)"),
                             NULL,
                             rel_heights = c(1,1,1,1,0.2),
                             ncol = 1, labels = c("(b)","(d)","(f)","(h)",NA),
                             hjust = -0.02, align = "v")


#-------------------------------------------------#
# New results to address the reviewers' comments
#-------------------------------------------------#
## DNA quantification
Fig_StdDNA <- readRDS("../../06_AdditionalAnalysis/03_StandardDNAR2/01_VisualizeStdCurveOut/StdDNA_QualityCheck.obj")
Fig_QuantusMiSeq <- readRDS("../../04_QuantifyDNA/01_VisualizeResultsOut/qMiSeq_Comparison_all.obj")
Fig_Metagenome <- readRDS("../../05_ShotogunMetagenome/15_VisualizeAllOut/qMiSeq_Shotgun_All.obj")

Fig_ValidateMiSeq1 <- plot_grid(Fig_StdDNA[[1]],
                                Fig_StdDNA[[2]],
                                ncol = 2,
                                labels = c("(a)","(b)"),
                                align = "hv", axis = "lrbt")
legend_plot <- get_legend(Fig_QuantusMiSeq[[1]])
Fig_ValidateMiSeq2 <- plot_grid(Fig_QuantusMiSeq[[1]] + theme(legend.position = "none"),
                                Fig_QuantusMiSeq[[2]] + theme(legend.position = "none"),
                                Fig_QuantusMiSeq[[3]] + theme(legend.position = "none"),
                                legend_plot,
                                align = "hv", axis = "lrtb",
                                ncol = 4, rel_widths = c(1,1,1,0.2),
                                labels = c("(c)","(d)","(e)", NA))
legend_meta <- get_legend(Fig_Metagenome[[1]])

fm1 <- Fig_Metagenome[[1]] +
  scale_x_discrete(labels = c("qMiSeq", "Shotgun")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin = margin(.8,.8,.8,.8,"cm"))
fm2 <- Fig_Metagenome[[2]] +
  scale_x_discrete(labels = c("qMiSeq", "Shotgun")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin = margin(.8,.8,.8,.8,"cm"))
fm3 <- Fig_Metagenome[[3]] +
  scale_x_discrete(labels = c("qMiSeq", "Shotgun")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin = margin(.8,.8,.8,.8,"cm"))
fm4 <- Fig_Metagenome[[4]] +
  scale_x_discrete(labels = c("qMiSeq", "Shotgun")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin = margin(.8,.8,.8,.8,"cm"))

Fig_ValidateMiSeq3 <- plot_grid(fm1, fm2, fm3, fm4, NULL, legend_meta,
                   ncol = 6, rel_widths = c(1,1,1,1,0.2,1),
                   labels = c("(f)","(g)","(h)","(i)",NA,NA),
                   align = "hv", axis = "lrbt")
Fig_ValidateMiSeq <- plot_grid(Fig_ValidateMiSeq1,
                               Fig_ValidateMiSeq2,
                               Fig_ValidateMiSeq3,
                               ncol = 1,
                               rel_heights = c(0.9, 0.7, 1.2))

## Subset stability
Fig_SubsetStab <- readRDS("../../06_AdditionalAnalysis/04_SubsetStability/01_SubsetStabilityOut/Subset_Stability2.obj")
Fig_Nint <- readRDS("../../06_AdditionalAnalysis/02_NegativePositiveInt/01_CompileNegPosOut/NintHistogram.obj")
Fig_StabGrid <- plot_grid(Fig_SubsetStab[[1]] + theme(axis.title = element_text(size = 9)),
                          Fig_SubsetStab[[2]] + theme(axis.title = element_text(size = 9)),
                          Fig_SubsetStab[[3]] + theme(axis.title = element_text(size = 9)),
                          Fig_SubsetStab[[4]] + theme(axis.title = element_text(size = 9)),
                          ncol = 4, align = "hv", axis = "lrtb",
                          labels = c("(d)","(e)","(f)","(g)"),
                          vjust = -1)
Fig_TimeVaryingProp2 <- plot_grid(Fig_Nint + xlab("No. of interactions\nfor eash ASV on each day"),
                                  Fig_TimeVaryingProp[[2]], ncol = 2, rel_widths = c(1, 2.3),
                                  align = "hv", axis = "lrbt", labels = c("(b)","(c)"))
Fig_TimeVaryingProp3 <- plot_grid(Fig_TimeVaryingProp[[1]],
                                  Fig_TimeVaryingProp2, rel_heights = c(1, 1), ncol = 1) + 
  theme(plot.margin = unit(c(0,0,0,1),"cm"))

Fig_TimeVaryingGrid <- plot_grid(Fig_TimeVaryingProp3,
                                 NULL,
                                 Fig_StabGrid,
                                 ncol = 1, rel_heights = c(3,0.2,1))


#-------------------------------------------------#
# Species-level patterns
#-------------------------------------------------#
Fig_Sp2 <- plot_grid(Fig_SpConnect,
                     Fig_SpICDNA, ncol = 1,
                     rel_heights = c(1,1.5))
Fig_SpIC_all <- plot_grid(Fig_SpIntCom, Fig_Sp2, ncol = 2,
                          rel_widths = c(1, 0.8))

#-------------------------------------------------#
# Meta-analysis
#-------------------------------------------------#
Fig_MetaGrid2 <- plot_grid(Fig_MetaAnalysis2[[1]],
                           Fig_MetaAnalysis2[[2]] + labs(fill = "Season"),
                           align = "hv", axis = "lrbt",
                           labels = c("(g)", "(h)"))
Fig_MetaAnalysis3 <- plot_grid(Fig_MetaAnalysis,
                               NULL,
                               Fig_MetaGrid2,
                               ncol = 1, rel_heights = c(2, 0.1, 1))


#-------------------------------------------------#
# Save figures
#-------------------------------------------------#
# Figure S1
pdf("../00_ReformatFigs/SI_Fig_01.pdf", width = 10, height = 8)
as.ggplot(Fig_PlotImage); dev.off()

# Figure S2
pdf("../00_ReformatFigs/SI_Fig_02.pdf", width = 24, height = 18)
Fig_AllNCPCRare; dev.off()

# Figure S3; to address the reviewers' comments
pdf("../00_ReformatFigs/SI_Fig_03.pdf", width = 12, height = 14)
Fig_ValidateMiSeq; dev.off()

# Figure S4
pdf("../00_ReformatFigs/SI_Fig_04.pdf", width = 18, height = 20)
plot_grid(Fig_DivTS_subset, Fig_BarTS_major,
          ncol = 2, rel_widths = c(1,1.8)) +
  theme(plot.margin = unit(c(0,0,0,1), "cm"))
dev.off()
 
# Figure S5
pdf("../00_ReformatFigs/SI_Fig_05.pdf", width = 10, height = 12)
Fig_TimeVaryingGrid; dev.off()

# Figure S6
# Additional network properties (e.g., the number of interactions)
pdf("../00_ReformatFigs/SI_Fig_06.pdf", width = 12, height = 8)
Fig_NetPropSI; dev.off()

# Figure S7
pdf("../00_ReformatFigs/SI_Fig_07.pdf", width = 12, height = 8)
Fig_RandomSurProp; dev.off()

# Figure S8
pdf("../00_ReformatFigs/SI_Fig_08.pdf", width = 12, height = 8)
Fig_SpIC_all; dev.off()

# Figure S9
pdf("../00_ReformatFigs/SI_Fig_09.pdf", width = 8, height = 7.5)
Fig_NetPropCCM; dev.off()

# Figure S10
pdf("../00_ReformatFigs/SI_Fig_10.pdf", width = 12, height = 12)
Fig_MetaAnalysis3; dev.off()

# Table S1
system("cp ../00_RawFigs/02_Fig_EDMnet/Table_GAMsummary.pdf ../00_ReformatFigs/SI_Table_01.pdf")
