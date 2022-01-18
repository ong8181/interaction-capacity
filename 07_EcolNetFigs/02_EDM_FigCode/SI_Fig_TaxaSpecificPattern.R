####
#### CERrice2017 All data analysis
#### Fig. Taxa-specific interaction strength
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
load("../../02_EcolComAnalysis/12_TaxaSpecificISOut/12_TaxaSpecificISOut.RData")

# Load library and functions
library(phyloseq); packageVersion("phyloseq") # 1.30.0, 2020.1.30
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.30
library(reshape2); packageVersion("reshape2") # 1.4.3, 2019.11.12
library(cowplot); packageVersion("cowplot") # 1.0.0, 2020.1.30
library(scales); packageVersion("scales") # 1.1.0, 2020.1.30
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.12
theme_set(theme_cowplot())

# Create output directory
fig_output <- "../00_RawFigs/02_Fig_EDMnet"
dir.create(fig_output)

aaa <- c()
for(i in 1:nrow(edna_tax3)) aaa <- c(aaa, sum(edna_all[,edna_var_coln[i]]>0, na.rm = T))

#-------------------- Visualize IS per species --------------------#
# Summarizing IS per species by taxa
edna_tax3$mean_int_per_sp <- colMeans(edna_ic, na.rm = T)
edna_tax3$converged_int_per_sp <- apply(edna_ic, 2, function(x) median(x, na.rm = T))
edna_tax3$max_int_per_sp <- apply(edna_ic, 2, function(x) max(x[is.finite(x)], na.rm = T))
edna_tax3$int_n <- colMeans(edna_intn, na.rm = T)
edna_tax3$superkingdom <- factor(edna_tax3$superkingdom, levels = c(levels(edna_tax3$superkingdom), "Undetermined"))
edna_tax3$superkingdom[edna_tax3$superkingdom == ""] <- "Undetermined"
edna_tax3$superkingdom <- drop(edna_tax3$superkingdom)
edna_tax3$category <- paste0("[", edna_tax3$superkingdom, "] ", edna_tax3$phylum)
prok_int_med <- median(edna_tax3[prok_id, "converged_int_per_sp"], na.rm = T)
euk_int_med <- median(edna_tax3[euk_id, "converged_int_per_sp"], na.rm = T)

g1 <- ggplot(edna_tax3, aes(x = category, y = converged_int_per_sp, color = superkingdom))
g1 <- g1 + geom_boxplot() + scale_y_log10(limits = c(1e-04, 1.1e01)) 
g1$layers <- g1$layers[-1]
g1 <- g1 + geom_boxplot(width = 0.5, outlier.shape = NULL, outlier.size = 0, outlier.colour = "white")
g1 <- g1 + geom_jitter(size = 1.5, width = 0.1, alpha = 0.5) + scale_color_startrek()
g1 <- g1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) #+ facet_grid(. ~ superkingdom, scales = "free")
g1 <- g1 + geom_segment(x = 2, xend = 16, y = log10(prok_int_med), yend = log10(prok_int_med), linetype = 2, color = "gray50", size = 0.5, alpha = 0.5)
g1 <- g1 + geom_segment(x = 17, xend = 30, y = log10(euk_int_med), yend = log10(euk_int_med), linetype = 2, color = "gray50", size = 0.5, alpha = 0.5)
g1 <- g1 + xlab("Phylum") + ylab("Interaction capacity")

g2 <- ggplot(edna_tax3, aes(x = converged_int_per_sp, fill = superkingdom))
g2 <- g2 + geom_histogram(alpha = 0.6, position = "identity") + scale_x_log10(limits = c(1e-04, 1.1e01)) + scale_fill_startrek()
g2 <- g2 + geom_vline(xintercept = prok_int_med, color = "darkblue", linetype = 2) +
  geom_vline(xintercept = euk_int_med, color = "darkgreen", linetype = 2) +
  xlab(NULL) + ylab(NULL) + theme(axis.text.y = element_blank()) + coord_flip()

g_all1 <- plot_grid(g1 + theme(legend.position = "none"), g2,
                    ncol = 2, rel_widths = c(1, 0.4), align = "hv", axis = "bt")


#-------------------- Visualize the number of ASV v.s. median IS per species (ASV) --------------------#
taxa_rank <- "kingdom"
med_int_phylum <- aggregate(edna_tax3$converged_int_per_sp, by = list(edna_tax3[,taxa_rank]),
                            function(x) median(x, na.rm = T))
n_phylum <- aggregate(edna_tax3[,taxa_rank], by = list(edna_tax3[,taxa_rank]), length)
int_phylum <- data.frame(taxa = med_int_phylum$Group.1,
                         converged_int = med_int_phylum$x,
                         phylum_n_sp = n_phylum$x)
int_phylum <- int_phylum[int_phylum$phylum_n_sp > 0, ]
int_phylum <- int_phylum[int_phylum$taxa != "",]
plot(int_phylum$converged_int, int_phylum$phylum_n_sp)


#-------------------- Visualize taxa-specific pattern --------------------#
# Evaluate all int_n v.s. mean IS curves
edna_ic_melt <- melt(cbind(1:610, clim_all$temp_mean[valid_row_id], edna_prop[pred_lib,"total_div"], edna_ic),
                     id.vars = 1:3)
edna_intn_melt <- melt(cbind(1:610, edna_intn), id.vars = 1)
edna_intmax_melt <- melt(cbind(1:610, edna_intmax), id.vars = 1)
edna_ic_melt$intn <- edna_intn_melt$value
edna_ic_melt$mean_int <- edna_ic_melt$value/edna_ic_melt$intn
edna_ic_melt$max_int <- edna_intmax_melt$value
edna_ic_melt$taxa <- edna_tax3[edna_ic_melt$variable,]$phylum
edna_ic_melt$taxa2 <- edna_tax3[edna_ic_melt$variable,]$class
edna_ic_melt$taxa3 <- edna_tax3[edna_ic_melt$variable,]$superkingdom
colnames(edna_ic_melt) <- c("time", "temp_mean", "total_div", "tax_id", "int_capacity",
                            "int_n", "mean_int", "max_int", "taxa", "taxa2", "taxa3")
edna_ic_melt$taxa3 <- factor(edna_ic_melt$taxa3, levels = c("Archaea", "Bacteria", "Eukaryota", "Undetermined", ""))
edna_ic_melt[edna_ic_melt$taxa3 == "","taxa3"] <- factor("Undetermined")
edna_ic_melt$taxa3 <- drop(edna_ic_melt$taxa3)

g3 <- ggplot(edna_ic_melt, aes(x = int_n, y = mean_int, group = tax_id, color = tax_id, facet = taxa))
g3 <- g3 + geom_point(size = 0.01, alpha = 0.5)
g3 <- g3 + geom_smooth(aes(group = tax_id), method = "loess", span = 1, se = F, size = 0.5)
g3 <- g3 + theme(legend.position = "none") + xlab("The number of interactions") + ylab("Mean interaction strength")
g3 <- g3 + geom_hline(yintercept = 0, linetype = 2)
g3 <- g3 + xlim(0,100) + ylim(0,0.3) + facet_wrap(.~taxa2, ncol = 5)

g4 <- ggplot(edna_ic_melt, aes(x = int_n, y = mean_int, group = tax_id, color = taxa3, facet = taxa3))
g4 <- g4 + geom_point(size = 0.02, alpha = 0.3)
g4 <- g4 + geom_smooth(aes(group = tax_id), method = "loess", span = 1, se = F, size = 0.2, alpha = 0.3)
g4 <- g4 + theme(legend.position = "none")
g4 <- g4 + geom_hline(yintercept = 0, linetype = 2) + scale_color_startrek()
g4 <- g4 + xlim(0,100) + ylim(0,0.3) + facet_wrap(.~taxa3) + panel_border()
g4 <- g4 + xlab("The number of interactions") + ylab("Mean interaction strength")

g5 <- ggplot(edna_ic_melt, aes(x = int_n, y = int_capacity, group = tax_id, color = taxa3, facet = taxa3))
g5 <- g5 + geom_point(size = 0.02, alpha = 0.3)
g5 <- g5 + geom_smooth(aes(group = tax_id), method = "loess", span = 1, se = F, size = 0.2, alpha = 0.3)
g5 <- g5 + theme(legend.position = "none")
g5 <- g5 + geom_hline(yintercept = 0, linetype = 2) + scale_color_startrek()
g5 <- g5 + facet_wrap(.~taxa3) + panel_border()
g5 <- g5 + xlab("The number of interactions") + ylab("(Realized) Interaction capacity")

g6 <- ggplot(edna_ic_melt, aes(x = total_div, y = int_capacity + 0.5, group = tax_id, color = taxa3, facet = taxa3))
g6 <- g6 + geom_point(size = 0.02, alpha = 0.3)
g6 <- g6 + geom_smooth(aes(group = tax_id), method = "loess", span = 1, se = F, size = 0.2, alpha = 0.3)
g6 <- g6 + theme(legend.position = "none")
g6 <- g6 + geom_hline(yintercept = 0, linetype = 2) + scale_color_startrek()
g6 <- g6 + facet_wrap(.~taxa3) + panel_border()
g6 <- g6 + xlab("ASV diversity") + ylab("(Realized) Interaction capacity")

g7 <- ggplot(edna_ic_melt, aes(x = temp_mean, y = max_int, group = tax_id, color = taxa3, facet = taxa3))
g7 <- g7 + geom_point(size = 0.02, alpha = 0.3)
g7 <- g7 + theme(legend.position = "none")
g7 <- g7 + geom_hline(yintercept = 0, linetype = 2) + scale_color_startrek()
g7 <- g7 + xlim(15,31) + ylim(0,1) + facet_wrap(.~taxa3) + panel_border()
g7 <- g7 + xlab(expression(paste("Mean air temperature (",degree,"C)"))) + ylab("Max interaction strength")

# Combine g_all1 and g4
g_all2 <- plot_grid(plot_grid(g4, NULL, rel_widths = c(1,0.3)),
                    g_all1,
                    ncol = 1, align = "hv",
                    rel_heights = c(1, 1), labels = c("a", "b"))


## Save figures
pdf(sprintf("%s/EDMFig_ISperSp_Phylum.pdf", fig_output), width = 12, height = 8)
g_all1; dev.off()
pdf(sprintf("%s/EDMFig_int_link_curve_sp.pdf", fig_output), width = 12, height = 12)
g3; dev.off()
pdf(sprintf("%s/EDMFig_int_link_curve_all.pdf", fig_output), width = 6, height = 6)
g4; dev.off()
ggsave(filename = sprintf("%s/EDMFig_ISpropertiesPerTaxa.jpg", fig_output),
       plot = g_all2, dpi = 200, width = 10, height = 12)
ggsave(filename = sprintf("%s/EDMFig_ISmaxPerTaxa.jpg", fig_output),
       plot = g7, dpi = 200, width = 6, height = 6)

saveRDS(g_all1, sprintf("%s/EDMFig_ISperSp_Phylum.obj", fig_output))
saveRDS(g_all2, sprintf("%s/EDMFig_ISpropertiesPerTaxa.obj", fig_output))
saveRDS(g3, sprintf("%s/EDMFig_int_link_curve_sp.obj", fig_output))
saveRDS(g4, sprintf("%s/EDMFig_int_link_curve_all.obj", fig_output))
saveRDS(g5, sprintf("%s/EDMFig_ISmaxPerTaxa.obj", fig_output))
