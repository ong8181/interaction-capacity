####
#### CERrice2017 All data analysis
#### Fig. Network Figure, version 2 (reference https://www.data-to-viz.com/)
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the RData file(s)!

# Load workspace
load("../../02_EcolComAnalysis/07_RegularizedSmapOut/07_RegularizedSmapOut.RData")

# Ceate output director
# Load library and functions
library(phyloseq); packageVersion("phyloseq") # 1.30.0, 2020.1.30
library(ggraph); packageVersion("ggraph") # 2.0.0, 2019.11.11
library(igraph); packageVersion("igraph") # 1.2.4.2, 2020.1.30
library(tidyverse); packageVersion("tidyverse") # 1.3.0, 2020.1.30
library(RColorBrewer); packageVersion("RColorBrewer") # 1.1.2, 2019.11.11
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.11
library(lubridate); packageVersion("lubridate") # 1.7.4, 2019.11.11

# Create output directory
fig_output <- "../00_RawFigs/02_Fig_EDMnet"
dir.create(fig_output)
network_output_folder <- "01_NetworkFigs_v2"
dir.create(network_output_folder)

#-------------------- Time-averaged network figures --------------------#
set.seed(8181)

# Rename taxa and collect additional information
taxa_cols <- colnames(edna_tax3)[c(2,3,5,6,10,14,17,19)]

# Assign edge attribute: Step 1
edna_tax3$superkingdom <- factor(edna_tax3$superkingdom, levels = c(levels(edna_tax3$superkingdom)[2:4], "Undetermined", levels(edna_tax3$superkingdom)[1]))
edna_tax3$phylum <- factor(edna_tax3$phylum, levels = c(levels(edna_tax3$phylum)[2:28], "Undetermined", levels(edna_tax3$phylum)[1]))
edna_tax3$superkingdom[edna_tax3$superkingdom == ""] <- "Undetermined"
edna_tax3$phylum[edna_tax3$phylum == ""] <- "Undetermined"

edna_tax3$taxa <- factor(edna_tax3$superkingdom, levels = c(levels(edna_tax3$superkingdom)[1:4], levels(edna_tax3$kingdom), "Other Eukaryota"))
edna_tax3$taxa[edna_tax3$superkingdom == "Eukaryota"] <- edna_tax3$kingdom[edna_tax3$superkingdom == "Eukaryota"]
edna_tax3$taxa[edna_tax3$superkingdom == "Eukaryota"][edna_tax3$taxa[edna_tax3$superkingdom == "Eukaryota"] == ""] <- "Other Eukaryota"
edna_tax3$taxa <- factor(edna_tax3$taxa, levels = c("Archaea", "Bacteria", "Fungi", "Metazoa", "Viridiplantae", "Other Eukaryota", "Undetermined"))

edna_tax3$tax_id <- rownames(edna_tax3)
edna_tax4 <- dplyr::arrange(edna_tax3, taxa, phylum)
rownames(edna_tax4) <- edna_tax4$tax_id

# Align interaction matrix
taxa_edge <- rbind(data.frame(from = "origin", to = as.character(levels(edna_tax4$taxa))),
                   data.frame(from = "Archaea", to = edna_tax4[edna_tax4$taxa == "Archaea", "tax_id"]),
                   data.frame(from = "Bacteria", to = edna_tax4[edna_tax4$taxa == "Bacteria", "tax_id"]),
                   data.frame(from = "Fungi", to = edna_tax4[edna_tax4$taxa == "Fungi", "tax_id"]),
                   data.frame(from = "Metazoa", to = edna_tax4[edna_tax4$taxa == "Metazoa", "tax_id"]),
                   data.frame(from = "Viridiplantae", to = edna_tax4[edna_tax4$taxa == "Viridiplantae", "tax_id"]),
                   data.frame(from = "Other Eukaryota", to = edna_tax4[edna_tax4$taxa == "Other Eukaryota", "tax_id"]),
                   data.frame(from = "Undetermined", to = edna_tax4[edna_tax4$taxa == "Undetermined", "tax_id"])
                   )


# Create a vertices data frame. One line per object of our hierarchy
taxa_vertices <- rbind(data.frame(name = c("origin", as.character(levels(edna_tax4$taxa))),
                                  log_abundance = NA,
                                  group = c(NA, rep("origin", length(unique(edna_tax4$taxa)))),
                                  subgroup = c(NA, rep("origin", length(unique(edna_tax4$taxa))))),
                       data.frame(name = rownames(edna_tax4),
                                  log_abundance = log(colSums(edna_all[,rownames(edna_tax4)], na.rm = T)),
                                  group = edna_tax4$taxa,
                                  subgroup = edna_tax4$phylum)
                       )
taxa_vertices$std_abundance <- (taxa_vertices$log_abundance - min(taxa_vertices$log_abundance, na.rm = T))/(max(taxa_vertices$log_abundance, na.rm = T) - min(taxa_vertices$log_abundance, na.rm = T))
taxa_vertices$std_abundance <- round(taxa_vertices$std_abundance, digits = 2)/20
  
# Add information concerning the label we are going to add: angle, horizontal adjustement and potential flip
taxa_vertices$id <- NA
my_leaves <- which(is.na(match(taxa_vertices$name, taxa_edge$from)))
n_leaves <- length(my_leaves)
taxa_vertices$id[my_leaves] <- seq(1:n_leaves)
taxa_vertices$angle <- 90 - 360 * (taxa_vertices$id + 300) / n_leaves


# Add subgroup2 names
taxa_vertices$subgroup2 <- NaN
for(i in 2:(nrow(taxa_vertices)-1)){
  if(as.character(taxa_vertices$subgroup[i]) == as.character(taxa_vertices$subgroup[i+1])){
    taxa_vertices$subgroup2[i] <- NaN
  }else{
    taxa_vertices$subgroup2[i] <- as.character(taxa_vertices$subgroup[i])
  }
}

taxa_vertices[taxa_vertices$subgroup2 == "NaN",]$subgroup2 <- ""
taxa_vertices$hjust <- 1

# Prepare time-averaged matrix
J1_matrix_sum <- matrix(0, ncol = nrow(edna_tax4), nrow = nrow(edna_tax4))
colnames(J1_matrix_sum) <- rownames(J1_matrix_sum) <- rownames(edna_tax4)
non_na_time <- 0; valid_times <- c(NULL)

for(time_i in 1:length(compiled_smap_all)){
  if(is.na(compiled_smap_all[[time_i]])[1]){
    J1_matrix_sum <- J1_matrix_sum + 0
  }else{
    J1_tmp <- compiled_smap_all[[time_i]]$J1_matrix_all[[1]]
    if(!sum(is.na(J1_tmp)) > 0){
      J1_matrix_sum[rownames(J1_tmp), colnames(J1_tmp)] <- J1_matrix_sum[rownames(J1_tmp), colnames(J1_tmp)] + J1_tmp
      non_na_time <- non_na_time + 1
      valid_times[time_i] <- time_i
      rm(J1_tmp)
    }else{
      J1_matrix_sum <- J1_matrix_sum + 0
      rm(J1_tmp)
    }
  }
}

J1_matrix_avr_nodiag <- J1_matrix_avr <- J1_matrix_sum/non_na_time
diag(J1_matrix_avr_nodiag) <- 0
J1_matrix_avr_adj <- J1_matrix_avr_nodiag
J1_matrix_avr_adj[J1_matrix_avr_adj != 0] <- 1

# Creat edge list
eco_connect <- as.data.frame(get.edgelist(graph.adjacency(J1_matrix_avr, weighted = TRUE, diag = FALSE, mode = "directed")))
colnames(eco_connect) <- c("from", "to")
eco_connect$value <- c(t(J1_matrix_avr_nodiag[J1_matrix_avr_nodiag != 0]))

# Create a plot
eco_graph <- igraph::graph_from_data_frame(taxa_edge, vertices=taxa_vertices)
eco_from <- match(eco_connect$from, taxa_vertices$name)
eco_to <- match(eco_connect$to, taxa_vertices$name)
eco_strength <- abs(eco_connect$value)

# Basic usual argument
gn1 <- ggraph(eco_graph, layout = 'dendrogram', circular = TRUE)
# Add phylum information
phylum_path <- data.frame(x = gn1$data$x, y = gn1$data$y)[9:1205,] * 1.1

for(id_i in 2:27){
  id_tmp <- which(taxa_vertices$subgroup == levels(taxa_vertices$subgroup)[id_i])
  if(length(id_tmp) > 15){
    id_tmp <- head(id_tmp - 10, length(id_tmp)-1)
    gn1 <- gn1 + annotate("path",
                          x = phylum_path$x[id_tmp],
                          y = phylum_path$y[id_tmp],
                          size = 5, color = "gray50")
    gn1  <- gn1 + annotate("text",
                          x = phylum_path$x[round(quantile(id_tmp, probs = 0.5))] * 1.15,
                          y = phylum_path$y[round(quantile(id_tmp, probs = 0.5))] * 1.15,
                          label = levels(taxa_vertices$subgroup)[id_i], size = 7)
  }
}
#gn1

gn2 <- gn1 + geom_node_point(aes(filter = leaf,
                                 x = x*1.03, y = y*1.03,
                                 size = std_abundance,
                                 fill = group),
                             color = "gray10",
                             shape = 21,
                             alpha = 0.5)
gn2 <- gn2 + geom_conn_bundle(data = get_con(from = eco_from, to = eco_to),
                              edge_alpha = 0.3,
                              width = 0.2,
                              aes(color = group),
                              arrow = grid::arrow(angle = 20, length = unit(0.01, "inches"), type = "closed"))

gn2 <- gn2 + theme_void() + theme(legend.position="bottom", plot.margin=unit(c(1.2,1.3,1.2,1.3),"cm")) #+ expand_limits(x = c(1.2, 1.2), y = c(1.2, 1.2))
gn2 <- gn2 + scale_fill_manual(values = c("gray80", "#5C88DA", "#CC0C00", "#00B5E2", "gray30", "#FFCD00", "#00AF66"))
gn2 <- gn2 + scale_edge_color_manual(values = c("gray80", "#5C88DA", "#CC0C00", "#00B5E2", "gray30", "#FFCD00", "#00AF66"))
gn2 <- gn2 + theme(legend.key.height = grid::unit(1, "cm"),
                   legend.key.width = grid::unit(1, "cm"),
                   legend.text = element_text(size = 22),
                   legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=30, alpha = 1)))

# Save network
ggsave(filename = sprintf("%s/EDMFig_NetworkTimeAvr_v2.jpg", fig_output),
       plot = gn2, dpi = 300, width = 26, height = 25)
saveRDS(gn2, sprintf("%s/EDMFig_NetworkTimeAvr_v2.obj", fig_output))




#-------------------- Network figures for animation --------------------#
# Compile date
network_date <- c(ymd("2017-04-30"), edna_all$date[1:764])

# Use J1_matrix_all to illustrate interaction matrix
for(time_i in c(37:145)){
  J1_matrix_tmp <- matrix(0, ncol = nrow(edna_tax4), nrow = nrow(edna_tax4))
  colnames(J1_matrix_tmp) <- rownames(J1_matrix_tmp) <- rownames(edna_tax4)
  
  J1_tmp <- compiled_smap_all[[time_i]]$J1_matrix_all[[1]]
  diag(J1_tmp) <- 0
  J1_matrix_tmp[rownames(J1_tmp), colnames(J1_tmp)] <- J1_tmp
  J1_matrix_posi <- J1_matrix_nega <- J1_matrix_tmp
  J1_matrix_posi[J1_matrix_posi < 0] <- 0
  J1_matrix_nega[J1_matrix_nega > 0] <- 0
  
  # Create edge list of positive interactions
  eco_connect_posi <- as.data.frame(get.edgelist(graph.adjacency(J1_matrix_posi, weighted = TRUE, diag = FALSE, mode = "directed")))
  colnames(eco_connect_posi) <- c("from", "to")
  eco_connect_posi$value <- c(t(J1_matrix_posi[J1_matrix_posi != 0]))
  
  eco_connect_nega <- as.data.frame(get.edgelist(graph.adjacency(J1_matrix_nega, weighted = TRUE, diag = FALSE, mode = "directed")))
  colnames(eco_connect_nega) <- c("from", "to")
  eco_connect_nega$value <- c(t(J1_matrix_nega[J1_matrix_nega != 0]))

  # Create a plot
  posi_from <- match(eco_connect_posi$from, taxa_vertices$name)
  posi_to <- match(eco_connect_posi$to, taxa_vertices$name)
  posi_stength <- abs(eco_connect_posi$value)
  
  nega_from <- match(eco_connect_nega$from, taxa_vertices$name)
  nega_to <- match(eco_connect_nega$to, taxa_vertices$name)
  nega_stength <- abs(eco_connect_nega$value)
  
  # Basic usual argument
  gn3 <- gn1 + geom_node_point(aes(filter = leaf,
                                   x = x*1.03, y = y*1.03,
                                   size = std_abundance,
                                   fill = group),
                               color = "gray10",
                               shape = 21,
                               alpha = 0.5)
  gn3 <- gn3 + geom_conn_bundle(data = get_con(from = posi_from, to = posi_to),
                                edge_alpha = 0.4,
                                width = 0.4,
                                color = "royalblue",
                                arrow = grid::arrow(angle = 20, length = unit(0.2, "inches"), type = "closed"))
  gn3 <- gn3 + geom_conn_bundle(data = get_con(from = nega_from, to = nega_to),
                                edge_alpha = 0.2,
                                width = 0.4,
                                color = "red3",
                                arrow = grid::arrow(angle = 20, length = unit(0.2, "inches"), type = "closed"))
  gn3 <- gn3 + theme_void() + theme(legend.position="none") #+ expand_limits(x = c(1.2, 1.2), y = c(1.2, 1.2))
  gn3 <- gn3 + scale_fill_manual(values = c("gray80", "#5C88DA", "#CC0C00", "#00B5E2", "gray30", "#FFCD00", "#00AF66"))
  gn3 <- gn3 + scale_edge_color_manual(values = c("gray80", "#5C88DA", "#CC0C00", "#00B5E2", "gray30", "#FFCD00", "#00AF66"))
  gn3 <- gn3 + theme(legend.key.height = grid::unit(1, "cm"),
                     legend.key.width = grid::unit(1, "cm"),
                     legend.text = element_text(size = 22),
                     legend.title = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=30, alpha = 1)))
  if(time_i < 150){
    gn3 <- gn3 + ggtitle(paste("Plot1:", network_date[time_i])) + theme(plot.title = element_text(size=42))
  }else{
    gn3 <- gn3 + ggtitle(paste("Plot2:", network_date[time_i])) + theme(plot.title = element_text(size=42))
  }
  
  # Save network
  ggsave(sprintf("%s/Network_%s.jpg", network_output_folder, time_i),
           plot = gn3, width = 26, height = 25, dpi = 100)
    
  cat(sprintf("Network t = %s drawn and saved\n\n", time_i))
}

