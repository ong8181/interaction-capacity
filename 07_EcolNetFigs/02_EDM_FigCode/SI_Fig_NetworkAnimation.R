####
#### CERrice2017 All data analysis
#### Fig. Network animation
####

# You need to re-run the analyses, or contact ong8181@gmail.com to obtain the jpg files!

# Load library and functions
library(ggplot2); packageVersion("ggplot2") # 3.2.1, 2020.1.30
library(network); packageVersion("network") # 1.16, 2020.1.30
library(ggnetwork); packageVersion("ggnetwork") # 0.5.6.9000, 2020.1.30
library(sna); packageVersion("sna") # 2.5, 2020.1.30
library(ggsci); packageVersion("ggsci") # 2.9, 2019.11.11

# Animation package
library(magick); packageVersion("magick") # 2.2, 2020.1.30

# Output directory
network_output_folder2 <- "01_NetworkFigs_v2"
fig_output <- "../00_RawFigs/02_Fig_EDMnet"
dir.create(fig_output)

# Creating animation version 2
img2 <- image_read(sprintf("%s/Network_37.jpg", network_output_folder2))
for(image_i in c(38:145)){
  image_file2 <- sprintf("%s/Network_%s.jpg", network_output_folder2, image_i)
  img_tmp2 <- image_read(image_file2)
  img2 <- c(img2, img_tmp2)
}

animation <- image_animate(image_scale(img2, "900x900"), fps = 10)
image_write(animation, sprintf("%s/Animation_InteractionNetwork_v2.gif", fig_output))
