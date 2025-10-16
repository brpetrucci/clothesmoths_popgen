###
## PCA and Adonis
###

###
# install packages
library(dplyr)
library(vegan)
library(ggplot2)
library(ggdendro)
library(pals)

###
# read data

# base directories
base_dir <- "/Users/petrucci/Documents/research/clothesmoths_popgen/"
pca_dir <- paste0(base_dir, "pca_adonis/")

# read metadata tables
bams_12p_loc <- read.delim(paste0(base_dir, "metadata/12p_loc.tsv"),
                           sep = "\t", header = TRUE)
bams_50p_loc <- read.delim(paste0(base_dir, "metadata/50p_loc.tsv"),
                           sep = "\t", header = TRUE)
bams_mass_loc <- read.delim(paste0(base_dir, "metadata/mass_loc.tsv"),
                            sep = "\t", header = TRUE)

# get donors list
donors_12p <- read.delim(paste0(base_dir, "metadata/donors_12p.tsv"))
donors_50p <- read.delim(paste0(base_dir, "metadata/donors_50p.tsv"))
donors_mass <- read.delim(paste0(base_dir, "metadata/donors_mass.tsv"))

# make columns for site_code
bams_12p_loc$site_code <- unlist(lapply(seq_along(bams_12p_loc$Donor), 
                                        function(i) {
                     donors_12p$site_code[donors_12p$site_name == bams_12p_loc$Donor[i]]
                                        }))
bams_50p_loc$site_code <- unlist(lapply(seq_along(bams_50p_loc$Donor), 
                                     function(i) {
                     donors_50p$site_code[donors_50p$site_name == bams_50p_loc$Donor[i]]
                                     }))
bams_mass_loc$site_code <- unlist(lapply(seq_along(bams_mass_loc$Donor), 
                                     function(i) {
                      donors_mass$site_code[donors_mass$site_name == bams_mass_loc$Donor[i]]
                                     }))

# make them factors (for better alphabetical legends)
bams_12p_loc$site_code <- as.factor(bams_12p_loc$site_code)
bams_50p_loc$site_code <- as.factor(bams_50p_loc$site_code)
bams_mass_loc$site_code <- as.factor(bams_mass_loc$site_code)

# read ibs matrices
ibs_12p <- as.matrix(read.table(paste0(pca_dir, "data/12p_ibs.ibsMat")))
ibs_50p <- as.matrix(read.table(paste0(pca_dir, "data/50p_ibs.ibsMat")))
ibs_mass <- as.matrix(read.table(paste0(pca_dir, "data/mass_ibs.ibsMat")))

# name columns and rows
rownames(ibs_12p) <- colnames(ibs_12p) <- bams_12p_loc$bam
rownames(ibs_50p) <- colnames(ibs_50p) <- bams_50p_loc$bam
rownames(ibs_mass) <- rownames(ibs_mass) <- bams_mass_loc$bam

# get colors 
donor_colors <- donors_12p$color
names(donor_colors) <- donors_12p$site_code

###
# make dendrograms

# function to make dendogram based on one thing at tips
dendo_id <- function(ibs, meta, pop_id, plot_name) {
  # make hclust
  hc <- hclust(as.dist(ibs), "ave")
  
  # make it a dendrogram
  dend <- as.dendrogram(hc)
  
  # make metadata rownames be the bams
  rownames(meta) <- meta[, 1]
  meta <- meta[, -1]
  
  # iterate through pop_id
  for (i in 1:length(pop_id)) {
    # prepping dendrogram to color it
    tree_labels <- dendro_data(dend, type = "rectangle")
    tree_labels$labels <- cbind(tree_labels$labels, 
                                col_by = meta[, pop_id[i]])
    
    tree_labels$segments$y[tree_labels$segments$yend == 0] <- 
      tree_labels$segments$y[tree_labels$segments$yend == 0]
    
    # name for legend
    legend_name <- ifelse(pop_id[i] == "site_code", "Collection site",
                          pop_id[i])
    
    # plot it with colored branches, no labels
    dend_plot <- 
      ggplot() +
      geom_segment(data = segment(tree_labels), aes(x = x, y = y, 
                                                    xend = xend, yend = yend)) +
      geom_segment(data = tree_labels$segments %>%
                     filter(yend == 0) %>%
                     left_join(tree_labels$labels, by = "x"), 
                   aes(x = x, y = y.x, xend = xend, yend = yend, 
                       color = col_by)) +
      scale_color_manual(name = legend_name, values = donor_colors) +
      guides(color = guide_legend(ncol = 1)) +
      coord_flip() +
      scale_y_reverse() +
      theme_dendro() +
      ggtitle(paste0(plot_name, " dendrogram clustered by ", legend_name))
    
    print(dend_plot)
  }
}

# make dendrogram for each
pdf(file = paste0(base_dir, "pca_adonis/plots/50p_dendro.pdf"),
    width = 10, height = 10)
dendo_id(ibs_50p, bams_50p_loc, c("site_code"), "Global")
dev.off()
pdf(file = paste0(base_dir, "pca_adonis/plots/mass_dendro.pdf"),
    width = 10, height = 10)
dendo_id(ibs_mass, bams_mass_loc, c("site_code"), "MA")
dev.off()

###
# make PCAs

# function to plot PCAs
make_pcas <- function(ibs, meta, pop_id, axes_to_plot, colors_vec, plot_name) {
  # get principal components
  pp <- capscale(ibs ~ 1)
  
  # get eigenvalues
  eigs <- as.data.frame(pp$CA$eig)
  
  # make rownames a column
  eigs$MDS <- rownames(eigs)
  
  # get variance explained by each PCA vector
  eig_var <- eigs[, 1] / sum(eigs[, 1])
  
  # first two pca axes
  pca_s <- as.data.frame(pp$CA$u[, axes_to_plot])
  
  # iterate through pop_id
  for (i in 1:length(pop_id)) {
    # get colors
    colors_pop <- colors_vec
    
    # get axes names
    x_axis <- colnames(pca_s[1])
    y_axis <- colnames(pca_s[2])
    
    # legend name
    legend_name <- ifelse(pop_id[i] == "site_code", "Collection site",
                          pop_id[i[]])
    
    # pca plot
    pca_plot <- ggplot(pca_s, aes_string(x = x_axis, y = y_axis)) +
      theme_bw() +
      geom_point(aes(colour = meta[, pop_id[i]]), size = 2, stroke = 1)  +
      stat_ellipse(type = "t",
                   aes(color = meta[, pop_id[i]]), show.legend = NA, lwd = 1) + 
      scale_color_manual(name = legend_name, values = colors_pop) +
      xlab(paste0(rownames(eigs)[axes_to_plot[1]], " (", 
                  round(eig_var[axes_to_plot[1]] * 100, 2), 
                  "% explained variance)")) +
      ylab(paste0(rownames(eigs)[axes_to_plot[2]], " (", 
                  round(eig_var[axes_to_plot[2]] * 100, 2), 
                  "% explained variance)")) +
      guides(color = guide_legend(ncol = 1)) + 
      ggtitle(paste0(plot_name, " PCA colored by ", legend_name))
    
    # print it
    print(pca_plot)
  }
}

# draw PCAs - axes 1-2
pdf(file = paste0(base_dir, "pca_adonis/plots/12p_pca.pdf"),
    width = 10, height = 7)
make_pcas(ibs_12p, bams_12p_loc, c("site_code"), c(1, 2), donor_colors, "Global (12p)")
dev.off()
pdf(file = paste0(base_dir, "pca_adonis/plots/50p_pca.pdf"),
    width = 10, height = 7)
make_pcas(ibs_50p, bams_50p_loc, c("site_code"), c(1, 2), donor_colors, "Global (50p)")
dev.off()
pdf(file = paste0(base_dir, "pca_adonis/plots/mass_pca.pdf"),
    width = 10, height = 7)
make_pcas(ibs_mass, bams_mass_loc, c("site_code"), c(1, 2), donor_colors, "MA")
dev.off()

# draw PCAs - axes 3-4
pdf(file = paste0(base_dir, "pca_adonis/plots/12p_pca_34.pdf"),
    width = 10, height = 7)
make_pcas(ibs_12p, bams_12p_loc, c("site_code"), c(3, 4), donor_colors, "Global (12p)")
dev.off()
pdf(file = paste0(base_dir, "pca_adonis/plots/50p_pca_34.pdf"),
    width = 10, height = 7)
make_pcas(ibs_50p, bams_50p_loc, c("site_code"), c(3, 4), donor_colors, "Global (50p)")
dev.off()
pdf(file = paste0(base_dir, "pca_adonis/plots/mass_pca_34.pdf"),
    width = 10, height = 7)
make_pcas(ibs_mass, bams_mass_loc, c("site_code"), c(3, 4), donor_colors, "MA")
dev.off()

# draw PCAs - axes 5-6
pdf(file = paste0(base_dir, "pca_adonis/plots/12p_pca_56.pdf"),
    width = 10, height = 7)
make_pcas(ibs_12p, bams_12p_loc, c("site_code"), c(5, 6), donor_colors, "Global (12p)")
dev.off()
pdf(file = paste0(base_dir, "pca_adonis/plots/50p_pca_56.pdf"),
    width = 10, height = 7)
make_pcas(ibs_50p, bams_50p_loc, c("site_code"), c(5, 6), donor_colors, "Global (50p)")
dev.off()
pdf(file = paste0(base_dir, "pca_adonis/plots/mass_pca_56.pdf"),
    width = 10, height = 7)
make_pcas(ibs_mass, bams_mass_loc, c("site_code"), c(5, 6), donor_colors, "MA")
dev.off()

###
# adonis analyses

adonis2(ibs_50p ~ site_code, bams_50p_loc)
adonis2(ibs_mass ~ site_code, bams_mass_loc)
