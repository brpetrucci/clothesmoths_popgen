########################
## Plotting NGSrelate ##
########################

###
# packages
library(plyr)
library(dplyr)
library(tidyverse)
library(vegan)
library(geosphere)

# base directory
base_dir <- "/Users/petrucci/Documents/research/clothesmoths_popgen/"

# this project's directory
ngs_dir <- paste0(base_dir, "ngsrelate/")

# get 50p and mass metadata table
bams_50p_loc <- read.delim(paste0(base_dir, "metadata/50p_loc.tsv"),
                           sep = "\t", header = TRUE)
bams_mass_loc <- read.delim(paste0(base_dir, "metadata/mass_loc.tsv"),
                            sep = "\t", header = TRUE)

# get donors list
donors_50p <- read.delim(paste0(base_dir, "metadata/donors_50p.tsv"))
donors_mass <- read.delim(paste0(base_dir, "metadata/donors_mass.tsv"))

# which donors have less than 5 samples
donors_remove_50p <- donors_50p$site_name[which(donors_50p$n_moths < 5)]
donors_remove_mass <- donors_mass$site_name[which(donors_mass$n_moths < 5)]

# find which bams to remove from bam locations
bams_remove_50p <- which(bams_50p_loc$Donor %in% donors_remove_50p)
bams_remove_mass <- which(bams_mass_loc$Donor %in% donors_remove_mass)

# remove them
bams_50p_loc <- bams_50p_loc[-bams_remove_50p, ]
bams_mass_loc <- bams_mass_loc[-bams_remove_mass, ]

# add merge id column
bams_50p_loc$merge_id <- 0:(nrow(bams_50p_loc) - 1)
bams_mass_loc$merge_id <- 0:(nrow(bams_mass_loc) - 1)

# levels for 50p, so we group same country and state together
levels_50p <- donors_50p$site_name

# make lineage column (each donor is a lineage)
bams_50p_loc$lineage <- as.integer(factor(bams_50p_loc$Donor, 
                                          levels = levels_50p))
bams_mass_loc$lineage <- bams_50p_loc$lineage[bams_50p_loc$Donor %in%
                                                bams_mass_loc$Donor]

# get lineage-donor names
donor_lin_50p <- bams_50p_loc[!(duplicated(bams_50p_loc$Donor)), 
                              c("Donor", "lineage")] %>% arrange(lineage)
donor_lin_mass <- bams_mass_loc[!(duplicated(bams_mass_loc$Donor)), 
                              c("Donor", "lineage")] %>% arrange(lineage)


# read in relatedness matrices
rel_50p <- read.table(paste0(ngs_dir, "data/newres_global_50p"), 
                      header = TRUE)[, c("a", "b", "rab")]
rel_mass <- read.table(paste0(ngs_dir, "data/newres_mass"), 
                       header = TRUE)[, c("a", "b", "rab")]

# function to remove all the necessary bams from a given rel
remove_bams <- function(rel, bams) {
  # make numerical list of bams
  bams_num <- as.numeric(gsub('([^0-9]+)', '', bams))
  
  # iterate through bams
  for (i in 1:length(bams_num)) {
    # bam number (0 indexed)
    bam <- bams_num[i] - 1
    
    # remove all mention of this bam
    rel <- rel[rel$a != bam & rel$b != bam, ]
    rel$a <- ifelse(rel$a > bam, rel$a - 1, rel$a)
    rel$b <- ifelse(rel$b > bam, rel$b - 1, rel$b)
    
    # subtract bams that are higher
    bams_num <- ifelse(bams_num > bam, bams_num - 1, bams_num)
  }
  
  # return final rel
  return(rel)
}

# remove bams
rel_50p <- remove_bams(rel_50p, bams_remove_50p)
rel_mass <- remove_bams(rel_mass, bams_remove_mass)

# function to return the data frame we want
final_rel <- function(df, metadata) {
  # merge with metadata
  df_merged <- left_join(df, metadata, by = c("a" = "merge_id"))
  
  # get just a, b, rab, and bam columns, naming bam as bam_a
  df_sel <- df_merged %>%
    select(a, b, rab, bam, lineage) %>%
    dplyr::rename(bam_a = bam) %>%
    dplyr::rename(lineage_a = lineage) %>%
    dplyr::mutate(bam_a = as.factor(bam_a))
  
  # merge with metadata again to get the sample_id_b
  df_sel2 <- left_join(df_sel, metadata, by = c("b" = "merge_id"))
  
  # now dataframe with sample_id_a and sample_id_b and relatedness
  df_final <- df_sel2 %>%
    select(bam_a, bam, rab, lineage_a, lineage) %>%
    dplyr::rename(bam_b = bam) %>%
    dplyr::rename(lineage_b = lineage) %>%
    dplyr::mutate(bam_b = as.factor(bam_b))
  
  # make lineage comparison column
  df_final$lineage_comparison <- as.factor(paste(df_final$lineage_a, 
                                                 df_final$lineage_b, sep = "_"))
  
  return(df_final)
}

# final relatedness matrices
rel_50p_final <- final_rel(rel_50p, bams_50p_loc)
rel_mass_final <- final_rel(rel_mass, bams_mass_loc)

# function to make heatmap based on relatedness data
rel_heatmap <- function(df, ggplot_heat = FALSE, threshold = 5) {
  # mean rab for each lineage comparison
  rel_means <- df %>%
    dplyr::group_by(lineage_a, lineage_b) %>%
    dplyr::summarise(mean_rab = mean(rab, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(lineage_a, lineage_b)
  
  # get all lineages
  all_lineages <- union(rel_means$lineage_a, rel_means$lineage_b)
  
  # make it symmetric
  sym_rel_means <- rel_means %>% 
    dplyr::bind_rows(
      rel_means %>%
        dplyr::rename(lineage_a = lineage_b, lineage_b = lineage_a)
    ) %>%
    dplyr::group_by(lineage_a, lineage_b) %>%
    dplyr::summarise(mean_rab = mean(mean_rab, na.rm = TRUE), .groups = "drop") %>%
    complete(lineage_a = all_lineages,
             lineage_b = all_lineages)
  
  # pivot into a matrix
  rel_mat <- sym_rel_means %>%
    pivot_wider(names_from = lineage_b, values_from = mean_rab) %>%
    dplyr::arrange(lineage_a)
  
  # convert to matrix, set rownames
  rel_mat_final <- as.matrix(rel_mat[, -1])
  rownames(rel_mat_final) <- rel_mat$lineage_a
  
  # remove rows and columns with na values
  na_rows <- apply(rel_mat_final, 1, function(x) any(is.na(x)))
  
  if (sum(na_rows) > 0) rel_mat_final <- rel_mat_final[!na_rows, !na_rows]
  
  # if we want a ggplot heatmap
  if (ggplot_heat) {
    # make it into a data frame
    df <- as.data.frame(rel_mat_final)
    df$row <- rownames(df)
    
    # Convert matrix to long format
    long_df <- df %>%
      pivot_longer(cols = -row,
                   names_to = "col",
                   values_to = "value") %>%
      filter(as.numeric(row) >= as.numeric(col)) %>%
      mutate(row = factor(row, levels = unique(c(as.character(df$row), 
                                                 as.character(df$col)))),
             col = factor(col, levels = unique(c(as.character(df$row), 
                                                 as.character(df$col)))))
    
    return(long_df)
  } else {
    return(rel_mat_final)
  }
}

# get a matrix for heatmap between lineages for mass
rel_50p_lins <- as.data.frame(rel_heatmap(rel_50p_final, ggplot_heat = TRUE))
rel_mass_lins <- as.data.frame(rel_heatmap(rel_mass_final, ggplot_heat = TRUE))

# log transform the values
rel_50p_lins$value <- log10(rel_50p_lins$value)
rel_mass_lins$value <- log10(rel_mass_lins$value)

# plot ggplot heatmap
ggplot(rel_50p_lins, aes(x = row, y = col, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.text.y = element_text(hjust = 1)
  ) +
  labs(
    x = NULL, y = NULL,
    fill = "Mean rab",
    title = "Mean log pairwise relatedness (Global)"
  )

ggplot(rel_mass_lins, aes(x = row, y = col, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.text.y = element_text(hjust = 1)
  ) +
  labs(
    x = NULL, y = NULL,
    fill = "Mean rab",
    title = "Mean log pairwise relatedness (MA)"
  )
