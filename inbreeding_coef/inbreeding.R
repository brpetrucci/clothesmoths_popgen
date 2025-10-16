##########################################################
## Calculating inbreeding coefficients for moth dataset ##
##########################################################

###
# packages

library(dplyr)
library(hierfstat)
library(vegan)
library(geosphere)

###
# read and manipulate data

# base directory
base_dir <- "/Users/petrucci/Documents/research/clothesmoths_popgen/"

# all bams at 50p
data_50p <- read.delim(paste0(base_dir, "inbreeding_coef/data/bams_50p.gz"), 
                       header = FALSE)
data_50p <- data_50p[, -ncol(data_50p)]
# removing NA at the end (whitespace)

# location information
bams_50p_loc <- read.delim(paste0(base_dir, "metadata/50p_loc.tsv"),
                           header = TRUE, sep = "\t")

# get donors lists
donors_50p <- read.delim(paste0(base_dir, "metadata/donors_50p.tsv"),
                         header = TRUE, sep = "\t")

# make names for data
colnames(data_50p) <- c("chr", "locus", 
                        paste(sort(rep(bams_50p_loc$bam, 3)), 
                              c("-0", "-1", "-2"), sep = ""))

# make loci name data frame for each dataset
loci <- data.frame(chr = data_50p$chr, locus_old = data_50p$locus,
                   chr_locus = paste0(data_50p$chr, "_", data_50p$locus))

# select rows for which the chr_locus row are unique
loci <- loci[!(duplicated(loci$chr_locus)), ]

# add row with new names for loci
loci$locus <- paste0("locus_", 1:nrow(loci))

# add locus number to loci data frame
loci$locus <- as.vector(loci %>%
  dplyr::filter(chr_locus %in% loci$chr_locus) %>%
  mutate(chr_locus = factor(chr_locus, levels = loci$chr_locus)) %>%
  arrange(chr_locus) %>%
  select(locus))[[1]]

# create final genotype data frame
gen_50p <- data_50p[, 3:ncol(data_50p)]

# row names based on loci
rownames(gen_50p) <- loci$locus

###
# functions to run the calculations

# calculating F for one locus
f_loc <- function(gen_locus) {
  # exclude all the NA values
  locus <- gen_locus[gen_locus != -1]
  
  # calculate p and q
  p <- (2 * sum(locus == 0) + sum(locus == 1)) / (2 * length(locus))
  q <- (2 * sum(locus == 2) + sum(locus == 1)) / (2 * length(locus))
  
  # calculate expected heterozygosity
  exp_h <- 2 * p * q
  
  # and observed heterozygosity
  obs_h <- sum(locus == 1) / length(locus)
  
  # calculate inbreeding coefficient
  # if p or q is 0, that's an obs_h of 0, so just return 1
  f <- if (exp_h == 0) 1 else 1 - obs_h / exp_h
}

# calculating F for a bunch of loci
f_total <- function(gen) {
  # check which rows have all missing data
  miss <- which(rowSums(gen == -1) == ncol(gen))
  
  # take those gens away
  if (length(miss) > 0) gen <- gen[-miss, ]
  
  # check which rows are all one allele
  all_one <- which(rowSums(gen == 1) == ncol(gen) |
                     rowSums(gen == 2) == ncol(gen))
  
  # take those away too
  if (length(all_one) > 0) gen <- gen[-all_one, ]
  
  # make final vector
  fs <- vector("numeric", nrow(gen))
  
  # iterate through loci
  for (i in 1:nrow(gen)) {
    # calculate f for this locus and add it to res
    fs[i] <- f_loc(gen[i, ])
  }
  
  # return average f
  sum(fs) / nrow(gen)
}

# calculating F for a number of populations
f_pops <- function(gen, locs_all, pop_id, threshold = 5) {
  # choose only the locs columns for the bam number and pop_id (e.g. Country)
  locs <- locs_all[, c(1, which(colnames(locs_all) == pop_id))]
  
  # start fs data frame
  fs <- data.frame(matrix(nrow = length(unique(locs[, pop_id])), ncol = 2))
  
  # iterate through populations
  for (i in 1:length(unique(locs[, pop_id]))) {
    # which pop_id is this
    pop_i <- unique(locs[, pop_id])[i]
    
    # get locs for that pop
    locs_pop <- locs[locs[, pop_id] == pop_i, ]
    
    # if there are less than threshold, skip this pop
    if (nrow(locs_pop) < threshold) next
    
    # select just the gen for that population
    gen_pop <- gen[, locs_pop$bam]
    
    # calculate the f for that pop
    fs[i, ] <- c(pop_i, f_total(gen_pop))
  }
  
  # remove NA
  fs <- fs[!is.na(fs[, 2]), ]
  
  # name fs
  colnames(fs) <- c(pop_id, "F")
  
  # return fs
  fs
}

# function to make a genotype matrix for basic.stats
make_new_gen <- function(gen, locs_all, pop_id, threshold = 5) {
  
  # changing genotype data frames to have 11, 12, 22, and NA
  gen_new <- gen
  gen_new[gen_new == 0] = 11
  gen_new[gen_new == 1] = 12
  gen_new[gen_new == 2] = 22
  gen_new[gen_new == -1] = NA
  
  # adding populations to genotype data frames
  gen_pop <- cbind(data.frame(pop =
                                locs_all[locs_all$bam %in% colnames(gen_new), ][, pop_id],
                              t(gen_new)))
  
  # check which populations have less than threshold
  nsamps_pop <- table(gen_pop$pop)
  pops_remove <- names(nsamps_pop)[which(nsamps_pop < threshold)]
  
  # remove those from gen_pop
  gen_pop <- gen_pop[!(gen_pop$pop %in% pops_remove), ]
  
  # return gen_pop and pops_remove
  return(list(GEN = gen_pop, POPS = pops_remove))
}

# make function to take a soft called genotype 
# and draw to make a hard called one
soft_to_hard <- function(gen) {
  # number of bams
  n_bams <- ncol(gen) / 3
  
  # make return data frame
  res <- data.frame(matrix(nrow = nrow(gen), ncol = n_bams))
  
  # make bam names vector
  bam_names <- c()
  
  # iterate through bams
  for (i in 1:n_bams) {
    # get the probabilities for 0, 1 and 2
    vals <- gen[, ((i - 1) * 3 + 1):(i * 3)]
    
    # sample from probabilities
    hard_gens <- unlist(lapply(1:nrow(vals), function(i) {
      sample(c(0, 1, 2), 1, prob = vals[i, ])
    }))
    
    # add that as a column to res
    res[, i] <- hard_gens
    
    # get bam name
    bam_names <- c(bam_names, sub("-.*", "", colnames(gen)[((i - 1) * 3 + 1)]))
  }
  
  # name bams
  colnames(res) <- bam_names
  
  # return res
  return(res)
}

# do a number of samplings of soft called genotypes
f_sampled <- function(gen, locs_all, pop_id,
                      threshold = 5, samples = 1000) {
  # get number of pop_id with more than threshold samples
  n_pops <- sum(table(locs_all[, pop_id]) >= threshold)
  
  # make fit and fstats data frames
  fit <- data.frame(matrix(nrow = samples, ncol = n_pops))
  fstats <- data.frame(matrix(nrow = samples, ncol = 10))
  
  # iterate through samples
  for (i in 1:samples) {

    # print sample
    print(paste0("Sample: ", i))
    
    # get a hard genotype matrix
    gen_hard <- soft_to_hard(gen)
    
    # get a hierfstate-friendly version
    gen_pop <- make_new_gen(gen_hard, locs_all, pop_id, threshold)$GEN
    
    # get fit
    fit_samp <- f_pops(gen_hard, locs_all, pop_id, threshold)
    fit[i, ] <- fit_samp[, 2]
    
    # and fstats
    fstats_samp <- basic.stats(gen_pop)$overall
    fstats[i, ] <- unname(fstats_samp)
  }
  
  # name columns
  colnames(fit) <- fit_samp[, 1]
  colnames(fstats) <- names(fstats_samp)
  
  # name rows
  rownames(fit) <- rownames(fstats) <- paste0("sample_", 1:samples)
  
  # return final data frames
  return(list(FIT = fit, FSTATS = fstats))
}

# get permutations for 50p
f_50p <- f_sampled(gen_50p, bams_50p_loc, "Donor", samples = 1000)
fit_50p <- f_50p$FIT
fstats_50p <- f_50p$FSTATS

# save them in files
write.table(fit_50p, paste0(base_dir, "inbreeding_coef/fit_50p.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(fstats_50p, paste0(base_dir, "inbreeding_coef/fstats_50p.tsv"), 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

###
# fixation index

# function to build fst matrix similar to those above
fst_matrix <- function(folder, donors, weighted = TRUE) {
  # make results matrix
  res <- data.frame(matrix(nrow = 0, ncol = 3))
  
  # iterate through donors twice
  for (i in 1:length(donors)) {
    # get first donor
    d1 <- donors[i]
    
    for (j in 1:length(donors)) {
      # ensure we get one unique i-j pairing
      if (j < i) next
      
      # get second donor
      d2 <- donors[j]
      
      if (j == i) {
        res <- rbind(res, c(d2, d1, NA))
        
        next
      }
      
      # get file for this comparison
      fst_file <- paste0(base_dir, "inbreeding_coef/data/",
                         folder, "p", d1, "_", d2, ".fst.stats")
      
      # check if it exists, if not do the opposite
      if (!file.exists(fst_file)) {
        fst_file <- paste0(base_dir, "inbreeding_coef/data/",
                           folder, "p", d2, "_", d1, ".fst.stats")
      }
      
      # get numbers
      fst <- unname(read.delim(fst_file, header = FALSE))
      
      # separate weighted and unweights fst
      fst_uw <- fst[1]
      fst_w <- fst[2]
      
      # add to res
      res <- rbind(res, c(d2, d1, fst_w))
    }
  }
  
  # name columns
  colnames(res) <- c("row", "col", "value")
  
  # make rows and columns factors
  res$row <- factor(res$row, levels = sort(unique(res$row)))
  res$col <- factor(res$col, levels = sort(unique(res$col)))
  
  # return res
  return(res)
}

# donors we want
donors_50p_kept <- donors_50p[donors_50p$n_moths > 4, ]
donors_50p_n <- as.numeric(sub('([a-zA-z]+)', '', donors_50p_kept$donor_name))
donors_mass_n <- c(9, 10, 12, 14)

# get data frames for each
fst_50p <- fst_matrix("fst_global/", donors_50p_n)
fst_mass <- fst_matrix("fst_mass/", donors_mass_n)

# plot heatmaps
ggplot(fst_50p, aes(x = row, y = col, fill = value)) +
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
    fill = "Fst",
    title = "Pairwise fixation index (Global)"
  )
ggplot(fst_mass, aes(x = row, y = col, fill = value)) +
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
    fill = "Fst",
    title = "Pairwise fixation index (MA)"
  )

# make symmetric fst matrix
symmetric_fst <- function(fst) {
  # make a matrix of the desired dimensions
  sym_fst <- as.data.frame(matrix(nrow = length(unique(fst$row)),
                                  ncol = length(unique(fst$col))))
  
  # iterate through rows and columns
  for (j in 1:length(unique(fst$col))) {
    for (i in j:length(unique(fst$row))) {
      # get row and column values
      c <- unique(fst$col)[j]
      r <- unique(fst$row)[i]
      
      # add to sym matrix
      sym_fst[i, j] <- fst[fst$row == r & fst$col == c, ]$value
      sym_fst[j, i] <- fst[fst$row == r & fst$col == c, ]$value
    }
  }
  
  # set NAs to 0
  sym_fst[is.na(sym_fst)] <- 0
  
  # name columns and rows
  colnames(sym_fst) <- rownames(sym_fst) <- paste0("donor", unique(fst$row))
  
  # return
  return(sym_fst)
}

# get symmetric fst
sym_fst_50p <- symmetric_fst(fst_50p)
sym_fst_mass <- symmetric_fst(fst_mass)

# rename matrices based on site code
colnames(sym_fst_50p) <- rownames(sym_fst_50p) <- 
  donors_50p_kept$site_code[donors_50p_kept$donor_name %in% colnames(sym_fst_50p)]
colnames(sym_fst_mass) <- rownames(sym_fst_mass) <- 
  donors_50p_kept$site_code[donors_50p_kept$donor_name %in% colnames(sym_fst_mass)]

# make it a dist object
fst_dist_50p <- as.dist(sym_fst_50p)
fst_dist_mass <- as.dist(sym_fst_mass)

# get latitude and longitude data frame
lat_long <- read.delim(paste0(base_dir, "inbreeding_coef/data/distances.csv"), 
                       sep = ",")[, c("site_code", "Latitude", "Longitude")]

# select only sites we're keeping
lat_long <- lat_long[lat_long$site_code %in% colnames(sym_fst_50p), ]

# make it into distances
geo_dist_50p <- distm(cbind(lat_long$Longitude, lat_long$Latitude), 
                      fun = distHaversine)
rownames(geo_dist_50p) <- colnames(geo_dist_50p) <- colnames(sym_fst_50p)

# make geo dist for mass
geo_dist_mass <- geo_dist_50p[4:7, 4:7]

# make them into dists
geo_dist_50p <- as.dist(geo_dist_50p)
geo_dist_mass <- as.dist(geo_dist_mass)

# run Mantel test
fst_mantel_50p <- mantel(fst_dist_50p, geo_dist_50p, method = "spearman",
                         permutations = 10000, na.rm = TRUE)
fst_mantel_mass <- mantel(fst_dist_mass, geo_dist_mass, method = "spearman",
                          permutations = 10000, na.rm = TRUE)
