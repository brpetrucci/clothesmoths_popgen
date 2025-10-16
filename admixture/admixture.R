###
## Admixture
###

###
# packages

# dplyr
library(dplyr)

# pals
library(pals)

# ggplot
library(ggplot2)

###
# read data

# base directory
base_dir <- "/Users/petrucci/Documents/research/clothesmoths_popgen/"

# metadata
bams_50p_loc <- read.delim(paste0(base_dir, "metadata/50p_loc.tsv"),
                           sep = "\t", header = TRUE)
bams_mass_loc <- read.delim(paste0(base_dir, "metadata/mass_loc.tsv"),
                            sep = "\t", header = TRUE)

# get donors list
donors_50p <- read.delim(paste0(base_dir, "metadata/donors_50p.tsv"))
donors_mass <- read.delim(paste0(base_dir, "metadata/donors_mass.tsv"))

# get colors 
donor_colors <- donors_50p$colors
names(donor_colors) <- donors_50p$site_code

###
# basic auxiliary functions

# finding the mode of each row
rowModes <- function(df) {
  # create results vector
  res <- vector("numeric", nrow(df))
  
  # iterate through rows
  for (i in 1:nrow(df)) {
    # get the table for the row
    tbl <- table(unlist(df[i, ]))
    
    # get the value for which there are more
    mode <- as.numeric(names(tbl)[which.max(tbl)])
    
    # add to results
    res[i] <- mode
  }
  
  # return res
  return(res)
}

###
# functions to find optimal K

# function to get the likelihoods for a dataset and the total ks we did
get_liks_k <- function(dataset, ks, nreps = 100) {
  # make a results vector
  liks <- data.frame(matrix(nrow = 0, ncol = 2))
  
  # iterate through ks
  for (i in 2:ks) {
    # iterate through reps
    for (rep in 1:nreps) {
      # read the log
      log <- readLines(paste0(base_dir, "admixture/", 
                              dataset, "/ngsadmix_k", i, ".", 
                              rep, ".log"))
      
      # get line with best likelihood
      idx_lik <- grep("best like", log)
      
      # get likelihood from it
      lik <- as.double(sub('best like=([-0-9.]+).*', '\\1', log[idx_lik]))
      
      # add it to liks
      liks <- rbind(liks, c(i, abs(lik)))
    }
  }
  
  # name columns
  colnames(liks) <- c("K", "lik")
  
  # return likelihoods
  return(liks)
}

# function to get likelihoods for a value of k
liks_k <- function(liks, k) {
  # return just the likelihoods for that value of k
  return(liks$lik[liks$K == k])
}

# get the deltaK values from Evanno et al 2005
deltaK <- function(liks) {
  # get number of ks
  ks <- unique(liks$K)
  
  # make return vector
  deltaK <- data.frame(matrix(nrow = 0, ncol = 5))
  
  # iterate through k
  for (k in ks[-c(1, length(ks))]) {
    # get mean of |L''(k)|
    mlk <- mean(abs(liks_k(liks, k + 1) - 2 * liks_k(liks, k) + 
                      liks_k(liks, k - 1)))
    
    # get stdev of L(K)
    sdlk <- sd(liks_k(liks, k))
    
    # highest L(K) rep
    max_rep <- which.max(liks_k(liks, k))
    
    # get deltaK
    deltaK <- rbind(deltaK, c(k, mlk, sdlk, mlk / sdlk, max_rep))
  }
  
  # name it
  colnames(deltaK) <- c("K", "mean_abs_lkpp", "sd_lk", "deltaK", "max_rep")
  return(deltaK)
}

# get likelihoods for each dataset
liks_50p <- get_liks_k("50p", length(unique(bams_50p_loc$Donor)) + 1)
liks_mass <- get_liks_k("mass", length(unique(bams_mass_loc$Donor)) + 1)

# get deltaK
deltaK_50p <- deltaK(liks_50p)
deltaK_mass <- deltaK(liks_mass)

###
# getting the lineages for each bam

# function to do the thing written two lines above this
get_lins <- function(dataset, meta, k, nreps = 100) {
  # make a data frame to hold all the possible lineages for each bam
  bam_lins <- data.frame(matrix(nrow = nrow(meta), ncol = nreps))
  
  # iterate through reps
  for (rep in 1:nreps) {
    # get the qopt file
    qopt <- read.table(paste0(base_dir, "admixture/", dataset, 
                              "/ngsadmix_k", k, ".", rep, ".qopt"))
    
    # get lineage assignments for this rep
    lins <- unlist(unname(lapply(1:nrow(qopt), 
                                 function(x) which.max(qopt[x, ]))))
    
    # copy it
    lins_reassigned <- lins
    
    # iterate through the lineages to ensure we have the same order for all reps
    for (i in 1:length(unique(lins))) {
      # make lin unique(lins)[i] into i, so as to make sure we always get
      # lineages in increasing order of appearance
      lins_reassigned[lins == unique(lins)[i]] <- i
    }
    
    # for each bam, find which lineage it is assigned to with highest probability
    bam_lins[, rep] <- lins_reassigned
  }
  
  # name bam_lins
  colnames(bam_lins) <- paste0("rep_", 1:nreps)
  rownames(bam_lins) <- meta$bam
  
  # return lineages
  return(bam_lins)
}

# function to get the rep closest to the mode, and the highest likelihood one
get_lin_mode <- function(lins, deltaK, k) {
  # get mode
  lin_mode <- rowModes(lins)
  
  # get lineages for the highest likelihood rep
  lin_max <- lins[, deltaK$max_rep[k - 2]]
  
  # get lineage that's closest to mode
  lin_closest_idx <- min(which.max(unlist(lapply(1:ncol(lins), 
                               function(x) sum(lins[, x] == lin_mode)))))
  lin_closest <- lins[, lin_closest_idx]
  
  # put them all in one data frame
  res <- cbind(lin_mode, lin_max, lin_closest)
  
  # name it
  colnames(res) <- c("mode", "max", "closest")
  
  # return it
  return(res)
}

###
# functions to plot admixture

# source plotAdmixture function
source(paste0(base_dir, "admixture/plot_admixture_function.R"))

# function to calculate the mean qopt for all the reps in an ngsadmix run
mean_qopt <- function(dataset, k, reps = 100) {
  # list for all qopt data frames
  qopt_df <- list()
  
  # iterate through reps
  for (rep in 1:reps) {
    # get the qopt file
    qopt <- read.table(paste0(base_dir, "admixture/", dataset, 
                              "/ngsadmix_k", k, ".", rep, ".qopt"))
    
    # get lineage assignments for this rep
    lins <- unlist(unname(lapply(1:nrow(qopt), 
                                 function(x) which.max(qopt[x, ]))))
    
    # reorder qopt for the order of lineages assigned to bams the most
    qopt_df[[rep]] <- qopt[, unique(lins)]
    
    # rename columns
    colnames(qopt_df[[rep]]) <- paste0("lin_", 1:k)
  }
  
  # get an average qopt file
  qopt_mean <- Reduce("+", qopt_df) / length(qopt_df)
  
  # return it
  return(qopt_mean)
}

# function to plot one admixture graph given a dataset and k
plot_admix_k <- function(dataset, meta, k, filename = NULL,
                         rep = 1, mean = FALSE, reps = 100) {
  # if mean is true, get mean qopt
  if (mean) {
    qopt <- mean_qopt(dataset, k, reps)
  } else {
    # if not, get it for rep
    qopt <- read.table(paste0(base_dir, "admixture/", 
                              dataset, "/ngsadmix_k", k, ".", rep, ".qopt"))
    
    # get lineage assignments for this rep
    lins <- unlist(unname(lapply(1:nrow(qopt), 
                                 function(x) which.max(qopt[x, ]))))
    
    # reorder qopt for the order of lineages assigned to bams the most
    qopt <- qopt[, unique(lins)]
    
    colnames(qopt) <- paste0("lin_", 1:ncol(qopt))
  }
  
  # get a number for each donor
  donors <- data.frame(old = sort(unique(meta$Donor)),
                       new = 1:length(unique(meta$Donor)))
  
  # get donor number instead of names
  meta$DonorN <- unlist(lapply(seq_along(meta$Donor), function(i) {
    donors$new[donors$old == meta$Donor[i]]
  }))
  
  # final dataset is the qopt bound to the metadata
  admix <- cbind(qopt, pop = meta$DonorN)
  
  # name it for the bams and remove bams column
  rownames(admix) <- meta$bam
  
  # putting populations in correct order
  admix$pop <- as.factor(admix$pop)
  
  # open quartz
  #quartz()
  
  # get colors
  lin_colors <- rainbow(k)
  
  # open pdf
  if (!is.null(filename)) {
    pdf(file = filename, width = 10, height = 7)
  }
  
  # plot admixture
  ords <- plotAdmixture(data = admix, npops = k, cex = 1,
                        grouping.method = "distance", colors = lin_colors)
  
  # close connection
  if (!is.null(filename)) {
    dev.off()
  }
}

# plot admixture for 50p
lapply(2:22, function(i) 
  plot_admix_k("50p", bams_50p_loc, i,
               paste0(base_dir, "admixture/plots/admix_50p_", i, ".pdf"),
               mean = TRUE, reps = 100))
plot_admix_k("50p", bams_50p_loc, 8, mean = TRUE, reps = 100)
#plot_admix_k("50p", bams_50p_loc, 22, mean = TRUE, reps = 100)
#plot_admix_k("50p", bams_50p_loc, 5, mean = TRUE, reps = 100)

# plot admixture for mass
lapply(2:5, function(i) 
  plot_admix_k("mass", bams_mass_loc, i,
               paste0(base_dir, "admixture/plots/admix_mass_", i, ".pdf"),
               mean = TRUE, reps = 100))
#plot_admix_k("mass", bams_mass_loc, 3, mean = TRUE, reps = 100)
#plot_admix_k("mass", bams_mass_loc, 5, mean = TRUE, reps = 100)
