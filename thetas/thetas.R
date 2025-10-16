####
## Calculating thetas and Tajimas D and allat
##

###
# packages

library(ggplot2)

###
# read data

# base directory
base_dir <- "/Users/petrucci/Documents/research/clothesmoths_popgen/"

# thetas directory
thetas_dir <- paste0(base_dir, "thetas/")

# location data per donor
meta <- read.delim(paste0(thetas_dir, "data/distances.csv"), sep = ",")

###
# extracting summary stats for each donor

# calculate weighted colMeans
weightedColMeans <- function(df, weights) {
  # results df
  res <- c()
  
  # iterate through columns
  for (i in 1:ncol(df)) {
    # get weighted mean
    res <- c(res, sum(df[, i] * weights) / sum(weights))
  }
  
  # name it
  names(res) <- colnames(df)
  
  # return res
  return(res)
}

# calculate stats for one donor
theta_stats_donor <- function(N, stats, window) {
  # get filename
  filename <- paste0("donor", N, ".thetas", 
                     ifelse(window, "Window.gz", ".idx"), ".pestPG")
  
  # read file
  thetas <- read.table(paste0(thetas_dir, "data/pestpg/", filename), 
                       header = FALSE)
  
  # name columns
  colnames(thetas) <- c("win_indices", "Chr", "WinCenter", 
                        "tW", "tP", "tF", "tH", "tL", "Tajima", 
                        "fuf", "fud", "fayh", "zeng", "nSites")
  
  # make extra column
  thetas$calc_Tajima <- (thetas$tP - thetas$tW) / sqrt(var(thetas$tP - thetas$tW))
  
  # filter out windows with no sites, if any
  thetas <- thetas[thetas$nSites > 0, ]
  
  # get summary statistics for columns we're interested in
  stat_cols <- weightedColMeans(thetas[, stats], thetas$nSites)
  
  # return means for each stats
  return(stat_cols)
}

# calculate stats for all donors
theta_stats <- function(Ns, stats, window) {
  # start results data frame
  res <- data.frame(matrix(nrow = 0, ncol = length(stats)))
  
  # iterate through Ns
  for (n in Ns) {
    # get theta stats for this donor
    stats_donor <- theta_stats_donor(n, stats, window)
    
    # add it to res
    res <- rbind(res, stats_donor)
  }
  
  # name columns
  colnames(res) <- stats
  
  # return res
  return(res)
}

# Ns (donors to try)
donor_ns <- as.numeric(gsub("([a-zA-Z]+)", "", 
                            meta[meta$n_moths > 4, ]$donor_name))

# stats
stats_to_do <- c("tW", "tP", "Tajima", "calc_Tajima", "nSites")

# get theta stats
chr_stats <- theta_stats(donor_ns, stats_to_do, FALSE)
win_stats <- theta_stats(donor_ns, stats_to_do, TRUE)

# rownames
rownames(chr_stats) <- rownames(win_stats) <- 
  meta$site_code[meta$donor_name %in% paste0("donor", donor_ns)]

