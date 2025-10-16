###
# packages

# dplyr
library(dplyr)

###
# read and manipulate data

# base directory
  base_dir <- "/Users/petrucci/Documents/research/clothesmoths_popgen/"

# read locations table
locs <- read.delim(paste0(base_dir, "bams_locations_12p.csv"), sep = ",",
                   header = TRUE)
locs$bam <- gsub("\\..*", "", locs$bam)

# read read counts table
counts <- read.delim(paste0(base_dir, "read_counts.tsv"), header = TRUE,
                     sep = "\t")

# read alignment rates
rates <- read.delim(paste0(base_dir, "alignment_rates.tsv"), header = TRUE,
                    sep = "\t")

# left join bams
rates_counts <- counts %>%
  left_join(rates, by = "bam") %>%
  inner_join(locs, by = "bam") %>%
  select(sample_name, read_counts, alignment_rates)

# save table
write.table(rates_counts, paste0(base_dir, "read_table.csv"), sep = ",",
            quote = FALSE, row.names = FALSE)
