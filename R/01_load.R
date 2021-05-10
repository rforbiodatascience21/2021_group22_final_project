# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")
library("readxl")

# Load data ---------------------------------------------------------------
my_data_raw_counts <- read_csv(file = "data/_raw/counts.csv")
my_data_raw_samples <- read_csv(file ="data/_raw/samples.csv")


# Write data
write_tsv(x = my_data_raw_counts,
          file = "data/01_data_counts.tsv")

write_tsv(x = my_data_raw_samples,
          file = "data/01_data_samples.tsv")