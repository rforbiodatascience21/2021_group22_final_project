# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")
library("magrittr")

# Load data ---------------------------------------------------------------
my_data_counts <- read_tsv(file = "data/01_my_data_counts.tsv")
my_data_samples <- read_tsv(file = "data/01_my_data_samples.tsv")
my_data_counts

# Wrangle data ------------------------------------------------------------


my_data_counts_wide <- rename(my_data_counts, genes = X1)
my_data_counts_wide

my_data_counts_wide <- my_data_counts_wide %>%
  pivot_longer(!genes, names_to = "experiment", values_to = "expression") %>%
  pivot_wider(names_from = "genes", values_from = "expression")
my_data_counts_wide




# Write data --------------------------------------------------------------
write_tsv(x = my_data_clean,
          file = "data/02_my_data_clean.tsv")
