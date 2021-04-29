# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")

# Load data ---------------------------------------------------------------
my_data_counts <- read_tsv(file = "data/01_my_data_counts.tsv")
my_data_samples <- read_tsv(file = "data/01_my_data_samples.tsv")
#my_data_counts
my_data_samples

# Wrangle data ------------------------------------------------------------

#Counts
my_data_counts_wide <- rename(my_data_counts, genes = X1) %>%
  pivot_longer(!genes, names_to = "experiment", values_to = "expression") %>%
  pivot_wider(names_from = "genes", values_from = "expression")
dim(my_data_counts_wide)

#Samples
my_data_samples <- rename(my_data_samples, experiment = X1) %>%
  select(experiment, treatment, time, replicate)
my_data_samples

my_data_clean <- inner_join(my_data_counts_wide, my_data_samples, by="experiment")
dim(my_data_clean)

# Write data --------------------------------------------------------------
write_tsv(x = my_data_clean,
          file = "data/02_my_data_clean.tsv")
