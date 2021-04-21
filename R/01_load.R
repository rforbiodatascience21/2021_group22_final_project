# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")


# Load data ---------------------------------------------------------------
my_data_raw <- read_tsv(file = "data/_raw/my_raw_data.tsv")
