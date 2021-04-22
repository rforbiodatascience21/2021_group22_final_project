# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")
library("readxl")


# Load data ---------------------------------------------------------------
my_data_raw <- read_xls("data/_raw/pone.0207943.s001.xls")
#my_data_raw

# Write data
write_tsv(x = my_data_raw,
          file = "data/01_my_data.tsv")
