# First some exploratory plots showing if the samples are equal: 

library("tidyverse")

data = read_tsv("data/02_my_data_clean.tsv")

data2 = read_tsv("data/01_my_data_counts.tsv")
