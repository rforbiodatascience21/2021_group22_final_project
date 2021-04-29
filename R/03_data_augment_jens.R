# Jens augment

library("tidyverse")
library("purrr")

data <- read_tsv("data/02_my_data_clean.tsv") %>% 
  mutate(treatment = case_when(experiment == ))
  
