# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")
library("magrittr")

# Load data ---------------------------------------------------------------
#Problem: the dataset has multiple headers
#Solution: leave second explanatory header out of dataset but keep as dictionary

#load only header
my_data_header <- read_tsv(file = "data/01_my_data.tsv", n_max = 0) %>%
  names()
#my_data_header

#load dataset except for the two first rows (headers) and then insert first row as header
my_data <- read_tsv(file = "data/01_my_data.tsv", skip = 2, col_names = my_data_header)
#my_data

#save headers as key value pairs in dictionary
headers_dict <- read_tsv(file = "data/01_my_data.tsv", n_max = 1) %>%
  gather(variable_name, variable_description)
#headers_dict

# Wrangle data ------------------------------------------------------------
#dim(my_data)
my_data_clean <- distinct(my_data) %>% 
  drop_na() 
dim(my_data_clean)




# Write data --------------------------------------------------------------
write_tsv(x = my_data_clean,
          file = "data/02_my_data_clean.tsv")
