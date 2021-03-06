# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")

# Load data ---------------------------------------------------------------
my_data_counts <- read_tsv(file = "data/01_data_counts.tsv")
my_data_samples <- read_tsv(file = "data/01_data_samples.tsv")

# Wrangle data ------------------------------------------------------------

#Prepare counts dataset for joining and replace zero values with low value 
my_data_counts_wide <- my_data_counts %>% 
  rename(genes = X1) %>%
  pivot_longer(cols = -genes,
               names_to = "experiment",
               values_to = "expression") %>%
  mutate(expression = case_when(expression == 0 ~ 0.00001, 
                                expression > 0.0 ~ expression)) %>%
  pivot_wider(names_from = "genes",
              values_from = "expression")


#Samples
my_data_samples <- my_data_samples %>%
  rename(experiment = X1) %>%
  select(experiment, 
         treatment,
         time,
         replicate)


#join dataframes
my_data_clean <- inner_join(my_data_samples,
                            my_data_counts_wide,
                            by="experiment") %>% 
  mutate(experiment = str_replace_all(experiment, c("\\dh$" = "2h_1",
                                                    "[ ]" = "_")))
  

# Write data --------------------------------------------------------------
write_tsv(x = my_data_clean,
          file = "data/02_data_clean.tsv")
