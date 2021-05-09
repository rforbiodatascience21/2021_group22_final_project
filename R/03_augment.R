# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")

# Load data ---------------------------------------------------------------
data_clean <- read_tsv(file = "data/02_data_clean.tsv")

# Wrangle data ------------------------------------------------------------
# Normalize data, change time variable to numeric
data_normalized <- data_clean %>% 
  pivot_longer(cols = c(-treatment,
                        -time,
                        -replicate,
                        -experiment),
               names_to = "genes",
               values_to = "counts") %>% 
  group_by(experiment) %>% 
  mutate(total_counts = sum(counts),
         normalized_counts = (6000000/total_counts)*counts, 
         time_as_numeric = as.numeric(str_extract(time, "\\d+"))) %>%   
  ungroup() 

# Calculate means of the normalized data
data_normalized_mean_across_replicates <- data_normalized %>% 
  select(c(-counts, -total_counts)) %>% 
  group_by(treatment, time, genes) %>% 
  mutate(mean_over_replicates = mean(normalized_counts)) %>% 
  ungroup() %>% 
  select(treatment, time, genes, replicate,
         normalized_counts, mean_over_replicates) %>% 
  distinct() 

# Convert back to tidy data with the normalized mean for each gene
data_mean <- data_normalized_mean_across_replicates %>%
  select(treatment, time, genes, mean_over_replicates) %>%
  distinct() %>%
  pivot_wider(names_from = "genes", values_from = "mean_over_replicates")

# Log2 transform all gene expression variables
data_log2 <- data_mean %>%
  mutate_at(vars(-c(treatment, time)), log2)

# Calculate log2 diff for each gene
data_log2_diff_long <- data_log2 %>%
  pivot_longer(cols = c(-time, -treatment), names_to = "genes", values_to = "log2") %>%
  pivot_wider(names_from = "treatment", values_from = "log2") %>%
  group_by(time) %>% 
  mutate(log2_diff = Virus-Control) %>%
  ungroup() %>%
  select(-Virus, -Control)

# Convert back to tidy data 
log2_diff <- data_log2_diff_long %>%
  pivot_wider(names_from = "genes", values_from = "log2_diff")

# Sort the genes by log2_diff and then time (high to low)
sorted_genes <- data_log2_diff_long %>%
  arrange(desc(log2_diff)) %>%
  arrange(desc(time)) %>%
  select(genes, time)

# Change the dataframe with means to fit the long format
data_mean_long <- data_mean %>% 
  pivot_longer(cols = c(-time, -treatment), names_to = "genes", values_to = "counts") 

# Change to order of the data mean long after highest logfold expression
sorted_means <- sorted_genes %>%
  full_join(x = ., y = data_mean_long, by = c("genes", "time"))

# Convert the sorted means back to tidy data format
sorted_means_wide <- sorted_means %>%
  pivot_wider(names_from = "genes", values_from = "counts")

# Write data --------------------------------------------------------------
write_tsv(x = data_normalized,
          file = "data/03_data_normalized_counts_and_raw_counts.tsv")

write_tsv(x = data_normalized_mean_across_replicates,
          file = "data/03_data_means.tsv")

write_tsv(x = data_log2,
          file = "data/03_data_mean_log2.tsv")

write_tsv(x = sorted_means_wide,
          file = "data/03_data_aug_sorted.tsv")

write_tsv(x = log2_diff,
          file = "data/03_data_mean_log2_diff.tsv")


