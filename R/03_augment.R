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
         time = as.numeric(str_extract(time, 
                                       "\\d+"))) %>%   
  ungroup %>% 
  select(c(-counts,
           -total_counts)) %>% 
  pivot_wider(names_from = genes, 
              values_from = normalized_counts)

# Calculate means of replicates of the normalized data
data_normalized_mean_across_replicates <- data_normalized %>% 
  pivot_longer(cols = c(-experiment,
                        -treatment,
                        -time,
                        -replicate),
               names_to = "genes", 
               values_to = "normalized_counts") %>% 
  group_by(treatment,
           time,
           genes) %>% 
  mutate(mean_over_replicates = mean(normalized_counts)) %>% 
  ungroup %>% 
  select(treatment,
         time,
         genes,
         replicate,
         normalized_counts,
         mean_over_replicates) %>% 
  distinct 

# Convert back to tidy data with the normalized mean for each gene
data_mean <- data_normalized_mean_across_replicates %>%
  select(treatment,
         time,
         genes,
         mean_over_replicates) %>%
  distinct %>%
  pivot_wider(names_from = "genes",
              values_from = "mean_over_replicates")

# Log2 transform all gene expression variables
data_log2 <- data_mean %>%
  mutate_at(vars(-c(treatment, 
                    time)), log2)

# Calculate log2 diff for each gene
data_log2_diff_long <- data_log2 %>%
  pivot_longer(cols = c(-time, 
                        -treatment),
               names_to = "genes",
               values_to = "log2") %>%
  pivot_wider(names_from = "treatment",
              values_from = "log2") %>%
  group_by(time) %>% 
  mutate(log2_diff = Virus-Control) %>%
  ungroup %>%
  select(-Virus, 
         -Control)

# Convert back to tidy data 
log2_diff <- data_log2_diff_long %>%
  pivot_wider(names_from = "genes",
              values_from = "log2_diff")

# Write data --------------------------------------------------------------
write_tsv(x = data_normalized,
          file = "data/03_data_normalized_counts.tsv")

write_tsv(x = data_normalized_mean_across_replicates,
          file = "data/03_data_normalized_count_mean_over_replicates.tsv")

write_tsv(x = log2_diff,
          file = "data/03_data_logFC_for_time.tsv")
