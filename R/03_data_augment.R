# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")

# Load data ---------------------------------------------------------------
data_clean <- read_tsv(file = "data/02_my_data_clean.tsv")

# Wrangle data ------------------------------------------------------------
# Take mean across of each set of replicates, change time variable to numeric
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

# Calculating mean expression across replicates for every  
data_normalized_mean_across_replicates <- data_normalized %>% 
  select(c(-time_as_numeric,
           -counts,
           -total_counts)) %>% 
  group_by(treatment,
           time,
           genes) %>% 
  mutate(mean_over_replicates = mean(normalized_counts)) %>% 
  ungroup() %>% 
  select(treatment,
         time, 
         genes,
         replicate,
         normalized_counts,
         mean_over_replicates) %>% 
  distinct() 

# Take mean across of each set of replicates, change time variable to numeric
data_mean <- data_clean %>%
  select(-replicate) %>%
  group_by(treatment, time) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  mutate(time = as.numeric(str_extract(time, "\\d+")))

# Log2 transform all gene expression variables
data_mean_log2 <- data_mean %>%
  mutate_at(vars(-c(treatment, time)), log2)

# Calculate the difference between virus and control expression for each time-point
# Move time variable to rownames
data_mean_log2_diff <- data_mean_log2 %>% 
  group_by(time) %>%
  summarise_if(is.numeric, diff) 

data_mean_log2_diff_2 <- data_mean_log2_diff %>% 
  column_to_rownames(var = "time")

# Extract order of highest differential expression based on 2, 6, 10 or 24 hours
rownum <- 4   # 4 = 24h
sort <- data_mean_log2_diff_2 %>%
  slice(rownum) %>% 
  order(decreasing = TRUE)
# Modify sorting key to not affect time and treatment columns
sort2 <- sort+2
sort3 <- prepend(sort2, c(1,2))

# Sort the unmodified mean data according to the key
data_sorted <- data_mean %>% relocate(all_of(sort3))

# Write data --------------------------------------------------------------
write_tsv(x = data_normalized,
          file = "data/03_data_normalized_counts_and_raw_counts.tsv")

write_tsv(x = data_normalized_mean_across_replicates,
          file = "data/03_data_normalized_mean_across_replicates.tsv")

write_tsv(x = data_mean_log2,
          file = "data/03_data_mean_log2.tsv")

write_tsv(x = data_sorted,
          file = "data/03_data_aug_sorted.tsv")

write_tsv(x = data_mean_log2_diff,
          file = "data/03_data_mean_log2_diff.tsv")



