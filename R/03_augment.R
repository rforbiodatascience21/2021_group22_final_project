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
         time_as_numeric = as.numeric(str_extract(time, "\\d+"))) %>%    #time as numeric
  ungroup() 

# Take mean across of each set of replicates, change time variable to numeric
data_mean <- data_clean %>%
  select(-replicate) %>%
  group_by(treatment, time) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  mutate(time = as.numeric(str_extract(time, "\\d+")))

# Log2 transform all gene expression variables
data_mean_log2 <- data_mean %>%
  mutate_at(vars(-c(treatment, time)), log2)

## Signes bud på log2 difference -----------
# Calculate log2 diff (Maybe change names)
new_data_mean_log2_diff_long <- data_mean_log2 %>%
  pivot_longer(cols = c(-time, -treatment), names_to = "genes", values_to = "log2") %>%
  pivot_wider(names_from = "treatment", values_from = "log2") %>%
  group_by(time) %>% 
  mutate(log2_diff = Virus-Control) %>%
  ungroup() %>%
  select(-Virus, -Control)

# Convert back to tidy data (maybe change names)
new_data_mean_log2_diff <- new_data_mean_log2_diff_long %>%
  pivot_wider(names_from = "genes", values_from = "log2_diff")

## Signes bud sort ---------------
# Sort the genes by log2_diff and then time (high to low)
sorted_genes <- new_data_mean_log2_diff_long %>%
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

## Original log2 ----------
# Calculate the difference between virus and control expression for each time-point
data_mean_log2_diff <- data_mean_log2 %>% 
  group_by(time) %>%
  summarise_if(is.numeric, diff)          #virus minus control eller omvendt?

## original sort -----------
# Extract order of highest differential expression based on 2, 6, 10 or 24 hours
# Put time in rows
data_mean_log2_diff_2 <- new_data_mean_log2_diff %>% 
  column_to_rownames(var = "time")

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

write_tsv(x = data_mean_log2,
          file = "data/03_data_mean_log2.tsv")

write_tsv(x = data_sorted,
          file = "data/03_data_aug_sorted.tsv")

write_tsv(x = new_data_mean_log2_diff,
          file = "data/03_data_mean_log2_diff.tsv")


