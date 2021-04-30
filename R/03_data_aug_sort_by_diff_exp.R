library("tidyverse")
library("forcats")

# load clean data and Move metadata variables to first columns
data_clean <- read_tsv(file = "data/02_my_data_clean.tsv") %>% 
  relocate(c(treatment, time, replicate), .after = experiment)


# Then take mean across of each set of replicates
# Finally change time variable to numeric
normalized_data <- data_clean %>% 
  pivot_longer(cols = c(-treatment, -time,-replicate,-experiment),
               names_to = "Genes",
               values_to = "Counts") %>% 
  group_by(experiment) %>% 
  mutate(total_counts = sum(Counts),
         normalized_counts = Counts/total_counts)


write_tsv(x = normalized_data,
          file = "data/3_normalized_counts_and_raw_counts.tsv")  
  

# Then take mean across of each set of replicates
# Finally change time variable to numeric

data_mean <- data_clean %>%
  select(-replicate) %>%
  group_by(treatment, time) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  mutate(time = as.numeric(str_extract(time, "\\d+")))



# Log2 transform all gene expression variables
# I think maybe you misunderstood - is log fold change not the log of the fold change? 
# not the fold change of the logs...

data_mean_log2 <- data_mean %>%
  mutate_at(vars(-c(treatment, time)), log2)

# Group by time to calculate the difference between virus and control expression for each time-point
# Move time variable to rownames


data_mean_log2_diff <- data_mean_log2 %>% 
  group_by(time) %>%
  summarise_if(is.numeric, diff) %>%
  column_to_rownames(var = "time")

# Extract order of highest differential expresison based 2, 6, 10 or 24 hours
rownum <- 4   # 4 = 24h
sort <- data_mean_log2_diff %>% slice(rownum) %>% order(decreasing = TRUE)
# Modify sorting key to not affect time and treatment columns
sort2 <- sort+2
sort3 <- prepend(sort2, c(1,2))

# Sort the unmodified mean data according to the key
data_sorted <- data_mean %>% relocate(all_of(sort3))

# Write data
write_tsv(x = data_sorted,
          file = "data/03_data_aug_sorted.tsv")
