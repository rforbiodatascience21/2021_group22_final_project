# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")

# Load data ---------------------------------------------------------------
clean_data = read_tsv("data/02_my_data_clean.tsv")
data_sorted <- read_tsv(file = "data/03_data_mean_log2.tsv")

# Wrangle data ------------------------------------------------------------
# Find mean of replicates
long_data <- clean_data %>%
  pivot_longer(cols = c(-experiment, -treatment, -time, -replicate), 
               names_to = "genes", 
               values_to = "counts") %>%
  select(treatment, time, genes, counts) %>%
  unite("experiment", treatment, time, sep = "_", remove = TRUE) %>%
  group_by(experiment, genes) %>%
  mutate(mean_counts = sum(counts)/3) %>%
  select(experiment, genes, mean_counts) %>%
  distinct()

# Calculate z-score
data_zscore <- long_data %>%
  group_by(genes) %>%
  mutate(mean_counts_for_gene = mean(mean_counts),
         sd_of_counts_for_gene = sd(mean_counts),
         count_minus_mean = mean_counts-mean_counts_for_gene,
         z_score = count_minus_mean/sd_of_counts_for_gene) %>%
  ungroup() %>%
  select(experiment, genes, z_score)

# Find the top n deferentially expressed genes to plot
num_genes <- 500

data_sorted_long <- data_sorted %>%
  select(1:all_of(num_genes+2)) %>%
  pivot_longer(!c(treatment, time),
               names_to = "genes",
               values_to = "count") %>%
  distinct(genes)

# Only keep z-score for the top n deferentially expressed genes chosen earlier
data_plot <- semi_join(data_zscore, data_sorted_long, by="genes")

# Plot and save heatmap
data_plot %>%
  mutate(experiment = as_factor(experiment)) %>%
  ggplot(aes(y=genes, x=experiment, fill=z_score)) + 
  geom_tile() +
  scale_fill_gradient2(low = "yellow", high = "red") + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


ggsave(path = "results",
       filename = "Heatmap.png")







