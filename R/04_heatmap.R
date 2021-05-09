# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")

# Load data ---------------------------------------------------------------
data_sorted <- read_tsv(file = "data/03_data_mean_log2.tsv")
data_normalized <- read_tsv("data/03_data_normalized_mean_across_replicates.tsv")
# Wrangle data ------------------------------------------------------------

data_normalized_long <- data_normalized %>%
  unite("experiment", treatment, time, sep = "_", remove = TRUE) %>%
  select(experiment, genes, mean_over_replicates) %>%
  distinct()

# Calculate z-score
data_zscore <- data_normalized_long %>%
  group_by(genes) %>%
  mutate(mean_counts_for_gene = mean(mean_over_replicates),
         sd_of_counts_for_gene = sd(mean_over_replicates),
         count_minus_mean = mean_over_replicates-mean_counts_for_gene,
         z_score = count_minus_mean/sd_of_counts_for_gene) %>%
  ungroup() %>%
  select(experiment, genes, z_score)

# Find the top 500 deferentially expressed genes to plot
data_sorted_long <- data_sorted %>%
  select(0:502) %>%
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







