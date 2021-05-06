# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")

# Load data ---------------------------------------------------------------
data_normalized = read_tsv("data/03_normalized_counts_and_raw_counts.tsv")
data_sorted <- read_tsv(file = "data/03_data_mean_log2.tsv")
# Wrangle data ------------------------------------------------------------
head(data_normalized)

num_genes <- 50

# Select the top n differentially expressed genes for plotting
# Then pivot longer to get 1 gene per row
data_sorted_long <- data_sorted %>%
  select(1:all_of(num_genes+2)) %>%
  pivot_longer(!c(treatment, time),
               names_to = "genes",
               values_to = "count") %>%
  distinct(genes)  
dim(data_sorted_long)
head(data_sorted_long)

# 
data_plot <- data_normalized %>%
  select(experiment, genes, normalized_counts) %>%
  group_by(genes) %>%
  mutate(z_score = normalized_counts-mean(normalized_counts)/sd(normalized_counts)) %>%
  ungroup()
head(data_plot)

data_plot2 <- data_plot %>%
  select(experiment, genes, z_score)

data_plot3 <- right_join(data_plot2, data_sorted_long, by="genes")

head(data_plot3)
dim(data_plot3)
head(data_plot3)

ggplot(data = data_plot3, aes(x = experiment, genes)) +
  geom_tile(aes(fill = z_score)) +
  scale_fill_gradient2(low = "yellow", high = "red")

ggsave(filename = "results/04_heatmap.png",
       plot = Heatmap,
       device = "png")
