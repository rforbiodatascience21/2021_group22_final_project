# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")

# Load data ---------------------------------------------------------------
clean_data = read_tsv("data/02_my_data_clean.tsv")
data_sorted <- read_tsv(file = "data/03_data_mean_log2.tsv")

# Wrangle data ------------------------------------------------------------
head(clean_data)

long_data <- clean_data %>%
  pivot_longer(cols = c(-experiment, -treatment, -time, -replicate), 
               names_to = "genes", values_to = "counts") %>%
  select(treatment, time, genes, counts) %>%
  unite("experiment", treatment, time, sep = "_", remove = TRUE) %>%
  group_by(experiment, genes) %>%
  mutate(mean_counts = sum(counts)/3) %>%
  select(experiment, genes, mean_counts) %>%
  distinct()


head(long_data)
dim(long_data)

calculate_zscore <- long_data %>%
  group_by(genes) %>%
  mutate(mean_counts_for_gene = mean(mean_counts)) %>%
  mutate(sd_of_counts_for_gene = sd(mean_counts)) %>%
  mutate(z_score_step1 = mean_counts-mean_counts_for_gene) %>%
  mutate(z_score = z_score_step1/sd_of_counts_for_gene) %>%
  ungroup() %>%
  arrange(mean_counts_for_gene)

head(calculate_zscore, n = 10)
dim(calculate_zscore)


# Wrangle data ------------------------------------------------------------

num_genes <- 500

# Select the top n differentially expressed genes for plotting
# Then pivot longer to get 1 gene per row
data_sorted_long <- data_sorted %>%
  select(1:all_of(num_genes+2)) %>%
  pivot_longer(!c(treatment, time),
               names_to = "genes",
               values_to = "count") %>%
  distinct(genes)

head(data_sorted_long)

data_plot3 <- semi_join(calculate_zscore, data_sorted_long, by="genes") %>%
  select(experiment) %>%
  rename(C1 = Control_2h)

  
  
dim(data_plot3)
head(data_plot3, n = 10)


ggplot(data = data_plot3, aes(x = rowid[1:8], y = genes)) +
  geom_tile(aes(fill = z_score)) +
  scale_fill_gradient2(low = "yellow", high = "red")
  #scale_x_discrete(labels=c("Control 2h", "Control 6h", "Control ","1h","2h","3h","6h","12h","24h","48h"))

ggsave(path = "results",
       filename = "Heatmap.png")
















