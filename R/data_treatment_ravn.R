# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")

# Load data ---------------------------------------------------------------
bojkova_counts <- read_tsv(file = "data/01_my_data_counts.tsv")
bojkova_samples <- read_tsv(file = "data/01_my_data_samples.tsv")

# Wrangle ---------------------------------------------------------------
# Clean the counts 
bojkova_counts_clean <- bojkova_counts %>%
  pivot_longer(-X1, names_to = "conditions", values_to = "expression") %>%
  pivot_wider(names_from = "X1", values_from = "expression")

# Clean the samples
bojkova_samples_clean <- bojkova_samples %>%
  rename(conditions = X1) %>%
  select(conditions, treatment, time, replicate)

# Join the two datasets and remove the confusing "conditions"-column
bojkova_clean <- bojkova_samples_clean %>% 
  full_join(bojkova_counts_clean, by = "conditions") 

# Get the time as a factor (Maybe remove the h?)
bojkova_clean <- bojkova_clean %>%
  mutate(time = as_factor(time))

# ----------------------------
# Get genes
bojkova_genes <- bojkova_counts %>%
  rename(genes = X1) %>%
  pivot_longer(-genes, names_to = "conditions", values_to = "expression") 

bojkova_long <- bojkova_genes %>%
  full_join(bojkova_samples_clean, by = "conditions") %>%
  mutate(time = as_factor(time))

# Mean expression
bojkova_means <- bojkova_long %>%
  group_by(treatment, genes, time) %>%
  summarise(average_expression = mean(expression))

bojkova_means_wide <- bojkova_means %>%
  pivot_wider(names_from = "genes", values_from = "average_expression")


# Trying plots ---------------------------------------------------------------
# Visualizing one gene
ggplot(data = bojkova_clean, mapping = aes(x = time, y = AARSD1)) +
  geom_point()




