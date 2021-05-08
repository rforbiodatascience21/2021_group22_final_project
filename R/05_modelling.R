# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")

# Load data ---------------------------------------------------------------
data_log2 <- read_tsv(file = "data/03_data_mean_log2_diff.tsv")
data <- read_tsv("data/03_data_normalized_mean_across_replicates.tsv")

# Wrangle data ------------------------------------------------------------
# General linear model -----------------------
# Fit general linear model to each gene 
model_nested <- data_log2 %>%
  pivot_longer(cols = -time, 
               names_to = "Gene", 
               values_to = "LogFC") %>%
  group_by(Gene) %>% 
  nest() %>%
  ungroup() %>%
  mutate(model = map(.x = data, 
                     .f = ~glm(formula = LogFC ~ time, data = .x)))

# Add more model statistics using broom
model_nested <- model_nested %>%
  mutate(tidymodel = map(.x = model, .f = ~tidy(.x))) %>%
  unnest(tidymodel) 

# Unnest the data again for later plotting
model_unnested <- model_nested %>%
  select(-model) %>%
  unnest(data) 

# Different DE expression analysis that uses all replicates ---------------

set.seed(934485)

data_DE <- data %>% 
  group_by(genes) %>% 
  nest() %>% 
  ungroup() %>% 
  unnest(cols = data)

data_DE_analysis_nested <- data_DE %>% 
  group_by(time, genes) %>% 
  nest() %>% 
  mutate(mdl = map(.x = data,
                   .f = ~glm(data = .x,
                   formula = normalized_counts ~ treatment))) 

data_tidy_model <- data_DE_analysis_nested %>% 
  mutate(mdl_tidy = map(.x = mdl, ~tidy(.x,conf.int=TRUE)))

unnested_tidy_model <- data_tidy_model %>% 
  unnest(mdl_tidy)

augmented_model_results <- unnested_tidy_model %>% 
  mutate(regulation = case_when(estimate > 0 ~ "Upregulated",
                                estimate < 0 ~ "Downregulated"),
         significance = case_when(p.value >= 0.05 ~ "Not significant",
                                  p.value < 0.05 ~ "Significant")) %>% 
  filter(term == "treatmentVirus") %>% 
  select(genes,
         time,
         data,
         estimate,
         p.value, 
         regulation,
         significance) %>% 
  unnest(data)

# Write data --------------------------------------------------------------
write_tsv(augmented_model_results, 
          file = "results/05_individual_times_ttest_and_data.tsv")
write_tsv(model_unnested, 
          file = "results/05_linear_model_time_results.tsv")
