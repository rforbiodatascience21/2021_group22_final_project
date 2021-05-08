# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")
library("purrr")  #Purrr er en del af tidyverse så vi behøver ikke kalde den

# Load data ---------------------------------------------------------------
data_log2 <- read_tsv(file = "data/03_data_mean_log2_diff.tsv")

# Wrangle data ------------------------------------------------------------
# Move the data around
data_log2_long <- data_log2 %>%
  select(-NFIC) %>% #Still trouble with Inf and this gene -> clean
  pivot_longer(-time, names_to = "gene", values_to = "log2_expr_level") 

# Converting to nested data
data_log2_nested <- data_log2_long %>%
  group_by(gene) %>%
  nest %>% 
  ungroup()

# Fitting linear model to each of the genes
data_log2_nested <- data_log2_nested  %>% 
  mutate(mdl = map(data, ~lm(log2_expr_level ~ time,
                              data = .x)))

# Extract more model data using the broom package
data_log2_nested <- data_log2_nested %>%
  mutate(mdl_tidy = map(mdl, ~tidy(.x, conf.int = TRUE))) %>% 
  unnest(mdl_tidy)

# Looking only at slopes
data_log2_nested <- data_log2_nested %>% 
  filter(str_detect(term, "time"))

# Adding significance
data_log2_nested <- data_log2_nested %>% 
  mutate(identified_as = case_when(p.value < 0.05 ~ "Significant",
                                   TRUE ~ "Non-significant"),
         gene_label = case_when(identified_as == "Significant" ~ gene,
                                identified_as == "Non-significant" ~ ""))


# Different DE expression analysis that uses all replicates ---------------

set.seed(934485)

data = read_tsv("data/03_data_normalized_mean_across_replicates.tsv") %>% 
  group_by(genes) %>% 
  nest() %>% 
  ungroup() %>% 
  unnest(cols = data)

data_DE_analysis_nested <- data %>% 
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