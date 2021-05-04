# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")
library("purrr")

# Load data ---------------------------------------------------------------
data_log2 <- read_tsv(file = "data/03_data_mean_log2_diff.tsv")

# Wrangle data ------------------------------------------------------------
## PCA (One script, 05_model_i.R --> Output 1 plot)
# Move the data around
data_log2_long <- data_log2 %>%
  pivot_longer(-time, names_to = "gene", values_to = "log2_expr_level") 

# Converting to nested data
data_log2_nested <- data_log2_long %>%
  group_by(gene) %>%
  nest %>% 
  ungroup()

#Select random genes
set.seed(934485)
data_log2_nested = data_log2_nested %>%
  sample_n(100)

# Fitting general linear model to each of the 100 genes
data_log2_nested <- data_log2_nested  %>%
  mutate(mdl = map(data, ~glm(time ~ log2_expr_level,
                              data = .x)))

# Add some more model data using broom
data_log2_nested <- data_log2_nested %>%
  mutate(mdl_tidy = map(mdl, ~tidy(.x, conf.int = TRUE))) %>% 
  unnest(mdl_tidy)

# Looking only at slopes
data_log2_nested <- data_log2_nested %>% 
  filter(str_detect(term, "level"))

# Adding significance
data_log2_nested <- data_log2_nested %>% 
  mutate(identified_as = case_when(p.value < 0.05 ~ "Significant",
                                   TRUE ~ "Non-significant"),
         gene_label = case_when(identified_as == "Significant" ~ gene,
                                identified_as == "Non-significant" ~ ""))

# Negative log p values
data_log2_nested <- data_log2_nested %>% 
  mutate(neg_log10_p = -log10(p.value))

#PCA plot
data_log2_wide <- data_log2_long %>%
  pivot_wider(names_from = "gene", values_from = "log2_expr_level")

data_log2_wide <- data_log2_wide %>%
  select(time, pull(data_log2_nested, gene))

pca_fit <- bojkova_data_wide %>%
  select(where(is.numeric)) %>%
  prcomp(scale = TRUE)

pca_fit %>%
  augment(bojkova_data_wide) %>% 
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, 
             color = as.numeric(time))) + 
  geom_point(size = 1.5)  + 
  theme_classic(base_family = "Avenir") + 
  theme(legend.position = "bottom", 
        panel.grid.major = element_line()) + 
  labs(title = "PCA coordinate plot", color = "Outcome", x = "fittedPC1",
       y = "fittedPC2")



# Different DE expression analysis that uses all replicates ---------------

set.seed(934485)

data = read_tsv("data/03_data_normalized_mean_across_replicates.tsv") %>% 
  group_by(genes) %>% 
  nest() %>% 
  ungroup() %>% 
  sample_n(500) %>% 
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
  
write_tsv(augmented_model_results, 
          file = "results/05_individual_times_ttest_and_data.tsv")


## K means clustering (One script, 05_model_ii.R --> Output 1 plot)

# Write data --------------------------------------------------------------