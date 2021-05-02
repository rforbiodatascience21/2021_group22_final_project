# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")

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
  mutate(time = as.factor(time)) %>%
  mutate(time = str_remove(time, "h"))

# Mean expression
bojkova_means <- bojkova_long %>%
  group_by(treatment, genes, time) %>%
  summarise(average_expression = mean(expression))

bojkova_means_wide <- bojkova_means %>%
  pivot_wider(names_from = "genes", values_from = "average_expression")

# ----------------------------
# Logfold changes
bojkova_log <- bojkova_means %>%
  group_by(genes, time) %>%
  summarise(log2_expr_level = diff(log2(average_expression))) 

# Trying plots ---------------------------------------------------------------
# Visualizing one gene
ggplot(data = bojkova_clean, mapping = aes(x = time, y = AARSD1)) +
  geom_point()

# Visualizing all genes and their logfold changes
ggplot(data = bojkova_clean, mapping = aes(x = time, y = AARSD1)) +
  geom_point()

# Trying models -------------------------------------------------------------
# Converting time to numeric to be able to model
bojkova_log2 <- bojkova_log %>%
  mutate(time = as.numeric(time))

# Converting to long nested data
bojkova_log_nested <- bojkova_log2 %>%
  group_by(genes) %>%
  nest %>% 
  ungroup

#Select random genes
set.seed(934485)
bojkova_log_nested = bojkova_log_nested %>%
  sample_n(100)

# Fitting general linear model to each of the 100 genes
bojkova_log_nested  <- bojkova_log_nested  %>%
  mutate(mdl = map(data, ~glm(time ~ log2_expr_level,
                              data = .x)))

# Add some more model data using broom
bojkova_log_nested <- bojkova_log_nested %>%
  mutate(mdl_tidy = map(mdl, ~tidy(.x, conf.int = TRUE))) %>% 
  unnest(mdl_tidy)

# Looking only at slopes
bojkova_log_nested = bojkova_log_nested %>% 
  filter(str_detect(term, "level"))

# Adding significance
bojkova_log_nested <- bojkova_log_nested %>% 
  mutate(identified_as = case_when(p.value < 0.05 ~ "Significant",
                                   TRUE ~ "Non-significant"),
         gene_label = case_when(identified_as == "Significant" ~ genes,
                                identified_as == "Non-significant" ~ ""))

#Negative log p values
bojkova_log_nested <- bojkova_log_nested %>% 
  mutate(neg_log10_p = -log10(p.value))

#PCA plot
bojkova_log2_wide <- bojkova_log2 %>%
  pivot_wider(names_from = "genes", values_from = "log2_expr_level")

bojkova_data_wide <- bojkova_log2_wide %>%
  select(time, pull(bojkova_log_nested, genes))

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
## Does not work
## Does it makes sense to reduce dimensions when we have already reduced by log

# ---------------------------------------
# Manhattan plot of 100 random genes
ggplot(data = bojkova_log_nested, 
       mapping = aes(x = genes,
             y = neg_log10_p,
             label = gene_label, # kan ikke få gene_label på
             colour = identified_as)) + 
  geom_point() +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed")  +
  theme_classic(base_family = "Avenir",
                base_size = 12) +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom") +
  labs(x = "Gene",
       y = "Minus log10(p)") +
  ggtitle("Manhattan plot of 100 random genes")

ggsave(filename = "results/04_Manhattan_plot_prøve.png", width = 16, height = 9, dpi = 72)

# CI plot effect
ggplot(data = bojkova_log_nested, mapping = aes(x = estimate,
             y = fct_reorder(genes, desc(estimate)),
             colour = identified_as,
             label = gene_label)) +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_errorbarh(aes(xmin = conf.low,
                     xmax = conf.high,
                     height = 0.2)) +
  geom_text(aes(x = conf.high),
            size = 3,
            colour = "black",
            nudge_x = 0.5) +
  theme_classic(base_family = "Avenir",
                base_size = 12) +
  theme(axis.text.y = element_blank(),
        legend.position = "bottom") +
  labs(y = "") + 
  ggtitle("CIplot of 100 random genes")

ggsave(filename = "results/04_CI_plot_prøve.png", width = 16, height = 9, dpi = 72)
