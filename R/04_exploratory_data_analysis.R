# First some exploratory plots showing if the samples are equal: 

library("tidyverse")
library("broom")
library("patchwork")


# Checking normalization --------------------------------------------------

data = read_tsv("data/03_normalized_counts_and_raw_counts.tsv") %>% 
  arrange(time_as_numeric) %>% 
  mutate(time = as_factor(time))

ggplot(data = data, mapping = aes(x = treatment, y = normalized_counts, fill = time)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10() + 
  facet_wrap(~replicate) + 
  scale_fill_viridis_d() +
  theme_minimal() + 
  labs(y = "Log10-scaled normalized counts", x = "Infection")



# Performing PCA on normalized data to see groupings ----------------------

PCA_data = read_tsv("data/03_data_normalized_counts_and_raw_counts.tsv") %>% 
  select(experiment,genes,counts,time) %>% 
  pivot_wider(values_from = counts,
              names_from = genes)

PCA_analysis <- PCA_data %>% 
  column_to_rownames("experiment") %>% 
  select(-time) %>% 
  prcomp(center = TRUE, scale = TRUE)
  
variance_PCA <- PCA_analysis %>% 
  tidy("pcs") 

augment_PCA <- PCA_analysis %>% 
  augment(PCA_data) %>% 
  mutate(treatment = case_when(str_detect(string = experiment, pattern = "^Virus") == TRUE ~"Corona",
                               str_detect(string = experiment, pattern = "^Control") == TRUE ~"Control"))

scree_plot <- variance_PCA %>% 
  ggplot(mapping = aes(x = PC,
                       y = percent)) +
  geom_col(fill = "lightblue") + 
  geom_line() + 
  labs(y = "Percent variance", 
       x = "PC number") +
  theme_minimal()
  
biplot_PCA <- augment_PCA %>% 
  ggplot(mapping = aes(x = .fittedPC1,
                       y = .fittedPC2,
                       color = time,
                       shape = treatment)) + 
  geom_point(size = 3) + 
  labs(x = "PC 1", 
       y = "PC 2") + 
  theme_minimal()

scree_plot / biplot_PCA

png(file="results/04_PCA_plot",
    width=600, height=350)

  
  
  