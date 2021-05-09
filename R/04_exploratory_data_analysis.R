# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")
library("patchwork")
library("ggrepel")

# Load data ---------------------------------------------------------------
data <- read_tsv("data/03_normalized_counts_and_raw_counts.tsv") 

# Checking normalization --------------------------------------------------

# Importing normalized count data 
data <- data %>% 
  arrange(time_as_numeric) %>% 
  mutate(time = as_factor(time))


# plotting boxplots to see if normalization worked
ggplot(data = data, mapping = aes(x = treatment, y = normalized_counts, fill = time)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10() + 
  facet_wrap(~replicate) + 
  scale_fill_viridis_d() +
  theme_minimal() + 
  labs(y = "Log10-scaled normalized counts", x = "Infection")



# Performing PCA on normalized data to see groupings ----------------------

# Importing normalized count data for PCA
PCA_data = read_tsv("data/03_data_normalized_counts_and_raw_counts.tsv") %>% 
  select(experiment,genes,counts,time) %>% 
  pivot_wider(values_from = counts,
              names_from = genes)

# Performing PCA analysis
PCA_analysis <- PCA_data %>% 
  column_to_rownames("experiment") %>% 
  select(-time) %>% 
  prcomp(center = TRUE, scale = TRUE)

# converting to tidy variance using broom 
variance_PCA <- PCA_analysis %>% 
  tidy("pcs") 

# getting PCA projected coordinates using broom
augment_PCA <- PCA_analysis %>% 
  augment(PCA_data) %>% 
  mutate(treatment = case_when(str_detect(string = experiment, pattern = "^Virus") == TRUE ~"Corona",
                               str_detect(string = experiment, pattern = "^Control") == TRUE ~"Control"))

# making scree plot
scree_plot <- variance_PCA %>% 
  ggplot(mapping = aes(x = PC,
                       y = percent)) +
  geom_col(fill = "lightblue") + 
  geom_line() + 
  labs(y = "Percent variance", 
       x = "PC number") +
  theme_minimal()

# making PCA biplot
biplot_PCA <- augment_PCA %>% 
  ggplot(mapping = aes(x = .fittedPC1,
                       y = .fittedPC2,
                       color = time,
                       shape = treatment)) + 
  geom_point(size = 3) + 
  labs(x = "PC 1", 
       y = "PC 2") + 
  theme_minimal()


# combining plots and saving
PCA_plot <- scree_plot / biplot_PCA

ggsave(filename = "results/04_PCA_plot.png",
       plot = PCA_plot,
       device = "png")


# Doing K-means on PCA transformed counts ---------------------------------

K_means_data <- augment_PCA %>% 
  select(.fittedPC1, .fittedPC2) %>% 
  kmeans(centers = 4) %>% 
  augment(augment_PCA) %>% 
  select(.fittedPC1, .fittedPC2, .cluster, experiment, time, treatment)

K_means_plot <- K_means_data %>% 
  ggplot(mapping = aes(x = .fittedPC1,
                       y = .fittedPC2,
                       color = .cluster)) + 
  geom_point(size = 2) + 
  theme_minimal() + 
  labs(x = "PC1",
       y = "PC2", 
       color = "Cluster") + 
  geom_label_repel(mapping = aes(label = experiment), size = 2.5)


