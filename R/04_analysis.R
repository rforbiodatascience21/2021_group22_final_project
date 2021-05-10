# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")
library("patchwork")
library("ggrepel")
library("viridis")          

# Load data ---------------------------------------------------------------
data_mean_over_replicates <- read_tsv(file = "data/03_data_normalized_count_mean_over_replicates.tsv")
data_log2 <- read_tsv(file = "data/03_data_logFC_for_time.tsv")
data_PCA_kmeans <- read_tsv(file = "data/03_data_normalized_counts.tsv") 

# Heat map ----------------------------------------------------------
data_normalized_long <- data_mean_over_replicates %>%
  unite("experiment", 
        treatment,
        time,
        sep = "_",
        remove = TRUE) %>%
  select(experiment, 
         genes,
         mean_over_replicates) %>%
  distinct

# Calculate z-score
data_zscore <- data_normalized_long %>%
  group_by(genes) %>%
  mutate(mean_counts_for_gene = mean(mean_over_replicates),
         sd_of_counts_for_gene = sd(mean_over_replicates),
         count_minus_mean = mean_over_replicates-mean_counts_for_gene,
         z_score = count_minus_mean/sd_of_counts_for_gene) %>%
  ungroup %>%
  select(experiment,
         genes, 
         z_score)

# Plot and save heatmap
heatmap_plot <- data_zscore %>%
  mutate(experiment = as_factor(experiment)) %>%
  ggplot(aes(y=genes,
             x=experiment,
             fill=z_score)) + 
  geom_tile() +
  scale_fill_gradient2(low = "yellow",
                       high = "red") + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle("Protein levels across all samples -standardized using z-scoring")

# Top expression ----------------------------------------------------------
# Creating a list of sorted genes
# Creating long version of log2
log2_diff_long <- data_log2 %>%
  pivot_longer(cols = -time,
               names_to = "genes",
              values_to = "log2_diff")

# Sort the genes by log2_diff and then time (high to low)
sorted_genes <- log2_diff_long %>%
  arrange(desc(log2_diff)) %>%
  arrange(desc(time)) %>%
  select(genes, time)

# Get the means in long format
data_mean <- data_mean_over_replicates %>%
  select(treatment,
         time,
         genes,
         mean_over_replicates) %>%
  distinct

# Change to order of the data mean long after highest logfold expression
sorted_means <- sorted_genes %>%
  full_join(x = .,
            y = data_mean,
            by = c("genes", "time"))

# Convert the sorted means back to tidy data format
sorted_means_wide <- sorted_means %>%
  pivot_wider(names_from = "genes",
              values_from = "mean_over_replicates")


# Select the top 20 differentially expressed genes for plotting
data_sorted_long <- sorted_means_wide %>%
  select(1:22) %>%
  pivot_longer(-c(treatment, time),
               names_to = "genes",
               values_to = "mean_over_replicates")

order_names <- data_sorted_long %>%
  ungroup %>%
  slice_head(n=20) %>%
  pull(genes) %>%
  factor

# Plot
top_20_genes <- ggplot(data = data_sorted_long,
       mapping = aes(factor(genes, 
                            level = as_factor(order_names)),
                     mean_over_replicates,
                     color = time,
                     shape = treatment)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = -45,
                                   vjust = -0.6,
                                   hjust = 0.4)) +
  ggtitle("Top 20 genes (virus, 24h)") +
  xlab("Top genes (by differential expression)") +
  ylab("count") +
  scale_y_log10()

# Now get just the top 8, as well as the bottom 8 (underexpressed)
data_sorted_long_top8 <- sorted_means_wide %>%
  select(1:10) %>%
  pivot_longer(cols = c(-treatment, 
                  -time),
               names_to = "genes",
               values_to = "mean_over_replicates")


data_sorted_long_bottom8 <- sorted_means_wide %>%
  select(cols = -c(3:last_col(8))) %>%
  pivot_longer(cols = c(-treatment, 
                  -time),
               names_to = "genes",
               values_to = "mean_over_replicates")

# Plot the expression of these selected genes over time
top_exp_plot <- ggplot(data = data_sorted_long_top8,
                  aes(x = time, 
                      y = mean_over_replicates,
                      shape = treatment)) +
  geom_point(size = 3,
             alpha = 0.5) +
  scale_y_log10() +
  theme_minimal() +
  ylab("count") +
  scale_x_discrete(limits = c(2, 6, 10, 24)) +
  scale_alpha(guide = 'none') +
  scale_size(guide = 'none') +
  ggtitle("Top 8 overexpressed genes in (virus, 24h)") +
  geom_line(aes(color = genes)) +
  geom_text_repel(aes(label = if_else(time == 24, 
                                   genes, '')),
                  hjust = 2,
                  size = 2.6,
                  xlim = c(24, 30)) +
  guides(color = FALSE) +
  theme(legend.position = c(0.15,0.85),
        legend.background = element_rect(),
        plot.margin = unit(c(10, 100, 10, 10), 
                           "points")) +
  coord_cartesian(clip = "off")

# Plot bottom
bottom_exp_plot <- ggplot(data = data_sorted_long_bottom8,
                     mapping = aes(x = time,
                         y = mean_over_replicates,
                         shape = treatment)) +
  geom_point(size = 3,
             alpha = 0.5) +
  scale_y_log10() +
  theme_minimal() +
  ylab("count") +
  scale_x_discrete(limits = c(2, 6, 10, 24)) +
  scale_alpha(guide = 'none') +
  scale_size(guide = 'none') +
  ggtitle("Bottom 8 underexpressed genes in (virus, 24h)") +
  geom_line(aes(color = genes)) +
  geom_text_repel(aes(label = if_else(time == 24, 
                                   genes, '')),
                  hjust = 2,
                  size = 2.6,
                  xlim = c(24, 30)) +
  guides(color = FALSE) +
  theme(legend.position = c(0.15, 0.15),
        legend.background = element_rect(),
        plot.margin=unit(c(10, 100, 10, 10), 
                         "points")) +
  coord_cartesian(clip = "off")

# Combine top and bottom
top_bottom <- top_exp_plot / bottom_exp_plot

# PCA and Kmeans ----------------------------------------------------------

data_time_as_factor_correct_order <- data_PCA_kmeans %>% 
  pivot_longer(cols = c(-time,
                        -treatment,
                        -replicate,
                        -experiment),
               names_to = "genes",
              values_to = "normalized_counts") %>% 
  arrange(time) %>% 
  mutate(time = as_factor(time))

# plotting boxplots to see if normalization worked
ggplot(data = data_time_as_factor_correct_order,
       mapping = aes(x = treatment,
                     y = normalized_counts,
                     fill = time)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10() + 
  facet_wrap(~replicate) + 
  scale_fill_viridis_d() +
  theme_minimal() + 
  labs(y = "Log10-scaled normalized counts", x = "Infection")

# Performing PCA on normalized data to see groupings ---

# Importing normalized count data for PCA
PCA_data <- data_time_as_factor_correct_order %>% 
  select(experiment,
         genes,
         normalized_counts,
         time) %>% 
  pivot_wider(values_from = normalized_counts,
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
  mutate(treatment = case_when(str_detect(string = experiment, 
                                          pattern = "^Virus") == TRUE ~"Corona",
                               str_detect(string = experiment, 
                                          pattern = "^Control") == TRUE ~"Control"))

# making scree plot
scree_plot <- variance_PCA %>% 
  ggplot(mapping = aes(x = PC,
                       y = percent)) +
  geom_col(fill = viridis(1)) + 
  geom_line() + 
  labs(y = "Percent variance", 
       x = "PC number") +
  theme_minimal() + 
  scale_fill_viridis_d()
  

# making PCA biplot
biplot_PCA <- augment_PCA %>% 
  ggplot(mapping = aes(x = .fittedPC1,
                       y = .fittedPC2,
                       color = time,
                       shape = treatment)) + 
  geom_point(size = 3) + 
  labs(x = "PC 1", 
       y = "PC 2") + 
  theme_minimal() + 
  scale_color_viridis_d()

# combining plots
PCA_plot <- scree_plot / biplot_PCA

# Doing K-means on PCA transformed counts -----
set.seed(1)

K_means_data <- augment_PCA %>% 
  select(matches("fittedPC[1-5]$")) %>% 
  kmeans(centers = 4) %>% 
  augment(augment_PCA)

K_means_plot <- K_means_data %>% 
  ggplot(mapping = aes(x = .fittedPC1,
                       y = .fittedPC2,
                       color = .cluster)) + 
  geom_point(size = 2) + 
  theme_minimal() + 
  labs(x = "PC1",
       y = "PC2", 
       color = "Cluster") + 
  geom_label_repel(mapping = aes(label = experiment), size = 2)


# Write data ------------------------------------------------------------
# Save heatmap
ggsave(path = "results",
       filename = "04_Heatmap.png", 
       plot = heatmap_plot)

# Save top expression
ggsave(path = "results",
       filename = str_c("04_diffexpGenes_top20.png"), 
       plot = top_20_genes)

ggsave(path = "results",
       filename = "04_top_and_bottom_8_over_time.png",
       width = 11,
       height = 10,
       plot = top_bottom)

ggsave(path = "results",
       filename = "04_top_8_over_time.png",
       plot = top_exp_plot)

ggsave(path = "results",
       filename = "04_bottom_8_over_time.png",
       plot = bottom_exp_plot)


# Save sorted means
write_tsv(x = sorted_means_wide,
          file = "results/04_genes_sorted_by_highest_logFC_per_time.tsv")

# Save PCA and K-means plots
ggsave(path = "results",
       filename = "04_PCA_plot.png",
       plot = PCA_plot)

ggsave(path = "results",
       filename = "04_Kmeans_plot.png",
       plot = K_means_plot)
