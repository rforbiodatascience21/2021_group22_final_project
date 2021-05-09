# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")
library("patchwork")
library("ggrepel")

# Load data ---------------------------------------------------------------
data_heatmap <- read_tsv("data/03_data_normalized_mean_across_replicates.tsv")
data_top_expr <- read_tsv(file = "data/03_data_mean_log2.tsv")
data_PCA_kmeans <- read_tsv("data/03_data_normalized_counts_and_raw_counts.tsv") 

# Heat map ----------------------------------------------------------
data_normalized_long <- data_heatmap %>%
  unite("experiment", treatment, time, sep = "_", remove = TRUE) %>%
  select(experiment, genes, mean_over_replicates) %>%
  distinct()

# Calculate z-score
data_zscore <- data_normalized_long %>%
  group_by(genes) %>%
  mutate(mean_counts_for_gene = mean(mean_over_replicates),
         sd_of_counts_for_gene = sd(mean_over_replicates),
         count_minus_mean = mean_over_replicates-mean_counts_for_gene,
         z_score = count_minus_mean/sd_of_counts_for_gene) %>%
  ungroup() %>%
  select(experiment, genes, z_score)

# Plot and save heatmap
heatmap_plot <- data_zscore %>%
  mutate(experiment = as_factor(experiment)) %>%
  ggplot(aes(y=genes, x=experiment, fill=z_score)) + 
  geom_tile() +
  scale_fill_gradient2(low = "yellow", high = "red") + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

# Top expression ----------------------------------------------------------
num_genes <- 20

# Select the top n differentially expressed genes for plotting
# Then pivot longer to get 1 gene per row
data_sorted_long <- data_top_expr %>%
  select(1:all_of(num_genes+2)) %>% # is this base R maybe change?
  pivot_longer(!c(treatment, time),
               names_to = "gene",
               values_to = "count")

order_names <- data_sorted_long %>%
  ungroup() %>%
  slice_head(n=num_genes) %>%
  pull(gene) %>%
  factor()

ggplot(data = data_sorted_long,
       mapping = aes(factor(gene, level = order_names),
                     count,
                     color = time,
                     shape=treatment)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = -45,
                                   vjust = -0.6,
                                   hjust=0.4)) +
  xlab("Top genes (by differential expression)") +
  ylab("log(count)") +
  scale_y_log10()

# Could we do this in tidy
gene_num <- 2
gene_name <- order_names[gene_num]

gene_plot_data <- data_sorted_long %>%
  filter(gene == gene_name)

top_expression_plot <- ggplot(data = gene_plot_data,
       aes(time, count,
           color = treatment,
           shape = treatment)) +
  geom_point(aes(size = 1.5,
                 alpha = 0.85)) +
  scale_y_log10() +
  theme_minimal() +
  ylab("log(count)") +
  scale_x_discrete(limits=c(2, 6, 10, 24)) +
  scale_alpha(guide = 'none') +
  scale_size(guide = 'none') +
  ggtitle(gene_name)

# PCA and Kmeans ----------------------------------------------------------
data_time_as_factor_correct_order <- data_PCA_kmeans %>% 
  arrange(time_as_numeric) %>% 
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
  geom_label_repel(mapping = aes(label = experiment), size = 2.5)

# Write data ------------------------------------------------------------
# Save heatmap
ggsave(path = "results",
       filename = "Heatmap.png", 
       plot = heatmap_plot )

# Save top expression
ggsave(path = "results",
       filename = paste("diffexpGenes_top",
                        as.character(num_genes),
                        ".png", sep=""), 
       plot = top_expression_plot)

ggsave(path = "results",
       filename = "04_PCA_plot.png",
       plot = PCA_plot)

ggsave(path = "results",
       filename = "04_Kmeans_plot.png",
       plot = K_means_plot)
