# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library("tidyverse")

# Load data ---------------------------------------------------------------
data_sorted <- read_tsv(file = "data/03_data_mean_log2.tsv")

# Wrangle data ------------------------------------------------------------
num_genes <- 20

# Select the top n differentially expressed genes for plotting
# Then pivot longer to get 1 gene per row
data_sorted_long <- data_sorted %>%
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

gene_num <- 2
gene_name <- order_names[gene_num]

gene_plot_data <- data_sorted_long %>%
  filter(gene == gene_name)

ggplot(data = gene_plot_data,
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
  
# Write data --------------------------------------------------------------
ggsave(path = "results",
       filename = paste("diffexpGenes_top",
                        as.character(num_genes),
                        ".png", sep=""))

