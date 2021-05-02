library("tidyverse")

data_sorted <- read_tsv(file = "data/02_my_data_clean.tsv")


num_genes <- 20

# Select the top n differentially expressed genes for plotting
# Then pivot longer to get 1 gene per row
data_sorted_long <- data_sorted %>% select(1:all_of(num_genes+2)) %>%
  pivot_longer(!c(treatment, time), names_to = "gene", values_to = "count")

order_names <- data_sorted_long %>% ungroup() %>% slice_head(n=num_genes) %>% pull(gene) %>% factor()

ggplot(data = data_sorted_long, aes(factor(gene, level = order_names), count, color = time, shape=treatment)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.6, hjust=0.4)) +
  xlab("Top genes (by differential expression)") +
  ylab("log(count)") +
  scale_y_log10()

ggsave(path = "results", filename = paste("diffexpGenes_top", as.character(num_genes), ".png", sep=""))

