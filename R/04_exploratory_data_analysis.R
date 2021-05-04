# First some exploratory plots showing if the samples are equal: 

library("tidyverse")

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

