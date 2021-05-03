# Jens augment

library("tidyverse")
library("purrr")

data <- read_tsv("data/02_my_data_clean.tsv") 

ggplot(data = data, mapping = aes(x = AAAS, y = AACS))

