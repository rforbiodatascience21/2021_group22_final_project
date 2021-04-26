# Subset a number of proteins. Random ?
set.seed(934485)
corona_data_long_nested = corona_data_long_nested %>%
  sample_n(100)

# Plot control vs corona-infected? Average over 3 replicates?
# 4 different times - 2 hours, 6 hours, 10 hours, 24 hours


# K means clustering (One script, 05_model_i.R --> Output 1 plot)

# PCA (One script, 05_model_ii.R --> Output 1 plot)

# Manhattan plot (One script, 05_model_iii.R --> Output 1 plot)