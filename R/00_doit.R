# Install packages --------------------------------------------------------

if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!requireNamespace("shiny", quietly = TRUE))
  install.packages("shiny")

if (!requireNamespace("shinythemes", quietly = TRUE))
  install.packages("shinythemes")

if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork")

if (!requireNamespace("readxl", quietly = TRUE))
  install.packages("readxl")

if (!requireNamespace("broom", quietly = TRUE))
  install.packages("broom")

if (!requireNamespace("ggrepel", quietly = TRUE))
  install.packages("ggrepel")


# Run all scripts ---------------------------------------------------------
source(file = "R/01_load.R")
source(file = "R/02_clean.R")
source(file = "R/03_augment.R")
source(file = "R/04_analysis.R")
source(file = "R/04_PCA_and_kmeans.R")
source(file = "R/04_heatmap.R")
source(file = "R/04_plot_top_expression.R")
source(file = "R/05_modelling.R")
