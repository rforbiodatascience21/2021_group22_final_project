# Project functions ------------------------------------------------

# Takes genes ordered from highest to lowest differential log expression
# Outputs top n (overexpressed) genes in long format for plotting
top_genes_wide_to_long <- function(data, num_genes){
  data_sorted_long <- data %>%
    select(1:all_of(num_genes+2)) %>%
    pivot_longer(!c(treatment, time),
                 names_to = "gene",
                 values_to = "count")
  return(data_sorted_long)
}

# Same as above but takes bottom n (underexpressed) genes
bottom_genes_wide_to_long <- function(data, num_genes){
  data_sorted_long <- data %>%
    select(-(3:last_col(num_genes))) %>%
    pivot_longer(!c(treatment, time),
                 names_to = "gene",
                 values_to = "count")
  return(data_sorted_long)
}

# Given a tibble as output in above functions, gives the gene names as factor for correct plotting order
top_gene_order <- function(data, num_genes){
  order_names <- data %>%
    ungroup() %>%
    slice_head(n=num_genes) %>%
    pull(gene) %>%
    factor()
  return(order_names)  
}


