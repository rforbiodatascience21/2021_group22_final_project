# this script contains the functions used in the projects: 

#CLR_transform <- function(tibble, ){
#  
#}


top_genes_wide_to_long <- function(data, num_genes){
  data_sorted_long <- data %>%
    select(1:all_of(num_genes+2)) %>%
    pivot_longer(!c(treatment, time),
                 names_to = "gene",
                 values_to = "count")
  return(data_sorted_long)
}


top_gene_order <- function(data, num_genes){
  order_names <- data %>%
    ungroup() %>%
    slice_head(n=num_genes) %>%
    pull(gene) %>%
    factor()
  return(order_names)  
}


