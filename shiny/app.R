# Load libraries ----------------------------------------------------------
library("tidyverse")
library("shiny")
library("shinythemes")
library("patchwork")

#setwd("/cloud/project/shiny")
source(file = "99_functions.R")

# Load data ---------------------------------------------------------------
data_for_plot_raw <- read_tsv("03_data_mean_log2.tsv")
data_for_plot_raw

# Wrangle data ------------------------------------------------------------
#prepare data for plot
data_for_plot <- data_for_plot_raw %>%
  pivot_longer(cols = c(-treatment, -time),
               names_to = "genes",
               values_to = "log_fold_change") 



# Loading data first panel ------------------------------------------------


data_times_seperated <- read_tsv("05_individual_times_ttest_and_data.tsv") %>% 
  mutate(time_as_numeric = as.numeric(str_extract(time, "\\d+"))) 

data_time_model <- read_tsv("05_linear_model_time_results.tsv")

# making list of genes to choose between: 
unique_gene_names <- data_times_seperated %>% 
  select(genes) %>% 
  distinct() %>% 
  arrange(genes)

#chosing only variables needed forplotting 
time_model_for_plotting <- data_time_model %>% 
  rename(time_stamp = time) %>% 
  pivot_wider(names_from = term,
              values_from = c(estimate, std.error),
              id_cols = c(Gene, time_stamp, LogFC))

# Loading data second panel -----------------------------------------------

data_tab2 <- read_tsv("03_data_aug_sorted.tsv")




# App creation ------------------------------------------------------------
# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage("SARS-CoV-19 protein expression",
                           tabPanel("Protein expression per time", 
                                    sidebarPanel(selectizeInput(inputId = "Selected_gene",
                                                             "Gene:",
                                                             choices = unique_gene_names,
                                                             selected = "SPIKE_WCPV",
                                                             options = list(maxItems = 1))),
                                    mainPanel(plotOutput("myplot"))
                                    ),
                           tabPanel("Plot top genes",
                                    sidebarPanel(sliderInput("numGenes",
                                                             "Number of genes: ",
                                                             min = 1,
                                                             max = 50,
                                                             value = 25),
                                                 checkboxInput("logScale",
                                                               "log-scale y", 
                                                               TRUE),
                                                 sliderInput("numGenes2",
                                                             "Number of genes: ",
                                                             min = 1,
                                                             max = 50,
                                                             value = 25)
                                    ),
                                    mainPanel(plotOutput("tab2_plot1"),
                                              plotOutput("tab2_plot2"))
                           ),
                           tabPanel("LogFC change over time", 
                                    sidebarPanel(selectizeInput(inputId = "Selected_gene2",
                                                                "Gene:",
                                                                choices = unique_gene_names,
                                                                selected = "SPIKE_WCPV",
                                                                options = list(maxItems = 1))),
                                    mainPanel(plotOutput("myplot2"))
                           )
                          )
                ) # fluidPage


# Define server function  
server <- function(input, output) {
  # tab 1 tibble for top panel 
  
  temp_tibble_for_plotting <- reactive(data_times_seperated  %>% 
    filter(genes == input$Selected_gene))
  
  tibble_for_time_model <- reactive(time_model_for_plotting %>% 
    filter(Gene == input$Selected_gene2))
  
  #tab2
  data_sorted_long <- reactive(
    top_genes_wide_to_long(data_tab2, input$numGenes)
  )
  order_names <- reactive(
    top_gene_order(data_sorted_long(), input$numGenes)
  )
  
  data_sorted_long_bottom <- reactive(
    bottom_genes_wide_to_long(data_tab2, input$numGenes2)
  )
  order_names_bottom <- reactive(
    top_gene_order(data_sorted_long_bottom(), input$numGenes2)
  )
  
  
  output$myplot <- renderPlot({
    ggplot(temp_tibble_for_plotting(),
           mapping = aes(x = fct_reorder(time,time_as_numeric),
                         y = normalized_counts,
                         fill = treatment,
                         color = significance)) +
    geom_boxplot() + 
    scale_colour_manual(values = c("black", "blue")) +
    labs(x = "Time",
         y = "Normalized Counts") + 
    theme_minimal()

  })
  output$myplot2 <- renderPlot({
      ggplot(tibble_for_time_model(),
             mapping = aes(x = time_stamp, y = LogFC)) + 
      geom_point() + 
      geom_abline(mapping = aes(intercept=mean(`estimate_(Intercept)`),
                                slope=mean(estimate_time))) +
      geom_abline(mapping = aes(slope=mean(estimate_time),
                                intercept=mean(`estimate_(Intercept)`)+mean(`std.error_(Intercept)`)),
                  lty = "dashed") +
      geom_abline(mapping = aes(slope=mean(estimate_time),
                                intercept=mean(`estimate_(Intercept)`)-mean(`std.error_(Intercept)`)),
                  lty = "dashed") +
      theme_minimal() +
      labs(x = "Time [Hours]",
           y = "Log2 Fold change")
  })
  
  # tab2 plot1 
  output$tab2_plot1 <- renderPlot({
      ggplot(data_sorted_long(),
             mapping = aes(factor(gene, level = order_names()),
                           count,
                           color=time,
                           shape=treatment,
                           size=1.2,
                           alpha=0.8)) +
        geom_point() +
        theme(axis.text.x = element_text(angle=-45,
                                         vjust=-0.6,
                                         hjust=0.4,
                                         size=10)) +
        xlab("Genes") +
        {if(input$logScale)ylab("log(count)")
        else ylab("count")} +
        {if(input$logScale)scale_y_log10()} +
        guides(size=FALSE, alpha=FALSE) +
        labs(title = "Top Genes",
             subtitle = "Ordered by differential expression at t=24h")
    })
  #tab2 plot2
  output$tab2_plot2 <- renderPlot({
    ggplot(data_sorted_long_bottom(),
           mapping = aes(factor(gene, level = order_names_bottom()),
                         count,
                         color=time,
                         shape=treatment,
                         size=1.2,
                         alpha=0.8)) +
      geom_point() +
      theme(axis.text.x = element_text(angle=-45,
                                       vjust=-0.6,
                                       hjust=0.4,
                                       size=10)) +
      xlab("Genes") +
      {if(input$logScale)ylab("log(count)")
        else ylab("count")} +
      {if(input$logScale)scale_y_log10()} +
      guides(size=FALSE, alpha=FALSE) +
      labs(title = "Top Genes",
           subtitle = "Ordered by differential expression at t=24h")
  })
}  #server


# Create Shiny object
shinyApp(ui = ui, server = server)
