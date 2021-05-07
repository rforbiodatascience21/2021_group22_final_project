library("tidyverse")
library("shiny")
library("shinythemes")

setwd("/cloud/project")
source(file = "R/99_functions.R")

#Load data
data_for_plot_raw <- read_tsv("data/03_data_mean_log2.tsv")

data_for_plot_raw

#prepare data for plot
data_for_plot <- data_for_plot_raw %>%
  pivot_longer(cols = c(-treatment, time),
               names_to = "genes",
               values_to = "log_fold_change") 

# Loading data for one of the subpanels: 
data_times_seperated <- read_tsv("results/05_individual_times_ttest_and_data.tsv") %>% 
  mutate(time_as_numeric = as.numeric(str_extract(time, "\\d+"))) 

# making list of genes to choose between: 
unique_gene_names <- data_times_seperated %>% 
  select(genes) %>% 
  distinct() %>% 
  arrange(genes)

data_tab2 <- read_tsv("data/03_data_aug_sorted.tsv")

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage("My first navbar page",
                           tabPanel("Protein expression per time", 
                                    sidebarPanel(selectizeInput(inputId = "Selected_gene",
                                                             "Gene:",
                                                             choices = unique_gene_names,
                                                             selected = "SPIKE_WCPV",
                                                             options = list(maxItems = 1))),
                                    mainPanel(plotOutput("myplot"))
                                    ),
                           tabPanel("Plot top genes",
                                    sidebarPanel(sliderInput("numGenes", "Number of genes: ",
                                                             min = 1, max = 50,
                                                             value = 25),
                                                 checkboxInput("logScale", "log-scale y", TRUE)
                                    ),
                                    mainPanel(plotOutput("plot2"))
                           ),
                           tabPanel("navbar 3")
                            )
                ) # fluidPage


# Define server function  
server <- function(input, output) {
  temp_tibble_for_plotting <- reactive(data_times_seperated  %>% 
    filter(genes == input$Selected_gene))
  
  #tab2
  data_sorted_long <- reactive(
    top_genes_wide_to_long(data_tab2, input$numGenes)
  )
  order_names <- reactive(
    top_gene_order(data_sorted_long(), input$numGenes)
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
  
  #plot2 
  output$plot2 <- renderPlot({
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
    
}  #server


# Create Shiny object
shinyApp(ui = ui, server = server)
