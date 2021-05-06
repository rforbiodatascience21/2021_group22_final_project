library("tidyverse")
library("shiny")
library("shinythemes")

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
                           tabPanel("navbar 2"),
                           tabPanel("navbar 3")
                            )
                ) # fluidPage


# Define server function  
server <- function(input, output) {
  temp_tibble_for_plotting <- reactive(data_times_seperated  %>% 
    filter(genes == input$Selected_gene))
  
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
}  #server


# Create Shiny object
shinyApp(ui = ui, server = server)
