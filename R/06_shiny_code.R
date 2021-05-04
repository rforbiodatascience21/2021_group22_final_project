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
data_times_seperated <- read_tsv("results/05_individual_times_ttest_and_data.tsv")

# making list of genes to choose between: 
unique_gene_names <- data_times_seperated %>% 
  select(genes) %>% 
  distinct() %>% 
  arrange(genes)

# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage("My first navbar page",
                           tabPanel("navbar 1", 
                                    sidebarPanel(selectInput(inputId = "bum",
                                                             "Gene:",
                                                             choices = unique_gene_names)),
                                    mainPanel('Boxplot', plotOutput("myplot"))
                                    ),
                           tabPanel("navbar 2"),
                           tabPanel("navbar 3")
                            )
                ) # fluidPage


# Define server function  
server <- function(input, output) {
  output$myplot <- renderPlot({
    boxplot((data$genes), data = "data")
  })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)
