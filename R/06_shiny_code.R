library("tidyverse")
library("shiny")
library("shinythemes")

#Load data
data_for_plot_raw <- read_tsv("../data/03_data_mean_log2.tsv")

data_for_plot_raw

#prepare data for plot
data_for_plot <- data_for_plot_raw %>%
  pivot_longer(cols = c(-treatment, time),
               names_to = "genes",
               values_to = "log_fold_change") 



# Define UI
ui <- fluidPage(theme = shinytheme("yeti"),
                titlePanel("Gene expression x hours after infection with coronavirus compared to control"),
                
                sidebarPanel(
                  selectInput("genes", "Gene:",
                              choices = c("genes"))
                ),
                mainPanel('Boxplot', plotOutput("myplot"))
                
                

) # fluidPage


# Define server function  
server <- function(input, output) {
  
  output$myplot <- renderPlot({
    
    
    boxplot((data$genes), data = "data")
    
  })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)
