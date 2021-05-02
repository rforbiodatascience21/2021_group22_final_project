library(shiny)
library(shinythemes)

#Load data
data <- read_tsv("data/03_data_aug_sorted.tsv")

# Define UI
ui <- fluidPage(theme = shinytheme("yeti"),
                titlePanel("Gene expression in coronavirus patient compared to control"),
                sidebarPanel(
                  selectInput("genes", "Gene:",
                                choices = c("ALL"))
                  ),
                mainPanel('Boxplot', plotOutput("myplot"))
                
                

) # fluidPage


# Define server function  
server <- function(input, output) {
  
  output$myplot <- renderPlot({
    boxplot(data ~ get(input$genes), data = "data")
  })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)
