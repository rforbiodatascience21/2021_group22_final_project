library(shiny)
library(shinythemes)

#Load data
data <- read_tsv("../data/03_normalized_counts_and_raw_counts.tsv")
data

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
