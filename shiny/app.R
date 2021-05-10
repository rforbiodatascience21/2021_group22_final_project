# Load libraries ----------------------------------------------------------
library("tidyverse")
library("shiny")
library("shinythemes")
library("patchwork")

setwd("/cloud/project/shiny")

# Load data ---------------------------------------------------------------
data_tab1 <- read_tsv("05_individual_times_ttest_and_data.tsv")
data_tab2 <- read_tsv("04_genes_sorted_by_highest_logFC_per_time.tsv")
data_tab3 <- read_tsv("05_linear_model_time_results.tsv")

# Wrangle data ------------------------------------------------------------
## Tab 1
data_times_seperated <- data_tab1 %>% 
  mutate(time_as_factor = as_factor(time)) 

# Making list of genes to choose between: 
unique_gene_names <- data_times_seperated %>% 
  select(genes) %>% 
  distinct() %>% 
  arrange(genes)

## Tab 3
# Chosing only variables needed forplotting 
time_model_for_plotting <- data_tab3 %>% 
  rename(time_stamp = time) %>% 
  pivot_wider(names_from = term,
              values_from = c(estimate, 
                              std.error),
              id_cols = c(Gene, 
                          time_stamp, 
                          LogFC))

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
                                                             "Number of upregulated genes: ",
                                                             min = 1,
                                                             max = 50,
                                                             value = 25),
                                                 sliderInput("numGenes2",
                                                             "Number of downregulated genes: ",
                                                             min = 1,
                                                             max = 50,
                                                             value = 25),
                                                 checkboxInput("logScale",
                                                               "log-scale y", 
                                                               TRUE)
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
  #tab1 -----
  # Tibble for top panel 
  temp_tibble_for_plotting <- reactive(data_times_seperated  %>% 
    filter(genes == input$Selected_gene))
  
  tibble_for_time_model <- reactive(time_model_for_plotting %>% 
    filter(Gene == input$Selected_gene2))
  
  #tab2 -----
  # data for top plot
  data_sorted_long <- reactive(data_tab2 %>%
                                 select(1:(input$numGenes+2)) %>%
                                 pivot_longer(-c(treatment, time),
                                            names_to = "gene",
                                            values_to = "count")
  )
  # plotting order for top plot
  order_names <- reactive(data_sorted_long() %>%
                            ungroup %>%
                            slice_head(n=input$numGenes) %>%
                            pull(gene) %>%
                            factor
  )
  
  # data for bottom plot
  data_sorted_long_bottom <- reactive(data_tab2 %>%
                                        select(-c(3:last_col(input$numGenes2))) %>%
                                        pivot_longer(-c(treatment, time),
                                                     names_to = "gene",
                                                     values_to = "count")
  )
  # plotting order for bottom plot
  order_names_bottom <- reactive(data_sorted_long_bottom() %>%
                                   ungroup %>%
                                   slice_head(n=input$numGenes2) %>%
                                   pull(gene) %>%
                                   factor
  )
  
  # tab2 plot1 overexpressed genes
  pal <- c("#FFCC00", "#FF3300", "#CC3399", "#660066")
  
  output$tab2_plot1 <- renderPlot({
    ggplot(data_sorted_long(),
           mapping = aes(factor(gene, level = order_names()),
                         count,
                         color=factor(time),
                         shape=treatment,
                         size=1.2,
                         alpha=0.8)) +
      geom_point() +
      theme(axis.text.x = element_text(angle=-45,
                                       vjust=-0.6,
                                       hjust=0.4,
                                       size=10)) +
      xlab("Genes") +
      {if(input$logScale)scale_y_log10()} +
      guides(size=FALSE, alpha=FALSE) +
      labs(title = "Top Genes",
           subtitle = "Genes over-expressed or upregulated in infected cells at t=24h") + 
      scale_colour_manual(values = pal)
  })
  
  #tab2 plot2 underexpressed genes
  output$tab2_plot2 <- renderPlot({
    ggplot(data_sorted_long_bottom(),
           mapping = aes(factor(gene, level = order_names_bottom()),
                         count,
                         color=factor(time),
                         shape=treatment,
                         size=1.2,
                         alpha=0.8)) +
      geom_point() +
      theme(axis.text.x = element_text(angle=-45,
                                       vjust=-0.6,
                                       hjust=0.4,
                                       size=10)) +
      xlab("Genes") +
      {if(input$logScale)scale_y_log10()} +
      guides(size=FALSE, alpha=FALSE) +
      labs(title = "Bottom Genes",
           subtitle = "Genes under-expressed or downregulated in infected cells at t=24h") + 
      scale_colour_manual(values = pal)
  })
  
  # Tab 3 -----
  output$myplot <- renderPlot({
    ggplot(temp_tibble_for_plotting(),
           mapping = aes(x = fct_reorder(time_as_factor,time),
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

}  #server end


# Create Shiny object
shinyApp(ui = ui, server = server)
