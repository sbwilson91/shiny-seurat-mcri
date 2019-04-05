

library(shiny)
library(tidyverse)
library(Seurat)

options(shiny.maxRequestSize = 3000*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Seurat analysis of Single Cell data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         fileInput(inputId = "dataset", label = NULL,  buttonLabel = ".rds data file", accept = ".rds")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         dataTableOutput("seurat"),
         plotOutput("tsne")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  SeuratObject <- reactive({
    
      # Load user-defined dataset
    load(input$dataset$datapath)
    
  })
   
   output$seurat <- renderDataTable({
      # generate bins based on input$bins from ui.R
     input$dataset
      
   })
   
   output$tsne <- renderPlot({
     load(file = input$dataset$datapath)
     TSNEPlot(get(SeuratObject()))
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

