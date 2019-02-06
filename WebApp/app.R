library(shiny)
library(Seurat)
library(dplyr)
library(Matrix)


ui <- fluidPage(
  titlePanel("Seurat Gallery"),
  fluidRow(
    column(3, 
           wellPanel(
             fileInput(inputId="file_name", label="Select Saved Seurat Object File"),
             radioButtons(inputId = "labelBoolean", label="Display labels?", choices = c("TRUE","FALSE"), selected="TRUE"),
             radioButtons(inputId = "legendBoolean", label="Hide legend?", choices = c("TRUE","FALSE"), selected="FALSE"),
             tags$hr(),
             sliderInput(inputId = "dotSize", label = "Set point size", value=1, min=0.01, max=10),
             sliderInput(inputId = "labelSize", label = "Set label size", value=4, min=0.5, max=10)
           ),
           wellPanel(
             textInput(inputId = "geneName", label="Enter Gene Name", value="Actb"),
             actionButton("goButton", "Plot Gene")
           )
    ),
    column(9,
           verbatimTextOutput("console1"),
           plotOutput("TSNE"),
           tags$hr(),
           tabsetPanel(
             tabPanel("TSNE", plotOutput("featureTSNE")),
             tabPanel("Violins", plotOutput("vlnPlot"))
           )
    )
  )
)


server <- function(input, output) {
  
  # increase the max upload file size
  options(shiny.maxRequestSize=1000*1024^2)
  
  # set demo dataset file here
  demo = "pancreas.Rda"
  
  SeuratObject <- reactive({
    if(is.null(input$file_name)){
      # Load a default demo dataset
      load(file=demo)
    } else {
      # Load user-defined dataset
      load(input$file_name$datapath)
    }
  })
  
  output$contents <- renderText({
    if (!is.null(input$file_name)){
      "Currently displaying 8.5k Pancreas dataset as demo"
    } else {
      load(file=inFile$name)
    }
  })
  
  output$console1 <- renderPrint({
    if (is.null(input$file_name)){
      get(load(demo))
    } else {
      # should be a more elegant way of doing this, but this works...
      get(load(input$file_name$name))
    }
  })
  
  # Display TSNE Plot as overview toggled options for display legends/labels and setting size. Default displays demo dataset
  output$TSNE <- renderPlot({
    if (is.null(input$file_name)){
      load(file=demo)
    } else {
      load(input$file_name$datapath)
    }
    TSNEPlot(get(SeuratObject()), do.label = input$labelBoolean, pt.size = input$dotSize, label.size = input$labelSize, no.legend=input$legendBoolean)
  })
  
  fPlot <- eventReactive(input$goButton, {
    if (is.null(input$file_name)){ 
      load(file=demo)
    } else {
      load(input$file_name$datapath)
    }
    FeaturePlot(get(SeuratObject()), input$geneName,pt.size = input$dotSize)
  })
  
  vPlot <- eventReactive(input$goButton, {
    if (is.null(input$file_name)){ 
      load(file=demo)
    } else {
      load(input$file_name$datapath)
    }
    VlnPlot(get(SeuratObject()), features.plot = input$geneName, size.x.use = NULLbo)
  })
  
  output$featureTSNE <- renderPlot(fPlot())
  output$vlnPlot <- renderPlot(vPlot())
  
}

shinyApp(ui = ui, server = server)