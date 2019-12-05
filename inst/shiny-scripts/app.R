library(shiny)
library(SAVERg)

# Define UI for random distribution app ----
ui <- fluidPage(
  # App title ----
  titlePanel("SAVERg Demo"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
      numericInput(inputId = "percent",
                   label = "Percent", min = 0, max = 1,
                   value = 0.1),
      numericInput(inputId = "ncores",
                   label = "Number of parallel cores for computing",
                   value = 4),
      br(),
      radioButtons("need", "Need imputation:",
                   c("Yes" = "Yes",
                     "No" = "No"),
                   selected = "No"),
      radioButtons("imputed", "Imputed data:",
                   c("Yes" = "Yes",
                     "No" = "No"),
                   selected = "Yes"),
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Cell Clustering", plotOutput("plot1")),
                  tabPanel("Trajectory Analysis", plotOutput("plot2"))
      )
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {
  # Generate a plot from cell clustering
  output$plot1 <- renderPlot({
    if (input$need == "Yes") {
      need = TRUE
    } else {
      need = FALSE
    }
    if (input$imputed == "Yes") {
      imputed = TRUE
    } else {
      imputed = FALSE
    }
    cell_clustering(ipsc_saver,
                    percent = input$percent,
                    ncores = input$ncores,
                    need.imputation = need,
                    imputed.data = imputed)
  })

  # Generate a plot from trajectory analysis
  output$plot2 <- renderPlot({
    if (input$need == "Yes") {
      need = TRUE
    } else {
      need = FALSE
    }
    if (input$imputed == "Yes") {
      imputed = TRUE
    } else {
      imputed = FALSE
    }
    deng_cellLabels <- factor(colnames(deng_saver),
                              levels=c('zygote', 'early 2-cell',
                                       'mid 2-cell', 'late 2-cell',
                                       '4-cell', '8-cell', '16-cell', 'early blastocyst',
                                       'mid blastocyst', 'late blastocyst'))
    trajectory_analysis(deng_saver, deng_cellLabels,
                        percent = input$percent,
                        ncores = input$ncores,
                        need.imputation = need,
                        imputed.data = imputed)
  })
}

# Create Shiny app ----
shinyApp(ui, server)
