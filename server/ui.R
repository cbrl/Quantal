library(shiny)
library(shinythemes)


shinyUI(fluidPage(theme = shinytheme("journal"),
    tags$style("#code {font-family: \"Lucida Console\", Monaco, monospace;}"),

    # Application title
    titlePanel("Quantum Computing"),
    
    sidebarLayout(
        #position = "right",

        sidebarPanel(
            width = 3,
            textAreaInput(
                "code",
                label = "Code Editor",
                value = "qubits 2\nH  1\nCX 1 2",
                rows = 30,
                resize = "vertical"
            ),
            actionButton("execute", "Execute")
        ),#sidebarPanel
        
        mainPanel(
            width = 9,
            
            fluidRow(style="height: 40vh;",
                h4("Circuit Diagram"),
                imageOutput("diagram", width = "100%", height = "90%")
            ),
            
            #hr(),
            p(style="width: 100%; border-top: 1px solid #eee;"),
            
            fluidRow(style="height: 50vh;",
                column(4, style="height: 100%;",
                    h4("State Probabilities"),
                    plotOutput("prob_plot", height="90%")
                ),
                column(8, style="height: 100%; border-left: 1px solid #eee;",
                    h4("Statevector"),
                    column(8, style="height: 100%;",
                        plotOutput("sv_plot", height="90%"),
                    ),
                    column(4,
                        verbatimTextOutput("states")
                    )
                )
            )
        )#mainPanel
    )#sidebarLayout
))
