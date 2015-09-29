library(shiny)
library(shinythemes)

chroms = c(22,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,'X','Y','M')

# Define the overall UI

shinyUI(navbarPage(
  "INDEL tools",
  tabPanel(
    "EIR",
    fluidPage(
      theme = shinytheme("flatly"), title = "INDEL Region Calculator",
      column(2),
      column(
        12,align = 'center',
        h1("INDEL Region Calculator",align =
             'center'),br(),
        
        # Create a new Row in the UI for selectInputs
        fluidRow(
          align = 'center',
          column(2,
                 selectInput(
                   "chr", "Chromosome:",unique(as.character(chroms))
                 )),
          column(2,
                 numericInput("loc", "Location (hg19):", value = 50885245)),
          column(4,
                 textInput("allele", "Alt allele", value = "TC")),
          column(4,
                 fileInput(
                   'file1', 'Choose CSV File',
                   accept =
                     c('text/csv',
                       'text/comma-separated-values,text/plain',
                       '.csv')
                 ))
        ),
        # Create a new row for the table.
        fluidRow(align = 'center',
                 plotOutput(outputId = 'plot',height = 150))
      ),
      column(2)
    )
  ),
  tabPanel(
    "Table",
    fluidPage(
      theme = shinytheme("flatly"), title = "INDEL Region Calculator",
      column(2),
      column(
        12,align = 'center',
        h1("INDEL Region Calculator",align =
             'center'),br(),
        
        # Create a new row for the table.
        fluidRow(
          align = 'center',
          h4("1000 Genome project variants within the EIR"),
          h4("(currently only for chr22)"),
          dataTableOutput(outputId = "table")
        )
      ),
      column(2)
    )
  )
))