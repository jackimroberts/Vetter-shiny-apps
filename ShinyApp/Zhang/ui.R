library(shiny)
library(shinydashboard)

shinyUI(dashboardPage(skin="purple",
  dashboardHeader(title="Jianmin's Data"),
  #sidebar ####
  dashboardSidebar(
    sidebarMenu(
      selectInput("dataset", label = "Select a dataset:",selected=4,
        choices = list(
          "PRC2 cko retina" = 4,
          "Jarid2 single cell" = 5)
        ),
      textInput("goi", label = "Ensembl gene name:", value = "Foxp1")
    )
  ),
  
    #body ####
    dashboardBody(
      h3(htmlOutput("dataheader")),
      htmlOutput("legend1"),
      uiOutput('geneboxes')
    )
))