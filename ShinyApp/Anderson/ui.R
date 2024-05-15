library(shiny)
library(shinydashboard)

shinyUI(dashboardPage(skin="purple",
  dashboardHeader(title="Sarah's data"),
  #sidebar ####
  dashboardSidebar(
    sidebarMenu(
      selectInput("dataset", label = "Select a dataset:",selected=2,
        choices = list(
          "retinal μglia timecourse" = 1,
          "Axl Mertk KO" =4,
          "retinal uglia single cell" = 2) 
        ),
      textInput("goi", label = "Ensembl gene name:", value = "Itgax")
    )
  ),
  
    #body ####
    dashboardBody(
      h3(htmlOutput("dataheader")),
      htmlOutput("legend1"),
      uiOutput('geneboxes')
    )
))