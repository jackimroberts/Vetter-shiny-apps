library(shiny)
library(shinydashboard)

shinyUI(dashboardPage(skin="purple",
  dashboardHeader(title="Data Explorer"),
  #sidebar ####
  dashboardSidebar(
    sidebarMenu(
      selectInput("dataset", label = "Select a dataset:",selected=4,
        choices = list(
          "retinal μglia timecourse" = 1,
          "Itgax lineage" =3,
          "Axl Mertk KO" =4,
          "retinal uglia single cell" = 2) 
        ),
      textInput("goi", label = "Ensembl gene name:", value = "Sparcl1")
    )
  ),
  
    #body ####
    dashboardBody(
      h3(htmlOutput("dataheader")),
      htmlOutput("legend1"),
      uiOutput('geneboxes')
    )
))