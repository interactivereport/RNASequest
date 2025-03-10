#This app enables users to create an online table listing all QuickOmics Projects with one click link to launch each project.
#Edit the project.csv file in this folder to add/remove projects. 
#The URLs for projects should be QuickOmics URL of your instance, followed by /?serverfile=projectID, or /?testfile=projectID. Use the config.csv file in QuickOmics directory to define server file and test file folders.

library(shiny)
library(DT)
library(dplyr)
library(stringr)

createLink <- function(val, name="View Project") {
  sprintf('<a href="%s" target="_blank" class="btn btn-primary">%s</a>',val, name)
}

ui <- fluidPage(
  titlePanel(
    fluidRow(
      column(4, img(height =75 , src = "ShinyOne_logo.png")), 
      column(8,  h2("Repository of RNASequest Results", align = 'left'))
    ),  windowTitle = "ShinyOne" ),
    tags$hr(style="border-color: RoyalBlue;"),
    dataTableOutput('table1')
) 

server <- function(input, output) {
  
  output$table1 <- renderDataTable({
    my_table=read.csv("projects.csv", check.names=F)
    my_table<-my_table%>%dplyr::mutate(Visualization=createLink(Visualization, name="Quickomics") )%>%
      dplyr::mutate(Bookdown=ifelse(str_detect(Bookdown, "http"), createLink(Bookdown, name="Bookdown"), Bookdown ), 
                    Slidedeck=ifelse(str_detect(Slidedeck, "http"), createLink(Slidedeck, name="Slidedeck"), Slidedeck ) )             
    #browser() #debug
    return(my_table) }
   , escape = FALSE)
}

shinyApp(ui, server)
