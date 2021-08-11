library(shiny)
library(leaflet)
library(spData)
library(dplyr)
library(rgdal)
library(jsonlite)
library(tidyverse)
library(yaml)



# Define UI for application that filters map points based on year and minimum population
ui <- fluidPage(
    
    # Application title
    titlePanel("World Population Over Time"),
    
    # Sidebar with a slider input for year, numeric input for population 
    sidebarLayout(
        sidebarPanel(
            
            sliderInput("year",
                        "Year",
                        min = 1950,
                        max = 2030,
                        step = 5,
                        sep = "",
                        value = 1950),
            
            numericInput("pop_min",
                         "Minimum Population (in millions)",
                         min = 1,
                         max = 20,
                         value = 10)
        ),
        
        # Show the map and table
        mainPanel(
            # plotOutput("distPlot"),
            leafletOutput("map"),
            dataTableOutput("table")
        )
    )
)
