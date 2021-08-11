
urban_agglomerations <- urban_agglomerations %>%
    mutate(longitude = map_dbl(.x=.$geometry,function(x) x[1]),
           latitude = map_dbl(.x=.$geometry,function(x) x[2])) %>%
    as.data.frame()
urban_agglomerations$geometry <- NULL

coordinates(urban_agglomerations) <- c("longitude", "latitude")

# Define server logic required to draw a map and table
server <- function(input, output) {
    
    
    output$map <- renderLeaflet({
        

        
        #geom <- urban_agglomerations$geometry
        pop_by_year <- filter(urban_agglomerations, 
                              year == input$year,
                              population_millions > input$pop_min)
        #geom -> urban_agglomerations$geometry
        
        leaflet(data = pop_by_year) %>%
            addTiles() %>%
            addMarkers()
    })
    
    output$table <- renderDataTable({
        
        pop_by_year <- filter(urban_agglomerations, 
                              year == input$year,
                              population_millions > input$pop_min)
        
        pop_by_year
        
    })
}
