library(shiny)
load("whales_and_possible_pop_indices.dat")

ui <- basicPage(
  plotOutput(
    "plot1", 
    click = "plot_click",
    height = 500
  ),
  # verbatimTextOutput("info")
  tableOutput("info")
)

server <- function(input, output) {
  # Set reactive for selected animals to all false
  vals <- reactiveValues(
    selected_animals = rep(F, nrow(whale_genotypes))
  )
  
  # Updated selected animal when clicked
  observeEvent(input$plot_click, {
    vals$selected_animals <- nearPoints(
      whale_genotypes, 
      input$plot_click, 
      maxpoints = 1,
      xvar = "x",
      yvar = "capture_year",
      allRows = TRUE
    )$selected_
  })
  
  output$plot1 <- renderPlot({
    par(mar = par()$mar - c(3.5, 0, 1, 0))
    par(cex = 0.95)
    
    plot.new()
    plot.window(
      xlim = c(0, 1),
      ylim = c(2018, 1995)
    )
    axis(
      2, 
      cex.axis = 0.75,
      at = 1995:2018
    )
    title(
      main = "New Zealand southern right whales and parent-offspring connections",
      ylab = "Capture year"
    )

    capture_locations <- whale_genotypes$capture_location
    selected_capture_location <- capture_locations[vals$selected_animals]
    
    selected_kin <- 
      (vals$selected_animals[first_possible_parent_offspring_pair_index] |
         vals$selected_animals[second_possible_parent_offspring_pair_index])

    segments(
      first_possible_parent_offspring_animal_x,
      first_possible_parent_offspring_animal_y,
      second_possible_parent_offspring_animal_x,
      second_possible_parent_offspring_animal_y,
      col = c(
        rgb(0.2, 0.2, 0.2, alpha = 0.1),
        4 - 3 * (selected_capture_location == "Auckland Islands") -
          2 * (selected_capture_location == "Campbell Island") -
          1 * (selected_capture_location == "Mainland NZ")
      )[1 + selected_kin],
      lwd = c(1, 2)[1 + selected_kin]
  )
    
    points(
      x = whale_genotypes$x,
      y = whale_genotypes$capture_year,
      col = (4 - 3 * (capture_locations == "Auckland Islands") -
               2 * (capture_locations == "Campbell Island") -
               1 * (capture_locations == "Mainland NZ")),
      cex = 0.95 + 0.7 * vals$selected_animals,
      pch = 20
    )
    
    legend(
      "bottomleft", 
      legend = c(
        "Auckland Islands", 
        "Campbell Island",
        "Mainland NZ",
        "Calf"
      ),
      col = 1:4,
      pch = 20
    )
  })
  
  # output$info <- renderPrint({
  output$info <- renderTable({
    selected_kin <- c(
      first_possible_parent_offspring_pair_index[
        vals$selected_animals[second_possible_parent_offspring_pair_index]
        ],
      second_possible_parent_offspring_pair_index[
        vals$selected_animals[first_possible_parent_offspring_pair_index]
        ]
    )
    
    whale_genotypes[
      c(which(vals$selected_animals), selected_kin), 
      c(1, 28, 29, 2:27)
      ]
  })
}

shinyApp(ui, server)