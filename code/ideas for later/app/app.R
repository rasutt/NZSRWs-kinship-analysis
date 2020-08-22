# App to display possible POPs.  Taken down because haven't distinguished which
# are likely POPs.  Also the year and location data was arbitrary from often
# multiple recaptures

# Load shiny and data for possible POPs
library(shiny)
load("whales_and_possible_pop_indices.dat")

# Define UI with a plot and a table
ui <- basicPage(
  plotOutput(
    "plot1", 
    click = "plot_click",
    height = 500
  ),
  # verbatimTextOutput("info")
  tableOutput("info")
)

# Define server logic
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
  
  # Display plot
  output$plot1 <- renderPlot({
    
    # Set margin widths and text size
    par(mar = par()$mar - c(3.5, 0, 1, 0))
    par(cex = 0.95)
    
    # Setup plot axes and title
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
    
    # Get capture locations and one for selected animal
    capture_locations <- whale_genotypes$capture_location
    selected_capture_location <- capture_locations[vals$selected_animals]
    
    # Find kin of selected animal
    selected_kin <- 
      (vals$selected_animals[first_possible_parent_offspring_pair_index] |
         vals$selected_animals[second_possible_parent_offspring_pair_index])
    
    # Draw lines for each pair, wider and coloured for selected ones
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
    
    # Draw points for animals, wider for selected
    points(
      x = whale_genotypes$x,
      y = whale_genotypes$capture_year,
      col = (4 - 3 * (capture_locations == "Auckland Islands") -
               2 * (capture_locations == "Campbell Island") -
               1 * (capture_locations == "Mainland NZ")),
      cex = 0.95 + 0.7 * vals$selected_animals,
      pch = 20
    )
    
    # Add legend
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
  
  # Print table for selected kin
  # output$info <- renderPrint({
  output$info <- renderTable({
    
    # Find selected kin
    selected_kin <- c(
      first_possible_parent_offspring_pair_index[
        vals$selected_animals[second_possible_parent_offspring_pair_index]
      ],
      second_possible_parent_offspring_pair_index[
        vals$selected_animals[first_possible_parent_offspring_pair_index]
      ]
    )
    
    # Select rows and columns to print
    whale_genotypes[
      c(which(vals$selected_animals), selected_kin), 
      c(1, 28, 29, 2:27)
    ]
  })
}

# Run app
shinyApp(ui, server)