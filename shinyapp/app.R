library(shiny)
library(viridisLite)

# ---- Load data ----
all_chrome <- read.csv("all_chrome.csv", stringsAsFactors = FALSE)
all_chrome$haploid <- as.numeric(all_chrome$haploid)

clade_choices <- sort(unique(all_chrome$clade))

# ---- UI ----
ui <- fluidPage(
  titlePanel("Chromosome Data Across the Tree of Life"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      checkboxGroupInput(
        "clades",
        "Select clades:",
        choices = clade_choices
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Distributions", plotOutput("distPlot", height = "600px")),
        tabPanel(
          "Data table",
          # Download button for the table
          downloadButton("downloadData", "Download CSV"),
          br(), br(),
          tableOutput("cladeTable")
        )
      )
    )
  )
)

# ---- SERVER ----
server <- function(input, output, session) {
  
  filtered_data <- reactive({
    req(input$clades)
    all_chrome[all_chrome$clade %in% input$clades, , drop = FALSE]
  })
  
  # ---- Plot densities ----
  output$distPlot <- renderPlot({
    
    df <- filtered_data()
    clades <- input$clades
    n <- length(clades)
    
    # All haploid values across selected clades
    all_vals <- df$haploid[!is.na(df$haploid)]
    xlim <- range(all_vals)
    
    # Precompute densities and global y max
    density_list <- list()
    max_y <- 0
    
    for (cl in clades) {
      x <- df$haploid[df$clade == cl & !is.na(df$haploid)]
      if (length(x) > 1) {
        d <- density(x)
        density_list[[cl]] <- d
        max_y <- max(max_y, max(d$y))
      }
    }
    
    if (max_y == 0) {
      plot.new()
      text(0.5, 0.5, "Not enough data to compute density curves")
      return()
    }
    
    ylim <- c(0, max_y * 1.05)
    
    # Base colors + transparent versions for fill
    base_cols <- viridis(n)
    fill_cols <- sapply(base_cols, function(z) adjustcolor(z, alpha.f = 0.3))
    
    # Empty canvas
    plot(0,
         type = "n",
         xlim = xlim,
         ylim = ylim,
         xlab = "Haploid chromosome number",
         ylab = "Density",
         main = "Haploid chromosome number distributions (density)"
    )
    
    # Draw filled polygons + outlines
    for (i in seq_along(clades)) {
      cl <- clades[i]
      d  <- density_list[[cl]]
      
      if (!is.null(d)) {
        # Filled polygon under the curve
        polygon(
          x = c(d$x, rev(d$x)),
          y = c(d$y, rep(0, length(d$y))),
          col = fill_cols[i],
          border = NA
        )
        
        # Outline on top
        lines(d, col = base_cols[i], lwd = 2)
      }
    }
    
    # Legend
    legend("topright",
           legend = clades,
           col = base_cols,
           lwd = 3,
           cex = 0.9,
           bty = "n")
  })
  
  # ---- Table ----
  output$cladeTable <- renderTable({
    filtered_data()
  })
  
  # ---- Download handler ----
  output$downloadData <- downloadHandler(
    filename = function() {
      clades <- input$clades
      if (!length(clades)) {
        "chromosome_data.csv"
      } else {
        # Safe-ish filename based on clades
        paste0("chromosome_data_", paste(clades, collapse = "_"), ".csv")
      }
    },
    content = function(file) {
      df <- filtered_data()
      write.csv(df, file, row.names = FALSE)
    }
  )
}

# ---- Run app ----
shinyApp(ui = ui, server = server)
