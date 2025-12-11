library(shiny)
library(viridisLite)

# ---- Load data ----
all_chrome <- read.csv("all_chrome.csv", stringsAsFactors = FALSE)
all_chrome$haploid <- as.numeric(all_chrome$haploid)

# Standardize clade names to lower case
all_chrome$clade <- tolower(all_chrome$clade)

# ---- Load higher classification ----
higher_class <- read.csv("higher_class.csv", stringsAsFactors = FALSE)

# Standardize clade names here too
higher_class$Clade <- tolower(higher_class$Clade)

# Keep only clades that are in all_chrome (safety)
higher_class <- higher_class[higher_class$Clade %in% all_chrome$clade, ]

# Unique higher-level groups (you can reorder this if you want a custom order)
higher_levels <- unique(higher_class$Higher.Classification)

# ---- UI ----
ui <- fluidPage(
  # Some light CSS to make the collapsible headers look nice
  tags$head(
    tags$style(HTML("
      .clade-scroll {
        max-height: 400px;
        overflow-y: auto;
        padding-right: 5px;
        border: 1px solid #ddd;
        border-radius: 4px;
      }
      details > summary {
        font-weight: bold;
        margin-top: 8px;
        cursor: pointer;
      }
      details > summary::-webkit-details-marker {
        display: none;
      }
      details summary::before {
        content: '▸ ';
      }
      details summary::before {
        content: '▾ ';
      }
    "))
  ),
  
  titlePanel("Chromosome Data Across the Tree of Life"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # Record counts
      textOutput("totalRecords"),
      textOutput("filteredRecords"),
      br(),
      
      h4("Select clades:"),
      
      # Scrollable container for the grouped clade selectors
      div(
        class = "clade-scroll",
        uiOutput("cladeSelector")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Distributions", plotOutput("distPlot", height = "600px")),
        tabPanel(
          "Data table",
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
  
  # ---- Dynamic grouped checkboxes in collapsible sections ----
  output$cladeSelector <- renderUI({
    tagList(
      lapply(seq_along(higher_levels), function(i) {
        hg <- higher_levels[i]
        
        # Clades that belong to this higher group
        clades_in_group <- sort(unique(higher_class$Clade[higher_class$Higher.Classification == hg]))
        
        # Pretty labels (capitalized) but values stay lowercase
        label_vec <- stats::setNames(clades_in_group, tools::toTitleCase(clades_in_group))
        
        # First group open by default, others collapsed
        tags$details(
          if (i == 1) list(open = "") else list(),
          tags$summary(hg),
          checkboxGroupInput(
            inputId = paste0("clades_", gsub("\\s+", "_", tolower(hg))),  # e.g., "clades_angiosperm"
            label  = NULL,   # label is handled by <summary>
            choices = label_vec
          )
        )
      })
    )
  })
  
  # Collect selected clades across *all* higher groups
  selected_clades <- reactive({
    unlist(
      lapply(higher_levels, function(hg) {
        input[[paste0("clades_", gsub("\\s+", "_", tolower(hg)))]]
      }),
      use.names = FALSE
    )
  })
  
  # Reactive filtered data
  filtered_data <- reactive({
    clades <- selected_clades()
    req(!is.null(clades), length(clades) > 0)
    all_chrome[all_chrome$clade %in% clades, , drop = FALSE]
  })
  
  # ---- Record counters ----
  output$totalRecords <- renderText({
    paste("Total records in dataset:", nrow(all_chrome))
  })
  
  output$filteredRecords <- renderText({
    clades <- selected_clades()
    if (is.null(clades) || !length(clades)) {
      return("Records shown: 0")
    }
    paste("Records shown:", nrow(filtered_data()))
  })
  
  # ---- Plot densities ----
  output$distPlot <- renderPlot({
    
    clades <- selected_clades()
    if (is.null(clades) || !length(clades)) {
      plot.new()
      text(0.5, 0.5, "Select at least one clade to show distributions")
      return()
    }
    
    df <- filtered_data()
    n  <- length(clades)
    
    # All haploid values across selected clades
    all_vals <- df$haploid[!is.na(df$haploid)]
    if (!length(all_vals)) {
      plot.new()
      text(0.5, 0.5, "No haploid data available for selected clades")
      return()
    }
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
    
    # Legend: show pretty labels (capitalized)
    legend_labels <- tools::toTitleCase(clades)
    legend("topright",
           legend = legend_labels,
           col = base_cols,
           lwd = 3,
           cex = 0.9,
           bty = "n")
  })
  
  # ---- Table ----
  output$cladeTable <- renderTable({
    clades <- selected_clades()
    if (is.null(clades) || !length(clades)) {
      return(NULL)
    }
    
    df <- filtered_data()
    
    # Show haploid as integers (no extra digits) in the table only
    if ("haploid" %in% names(df)) {
      df$haploid <- as.integer(df$haploid)
    }
    
    df
  })
  
  # ---- Download handler ----
  output$downloadData <- downloadHandler(
    filename = function() {
      clades <- selected_clades()
      if (is.null(clades) || !length(clades)) {
        "chromosome_data.csv"
      } else {
        paste0("chromosome_data_", paste(clades, collapse = "_"), ".csv")
      }
    },
    content = function(file) {
      clades <- selected_clades()
      if (is.null(clades) || !length(clades)) {
        # Write empty file if nothing selected
        write.csv(all_chrome[0, ], file, row.names = FALSE)
      } else {
        df <- filtered_data()
        write.csv(df, file, row.names = FALSE)
      }
    }
  )
}

# ---- Run app ----
shinyApp(ui = ui, server = server)
