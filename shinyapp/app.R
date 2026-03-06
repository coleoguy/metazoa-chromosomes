## ============================
##  Libraries
## ============================
library(shiny)
library(viridisLite)
library(bslib)
library(DT)
library(shinycssloaders)

## ============================
##  Load data
## ============================
all_chrome <- read.csv("all_chrome.csv", stringsAsFactors = FALSE)
all_chrome$haploid <- as.numeric(all_chrome$haploid)
all_chrome$clade   <- tolower(all_chrome$clade)

higher_class <- read.csv("higher_class.csv", stringsAsFactors = FALSE)
higher_class$Clade <- tolower(higher_class$Clade)

# Keep only valid clades
higher_class <- higher_class[higher_class$Clade %in% all_chrome$clade, ]
higher_levels <- unique(higher_class$Higher.Classification)

## ============================
##  UI
## ============================
ui <- bslib::page_sidebar(
  theme = bslib::bs_theme(
    version    = 5,
    bootswatch = "flatly",
    primary    = "#7FA58C",   # sage green
    secondary  = "#E6EFEA"
  ),
  title = "Chromosome Data Across the Tree of Life",
  
  sidebar = bslib::sidebar(
    width = 350,
    
    # Record counters
    bslib::card(
      bslib::card_body(
        strong(textOutput("totalRecords")),
        strong(textOutput("filteredRecords"))
      )
    ),
    br(),
    
    # Clade selectors
    h5("Select clades"),
    div(class = "clade-scroll", uiOutput("cladeSelector"))
  ),
  
  bslib::navset_card_tab(
    bslib::nav_panel(
      "Distributions",
      bslib::card_body(
        fluidRow(
          column(
            5,
            radioButtons(
              "plotType",
              NULL,
              choices = c("Density" = "density",
                          "Histogram" = "hist"),
              inline = TRUE
            )
          ),
          column(
            7,
            fluidRow(
              column(6, checkboxInput("showMedian", "Show median lines", TRUE)),
              column(6, checkboxInput("showN", "Show n in legend", TRUE))
            )
          )
        ),
        shinycssloaders::withSpinner(
          plotOutput("distPlot", height = "65vh")
        )
      )
    ),
    bslib::nav_panel(
      "Data table",
      bslib::card_body(
        downloadButton("downloadData", "Download CSV"),
        br(), br(),
        shinycssloaders::withSpinner(DT::DTOutput("cladeTable"))
      )
    )
  ),
  
  ## ---- CSS polish ----
  tags$head(
    tags$style(HTML("
      .clade-scroll {
        max-height: 520px;
        overflow-y: auto;
        padding-right: 6px;
        border: 1px solid #e5e5e5;
        border-radius: 8px;
        padding: 8px;
        background: #fff;
      }

      details > summary {
        font-weight: 600;
        margin-top: 8px;
        cursor: pointer;
        list-style: none;
      }

      details > summary::-webkit-details-marker {
        display: none;
      }

      details > summary::before {
        content: '▸ ';
      }

      details[open] > summary::before {
        content: '▾ ';
      }

      .card {
        border-color: #D6E2DA;
      }

      .card-header {
        background-color: #EEF4F1;
      }

      .card-body {
        overflow: hidden;
      }

      /* spacing between Density / Histogram radio buttons */
      .shiny-options-group label {
        margin-right: 18px;
      }
    "))
  )
)

## ============================
##  SERVER
## ============================
server <- function(input, output, session) {
  
  ## ---- Dynamic grouped clade selectors ----
  output$cladeSelector <- renderUI({
    tagList(
      lapply(seq_along(higher_levels), function(i) {
        hg <- higher_levels[i]
        
        clades_in_group <- sort(
          unique(higher_class$Clade[higher_class$Higher.Classification == hg])
        )
        
        label_vec <- stats::setNames(
          clades_in_group,
          tools::toTitleCase(clades_in_group)
        )
        
        tags$details(
          if (i == 1) list(open = "") else list(),
          tags$summary(hg),
          checkboxGroupInput(
            inputId = paste0("clades_", gsub("\\s+", "_", tolower(hg))),
            label   = NULL,
            choices = label_vec
          )
        )
      })
    )
  })
  
  ## ---- Collect selected clades ----
  selected_clades <- reactive({
    out <- unlist(
      lapply(higher_levels, function(hg) {
        input[[paste0("clades_", gsub("\\s+", "_", tolower(hg)))]]
      }),
      use.names = FALSE
    )
    out[!is.na(out)]
  })
  
  ## ---- Filtered data ----
  filtered_data <- reactive({
    clades <- selected_clades()
    req(length(clades) > 0)
    all_chrome[all_chrome$clade %in% clades, , drop = FALSE]
  })
  
  ## ---- Record counters ----
  output$totalRecords <- renderText({
    paste("Total records in dataset:", nrow(all_chrome))
  })
  
  output$filteredRecords <- renderText({
    clades <- selected_clades()
    if (!length(clades)) return("Records shown: 0")
    paste("Records shown:", nrow(filtered_data()))
  })
  
  ## ---- Plot ----
  output$distPlot <- renderPlot({
    par(mar = c(4, 4, 3, 1))
    
    clades <- selected_clades()
    if (!length(clades)) {
      plot.new()
      text(0.5, 0.5, "Select at least one clade to show distributions")
      return()
    }
    
    df <- filtered_data()
    all_vals <- df$haploid[!is.na(df$haploid)]
    if (!length(all_vals)) {
      plot.new()
      text(0.5, 0.5, "No haploid data available")
      return()
    }
    
    if (input$plotType == "hist") {
      hist(all_vals, breaks = "FD",
           main = "Haploid chromosome number",
           xlab = "Haploid chromosome number")
      return()
    }
    
    xlim <- range(all_vals)
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
    base_cols <- viridisLite::viridis(length(clades), option = "C")
    fill_cols <- sapply(base_cols, adjustcolor, alpha.f = 0.25)
    
    plot(0, type = "n", xlim = xlim, ylim = ylim,
         xlab = "Haploid chromosome number",
         ylab = "Density",
         main = "Haploid chromosome number distributions")
    
    for (i in seq_along(clades)) {
      d <- density_list[[clades[i]]]
      if (!is.null(d)) {
        polygon(c(d$x, rev(d$x)),
                c(d$y, rep(0, length(d$y))),
                col = fill_cols[i], border = NA)
        lines(d, col = base_cols[i], lwd = 2)
      }
    }
    
    if (isTRUE(input$showMedian)) {
      for (i in seq_along(clades)) {
        x <- df$haploid[df$clade == clades[i] & !is.na(df$haploid)]
        if (length(x)) {
          abline(v = median(x), col = base_cols[i], lty = 2, lwd = 2)
        }
      }
    }
    
    legend_labels <- tools::toTitleCase(clades)
    if (isTRUE(input$showN)) {
      ns <- sapply(clades, function(cl)
        sum(df$clade == cl & !is.na(df$haploid)))
      legend_labels <- paste0(legend_labels, " (n=", ns, ")")
    }
    
    legend("topright",
           legend = legend_labels,
           col = base_cols,
           lwd = 3,
           bty = "n")
  })
  
  ## ---- Data table ----
  output$cladeTable <- DT::renderDT({
    clades <- selected_clades()
    if (!length(clades)) return(NULL)
    
    df <- filtered_data()
    df$haploid <- as.integer(df$haploid)
    
    DT::datatable(
      df,
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  ## ---- Download ----
  output$downloadData <- downloadHandler(
    filename = function() {
      clades <- selected_clades()
      if (!length(clades)) {
        "chromosome_data.csv"
      } else {
        paste0("chromosome_data_", paste(clades, collapse = "_"), ".csv")
      }
    },
    content = function(file) {
      clades <- selected_clades()
      if (!length(clades)) {
        write.csv(all_chrome[0, ], file, row.names = FALSE)
      } else {
        write.csv(filtered_data(), file, row.names = FALSE)
      }
    }
  )
}

## ============================
##  Run app
## ============================
shinyApp(ui = ui, server = server)
