# ============================================================================
# app.R — Shiny app wrapping the community cell-cell communication pipeline
# ============================================================================
# Run with: shiny::runApp("community_shiny_app")
# ============================================================================

library(shiny)
library(DT)

# Remove upload size limit (default is 5MB, too small for scRNA-seq data)
options(shiny.maxRequestSize = Inf)

source(file.path("R", "helpers.R"))

# ============================================================================
# UI
# ============================================================================

ui <- navbarPage(
  title = "community — Cell-Cell Communication",
  id    = "main_tabs",
  theme = NULL,
  
  # ==========================================================================
  # TAB 1 — Inputs
  # ==========================================================================
  tabPanel("1. Inputs",
    fluidRow(
      column(4,
        wellPanel(
          h4("Upload Input Files"),
          fileInput("file_counts", "Counts matrix (.csv, .csv.gz)",
                    accept = c(".csv", ".csv.gz", ".gz")),
          fileInput("file_anno_cells", "Cell annotation (.txt, .tsv, .csv)",
                    accept = c(".txt", ".tsv", ".csv")),
          fileInput("file_anno_samples", "Sample annotation (.txt, .tsv, .csv)",
                    accept = c(".txt", ".tsv", ".csv")),
          hr(),
          h4("Ligand-Receptor Database"),
          radioButtons("db_choice", NULL,
                       choices = c("Use default LR_database" = "default",
                                   "Upload custom database"  = "custom"),
                       selected = "default"),
          conditionalPanel(
            condition = "input.db_choice == 'custom'",
            fileInput("file_custom_db",
                      "Custom DB (.rds, .RData, .csv)",
                      accept = c(".rds", ".RData", ".rda", ".csv")),
            textInput("db_object_name",
                      "Object name in .RData (leave blank for auto-detect)",
                      value = "")
          )
        )
      ),
      column(8,
        h4("File Previews"),
        tabsetPanel(
          tabPanel("Counts",
            verbatimTextOutput("counts_dim"),
            DT::dataTableOutput("preview_counts")
          ),
          tabPanel("Cell Annotation",
            DT::dataTableOutput("preview_anno_cells")
          ),
          tabPanel("Sample Annotation",
            DT::dataTableOutput("preview_anno_samples")
          ),
          tabPanel("LR Database",
            DT::dataTableOutput("preview_database")
          )
        )
      )
    )
  ),
  
  # ==========================================================================
  # TAB 2 — Parameters
  # ==========================================================================
  tabPanel("2. Parameters",
    fluidRow(
      column(4,
        wellPanel(
          h4("Cell Barcode Column"),
          uiOutput("barcode_col_ui"),
          helpText("Select the column in cell annotation that contains cell barcodes.",
                   "If 'cell_ID.1' is present it will be auto-renamed to 'cell_ID'.")
        )
      ),
      column(4,
        wellPanel(
          h4("Communication Thresholds"),
          numericInput("threshold_celltype_size",
                       "Min cell-type size per sample",
                       value = 6, min = 1, step = 1),
          numericInput("threshold_nr_active_cells",
                       "Min active cells",
                       value = 6, min = 1, step = 1),
          numericInput("threshold_expr",
                       "Expression threshold",
                       value = 0.1, min = 0, step = 0.01)
        )
      ),
      column(4,
        wellPanel(
          h4("Filtering & Differential Thresholds"),
          numericInput("threshold_log10_cum_weight",
                       "log10 cumulative weight",
                       value = 0.01, min = 0, step = 0.001),
          numericInput("threshold_frac_samples",
                       "Fraction samples per condition",
                       value = 0.6, min = 0, max = 1, step = 0.05),
          numericInput("threshold_log10_meanexpr",
                       "log10 mean expression per condition",
                       value = 0.02, min = 0, step = 0.001),
          numericInput("threshold_log2FC",
                       "log2 fold-change",
                       value = 1, min = 0, step = 0.1),
          numericInput("threshold_fdr",
                       "FDR threshold",
                       value = 0.1, min = 0, max = 1, step = 0.01)
        )
      )
    )
  ),
  
  # ==========================================================================
  # TAB 3 — Run
  # ==========================================================================
  tabPanel("3. Run",
    fluidRow(
      column(6, offset = 3,
        wellPanel(
          h4("Run Full Pipeline"),
          helpText("This will execute the complete community workflow:",
                   "calculate_communication → general_stat → filter_interactions",
                   " → test_diff → interaction_classes"),
          actionButton("btn_run", "Run Analysis",
                       class = "btn-primary btn-lg",
                       width = "100%"),
          hr(),
          verbatimTextOutput("run_status")
        )
      )
    )
  ),
  
  # ==========================================================================
  # TAB 4 — Results
  # ==========================================================================
  tabPanel("4. Results",
    fluidRow(
      column(12,
        conditionalPanel(
          condition = "output.has_results",
          h4("Thresholds Used"),
          verbatimTextOutput("result_thresholds"),
          hr(),
          tabsetPanel(
            tabPanel("anno_interactions",
              DT::dataTableOutput("tbl_anno_interactions"),
              downloadButton("dl_anno_interactions", "Download CSV")
            ),
            tabPanel("weights",
              DT::dataTableOutput("tbl_weights"),
              downloadButton("dl_weights", "Download CSV")
            ),
            tabPanel("Full Object",
              helpText("Download the complete interactions object as .rds"),
              downloadButton("dl_rds", "Download interactions.rds")
            )
          )
        ),
        conditionalPanel(
          condition = "!output.has_results",
          div(style = "padding:40px; text-align:center;",
            h4("No results yet. Go to the 'Run' tab to execute the pipeline.")
          )
        )
      )
    )
  ),
  
  # ==========================================================================
  # TAB 5 — Filter + Plot
  # ==========================================================================
  tabPanel("5. Filter & Plot",
    fluidRow(
      column(4,
        wellPanel(
          h4("Re-filter Parameters"),
          helpText("Adjust filters and re-run only filter_interactions →",
                   "test_diff → interaction_classes (skips the expensive",
                   "calculate_communication step)."),
          numericInput("refilter_log10_cum_weight",
                       "log10 cumulative weight",
                       value = 0.01, min = 0, step = 0.001),
          numericInput("refilter_frac_samples",
                       "Fraction samples per condition",
                       value = 0.6, min = 0, max = 1, step = 0.05),
          numericInput("refilter_log10_meanexpr",
                       "log10 mean expression per condition",
                       value = 0.02, min = 0, step = 0.001),
          numericInput("refilter_log2FC",
                       "log2 fold-change",
                       value = 1, min = 0, step = 0.1),
          numericInput("refilter_fdr",
                       "FDR threshold",
                       value = 0.1, min = 0, max = 1, step = 0.01),
          actionButton("btn_refilter", "Re-filter & Plot",
                       class = "btn-info btn-lg", width = "100%")
        )
      ),
      column(8,
        conditionalPanel(
          condition = "output.has_results",
          h4("Interaction Summary Plot"),
          plotOutput("interaction_plot", height = "500px"),
          hr(),
          fluidRow(
            column(3, downloadButton("dl_plot_png", "Download PNG")),
            column(3, downloadButton("dl_plot_pdf", "Download PDF"))
          ),
          hr(),
          h4("Filtered anno_interactions"),
          DT::dataTableOutput("tbl_filtered_anno")
        ),
        conditionalPanel(
          condition = "!output.has_results",
          div(style = "padding:40px; text-align:center;",
            h4("Run the pipeline first (Tab 3).")
          )
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {
  
  # Reactive values ----------------------------------------------------------
  rv <- reactiveValues(
    counts        = NULL,
    anno_cells    = NULL,
    anno_samples  = NULL,
    lr_database   = NULL,
    interactions  = NULL,       # full result after pipeline
    interactions_base = NULL,   # result after general_stat (for re-filtering)
    run_log       = "",
    current_plot  = NULL
  )
  
  # --------------------------------------------------------------------------
  # File loading reactives
  # --------------------------------------------------------------------------
  
  observeEvent(input$file_counts, {
    req(input$file_counts)
    tryCatch({
      rv$counts <- read_counts(input$file_counts$datapath)
      rv$run_log <- paste0(rv$run_log, "\n✓ Counts loaded: ",
                           nrow(rv$counts), " genes × ",
                           ncol(rv$counts), " cells")
    }, error = function(e) {
      showNotification(paste("Error loading counts:", e$message),
                       type = "error", duration = 10)
    })
  })
  
  observeEvent(input$file_anno_cells, {
    req(input$file_anno_cells)
    tryCatch({
      # Read with default barcode col first; user can re-select later
      rv$anno_cells <- read_anno_cells(input$file_anno_cells$datapath,
                                        barcode_col = "cell_ID")
      rv$run_log <- paste0(rv$run_log, "\n✓ Cell annotation loaded: ",
                           nrow(rv$anno_cells), " cells, ",
                           ncol(rv$anno_cells), " columns")
    }, error = function(e) {
      # Try reading without barcode enforcement to at least get columns
      tryCatch({
        ext <- tolower(tools::file_ext(input$file_anno_cells$name))
        sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
        rv$anno_cells <- read.table(input$file_anno_cells$datapath,
                                     sep = sep, header = TRUE,
                                     stringsAsFactors = FALSE,
                                     check.names = FALSE)
        showNotification(
          paste("Loaded cell annotation but barcode column needs attention:",
                e$message),
          type = "warning", duration = 10)
      }, error = function(e2) {
        showNotification(paste("Error loading cell annotation:", e2$message),
                         type = "error", duration = 10)
      })
    })
  })
  
  observeEvent(input$file_anno_samples, {
    req(input$file_anno_samples)
    tryCatch({
      rv$anno_samples <- read_anno_samples(input$file_anno_samples$datapath)
      rv$run_log <- paste0(rv$run_log, "\n✓ Sample annotation loaded: ",
                           nrow(rv$anno_samples), " samples, ",
                           ncol(rv$anno_samples), " columns")
    }, error = function(e) {
      showNotification(paste("Error loading sample annotation:", e$message),
                       type = "error", duration = 10)
    })
  })
  
  # LR database
  observe({
    if (input$db_choice == "default") {
      tryCatch({
        rv$lr_database <- load_lr_database(use_default = TRUE)
      }, error = function(e) {
        showNotification(paste("Error loading default database:", e$message),
                         type = "error", duration = 10)
      })
    }
  })
  
  observeEvent(input$file_custom_db, {
    req(input$file_custom_db, input$db_choice == "custom")
    tryCatch({
      obj_name <- if (nzchar(input$db_object_name)) input$db_object_name else NULL
      rv$lr_database <- load_lr_database(use_default = FALSE,
                                          filepath = input$file_custom_db$datapath,
                                          object_name = obj_name)
      showNotification("Custom database loaded successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading custom database:", e$message),
                       type = "error", duration = 10)
    })
  })
  
  # --------------------------------------------------------------------------
  # Dynamic UI: barcode column selector
  # --------------------------------------------------------------------------
  
  output$barcode_col_ui <- renderUI({
    cols <- if (!is.null(rv$anno_cells)) colnames(rv$anno_cells) else "cell_ID"
    selected <- if ("cell_ID" %in% cols) "cell_ID"
                else if ("cell_ID.1" %in% cols) "cell_ID.1"
                else cols[1]
    selectInput("barcode_col", "Barcode column", choices = cols,
                selected = selected)
  })
  
  # --------------------------------------------------------------------------
  # File previews
  # --------------------------------------------------------------------------
  
  output$counts_dim <- renderText({
    req(rv$counts)
    paste0("Dimensions: ", nrow(rv$counts), " genes × ", ncol(rv$counts), " cells")
  })
  
  output$preview_counts <- DT::renderDataTable({
    req(rv$counts)
    # Show first 10 rows and first 10 columns
    nr <- min(10, nrow(rv$counts))
    nc <- min(10, ncol(rv$counts))
    df <- rv$counts[1:nr, 1:nc, drop = FALSE]
    df <- cbind(gene_symbol = rownames(df), df)
    DT::datatable(df, options = list(scrollX = TRUE, pageLength = 10),
                  rownames = FALSE)
  })
  
  output$preview_anno_cells <- DT::renderDataTable({
    req(rv$anno_cells)
    DT::datatable(head(rv$anno_cells, 10),
                  options = list(scrollX = TRUE, pageLength = 10),
                  rownames = FALSE)
  })
  
  output$preview_anno_samples <- DT::renderDataTable({
    req(rv$anno_samples)
    DT::datatable(head(rv$anno_samples, 10),
                  options = list(scrollX = TRUE, pageLength = 10),
                  rownames = FALSE)
  })
  
  output$preview_database <- DT::renderDataTable({
    req(rv$lr_database)
    DT::datatable(head(rv$lr_database, 20),
                  options = list(scrollX = TRUE, pageLength = 10),
                  rownames = FALSE)
  })
  
  # --------------------------------------------------------------------------
  # Run pipeline
  # --------------------------------------------------------------------------
  
  observeEvent(input$btn_run, {
    
    # ---- Validation ----
    errors <- c()
    
    if (is.null(rv$counts))
      errors <- c(errors, "Counts matrix not uploaded.")
    if (is.null(rv$anno_cells))
      errors <- c(errors, "Cell annotation not uploaded.")
    if (is.null(rv$anno_samples))
      errors <- c(errors, "Sample annotation not uploaded.")
    if (is.null(rv$lr_database))
      errors <- c(errors, "LR database not loaded.")
    
    if (length(errors) > 0) {
      showNotification(paste(errors, collapse = "\n"),
                       type = "error", duration = 10)
      return()
    }
    
    # Re-read anno_cells with selected barcode column
    barcode_col <- if (!is.null(input$barcode_col)) input$barcode_col else "cell_ID"
    tryCatch({
      anno_cells_prep <- read_anno_cells(input$file_anno_cells$datapath,
                                          barcode_col = barcode_col)
    }, error = function(e) {
      showNotification(paste("Barcode column error:", e$message),
                       type = "error", duration = 10)
      return()
    })
    
    # Dimension check
    n_cells_counts <- ncol(rv$counts)
    n_cells_anno   <- nrow(anno_cells_prep)
    if (n_cells_counts != n_cells_anno) {
      showNotification(
        paste0("Dimension mismatch! Counts has ", n_cells_counts,
               " columns but cell annotation has ", n_cells_anno, " rows. ",
               "These must match."),
        type = "error", duration = 15)
      return()
    }
    
    # ---- Run ----
    rv$run_log <- "Starting pipeline..."
    
    withProgress(message = "Running community pipeline", value = 0, {
      tryCatch({
        result <- run_pipeline(
          counts       = rv$counts,
          anno_cells   = anno_cells_prep,
          anno_samples = rv$anno_samples,
          lr_database  = rv$lr_database,
          threshold_celltype_size   = input$threshold_celltype_size,
          threshold_nr_active_cells = input$threshold_nr_active_cells,
          threshold_expr            = input$threshold_expr,
          threshold_frac_samples_per_condition  = input$threshold_frac_samples,
          threshold_log10_cum_weight            = input$threshold_log10_cum_weight,
          threshold_log10_meanexpr_per_condition = input$threshold_log10_meanexpr,
          threshold_log2FC = input$threshold_log2FC,
          threshold_fdr    = input$threshold_fdr,
          progress_fn = function(val, msg) {
            setProgress(value = val, message = msg)
          }
        )
        
        rv$interactions      <- result
        rv$interactions_base <- result  # snapshot for re-filtering
        rv$run_log <- paste0("Pipeline completed successfully!\n",
                             "anno_interactions: ",
                             nrow(result$anno_interactions), " rows\n",
                             "weights matrix: ",
                             nrow(result$weights), " × ",
                             ncol(result$weights))
        
        # Build initial plot
        tryCatch({
          rv$current_plot <- plot_interaction_summary(result$anno_interactions)
        }, error = function(e) {
          rv$run_log <- paste0(rv$run_log, "\nPlot warning: ", e$message)
        })
        
        # Sync refilter inputs
        updateNumericInput(session, "refilter_log10_cum_weight",
                           value = input$threshold_log10_cum_weight)
        updateNumericInput(session, "refilter_frac_samples",
                           value = input$threshold_frac_samples)
        updateNumericInput(session, "refilter_log10_meanexpr",
                           value = input$threshold_log10_meanexpr)
        updateNumericInput(session, "refilter_log2FC",
                           value = input$threshold_log2FC)
        updateNumericInput(session, "refilter_fdr",
                           value = input$threshold_fdr)
        
        showNotification("Pipeline completed!", type = "message")
        
      }, error = function(e) {
        rv$run_log <- paste0("ERROR: ", e$message)
        showNotification(paste("Pipeline error:", e$message),
                         type = "error", duration = 15)
      })
    })
  })
  
  output$run_status <- renderText({ rv$run_log })
  
  # --------------------------------------------------------------------------
  # Results display
  # --------------------------------------------------------------------------
  
  output$has_results <- reactive({ !is.null(rv$interactions) })
  outputOptions(output, "has_results", suspendWhenHidden = FALSE)
  
  output$result_thresholds <- renderText({
    req(rv$interactions)
    th <- rv$interactions$thresholds
    if (is.null(th)) return("(thresholds slot not found)")
    paste(capture.output(str(th)), collapse = "\n")
  })
  
  output$tbl_anno_interactions <- DT::renderDataTable({
    req(rv$interactions)
    DT::datatable(rv$interactions$anno_interactions,
                  options = list(scrollX = TRUE, pageLength = 15),
                  rownames = FALSE, filter = "top")
  })
  
  output$tbl_weights <- DT::renderDataTable({
    req(rv$interactions)
    w <- rv$interactions$weights
    if (!is.null(w)) {
      # Weights can be very wide; show first 50 cols
      nc <- min(50, ncol(w))
      df <- as.data.frame(w[, 1:nc, drop = FALSE])
      df <- cbind(row = rownames(df), df)
      DT::datatable(df, options = list(scrollX = TRUE, pageLength = 15),
                    rownames = FALSE)
    }
  })
  
  # Downloads ----------------------------------------------------------------
  
  output$dl_anno_interactions <- downloadHandler(
    filename = function() "community_anno_interactions.csv",
    content  = function(file) {
      write.csv(rv$interactions$anno_interactions, file, row.names = FALSE)
    }
  )
  
  output$dl_weights <- downloadHandler(
    filename = function() "community_weights.csv",
    content  = function(file) {
      write.csv(rv$interactions$weights, file)
    }
  )
  
  output$dl_rds <- downloadHandler(
    filename = function() "interactions.rds",
    content  = function(file) {
      saveRDS(rv$interactions, file)
    }
  )
  
  # --------------------------------------------------------------------------
  # Filter + Plot tab
  # --------------------------------------------------------------------------
  
  observeEvent(input$btn_refilter, {
    req(rv$interactions_base)
    
    withProgress(message = "Re-filtering...", value = 0, {
      tryCatch({
        result <- refilter_pipeline(
          interactions = rv$interactions_base,
          threshold_frac_samples_per_condition  = input$refilter_frac_samples,
          threshold_log10_cum_weight            = input$refilter_log10_cum_weight,
          threshold_log10_meanexpr_per_condition = input$refilter_log10_meanexpr,
          threshold_log2FC = input$refilter_log2FC,
          threshold_fdr    = input$refilter_fdr,
          progress_fn = function(val, msg) setProgress(value = val, message = msg)
        )
        rv$interactions <- result
        
        tryCatch({
          rv$current_plot <- plot_interaction_summary(result$anno_interactions)
        }, error = function(e) {
          showNotification(paste("Plot error:", e$message),
                           type = "warning", duration = 8)
        })
        
        showNotification(
          paste0("Re-filtered: ", nrow(result$anno_interactions),
                 " interactions remaining."),
          type = "message")
        
      }, error = function(e) {
        showNotification(paste("Re-filter error:", e$message),
                         type = "error", duration = 15)
      })
    })
  })
  
  output$interaction_plot <- renderPlot({
    req(rv$current_plot)
    rv$current_plot
  })
  
  output$tbl_filtered_anno <- DT::renderDataTable({
    req(rv$interactions)
    DT::datatable(rv$interactions$anno_interactions,
                  options = list(scrollX = TRUE, pageLength = 15),
                  rownames = FALSE, filter = "top")
  })
  
  # Plot downloads
  output$dl_plot_png <- downloadHandler(
    filename = function() "interaction_summary.png",
    content  = function(file) {
      req(rv$current_plot)
      ggplot2::ggsave(file, plot = rv$current_plot,
                      width = 10, height = 7, dpi = 150)
    }
  )
  
  output$dl_plot_pdf <- downloadHandler(
    filename = function() "interaction_summary.pdf",
    content  = function(file) {
      req(rv$current_plot)
      ggplot2::ggsave(file, plot = rv$current_plot,
                      width = 10, height = 7)
    }
  )
}

# ============================================================================
# Launch
# ============================================================================
shinyApp(ui, server)
