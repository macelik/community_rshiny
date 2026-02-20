# ============================================================================
# R/helpers.R — Helper functions for the community Shiny app
# ============================================================================

# Required packages ---------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(community)
})

# ---------------------------------------------------------------------------
# 1. File readers
# ---------------------------------------------------------------------------

#' Read counts matrix (CSV or CSV.GZ)
#' Expects first column = gene_symbol; remaining columns = cell barcodes.
#' Returns a data.frame with gene_symbol as row names, first column dropped.
read_counts <- function(filepath) {
  counts <- data.table::fread(filepath, header = TRUE)
  counts <- as.data.frame(counts)
  
  # Validate gene_symbol column
  if (!"gene_symbol" %in% colnames(counts)) {
    if ("V1" %in% colnames(counts)) {
      colnames(counts)[1] <- "gene_symbol"
    } else {
      stop(
        "Counts file must have a 'gene_symbol' column (first column). ",
        "Found columns: ", paste(head(colnames(counts), 5), collapse = ", ")
      )
    }
  }
  
  rownames(counts) <- counts$gene_symbol
  counts <- counts[, -which(colnames(counts) == "gene_symbol"), drop = FALSE]
  counts
}

#' Read cell annotation table (tab/comma separated)
#' Returns a data.frame. Handles cell_ID.1 → cell_ID rename.
read_anno_cells <- function(filepath, barcode_col = "cell_ID") {
  ext <- tolower(tools::file_ext(filepath))
  
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  
  anno <- read.table(filepath, sep = sep, header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
  
  # Handle cell_ID.1 → cell_ID renaming (as in the tutorial)
  if ("cell_ID.1" %in% colnames(anno) && !"cell_ID" %in% colnames(anno)) {
    colnames(anno)[colnames(anno) == "cell_ID.1"] <- "cell_ID"
  }
  
  # If user-selected barcode column differs, rename it
 if (barcode_col != "cell_ID" && barcode_col %in% colnames(anno)) {
    # Keep original cell_ID if it exists
    colnames(anno)[colnames(anno) == barcode_col] <- "cell_ID"
  }
  
  if (!"cell_ID" %in% colnames(anno)) {
    stop(
      "Cell annotation must contain a cell barcode column. ",
      "Expected '", barcode_col, "' or 'cell_ID'. ",
      "Found: ", paste(head(colnames(anno), 10), collapse = ", ")
    )
  }
  
  rownames(anno) <- anno$cell_ID
  anno
}

#' Read sample annotation table (tab/comma separated)
read_anno_samples <- function(filepath) {
  ext <- tolower(tools::file_ext(filepath))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  
  anno <- read.table(filepath, sep = sep, header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
  anno
}

# ---------------------------------------------------------------------------
# 2. Database loader
# ---------------------------------------------------------------------------

#' Load LR database: either the package default or a user-uploaded file.
#' For .rds: expects the object itself to be a data.frame.
#' For .RData: loads all objects and picks the first data.frame found
#'   (or the one matching `object_name` if supplied).
#' For .csv: reads as data.frame.
load_lr_database <- function(use_default = TRUE, filepath = NULL,
                             object_name = NULL) {
  if (use_default) {
    data("LR_database", package = "community", envir = environment())
    return(LR_database)
  }
  
  if (is.null(filepath) || !file.exists(filepath)) {
    stop("Custom database file not found.")
  }
  
  ext <- tolower(tools::file_ext(filepath))
  
  if (ext == "rds") {
    db <- readRDS(filepath)
  } else if (ext == "rdata" || ext == "rda") {
    e <- new.env()
    load(filepath, envir = e)
    objs <- ls(e)
    if (!is.null(object_name) && object_name %in% objs) {
      db <- get(object_name, envir = e)
    } else {
      # Pick first data.frame-like object
      for (nm in objs) {
        obj <- get(nm, envir = e)
        if (is.data.frame(obj) || is.list(obj)) {
          db <- obj
          break
        }
      }
      if (!exists("db")) stop("No data.frame found in uploaded .RData file.")
    }
  } else if (ext == "csv") {
    db <- read.csv(filepath, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported database format: ", ext,
         ". Use .rds, .RData, or .csv")
  }
  
  # Basic validation: community expects Pair.Name, Ligand, Receptor
  required <- c("Pair.Name", "Ligand", "Receptor")
  missing <- setdiff(required, colnames(db))
  if (length(missing) > 0) {
    stop("Custom database is missing required columns: ",
         paste(missing, collapse = ", "),
         ". Required: Pair.Name, Ligand, Receptor")
  }
  
  db
}

# ---------------------------------------------------------------------------
# 3. Pipeline runner
# ---------------------------------------------------------------------------

#' Run the full community pipeline.
#' Returns the `interactions` object.
#' `progress_fn` is an optional callback: progress_fn(value, message)
run_pipeline <- function(counts, anno_cells, anno_samples, lr_database,
                         threshold_celltype_size = 6,
                         threshold_nr_active_cells = 6,
                         threshold_expr = 0.1,
                         threshold_frac_samples_per_condition = 0.6,
                         threshold_log10_cum_weight = 0.01,
                         threshold_log10_meanexpr_per_condition = 0.02,
                         threshold_log2FC = 1,
                         threshold_fdr = 0.1,
                         progress_fn = NULL) {
  
  pfn <- if (!is.null(progress_fn)) progress_fn else function(v, m) invisible(NULL)
  
  # Align counts columns with anno_cells
  colnames(counts) <- anno_cells$cell_ID
  
  # Step 1: calculate_communication
  pfn(0.15, "Running calculate_communication()...")
  interactions <- community::calculate_communication(
    counts                    = counts,
    anno_samples              = anno_samples,
    anno_cells                = anno_cells,
    threshold_celltype_size   = threshold_celltype_size,
    threshold_nr_active_cells = threshold_nr_active_cells,
    threshold_expr            = threshold_expr,
    lrp_database              = lr_database
  )
  
  # Step 2: general_stat
  pfn(0.40, "Computing general statistics...")
  interactions <- community::general_stat(
    comm_result = interactions,
    verbose     = FALSE
  )
  
  # Step 3: filter_interactions
  pfn(0.60, "Filtering weak interactions...")
  interactions <- community::filter_interactions(
    comm_result                           = interactions,
    threshold_frac_samples_per_condition  = threshold_frac_samples_per_condition,
    threshold_log10_cum_weight            = threshold_log10_cum_weight,
    threshold_log10_meanexpr_per_condition = threshold_log10_meanexpr_per_condition
  )
  
  # Step 4: test_diff
  pfn(0.80, "Testing differential communication...")
  interactions <- community::test_diff(
    comm_result    = interactions,
    threshold_fdr  = threshold_fdr,
    which_test     = "t-test",
    threshold_log2FC = threshold_log2FC
  )
  
  # Step 5: interaction_classes
  pfn(0.95, "Computing interaction classes...")
  interactions <- community::interaction_classes(
    interactions,
    threshold = threshold_log2FC
  )
  
  pfn(1.0, "Done!")
  interactions
}

# ---------------------------------------------------------------------------
# 4. Re-filter only (skip steps 1-2)
# ---------------------------------------------------------------------------

#' Re-run filter + test_diff + interaction_classes on an existing result.
#' Requires the interactions object from after general_stat().
refilter_pipeline <- function(interactions,
                              threshold_frac_samples_per_condition = 0.6,
                              threshold_log10_cum_weight = 0.01,
                              threshold_log10_meanexpr_per_condition = 0.02,
                              threshold_log2FC = 1,
                              threshold_fdr = 0.1,
                              progress_fn = NULL) {
  pfn <- if (!is.null(progress_fn)) progress_fn else function(v, m) invisible(NULL)
  
  pfn(0.3, "Re-filtering interactions...")
  interactions <- community::filter_interactions(
    comm_result                           = interactions,
    threshold_frac_samples_per_condition  = threshold_frac_samples_per_condition,
    threshold_log10_cum_weight            = threshold_log10_cum_weight,
    threshold_log10_meanexpr_per_condition = threshold_log10_meanexpr_per_condition
  )
  
  pfn(0.6, "Re-testing differential communication...")
  interactions <- community::test_diff(
    comm_result    = interactions,
    threshold_fdr  = threshold_fdr,
    which_test     = "t-test",
    threshold_log2FC = threshold_log2FC
  )
  
  pfn(0.9, "Re-computing interaction classes...")
  interactions <- community::interaction_classes(
    interactions,
    threshold = threshold_log2FC
  )
  
  pfn(1.0, "Done!")
  interactions
}

# ---------------------------------------------------------------------------
# 5. Summary plotting (ggplot2-based, no package dependency beyond ggplot2)
# ---------------------------------------------------------------------------

#' Build a summary bar plot of interaction classes from anno_interactions.
#' Falls back gracefully if expected columns are missing.
plot_interaction_summary <- function(anno_interactions) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting.")
  }
  
  library(ggplot2)
  
  # Try to detect a useful categorical column for coloring
  class_col <- NULL
  for (candidate in c("interaction_class", "class", "diff_class",
                       "interaction_type", "type")) {
    if (candidate %in% colnames(anno_interactions)) {
      class_col <- candidate
      break
    }
  }
  
  # --- Plot A: If we have a classification column, bar chart of counts ---
  if (!is.null(class_col)) {
    tbl <- as.data.frame(table(anno_interactions[[class_col]]))
    colnames(tbl) <- c("Class", "Count")
    tbl <- tbl[order(-tbl$Count), ]
    tbl$Class <- factor(tbl$Class, levels = tbl$Class)
    
    p <- ggplot(tbl, aes(x = Class, y = Count, fill = Class)) +
      geom_col(show.legend = FALSE) +
      theme_minimal(base_size = 14) +
      labs(title = "Interaction Classes", x = NULL, y = "Count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return(p)
  }
  
  # --- Plot B: If we have log2FC and FDR, volcano-like ---
  fc_col <- NULL
  fdr_col <- NULL
  for (c in c("log2FC", "log2fc", "logFC")) {
    if (c %in% colnames(anno_interactions)) { fc_col <- c; break }
  }
  for (c in c("fdr", "FDR", "p_adj", "padj", "adj.P.Val")) {
    if (c %in% colnames(anno_interactions)) { fdr_col <- c; break }
  }
  
  if (!is.null(fc_col) && !is.null(fdr_col)) {
    df <- anno_interactions
    df$neg_log10_fdr <- -log10(df[[fdr_col]] + 1e-300)
    
    p <- ggplot(df, aes(x = .data[[fc_col]], y = neg_log10_fdr)) +
      geom_point(alpha = 0.5, size = 1.2, color = "steelblue") +
      theme_minimal(base_size = 14) +
      labs(title = "Volcano Plot of Interactions",
           x = fc_col, y = paste0("-log10(", fdr_col, ")"))
    return(p)
  }
  
  # --- Plot C: Fallback – top interactions by first numeric column ---
  num_cols <- names(which(sapply(anno_interactions, is.numeric)))
  if (length(num_cols) == 0) {
    stop("No numeric columns found in anno_interactions for plotting.")
  }
  
  # Pick a weight-like column
  weight_col <- NULL
  for (c in c("log10_cum_weight", "cum_weight", "weight", num_cols[1])) {
    if (c %in% colnames(anno_interactions)) { weight_col <- c; break }
  }
  if (is.null(weight_col)) weight_col <- num_cols[1]
  
  # Label column
  label_col <- NULL
  for (c in c("interaction", "Pair.Name", "pair_name", "lrp")) {
    if (c %in% colnames(anno_interactions)) { label_col <- c; break }
  }
  
  df <- anno_interactions
  df <- df[order(-abs(df[[weight_col]])), ]
  df <- head(df, 30)
  
  if (!is.null(label_col)) {
    df$label <- df[[label_col]]
  } else {
    df$label <- paste0("int_", seq_len(nrow(df)))
  }
  df$label <- factor(df$label, levels = rev(df$label))
  
  p <- ggplot(df, aes(x = .data[[weight_col]], y = label)) +
    geom_col(fill = "steelblue") +
    theme_minimal(base_size = 12) +
    labs(title = paste0("Top Interactions by ", weight_col),
         x = weight_col, y = NULL)
  
  return(p)
}
