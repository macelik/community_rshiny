# community Shiny App

A local-only RShiny application that wraps the **community** R package workflow for cell-cell communication analysis from single-cell RNA-seq data.

The app mirrors the exact pipeline from the [Lasry tutorial script](https://github.com/SoloveyMaria/community/blob/main/docs/showcase_notebooks/Lasry/calculate_communication.r).

---

## Prerequisites

- **R ≥ 4.2**
- **RStudio** (recommended) or a terminal with `Rscript`

## Installation

### 1. Install the community package from GitHub

```r
install.packages("devtools")
devtools::install_github("SoloveyMaria/community")
```

### 2. Install CRAN dependencies

```r
install.packages(c(
  "shiny",
  "DT",
  "data.table",
  "ggplot2"
))
```

> **Note:** The `community` package itself may pull in additional dependencies (e.g., `tidyverse` components). Follow any install prompts.

### 3. Verify installation

```r
library(community)
library(shiny)
library(DT)
library(data.table)
library(ggplot2)
```

---

## How to Run

### Option A: From RStudio

1. Open RStudio
2. Set your working directory to the app folder:  
   `setwd("/path/to/community_shiny_app")`
3. Run:  
   ```r
   shiny::runApp()
   ```

### Option B: From terminal

```bash
cd /path/to/community_shiny_app
Rscript -e "shiny::runApp('.', port = 3838)"
```

Then open `http://localhost:3838` in your browser.

---

## App Walkthrough

### Tab 1 — Inputs

Upload the three required files:

| File | Format | Description |
|------|--------|-------------|
| **Counts matrix** | `.csv` or `.csv.gz` | First column = `gene_symbol`; remaining columns = cell barcodes with expression values |
| **Cell annotation** | `.txt`, `.tsv`, or `.csv` | Tab-separated (for `.txt`/`.tsv`) or comma-separated. Must contain a cell barcode column (default: `cell_ID`) |
| **Sample annotation** | `.txt`, `.tsv`, or `.csv` | Tab-separated or comma-separated. Must contain sample-level metadata expected by the community package |

**Demo data** can be downloaded from:  
<https://github.com/SoloveyMaria/community/tree/main/docs/showcase_notebooks/Lasry/input_data>

Optionally upload a **custom LR database** (`.rds`, `.RData`, or `.csv`). It must contain columns: `Pair.Name`, `Ligand`, `Receptor`. If not provided, the package default `LR_database` is used.

### Tab 2 — Parameters

All parameters from the tutorial are exposed with their defaults:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `threshold_celltype_size` | 6 | Min cells of a type per sample |
| `threshold_nr_active_cells` | 6 | Min number of active cells |
| `threshold_expr` | 0.1 | Expression threshold |
| `threshold_log10_cum_weight` | 0.01 | Filter: cumulative weight |
| `threshold_frac_samples_per_condition` | 0.6 | Filter: fraction of samples |
| `threshold_log10_meanexpr_per_condition` | 0.02 | Filter: mean expression |
| `threshold_log2FC` | 1 | Differential: fold-change |
| `threshold_fdr` | 0.1 | Differential: FDR |

You can also select which column in the cell annotation contains cell barcodes.

### Tab 3 — Run

Click **"Run Analysis"** to execute the full pipeline:

1. `calculate_communication()`
2. `general_stat()`
3. `filter_interactions()`
4. `test_diff()` (t-test)
5. `interaction_classes()`

Progress is shown in real time. Errors are displayed as notifications.

### Tab 4 — Results

Browse and download:

- `anno_interactions` table (interactive DT with column filters)
- `weights` matrix
- Full `interactions` object as `.rds`

### Tab 5 — Filter & Plot

Adjust filtering thresholds and **re-run only** `filter_interactions` → `test_diff` → `interaction_classes` without recomputing the expensive `calculate_communication` step.

A summary plot is generated automatically (interaction class bar chart, volcano plot, or ranked bar chart depending on available columns). Download as PNG or PDF.

---

## File Structure

```
community_shiny_app/
├── app.R              # Main Shiny app (UI + server)
├── R/
│   └── helpers.R      # Pipeline functions, file readers, plotting
└── README.md          # This file
```

---

## Notes

- The app runs **entirely locally** — no cloud services or external APIs.
- Outputs are written to R temp directories by default (via download handlers).
- The `cell_ID.1` → `cell_ID` rename from the tutorial is handled automatically.
- If the community package's functions throw errors about missing columns, the app displays the error message to help you diagnose the issue.

---

## Troubleshooting

**"Error loading counts: Counts file must have a 'gene_symbol' column"**  
→ Ensure the first column of your counts CSV is named `gene_symbol`.

**"Dimension mismatch"**  
→ The number of columns in the counts matrix (after removing `gene_symbol`) must equal the number of rows in cell annotation.

**Package installation fails**  
→ Try: `devtools::install_github("SoloveyMaria/community", dependencies = TRUE, force = TRUE)`

---

## License

This Shiny app is a wrapper around the [community](https://github.com/SoloveyMaria/community) package (MIT license).
