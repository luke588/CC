rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
options(stringsAsFactors = FALSE)

# =========================================================
# 00.A3_config.R
# A3 supporting evidence for HMMR virtual KO
# Fixed paths, theme, helper functions
# Final corrected version: rank-based tertile grouping
# =========================================================

# -----------------------------
# 1. Fixed base paths
# -----------------------------
base_dir <- "/Users/donglinlu/Desktop/CC1/9.scTenifoldKnk/05B.HMMR_KO_supporting_evidence"

epi_rdata_file <- "/Users/donglinlu/Desktop/CC1/9.scTenifoldKnk/Epithelial_subset.Rdata"

ko_base_dir <- "/Users/donglinlu/Desktop/CC1/9.scTenifoldKnk/05.HMMR_virtual_KO_functional_summary/tables"
ko_all_diff_file   <- file.path(ko_base_dir, "01.HMMR_KO_all_diffRegulation_cleaned.csv")
ko_sig_diff_file   <- file.path(ko_base_dir, "02.HMMR_KO_sigDiff_cleaned.csv")
ko_go_file         <- file.path(ko_base_dir, "03.HMMR_KO_GO_BP_results.csv")
ko_kegg_file       <- file.path(ko_base_dir, "04.HMMR_KO_KEGG_results.csv")
ko_focus_go_file   <- file.path(ko_base_dir, "05.HMMR_KO_GO_BP_focus_terms.csv")
ko_focus_kegg_file <- file.path(ko_base_dir, "06.HMMR_KO_KEGG_focus_terms.csv")

input_manifest_dir <- file.path(base_dir, "input_manifest")
intermediate_dir   <- file.path(base_dir, "intermediate")
table_dir          <- file.path(base_dir, "tables")
plot_dir           <- file.path(base_dir, "plots")

dir.create(base_dir,           showWarnings = FALSE, recursive = TRUE)
dir.create(input_manifest_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(intermediate_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir,          showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir,           showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Analysis parameters
# -----------------------------
hmmr_gene <- "HMMR"
target_assay <- "RNA"

target_clusters  <- c("6", "9", "10")
target_condition <- "CA_HPV"
state_method     <- "rank_tertile_all_cells"

min_pct_de       <- 0.10
logfc_threshold  <- 0.25
padj_threshold   <- 0.05

focus_modules <- c(
  "ribosome",
  "cytoplasmic translation",
  "translation",
  "peptide biosynthetic process",
  "ribosome biogenesis",
  "cell cycle"
)

# -----------------------------
# 3. Packages
# -----------------------------
cran_pkgs <- c(
  "ggplot2", "dplyr", "readr", "stringr", "forcats",
  "scales", "tibble", "patchwork", "tidyr"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat", dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(scales)
  library(tibble)
  library(patchwork)
  library(tidyr)
  library(Seurat)
})

# -----------------------------
# 4. Unified publication theme
# -----------------------------
col_text   <- "#222222"
col_border <- "#4A4A4A"
col_grid   <- "#EDEDED"
col_light  <- "#F6ECEC"
col_light2 <- "#EBCFCF"
col_mid    <- "#C85C5C"
col_mid2   <- "#D98888"
col_dark   <- "#8F1D21"
col_grey   <- "#D9D9D9"
col_grey2  <- "#BFBFBF"

theme_pub_cc <- function(base_size = 13) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = col_text),
      plot.subtitle = element_text(hjust = 0.5, size = 11.2, color = "#4F4F4F", margin = margin(b = 8)),
      axis.title = element_text(face = "bold", size = 12.5, color = col_text),
      axis.text = element_text(size = 11, color = col_text),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
      panel.grid.major.x = element_line(color = col_grid, linewidth = 0.55),
      panel.grid.major.y = element_line(color = col_grid, linewidth = 0.45),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = col_border, fill = NA, linewidth = 0.75),
      strip.background = element_rect(fill = col_light, color = NA),
      strip.text = element_text(face = "bold", color = col_text, size = 11),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11.2, color = col_text),
      legend.text = element_text(size = 10.4, color = col_text),
      plot.margin = margin(10, 14, 10, 14),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

save_plot_dual <- function(plot_obj, filename, width, height) {
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".pdf")),
    plot = plot_obj, width = width, height = height, units = "in"
  )
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".png")),
    plot = plot_obj, width = width, height = height, units = "in", dpi = 320
  )
}

# -----------------------------
# 5. Helper functions
# -----------------------------
check_required_files <- function(files) {
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop("These required files are missing:\n", paste(missing_files, collapse = "\n"))
  }
}

load_first_seurat_object <- function(rdata_file) {
  env <- new.env()
  loaded_names <- load(rdata_file, envir = env)

  if (length(loaded_names) == 0) {
    stop("No objects were loaded from: ", rdata_file)
  }

  obj_list <- mget(loaded_names, envir = env)
  seurat_hits <- names(obj_list)[sapply(obj_list, function(x) inherits(x, "Seurat"))]

  if (length(seurat_hits) == 0) {
    stop("No Seurat object found in: ", rdata_file)
  }

  obj_list[[seurat_hits[1]]]
}

detect_condition_col <- function(meta, target_value) {
  char_or_factor_cols <- names(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))]
  preferred_order <- c("Type", "type", "group", "Group", "condition", "Condition", "stage", "Stage", "sample_group")

  preferred_hits <- preferred_order[preferred_order %in% char_or_factor_cols]
  for (cc in preferred_hits) {
    vals <- as.character(meta[[cc]])
    if (target_value %in% vals) {
      return(cc)
    }
  }

  for (cc in char_or_factor_cols) {
    vals <- as.character(meta[[cc]])
    if (target_value %in% vals) {
      return(cc)
    }
  }

  stop("Could not detect a metadata column containing target condition: ", target_value)
}

detect_cluster_col <- function(meta, target_clusters) {
  candidate_cols <- c("seurat_clusters", "cluster", "Cluster", "clusters", "cluster_id", "clean_cluster")
  candidate_hits <- candidate_cols[candidate_cols %in% names(meta)]

  for (cc in candidate_hits) {
    vals <- as.character(meta[[cc]])
    if (all(target_clusters %in% unique(vals))) {
      return(cc)
    }
  }

  possible_cols <- names(meta)[sapply(meta, function(x) {
    is.character(x) || is.factor(x) || is.numeric(x) || is.integer(x)
  })]

  for (cc in possible_cols) {
    vals <- as.character(meta[[cc]])
    if (all(target_clusters %in% unique(vals))) {
      return(cc)
    }
  }

  stop("Could not detect a metadata column containing target clusters: ",
       paste(target_clusters, collapse = ", "))
}

extract_gene_expression <- function(seu, gene_symbol, assay_use = target_assay) {
  if (!(assay_use %in% Assays(seu))) {
    stop("Assay not found in Seurat object: ", assay_use)
  }

  genes <- rownames(seu[[assay_use]])
  if (is.null(genes) || length(genes) == 0) {
    stop("No feature names found in assay: ", assay_use)
  }

  idx <- match(toupper(gene_symbol), toupper(genes))
  if (is.na(idx)) {
    stop("Gene not found in assay ", assay_use, ": ", gene_symbol)
  }

  gene_use <- genes[idx]

  expr <- NULL

  expr <- tryCatch({
    old_assay <- DefaultAssay(seu)
    on.exit(DefaultAssay(seu) <- old_assay, add = TRUE)
    DefaultAssay(seu) <- assay_use
    df <- FetchData(seu, vars = gene_use, cells = colnames(seu))
    val <- as.numeric(df[[1]])
    names(val) <- rownames(df)
    val
  }, error = function(e) NULL)

  if (!is.null(expr) && length(expr) == ncol(seu)) {
    return(expr)
  }

  expr <- tryCatch({
    mat <- LayerData(seu, assay = assay_use, layer = "data")
    val <- as.numeric(mat[gene_use, colnames(seu)])
    names(val) <- colnames(seu)
    val
  }, error = function(e) NULL)

  if (!is.null(expr) && length(expr) == ncol(seu)) {
    return(expr)
  }

  expr <- tryCatch({
    mat <- GetAssayData(seu, assay = assay_use, slot = "data")
    val <- as.numeric(mat[gene_use, colnames(seu)])
    names(val) <- colnames(seu)
    val
  }, error = function(e) NULL)

  if (!is.null(expr) && length(expr) == ncol(seu)) {
    return(expr)
  }

  expr <- tryCatch({
    mat <- LayerData(seu, assay = assay_use, layer = "counts")
    val <- as.numeric(mat[gene_use, colnames(seu)])
    names(val) <- colnames(seu)
    val
  }, error = function(e) NULL)

  if (!is.null(expr) && length(expr) == ncol(seu)) {
    return(expr)
  }

  expr <- tryCatch({
    mat <- GetAssayData(seu, assay = assay_use, slot = "counts")
    val <- as.numeric(mat[gene_use, colnames(seu)])
    names(val) <- colnames(seu)
    val
  }, error = function(e) NULL)

  if (!is.null(expr) && length(expr) == ncol(seu)) {
    return(expr)
  }

  stop("Could not extract expression for gene: ", gene_symbol,
       " | assay fixed to: ", assay_use,
       " | rowname used: ", gene_use)
}

make_state_group <- function(expr_vector, method = "rank_tertile_all_cells") {
  expr_vector <- as.numeric(expr_vector)

  if (method != "rank_tertile_all_cells") {
    stop("Currently only 'rank_tertile_all_cells' method is implemented.")
  }

  if (all(is.na(expr_vector))) {
    stop("Expression vector is entirely NA.")
  }

  non_na_idx <- which(!is.na(expr_vector))
  if (length(non_na_idx) < 30) {
    stop("Too few non-missing cells to define state groups robustly.")
  }

  rank_all <- rank(expr_vector, ties.method = "first", na.last = "keep")
  tile_all <- rep(NA_integer_, length(expr_vector))
  tile_all[non_na_idx] <- dplyr::ntile(rank_all[non_na_idx], 3)

  group <- rep(NA_character_, length(expr_vector))
  group[tile_all == 1] <- "HMMR_low"
  group[tile_all == 2] <- "middle"
  group[tile_all == 3] <- "HMMR_high"

  low_upper_value   <- max(expr_vector[group == "HMMR_low"], na.rm = TRUE)
  high_lower_value  <- min(expr_vector[group == "HMMR_high"], na.rm = TRUE)
  middle_lower_value <- min(expr_vector[group == "middle"], na.rm = TRUE)
  middle_upper_value <- max(expr_vector[group == "middle"], na.rm = TRUE)

  list(
    group = group,
    low_upper_value = low_upper_value,
    high_lower_value = high_lower_value,
    middle_lower_value = middle_lower_value,
    middle_upper_value = middle_upper_value,
    n_low = sum(group == "HMMR_low", na.rm = TRUE),
    n_middle = sum(group == "middle", na.rm = TRUE),
    n_high = sum(group == "HMMR_high", na.rm = TRUE)
  )
}

write_input_manifest <- function() {
  manifest_df <- tibble(
    item = c(
      "epi_rdata_file",
      "ko_all_diff_file",
      "ko_sig_diff_file",
      "ko_go_file",
      "ko_kegg_file",
      "ko_focus_go_file",
      "ko_focus_kegg_file"
    ),
    path = c(
      epi_rdata_file,
      ko_all_diff_file,
      ko_sig_diff_file,
      ko_go_file,
      ko_kegg_file,
      ko_focus_go_file,
      ko_focus_kegg_file
    ),
    exists = file.exists(c(
      epi_rdata_file,
      ko_all_diff_file,
      ko_sig_diff_file,
      ko_go_file,
      ko_kegg_file,
      ko_focus_go_file,
      ko_focus_kegg_file
    ))
  )

  write.csv(
    manifest_df,
    file.path(input_manifest_dir, "00.A3_input_manifest.csv"),
    row.names = FALSE
  )
}

# -----------------------------
# 6. Save config snapshot
# -----------------------------
write_input_manifest()

config_snapshot <- tibble(
  parameter = c(
    "base_dir",
    "hmmr_gene",
    "target_assay",
    "target_condition",
    "target_clusters",
    "state_method",
    "min_pct_de",
    "logfc_threshold",
    "padj_threshold"
  ),
  value = c(
    base_dir,
    hmmr_gene,
    target_assay,
    target_condition,
    paste(target_clusters, collapse = ","),
    state_method,
    as.character(min_pct_de),
    as.character(logfc_threshold),
    as.character(padj_threshold)
  )
)

write.csv(
  config_snapshot,
  file.path(input_manifest_dir, "01.A3_config_snapshot.csv"),
  row.names = FALSE
)

cat("\nA3 config loaded successfully.\n")
cat("Base directory:\n", base_dir, "\n")
cat("Target assay fixed to:\n", target_assay, "\n")
cat("Intermediate directory:\n", intermediate_dir, "\n")
cat("Table directory:\n", table_dir, "\n")
cat("Plot directory:\n", plot_dir, "\n")