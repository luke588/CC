rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
options(stringsAsFactors = FALSE)

# =========================================================
# 00.A1_config.R
# A1 signature robustness - fixed config, paths, theme, helpers
# Self-contained version for current folder
# =========================================================

# -----------------------------
# 1. Fixed base paths
# -----------------------------
base_dir <- "/Users/donglinlu/Desktop/CC1/10.signature_robustness"

# marker files are now stored directly in current A1 folder
cluster6_marker_file  <- file.path(base_dir, "step1_5.cluster_6.markers_vs_rest.csv")
cluster9_marker_file  <- file.path(base_dir, "step1_5.cluster_9.markers_vs_rest.csv")
cluster10_marker_file <- file.path(base_dir, "step1_5.cluster_10.markers_vs_rest.csv")

# downstream cohort files remain in original project folders
gse63514_expr_file <- "/Users/donglinlu/Desktop/CC1/8.Hub Gene/GSE63514_geneMatrix.txt"
gse63514_meta_file <- "/Users/donglinlu/Desktop/CC1/8.Hub Gene/GSE63514_metadata.txt"

gse44001_base_dir  <- "/Users/donglinlu/Desktop/CC1/6.bulk_validation/GSE44001"
gse44001_expr_file <- file.path(gse44001_base_dir, "geneMatrix.txt")
gse44001_clin_file <- file.path(gse44001_base_dir, "04.GSE44001_clinical_clean_matched.csv")

input_manifest_dir <- file.path(base_dir, "input_manifest")
signature_dir      <- file.path(base_dir, "signatures")
intermediate_dir   <- file.path(base_dir, "intermediate")
table_dir          <- file.path(base_dir, "tables")
plot_dir           <- file.path(base_dir, "plots")

dir.create(base_dir,           showWarnings = FALSE, recursive = TRUE)
dir.create(input_manifest_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(signature_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(intermediate_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir,          showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir,           showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Analysis parameters
# -----------------------------
signature_sizes <- c(30, 50, 100)

progression_group_levels <- c("Normal", "LSIL", "HSIL", "Cancer")
gse44001_stage_levels    <- c("IA2", "IB1", "IB2", "IIA")

min_pct1  <- 0.20
max_pct2  <- 0.50
min_logfc <- 0.25

program_info <- data.frame(
  cluster_id    = c("6", "9", "10"),
  program_code  = c("cluster6", "cluster9", "cluster10"),
  program_label = c(
    "proliferation_core_epithelial",
    "malignant_squamous_epithelial",
    "EMT_like_epithelial"
  ),
  marker_file = c(
    cluster6_marker_file,
    cluster9_marker_file,
    cluster10_marker_file
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 3. Packages
# -----------------------------
cran_pkgs <- c(
  "ggplot2", "dplyr", "readr", "stringr", "forcats",
  "scales", "tibble", "purrr", "tidyr"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(scales)
  library(tibble)
  library(purrr)
  library(tidyr)
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
      plot.subtitle = element_text(hjust = 0.5, size = 11.2, color = "#4F4F4F", margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = 12.5, color = col_text),
      axis.text = element_text(size = 11, color = col_text),
      axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
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

read_marker_file <- function(file) {
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)

  if (ncol(df) == 0) {
    stop("No columns detected in marker file: ", file)
  }

  if (!("gene" %in% colnames(df))) {
    if (!is.null(rownames(df)) &&
        !all(rownames(df) %in% as.character(seq_len(nrow(df))))) {
      df$gene <- rownames(df)
    }
  }

  df
}

resolve_colname <- function(df, candidates, label_for_error) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    stop(
      "Cannot find required column for ", label_for_error, ". Checked: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[1]
}

prepare_marker_table <- function(df) {
  gene_col  <- resolve_colname(df, c("gene", "Gene", "SYMBOL", "symbol", "feature"), "gene")
  logfc_col <- resolve_colname(df, c("avg_log2FC", "avg_logFC", "log2FC", "logFC"), "logFC")
  pct1_col  <- resolve_colname(df, c("pct.1", "pct_1", "pct1"), "pct.1")
  pct2_col  <- resolve_colname(df, c("pct.2", "pct_2", "pct2"), "pct.2")

  p_adj_col <- NULL
  p_adj_candidates <- c("p_val_adj", "p.adjust", "padj", "adj.P.Val")
  p_adj_hit <- p_adj_candidates[p_adj_candidates %in% colnames(df)]
  if (length(p_adj_hit) > 0) {
    p_adj_col <- p_adj_hit[1]
  }

  out <- tibble(
    gene       = as.character(df[[gene_col]]),
    avg_log2FC = suppressWarnings(as.numeric(df[[logfc_col]])),
    pct.1      = suppressWarnings(as.numeric(df[[pct1_col]])),
    pct.2      = suppressWarnings(as.numeric(df[[pct2_col]]))
  )

  if (!is.null(p_adj_col)) {
    out$p_val_adj <- suppressWarnings(as.numeric(df[[p_adj_col]]))
  } else {
    out$p_val_adj <- NA_real_
  }

  out %>%
    filter(!is.na(gene), gene != "") %>%
    group_by(gene) %>%
    slice_max(order_by = avg_log2FC, n = 1, with_ties = FALSE) %>%
    ungroup()
}

filter_marker_candidates <- function(marker_tbl) {
  marker_tbl %>%
    filter(
      !is.na(avg_log2FC),
      !is.na(pct.1),
      !is.na(pct.2),
      pct.1 >= min_pct1,
      pct.2 <= max_pct2,
      avg_log2FC >= min_logfc
    ) %>%
    arrange(desc(avg_log2FC), desc(pct.1), pct.2, p_val_adj, gene)
}

make_signature_filename <- function(program_code, program_label, size) {
  paste0(program_code, "_", program_label, "_top", size, ".txt")
}

make_signature_genes <- function(candidate_tbl, size) {
  candidate_tbl %>%
    slice_head(n = min(size, nrow(candidate_tbl))) %>%
    pull(gene) %>%
    unique()
}

pairwise_overlap_summary <- function(signature_list) {
  nm <- names(signature_list)

  if (length(nm) < 2) {
    return(data.frame())
  }

  purrr::map_dfr(combn(nm, 2, simplify = FALSE), function(x) {
    g1 <- unique(signature_list[[x[1]]])
    g2 <- unique(signature_list[[x[2]]])

    inter <- intersect(g1, g2)
    uni   <- union(g1, g2)

    tibble(
      set1 = x[1],
      set2 = x[2],
      n_set1 = length(g1),
      n_set2 = length(g2),
      n_intersection = length(inter),
      jaccard = ifelse(length(uni) == 0, NA_real_, length(inter) / length(uni)),
      overlap_genes = paste(inter, collapse = ";")
    )
  })
}

write_input_manifest <- function() {
  manifest_df <- tibble(
    item = c(
      "cluster6_marker_file",
      "cluster9_marker_file",
      "cluster10_marker_file",
      "gse63514_expr_file",
      "gse63514_meta_file",
      "gse44001_expr_file",
      "gse44001_clin_file"
    ),
    path = c(
      cluster6_marker_file,
      cluster9_marker_file,
      cluster10_marker_file,
      gse63514_expr_file,
      gse63514_meta_file,
      gse44001_expr_file,
      gse44001_clin_file
    ),
    exists = file.exists(c(
      cluster6_marker_file,
      cluster9_marker_file,
      cluster10_marker_file,
      gse63514_expr_file,
      gse63514_meta_file,
      gse44001_expr_file,
      gse44001_clin_file
    ))
  )

  write.csv(
    manifest_df,
    file.path(input_manifest_dir, "00.A1_input_manifest.csv"),
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
    "signature_dir",
    "table_dir",
    "plot_dir",
    "signature_sizes",
    "min_pct1",
    "max_pct2",
    "min_logfc"
  ),
  value = c(
    base_dir,
    signature_dir,
    table_dir,
    plot_dir,
    paste(signature_sizes, collapse = ","),
    as.character(min_pct1),
    as.character(max_pct2),
    as.character(min_logfc)
  )
)

write.csv(
  config_snapshot,
  file.path(input_manifest_dir, "01.A1_config_snapshot.csv"),
  row.names = FALSE
)

cat("\nA1 config loaded successfully.\n")
cat("Base directory:\n", base_dir, "\n")
cat("Signature output directory:\n", signature_dir, "\n")
cat("Table output directory:\n", table_dir, "\n")
cat("Plot output directory:\n", plot_dir, "\n")