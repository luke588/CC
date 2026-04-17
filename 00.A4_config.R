rm(list = ls())
options(stringsAsFactors = FALSE)

# =========================================================
# 00.A4_config.R
# A4 orthogonal public support
# HPA support for HMMR and selected supporting genes
# =========================================================

# -----------------------------
# 1. Fixed base paths
# -----------------------------
base_dir <- "/Users/donglinlu/Desktop/CC1/13.public_orthogonal_support/HPA_HMMR_support"

a3_table_dir <- "/Users/donglinlu/Desktop/CC1/9.scTenifoldKnk/05B.HMMR_KO_supporting_evidence/tables"

a3_overlap_file     <- file.path(a3_table_dir, "11.A3_representative_gene_overlap_summary.csv")
a3_module_file      <- file.path(a3_table_dir, "10.A3_KO_vs_HMMR_state_module_consistency.csv")
a3_de_summary_file  <- file.path(a3_table_dir, "06A.A3_HMMR_high_low_DE_summary.csv")
a3_focus_terms_file <- file.path(a3_table_dir, "09.A3_HMMR_high_focus_terms.csv")

input_manifest_dir <- file.path(base_dir, "input_manifest")
table_dir          <- file.path(base_dir, "tables")
plot_dir           <- file.path(base_dir, "plots")
notes_dir          <- file.path(base_dir, "notes")
template_dir       <- file.path(base_dir, "templates")

dir.create(base_dir,           showWarnings = FALSE, recursive = TRUE)
dir.create(input_manifest_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir,          showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir,           showWarnings = FALSE, recursive = TRUE)
dir.create(notes_dir,          showWarnings = FALSE, recursive = TRUE)
dir.create(template_dir,       showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Analysis parameters
# -----------------------------
primary_gene <- "HMMR"

max_shared_support_genes   <- 6
max_hmmr_high_only_genes   <- 4
max_ko_only_genes          <- 4
max_focus_terms_for_notes  <- 6

support_class_levels <- c("Primary", "Shared", "HMMR-high only", "KO only")
module_levels <- c("Proliferation-support", "Translation-support")

# -----------------------------
# 3. Packages
# -----------------------------
cran_pkgs <- c(
  "ggplot2", "dplyr", "readr", "stringr", "forcats",
  "scales", "tibble", "tidyr"
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

clean_term_text <- function(x) {
  y <- as.character(x)
  y <- gsub("^KEGG_", "", y)
  y <- gsub("_", " ", y)
  y <- gsub("\\s+", " ", y)
  y <- trimws(y)
  y <- tolower(y)
  y <- tools::toTitleCase(y)
  y
}

write_input_manifest <- function() {
  manifest_df <- tibble(
    item = c(
      "a3_overlap_file",
      "a3_module_file",
      "a3_de_summary_file",
      "a3_focus_terms_file"
    ),
    path = c(
      a3_overlap_file,
      a3_module_file,
      a3_de_summary_file,
      a3_focus_terms_file
    ),
    exists = file.exists(c(
      a3_overlap_file,
      a3_module_file,
      a3_de_summary_file,
      a3_focus_terms_file
    ))
  )

  write.csv(
    manifest_df,
    file.path(input_manifest_dir, "00.A4_input_manifest.csv"),
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
    "primary_gene",
    "max_shared_support_genes",
    "max_hmmr_high_only_genes",
    "max_ko_only_genes",
    "max_focus_terms_for_notes"
  ),
  value = c(
    base_dir,
    primary_gene,
    as.character(max_shared_support_genes),
    as.character(max_hmmr_high_only_genes),
    as.character(max_ko_only_genes),
    as.character(max_focus_terms_for_notes)
  )
)

write.csv(
  config_snapshot,
  file.path(input_manifest_dir, "01.A4_config_snapshot.csv"),
  row.names = FALSE
)

cat("\nA4 config loaded successfully.\n")
cat("Base directory:\n", base_dir, "\n")
cat("Table directory:\n", table_dir, "\n")
cat("Plot directory:\n", plot_dir, "\n")
cat("Notes directory:\n", notes_dir, "\n")
cat("Template directory:\n", template_dir, "\n")