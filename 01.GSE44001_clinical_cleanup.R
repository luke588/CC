source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
rm(list = ls())
graphics.off()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
})

# =========================================================
# GSE44001 clinical cleanup
# Unified and robust version
# Input and output under:
# /Users/donglinlu/Desktop/CC1/6.bulk_validation/GSE44001
# =========================================================

# =========================
# USER SETTINGS
# =========================
base_dir <- "/Users/donglinlu/Desktop/CC1/6.bulk_validation/GSE44001"

clinical_file   <- file.path(base_dir, "gse44001_clinical_raw.txt")
expression_file <- file.path(base_dir, "geneMatrix.txt")

out_dir <- file.path(base_dir, "01.clinical_cleanup")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# HELPERS
# =========================
find_first_col <- function(df, candidates, required = TRUE, label = NULL) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])

  if (required) {
    stop(
      "Cannot find required column",
      if (!is.null(label)) paste0(" for ", label) else "",
      ". Available columns:\n",
      paste(colnames(df), collapse = ", ")
    )
  } else {
    return(NA_character_)
  }
}

clean_stage_text <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- str_replace_all(x, "\\s+", " ")
  x[x == ""] <- NA_character_
  x
}

clean_status_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))

  x0 <- as.character(x)
  x1 <- str_to_lower(trimws(x0))

  out <- case_when(
    x1 %in% c("0", "censored", "disease free", "disease-free", "alive without disease",
              "no event", "event-free", "negative", "none") ~ 0,
    x1 %in% c("1", "recurred", "recurrence", "event", "dead of disease",
              "progressed", "positive") ~ 1,
    TRUE ~ suppressWarnings(as.numeric(x0))
  )

  as.numeric(out)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

# =========================
# CHECK INPUT FILES
# =========================
if (!file.exists(clinical_file)) {
  stop("Cannot find clinical file: ", clinical_file)
}
if (!file.exists(expression_file)) {
  stop("Cannot find expression file: ", expression_file)
}

# =========================
# 1. READ CLINICAL TABLE
# =========================
clin <- read.delim(
  clinical_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(clin) <- trimws(colnames(clin))

writeLines(
  colnames(clin),
  file.path(out_dir, "00.original_clinical_colnames.txt")
)

write.csv(
  clin,
  file.path(out_dir, "01.original_clinical_table_backup.csv"),
  row.names = FALSE
)

# =========================
# 2. DETECT AND CLEAN KEY COLUMNS
# =========================
sample_col <- find_first_col(
  clin,
  candidates = c("ID", "Sample", "sample", "sample_id", "gsm", "GSM"),
  required = TRUE,
  label = "sample"
)

stage_col <- find_first_col(
  clin,
  candidates = c("Stage", "stage", "FIGO_stage", "figo_stage", "Tumor stage"),
  required = TRUE,
  label = "stage"
)

largest_diameter_col <- find_first_col(
  clin,
  candidates = c("largest diameter", "largest_diameter", "Largest diameter", "tumor_size"),
  required = FALSE,
  label = "largest diameter"
)

dfs_months_col <- find_first_col(
  clin,
  candidates = c(
    "disease_free_survival_(dfs)_(months)",
    "DFS_months",
    "dfs_months",
    "disease free survival (months)",
    "disease_free_survival"
  ),
  required = TRUE,
  label = "DFS months"
)

dfs_status_col <- find_first_col(
  clin,
  candidates = c(
    "status_of_dfs",
    "DFS_status",
    "dfs_status",
    "disease_free_survival_status",
    "status"
  ),
  required = TRUE,
  label = "DFS status"
)

clin_clean <- clin %>%
  transmute(
    sample = trimws(as.character(.data[[sample_col]])),
    stage = clean_stage_text(.data[[stage_col]]),
    largest_diameter = if (!is.na(largest_diameter_col)) safe_numeric(.data[[largest_diameter_col]]) else NA_real_,
    dfs_months = safe_numeric(.data[[dfs_months_col]]),
    dfs_status = clean_status_numeric(.data[[dfs_status_col]])
  ) %>%
  mutate(
    sample = str_replace_all(sample, "\\s+", ""),
    dfs_event = dfs_status
  )

# 去重：同一样本如果重复，优先保留 DFS 信息更完整的一行
clin_clean <- clin_clean %>%
  mutate(
    non_missing_score =
      (!is.na(stage)) +
      (!is.na(largest_diameter)) +
      (!is.na(dfs_months)) +
      (!is.na(dfs_status))
  ) %>%
  arrange(sample, desc(non_missing_score)) %>%
  distinct(sample, .keep_all = TRUE) %>%
  select(-non_missing_score)

write.csv(
  clin_clean,
  file.path(out_dir, "02.GSE44001_clinical_clean.csv"),
  row.names = FALSE
)

# =========================
# 3. READ EXPRESSION HEADER
# =========================
expr_header <- read.delim(
  expression_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  nrows = 5
)

expr_gene_col <- colnames(expr_header)[1]
expr_samples  <- colnames(expr_header)[-1]
expr_samples  <- trimws(expr_samples)
expr_samples  <- str_replace_all(expr_samples, "\\s+", "")

if (expr_gene_col != "geneNames") {
  warning(
    "First expression column is '", expr_gene_col,
    "', not 'geneNames'. Please double-check."
  )
}

# =========================
# 4. MATCH CLINICAL WITH EXPRESSION
# =========================
match_df <- data.frame(
  expr_sample = expr_samples,
  in_clinical = expr_samples %in% clin_clean$sample,
  stringsAsFactors = FALSE
)

write.csv(
  match_df,
  file.path(out_dir, "03.expression_clinical_match.csv"),
  row.names = FALSE
)

clin_matched <- data.frame(sample = expr_samples, stringsAsFactors = FALSE) %>%
  left_join(clin_clean, by = "sample") %>%
  filter(!is.na(stage) | !is.na(dfs_months) | !is.na(dfs_status))

write.csv(
  clin_matched,
  file.path(out_dir, "04.GSE44001_clinical_clean_matched.csv"),
  row.names = FALSE
)

only_in_expr <- setdiff(expr_samples, clin_clean$sample)
only_in_clin <- setdiff(clin_clean$sample, expr_samples)

writeLines(
  only_in_expr,
  file.path(out_dir, "05.samples_only_in_expression.txt")
)

writeLines(
  only_in_clin,
  file.path(out_dir, "06.samples_only_in_clinical.txt")
)

# =========================
# 5. QC SUMMARY
# =========================
stage_tab <- table(clin_clean$stage, useNA = "ifany")
dfs_tab   <- table(clin_clean$dfs_status, useNA = "ifany")

qc_lines <- c(
  paste0("Clinical file: ", clinical_file),
  paste0("Expression file: ", expression_file),
  paste0("Expression first column: ", expr_gene_col),
  "",
  "Detected clinical columns:",
  paste0("sample_col = ", sample_col),
  paste0("stage_col = ", stage_col),
  paste0("largest_diameter_col = ", ifelse(is.na(largest_diameter_col), "NOT FOUND", largest_diameter_col)),
  paste0("dfs_months_col = ", dfs_months_col),
  paste0("dfs_status_col = ", dfs_status_col),
  "",
  paste0("Clinical sample count (raw): ", nrow(clin)),
  paste0("Clinical sample count (clean unique): ", nrow(clin_clean)),
  paste0("Expression sample count: ", length(expr_samples)),
  paste0("Matched sample count: ", nrow(clin_matched)),
  paste0("Samples only in expression: ", length(only_in_expr)),
  paste0("Samples only in clinical: ", length(only_in_clin)),
  "",
  paste0("Missing stage: ", sum(is.na(clin_clean$stage))),
  paste0("Missing largest_diameter: ", sum(is.na(clin_clean$largest_diameter))),
  paste0("Missing dfs_months: ", sum(is.na(clin_clean$dfs_months))),
  paste0("Missing dfs_status: ", sum(is.na(clin_clean$dfs_status))),
  "",
  paste0("Stage levels: ", paste(names(stage_tab), as.integer(stage_tab), collapse = "; ")),
  paste0("DFS status table: ", paste(names(dfs_tab), as.integer(dfs_tab), collapse = "; "))
)

writeLines(
  qc_lines,
  file.path(out_dir, "07.cleanup_summary.txt")
)

# =========================
# 6. OPTIONAL DFS-READY TABLE
# =========================
clin_surv_ready <- clin_matched %>%
  filter(!is.na(dfs_months), !is.na(dfs_status)) %>%
  mutate(
    dfs_months = as.numeric(dfs_months),
    dfs_status = as.numeric(dfs_status)
  )

write.csv(
  clin_surv_ready,
  file.path(out_dir, "08.GSE44001_clinical_survival_ready.csv"),
  row.names = FALSE
)

# =========================
# 7. FINISH
# =========================
cat("Done.\n")
cat("Output directory:\n", out_dir, "\n\n")
cat("Key output files:\n")
cat("02.GSE44001_clinical_clean.csv\n")
cat("04.GSE44001_clinical_clean_matched.csv\n")
cat("07.cleanup_summary.txt\n")
cat("08.GSE44001_clinical_survival_ready.csv\n")