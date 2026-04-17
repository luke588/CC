rm(list = ls())
options(stringsAsFactors = FALSE)

# =========================================================
# 02.A4_HPA_support_summary.R
# Summarize manually curated HPA orthogonal support
# =========================================================

config_file <- "/Users/donglinlu/Desktop/CC1/13.public_orthogonal_support/HPA_HMMR_support/00.A4_config.R"

if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}

source(config_file)

# -----------------------------
# Extra packages
# -----------------------------
extra_pkgs <- c("patchwork")
for (pkg in extra_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(patchwork)
})

# -----------------------------
# Input files
# -----------------------------
manual_file  <- file.path(template_dir, "04.A4_HPA_manual_curation_template.csv")
target_file  <- file.path(table_dir, "03.A4_HPA_target_manifest.csv")
module_file  <- file.path(table_dir, "05.A4_support_module_context.csv")
focus_file   <- file.path(table_dir, "06.A4_focus_terms_for_notes.csv")
de_file      <- file.path(table_dir, "07.A4_A3_DE_context_summary.csv")

check_required_files(c(
  manual_file,
  target_file,
  module_file,
  focus_file,
  de_file
))

manual_df <- read.csv(manual_file, check.names = FALSE, stringsAsFactors = FALSE)
target_df <- read.csv(target_file, check.names = FALSE, stringsAsFactors = FALSE)
module_df <- read.csv(module_file, check.names = FALSE, stringsAsFactors = FALSE)
focus_df  <- read.csv(focus_file, check.names = FALSE, stringsAsFactors = FALSE)
de_df     <- read.csv(de_file, check.names = FALSE, stringsAsFactors = FALSE)

# -----------------------------
# Helper functions
# -----------------------------
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

clean_text <- function(x) {
  y <- as.character(x)
  y[is.na(y)] <- ""
  y <- trimws(y)
  y
}

is_filled_text <- function(x) {
  y <- clean_text(x)
  y != ""
}

classify_support_level <- function(x) {
  y <- tolower(clean_text(x))

  dplyr::case_when(
    y == "" ~ "Not reviewed",
    str_detect(y, "strong|high|clear|robust|main orthogonal") ~ "Strong",
    str_detect(y, "moderate|medium|detectable|available|cell-cycle-associated") ~ "Moderate",
    str_detect(y, "limited|low|weak|not detected|not prognostic|broad|non-specific") ~ "Limited",
    TRUE ~ "Moderate"
  )
}

classify_context_level <- function(x) {
  y <- tolower(clean_text(x))

  dplyr::case_when(
    y == "" ~ "Not reviewed",
    str_detect(y, "cervical-specific|direct cervical") ~ "Strong",
    str_detect(y, "cervical cancer|cervical category available|cervical context available") ~ "Moderate",
    str_detect(y, "no cervical-specific|not cervical-specific|limited") ~ "Limited",
    TRUE ~ "Moderate"
  )
}

level_to_score <- function(level_vec) {
  dplyr::case_when(
    level_vec == "Strong" ~ 3,
    level_vec == "Moderate" ~ 2,
    level_vec == "Limited" ~ 1,
    TRUE ~ 0
  )
}

short_label <- function(x, width = 42) {
  out <- clean_text(x)
  out <- stringr::str_wrap(out, width = width)
  out
}

# -----------------------------
# Standardize manual curation table
# -----------------------------
gene_col <- resolve_colname(manual_df, c("gene", "Gene", "hpa_gene_symbol"), "gene")
support_class_col <- resolve_colname(manual_df, c("support_class", "SupportClass"), "support_class")
tier_col <- resolve_colname(manual_df, c("target_tier", "tier", "TargetTier"), "target_tier")

path_col <- resolve_colname(manual_df, c("pathology_support"), "pathology_support")
cell_col <- resolve_colname(manual_df, c("cell_line_support"), "cell_line_support")
context_col <- resolve_colname(manual_df, c("cervical_context_note"), "cervical_context_note")
interp_col <- resolve_colname(manual_df, c("summary_interpretation"), "summary_interpretation")
final_use_col <- resolve_colname(manual_df, c("final_use"), "final_use")

hpa_tbl <- tibble(
  gene = as.character(manual_df[[gene_col]]),
  support_class = as.character(manual_df[[support_class_col]]),
  target_tier = as.character(manual_df[[tier_col]]),
  pathology_support = clean_text(manual_df[[path_col]]),
  cell_line_support = clean_text(manual_df[[cell_col]]),
  cervical_context_note = clean_text(manual_df[[context_col]]),
  summary_interpretation = clean_text(manual_df[[interp_col]]),
  final_use = clean_text(manual_df[[final_use_col]])
) %>%
  mutate(
    support_class = dplyr::case_when(
      support_class %in% c("Primary") ~ "Primary",
      support_class %in% c("Shared") ~ "Shared",
      support_class %in% c("State only", "HMMR-high only") ~ "State only",
      support_class %in% c("KO only") ~ "KO only",
      TRUE ~ support_class
    ),
    target_tier = dplyr::case_when(
      target_tier %in% c("Primary") ~ "Primary",
      TRUE ~ "Secondary"
    ),
    pathology_level = classify_support_level(pathology_support),
    cell_line_level = classify_support_level(cell_line_support),
    context_level = classify_context_level(cervical_context_note),
    reviewed_any = is_filled_text(pathology_support) |
      is_filled_text(cell_line_support) |
      is_filled_text(cervical_context_note) |
      is_filled_text(summary_interpretation)
  ) %>%
  distinct(gene, .keep_all = TRUE)

# -----------------------------
# Export cleaned tables
# -----------------------------
clean_hpa_tbl <- hpa_tbl %>%
  select(
    gene, target_tier, support_class,
    pathology_support, pathology_level,
    cell_line_support, cell_line_level,
    cervical_context_note, context_level,
    summary_interpretation, final_use, reviewed_any
  ) %>%
  arrange(
    factor(target_tier, levels = c("Primary", "Secondary")),
    factor(support_class, levels = c("Primary", "Shared", "State only", "KO only")),
    gene
  )

write.csv(
  clean_hpa_tbl,
  file.path(table_dir, "08.A4_HPA_curated_support_summary.csv"),
  row.names = FALSE
)

review_summary_tbl <- hpa_tbl %>%
  group_by(support_class) %>%
  summarise(
    n_targets = n(),
    n_reviewed = sum(reviewed_any, na.rm = TRUE),
    n_pathology_reviewed = sum(pathology_level != "Not reviewed", na.rm = TRUE),
    n_cell_line_reviewed = sum(cell_line_level != "Not reviewed", na.rm = TRUE),
    n_context_reviewed = sum(context_level != "Not reviewed", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    support_class = factor(support_class, levels = c("Primary", "Shared", "State only", "KO only"))
  ) %>%
  arrange(support_class)

write.csv(
  review_summary_tbl,
  file.path(table_dir, "09.A4_HPA_review_progress_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# Plot 1: HPA support summary matrix
# -----------------------------
matrix_df <- hpa_tbl %>%
  transmute(
    gene = gene,
    support_class = factor(support_class, levels = c("Primary", "Shared", "State only", "KO only")),
    Pathology = pathology_level,
    `Cell line` = cell_line_level,
    `Cervical context` = context_level
  ) %>%
  tidyr::pivot_longer(
    cols = c("Pathology", "Cell line", "Cervical context"),
    names_to = "Evidence",
    values_to = "Level"
  ) %>%
  mutate(
    Evidence = factor(Evidence, levels = c("Pathology", "Cell line", "Cervical context")),
    Level = factor(Level, levels = c("Not reviewed", "Limited", "Moderate", "Strong"))
  )

gene_order <- hpa_tbl %>%
  arrange(
    factor(target_tier, levels = c("Primary", "Secondary")),
    factor(support_class, levels = c("Primary", "Shared", "State only", "KO only")),
    gene
  ) %>%
  pull(gene)

matrix_df$gene <- factor(matrix_df$gene, levels = rev(gene_order))

p_matrix <- ggplot(matrix_df, aes(x = Evidence, y = gene)) +
  geom_tile(aes(fill = Level), color = "white", linewidth = 0.95) +
  facet_grid(support_class ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(
    values = c(
      "Not reviewed" = col_grey,
      "Limited" = col_light2,
      "Moderate" = col_mid2,
      "Strong" = col_dark
    ),
    drop = FALSE
  ) +
  labs(
    title = "HPA support summary",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_pub_cc(base_size = 15) +
  theme(
    axis.text.x = element_text(face = "bold", size = 12.5),
    axis.text.y = element_text(size = 12, color = col_text),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12.5, color = col_text),
    legend.text = element_text(size = 11.5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20, color = col_text),
    panel.spacing.y = grid::unit(0.55, "lines")
  )

save_plot_dual(
  p_matrix,
  filename = "02.A4_HPA_support_summary",
  width = 9.8,
  height = 8.6
)

# -----------------------------
# Plot 2: HMMR orthogonal support
# -----------------------------
hmmr_row <- hpa_tbl %>%
  filter(gene == primary_gene)

if (nrow(hmmr_row) == 0) {
  hmmr_row <- tibble(
    gene = primary_gene,
    pathology_support = "",
    cell_line_support = "",
    cervical_context_note = "",
    summary_interpretation = "",
    pathology_level = "Not reviewed",
    cell_line_level = "Not reviewed",
    context_level = "Not reviewed"
  )
}

hmmr_card_df <- tibble(
  Domain = factor(
    c("Pathology", "Cell line", "Cervical context"),
    levels = c("Pathology", "Cell line", "Cervical context")
  ),
  Level = factor(
    c(hmmr_row$pathology_level[1], hmmr_row$cell_line_level[1], hmmr_row$context_level[1]),
    levels = c("Not reviewed", "Limited", "Moderate", "Strong")
  ),
  Score = level_to_score(c(
    hmmr_row$pathology_level[1],
    hmmr_row$cell_line_level[1],
    hmmr_row$context_level[1]
  )),
  Note = c(
    short_label(hmmr_row$pathology_support[1], width = 52),
    short_label(hmmr_row$cell_line_support[1], width = 52),
    short_label(hmmr_row$cervical_context_note[1], width = 52)
  )
)

p_card_left <- ggplot(hmmr_card_df, aes(x = Score, y = Domain)) +
  geom_segment(
    aes(x = 0, xend = Score, y = Domain, yend = Domain),
    linewidth = 1.2,
    color = col_grey2
  ) +
  geom_point(
    aes(fill = Level),
    shape = 21,
    size = 5.2,
    color = col_border,
    stroke = 0.45
  ) +
  scale_fill_manual(
    values = c(
      "Not reviewed" = col_grey,
      "Limited" = col_light2,
      "Moderate" = col_mid2,
      "Strong" = col_dark
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 3.2),
    breaks = c(0, 1, 2, 3),
    labels = c("Not reviewed", "Limited", "Moderate", "Strong"),
    expand = c(0.02, 0.02)
  ) +
  labs(
    title = "HMMR orthogonal support",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_pub_cc(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 11.5),
    axis.text.y = element_text(face = "bold", size = 12.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20, color = col_text)
  )

note_plot_df <- hmmr_card_df %>%
  mutate(
    y_pos = c(3, 2, 1)
  )

p_card_right <- ggplot(note_plot_df, aes(x = 1, y = y_pos)) +
  geom_text(
    aes(label = Note),
    hjust = 0,
    vjust = 0.5,
    size = 4.3,
    color = col_text,
    lineheight = 1.10
  ) +
  xlim(1, 1.02) +
  ylim(0.5, 3.5) +
  theme_void() +
  theme(
    plot.margin = margin(44, 6, 10, 6)
  )

p_card <- p_card_left + p_card_right + plot_layout(widths = c(1.25, 1.05))

save_plot_dual(
  p_card,
  filename = "03.A4_HMMR_orthogonal_support",
  width = 12.4,
  height = 5.6
)

# -----------------------------
# Integrated notes export
# -----------------------------
module_short_tbl <- module_df %>%
  mutate(
    Analysis = if ("Analysis" %in% colnames(module_df)) as.character(Analysis) else "",
    module = if ("module" %in% colnames(module_df)) as.character(module) else "",
    representative_terms = if ("representative_terms" %in% colnames(module_df)) as.character(representative_terms) else ""
  ) %>%
  select(Analysis, module, representative_terms)

write.csv(
  module_short_tbl,
  file.path(notes_dir, "02.A4_module_context_for_summary.csv"),
  row.names = FALSE
)

write.csv(
  focus_df,
  file.path(notes_dir, "03.A4_focus_terms_snapshot.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "A4 HPA SUPPORT SUMMARY",
    "",
    paste0("Primary gene: ", primary_gene),
    paste0("Reviewed targets: ", sum(hpa_tbl$reviewed_any, na.rm = TRUE), "/", nrow(hpa_tbl)),
    "",
    "Interpretation reminder:",
    "- Use HPA as orthogonal support rather than as the main discovery axis.",
    "- Keep the final narrative centered on HMMR as a proliferation/cell-cycle-associated lead gene.",
    "- Avoid overstating cervical-specificity when HPA indicates broad or low cancer specificity."
  ),
  con = file.path(notes_dir, "04.A4_summary_notes.txt"),
  useBytes = TRUE
)

cat("\nA4 HPA support summary finished successfully.\n")
cat("Reviewed targets: ", sum(hpa_tbl$reviewed_any, na.rm = TRUE), "/", nrow(hpa_tbl), "\n")
cat("Tables saved to:\n", table_dir, "\n")
cat("Plots saved to:\n", plot_dir, "\n")
cat("Notes saved to:\n", notes_dir, "\n")