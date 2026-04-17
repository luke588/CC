rm(list = ls())
options(stringsAsFactors = FALSE)

# =========================================================
# 01.A4_prepare_HPA_support_targets.R
# Prepare HPA orthogonal-support targets and manual curation sheets
# Refined manuscript-style version
# =========================================================

config_file <- "/Users/donglinlu/Desktop/CC1/13.public_orthogonal_support/HPA_HMMR_support/00.A4_config.R"

if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}

source(config_file)

# -----------------------------
# Input
# -----------------------------
check_required_files(c(
  a3_overlap_file,
  a3_module_file,
  a3_de_summary_file,
  a3_focus_terms_file
))

overlap_df <- read.csv(a3_overlap_file, check.names = FALSE, stringsAsFactors = FALSE)
module_df  <- read.csv(a3_module_file, check.names = FALSE, stringsAsFactors = FALSE)
de_sum_df  <- read.csv(a3_de_summary_file, check.names = FALSE, stringsAsFactors = FALSE)
focus_df   <- read.csv(a3_focus_terms_file, check.names = FALSE, stringsAsFactors = FALSE)

# -----------------------------
# Standardize overlap table
# -----------------------------
gene_col <- resolve_colname(
  overlap_df,
  c("gene", "Gene", "SYMBOL", "symbol"),
  "gene"
)

support_class_col <- resolve_colname(
  overlap_df,
  c("support_class", "SupportClass", "class"),
  "support_class"
)

a3_metric_col <- NULL
ko_metric_col <- NULL

a3_metric_hit <- c("A3_metric", "avg_log2FC", "A3_full_metric")[c("A3_metric", "avg_log2FC", "A3_full_metric") %in% colnames(overlap_df)]
ko_metric_hit <- c("KO_metric", "diffRegulation", "avg_log2FC")[c("KO_metric", "diffRegulation", "avg_log2FC") %in% colnames(overlap_df)]

if (length(a3_metric_hit) > 0) a3_metric_col <- a3_metric_hit[1]
if (length(ko_metric_hit) > 0) ko_metric_col <- ko_metric_hit[1]

support_tbl <- tibble(
  gene = as.character(overlap_df[[gene_col]]),
  support_class = as.character(overlap_df[[support_class_col]]),
  A3_metric = if (!is.null(a3_metric_col)) suppressWarnings(as.numeric(overlap_df[[a3_metric_col]])) else NA_real_,
  KO_metric = if (!is.null(ko_metric_col)) suppressWarnings(as.numeric(overlap_df[[ko_metric_col]])) else NA_real_
) %>%
  filter(!is.na(gene), gene != "") %>%
  mutate(
    support_class = dplyr::case_when(
      support_class %in% c("Shared", "shared") ~ "Shared",
      support_class %in% c("HMMR-high only", "A3 only", "HMMR-high only ") ~ "State only",
      support_class %in% c("KO only", "KO_only") ~ "KO only",
      TRUE ~ support_class
    ),
    priority_score = pmax(abs(A3_metric), abs(KO_metric), na.rm = TRUE)
  ) %>%
  distinct(gene, .keep_all = TRUE)

# -----------------------------
# Primary / secondary target selection
# -----------------------------
primary_tbl <- tibble(
  gene = primary_gene,
  target_tier = "Primary",
  support_class = "Primary",
  selection_reason = "Lead gene prioritized by integrated epithelial-program, hub-gene, KO, and state-support evidence"
)

shared_tbl <- support_tbl %>%
  filter(support_class == "Shared", gene != primary_gene) %>%
  arrange(desc(priority_score), gene) %>%
  slice_head(n = max_shared_support_genes) %>%
  transmute(
    gene = gene,
    target_tier = "Secondary",
    support_class = "Shared",
    selection_reason = "Shared support between HMMR-high state and HMMR virtual knockout"
  )

state_only_tbl <- support_tbl %>%
  filter(support_class == "State only", gene != primary_gene) %>%
  arrange(desc(priority_score), gene) %>%
  slice_head(n = max_hmmr_high_only_genes) %>%
  transmute(
    gene = gene,
    target_tier = "Secondary",
    support_class = "State only",
    selection_reason = "State-associated support gene enriched in HMMR-high epithelial cells"
  )

ko_only_tbl <- support_tbl %>%
  filter(support_class == "KO only", gene != primary_gene) %>%
  arrange(desc(priority_score), gene) %>%
  slice_head(n = max_ko_only_genes) %>%
  transmute(
    gene = gene,
    target_tier = "Secondary",
    support_class = "KO only",
    selection_reason = "KO-associated support gene perturbed in HMMR virtual knockout"
  )

target_manifest <- bind_rows(
  primary_tbl,
  shared_tbl,
  state_only_tbl,
  ko_only_tbl
) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(
    support_class = factor(
      support_class,
      levels = c("Primary", "Shared", "State only", "KO only")
    ),
    target_tier = factor(
      target_tier,
      levels = c("Primary", "Secondary")
    )
  ) %>%
  arrange(target_tier, support_class, gene)

# -----------------------------
# Context summary
# -----------------------------
module_context_tbl <- module_df %>%
  mutate(
    module = if ("module" %in% colnames(module_df)) as.character(module_df$module) else NA_character_,
    Analysis = if ("Analysis" %in% colnames(module_df)) as.character(module_df$Analysis) else NA_character_,
    representative_terms = if ("representative_terms" %in% colnames(module_df)) as.character(module_df$representative_terms) else NA_character_
  ) %>%
  filter(module %in% module_levels) %>%
  select(Analysis, module, representative_terms)

focus_terms_tbl <- focus_df %>%
  mutate(
    Description = if ("Description" %in% colnames(focus_df)) clean_term_text(Description) else NA_character_,
    module = if ("module" %in% colnames(focus_df)) as.character(module) else NA_character_,
    score = if ("score" %in% colnames(focus_df)) suppressWarnings(as.numeric(score)) else NA_real_
  ) %>%
  filter(module %in% module_levels) %>%
  arrange(module, desc(score)) %>%
  group_by(module) %>%
  slice_head(n = max_focus_terms_for_notes) %>%
  ungroup()

# -----------------------------
# Manual HPA curation template
# -----------------------------
hpa_template <- target_manifest %>%
  mutate(
    hpa_gene_symbol = gene,
    hpa_main_url = NA_character_,
    hpa_cancer_url = NA_character_,
    hpa_cell_line_url = NA_character_,
    pathology_support = NA_character_,
    cell_line_support = NA_character_,
    cervical_context_note = NA_character_,
    summary_interpretation = NA_character_,
    final_use = NA_character_
  )

hpa_template <- hpa_template %>%
  mutate(
    hpa_main_url = ifelse(
      gene == "HMMR",
      "https://www.proteinatlas.org/ENSG00000072571-HMMR",
      hpa_main_url
    ),
    hpa_cancer_url = ifelse(
      gene == "HMMR",
      "https://www.proteinatlas.org/ENSG00000072571-HMMR/cancer",
      hpa_cancer_url
    ),
    hpa_cell_line_url = ifelse(
      gene == "HMMR",
      "https://www.proteinatlas.org/ENSG00000072571-HMMR/cell%2Bline",
      hpa_cell_line_url
    )
  )

# -----------------------------
# Export tables
# -----------------------------
write.csv(
  primary_tbl,
  file.path(table_dir, "01.A4_HPA_primary_target_gene.csv"),
  row.names = FALSE
)

write.csv(
  bind_rows(shared_tbl, state_only_tbl, ko_only_tbl),
  file.path(table_dir, "02.A4_HPA_secondary_support_genes.csv"),
  row.names = FALSE
)

write.csv(
  target_manifest,
  file.path(table_dir, "03.A4_HPA_target_manifest.csv"),
  row.names = FALSE
)

write.csv(
  hpa_template,
  file.path(template_dir, "04.A4_HPA_manual_curation_template.csv"),
  row.names = FALSE
)

write.csv(
  module_context_tbl,
  file.path(table_dir, "05.A4_support_module_context.csv"),
  row.names = FALSE
)

write.csv(
  focus_terms_tbl,
  file.path(table_dir, "06.A4_focus_terms_for_notes.csv"),
  row.names = FALSE
)

write.csv(
  de_sum_df,
  file.path(table_dir, "07.A4_A3_DE_context_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# Notes file
# -----------------------------
notes_lines <- c(
  "A4 HPA ORTHOGONAL SUPPORT NOTES",
  "",
  paste0("Primary target gene: ", primary_gene),
  "",
  "Recommended interpretation scope:",
  "1. Use HPA as orthogonal support for HMMR-centered relevance.",
  "2. Do not expand A4 into a new independent biological storyline.",
  "3. Prefer concise pathology + cell-line support summaries.",
  "",
  "Official HMMR HPA links:",
  "Main page: https://www.proteinatlas.org/ENSG00000072571-HMMR",
  "Cancer page: https://www.proteinatlas.org/ENSG00000072571-HMMR/cancer",
  "Cell line page: https://www.proteinatlas.org/ENSG00000072571-HMMR/cell%2Bline",
  "",
  "Suggested manuscript interpretation frame:",
  "- HMMR has orthogonal support from a public protein/pathology resource.",
  "- This support is used to reinforce, not redefine, the HMMR-centered epithelial program model."
)

writeLines(
  notes_lines,
  con = file.path(notes_dir, "01.A4_HPA_notes.txt"),
  useBytes = TRUE
)

# -----------------------------
# Plot: target overview
# -----------------------------
plot_df <- target_manifest %>%
  mutate(
    support_class = factor(
      as.character(support_class),
      levels = c("Primary", "Shared", "State only", "KO only")
    ),
    plot_order = dplyr::case_when(
      support_class == "Primary" ~ 1,
      support_class == "Shared" ~ 2,
      support_class == "State only" ~ 3,
      support_class == "KO only" ~ 4,
      TRUE ~ 5
    )
  ) %>%
  arrange(plot_order, gene) %>%
  mutate(
    gene = factor(gene, levels = rev(unique(gene)))
  )

p_targets <- ggplot(plot_df, aes(x = support_class, y = gene)) +
  geom_point(
    aes(size = target_tier, fill = support_class),
    shape = 21,
    color = col_border,
    stroke = 0.50
  ) +
  scale_fill_manual(
    values = c(
      "Primary" = col_dark,
      "Shared" = col_mid,
      "State only" = col_mid2,
      "KO only" = col_light2
    )
  ) +
  scale_size_manual(values = c("Primary" = 7.2, "Secondary" = 4.8), guide = "none") +
  labs(
    title = "Orthogonal support targets",
    x = "Support category",
    y = NULL,
    fill = NULL
  ) +
  theme_pub_cc(base_size = 15) +
  theme(
    axis.text.x = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 12.5, color = col_text),
    axis.title.x = element_text(face = "bold", size = 13, color = col_text),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20, color = col_text),
    legend.position = "none",
    plot.margin = margin(12, 18, 12, 20)
  )

save_plot_dual(
  p_targets,
  filename = "01.A4_HPA_target_overview",
  width = 9.2,
  height = 7.4
)

cat("\nA4 HPA target preparation finished successfully.\n")
cat("Primary target gene: ", primary_gene, "\n")
cat("Total selected targets: ", nrow(target_manifest), "\n")
cat("Shared support genes selected: ", nrow(shared_tbl), "\n")
cat("State-only genes selected: ", nrow(state_only_tbl), "\n")
cat("KO-only genes selected: ", nrow(ko_only_tbl), "\n")
cat("Tables saved to:\n", table_dir, "\n")
cat("Templates saved to:\n", template_dir, "\n")
cat("Notes saved to:\n", notes_dir, "\n")
cat("Plots saved to:\n", plot_dir, "\n")