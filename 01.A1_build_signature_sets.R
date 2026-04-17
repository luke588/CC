rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
options(stringsAsFactors = FALSE)

# =========================================================
# 01.A1_build_signature_sets.R
# Build top30 / top50 / top100 signatures for clusters 6 / 9 / 10
# using fixed marker files and fixed filtering logic
# Self-contained with config loading check
# =========================================================

config_file <- "/Users/donglinlu/Desktop/CC1/10.signature_robustness/00.A1_config.R"

if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}

source(config_file)

# -----------------------------
# 1. Check required marker files
# -----------------------------
check_required_files(c(
  cluster6_marker_file,
  cluster9_marker_file,
  cluster10_marker_file
))

# -----------------------------
# 2. Build signatures
# -----------------------------
signature_store <- list()
manifest_rows   <- list()
candidate_rows  <- list()

for (i in seq_len(nrow(program_info))) {
  program_code  <- program_info$program_code[i]
  program_label <- program_info$program_label[i]
  marker_file   <- program_info$marker_file[i]

  raw_df  <- read_marker_file(marker_file)
  prep_df <- prepare_marker_table(raw_df)
  cand_df <- filter_marker_candidates(prep_df)

  candidate_rows[[program_code]] <- cand_df %>%
    mutate(
      program_code  = program_code,
      program_label = program_label,
      marker_file   = marker_file
    )

  for (size in signature_sizes) {
    sig_genes <- make_signature_genes(cand_df, size = size)

    sig_name <- paste0(program_code, "_top", size)
    signature_store[[sig_name]] <- sig_genes

    out_file <- file.path(
      signature_dir,
      make_signature_filename(program_code, program_label, size)
    )

    writeLines(sig_genes, out_file, useBytes = TRUE)

    preview_genes <- paste(utils::head(sig_genes, 10), collapse = ";")

    manifest_rows[[length(manifest_rows) + 1]] <- tibble(
      program_code       = program_code,
      program_label      = program_label,
      size               = size,
      marker_file        = marker_file,
      n_raw_markers      = nrow(prep_df),
      n_filtered_marker  = nrow(cand_df),
      n_exported         = length(sig_genes),
      output_file        = out_file,
      preview_top10      = preview_genes
    )
  }
}

# -----------------------------
# 3. Export manifests and candidate tables
# -----------------------------
manifest_df <- bind_rows(manifest_rows) %>%
  arrange(program_code, size)

candidate_df <- bind_rows(candidate_rows) %>%
  dplyr::select(
    program_code, program_label, marker_file,
    gene, avg_log2FC, pct.1, pct.2, p_val_adj
  ) %>%
  arrange(program_code, desc(avg_log2FC), desc(pct.1), pct.2, p_val_adj, gene)

write.csv(
  manifest_df,
  file.path(table_dir, "01.signature_export_manifest.csv"),
  row.names = FALSE
)

write.csv(
  candidate_df,
  file.path(table_dir, "02.signature_candidate_genes_after_filter.csv"),
  row.names = FALSE
)

# -----------------------------
# 4. Export overlap summaries
# -----------------------------
overlap_df <- pairwise_overlap_summary(signature_store)

if (nrow(overlap_df) > 0) {
  overlap_df <- overlap_df %>%
    tidyr::separate(set1, into = c("program1", "size1"), sep = "_top", remove = FALSE) %>%
    tidyr::separate(set2, into = c("program2", "size2"), sep = "_top", remove = FALSE) %>%
    mutate(
      size1 = as.integer(size1),
      size2 = as.integer(size2)
    ) %>%
    arrange(program1, size1, program2, size2)
}

write.csv(
  overlap_df,
  file.path(table_dir, "03.signature_overlap_manifest.csv"),
  row.names = FALSE
)

# -----------------------------
# 5. Export simple wide overlap matrix
# -----------------------------
sig_names <- names(signature_store)

overlap_mat <- matrix(
  NA_real_,
  nrow = length(sig_names),
  ncol = length(sig_names),
  dimnames = list(sig_names, sig_names)
)

for (r in sig_names) {
  for (c in sig_names) {
    g1 <- unique(signature_store[[r]])
    g2 <- unique(signature_store[[c]])
    inter <- intersect(g1, g2)
    uni   <- union(g1, g2)

    overlap_mat[r, c] <- ifelse(length(uni) == 0, NA_real_, length(inter) / length(uni))
  }
}

write.csv(
  as.data.frame(overlap_mat),
  file.path(table_dir, "04.signature_overlap_jaccard_matrix.csv"),
  row.names = TRUE
)

# -----------------------------
# 6. Plot 1: filtered marker counts by program
# -----------------------------
count_df <- manifest_df %>%
  distinct(program_code, program_label, n_raw_markers, n_filtered_marker) %>%
  mutate(
    display_label = factor(
      program_label,
      levels = c(
        "malignant_squamous_epithelial",
        "EMT_like_epithelial",
        "proliferation_core_epithelial"
      )
    )
  ) %>%
  tidyr::pivot_longer(
    cols = c(n_raw_markers, n_filtered_marker),
    names_to = "count_type",
    values_to = "n_genes"
  ) %>%
  mutate(
    count_type = factor(
      count_type,
      levels = c("n_raw_markers", "n_filtered_marker"),
      labels = c("Raw markers", "Filtered candidates")
    )
  )

p_count <- ggplot(count_df, aes(x = display_label, y = n_genes, fill = count_type)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62,
    color = col_border,
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = comma(n_genes)),
    position = position_dodge(width = 0.72),
    vjust = -0.28,
    size = 3.8,
    color = col_text
  ) +
  scale_fill_manual(values = c(col_grey, col_mid)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Marker filtering summary for prioritized epithelial programs",
    x = NULL,
    y = "Gene count",
    fill = NULL
  ) +
  theme_pub_cc() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

save_plot_dual(
  p_count,
  filename = "01.A1_marker_filtering_summary",
  width = 9.6,
  height = 6.2
)

# -----------------------------
# 7. Plot 2: exported signature sizes
# -----------------------------
size_df <- manifest_df %>%
  mutate(
    display_label = factor(
      program_label,
      levels = c(
        "malignant_squamous_epithelial",
        "EMT_like_epithelial",
        "proliferation_core_epithelial"
      )
    ),
    size = factor(size, levels = c(30, 50, 100))
  )

p_size <- ggplot(size_df, aes(x = size, y = n_exported, fill = size)) +
  geom_col(width = 0.62, color = col_border, linewidth = 0.35) +
  geom_text(
    aes(label = n_exported),
    vjust = -0.28,
    size = 3.8,
    color = col_text
  ) +
  facet_grid(. ~ display_label) +
  scale_fill_manual(values = c(col_light2, col_mid2, col_dark), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Exported signature sizes for prioritized epithelial programs",
    x = "Signature size",
    y = "Exported gene count"
  ) +
  theme_pub_cc()

save_plot_dual(
  p_size,
  filename = "02.A1_signature_size_export_summary",
  width = 10.2,
  height = 5.8
)

# -----------------------------
# 8. Save RDS objects for downstream scripts
# -----------------------------
saveRDS(
  signature_store,
  file = file.path(intermediate_dir, "01.signature_store.rds")
)

saveRDS(
  manifest_df,
  file = file.path(intermediate_dir, "02.signature_manifest.rds")
)

saveRDS(
  candidate_df,
  file = file.path(intermediate_dir, "03.signature_candidates_after_filter.rds")
)

# -----------------------------
# 9. Done
# -----------------------------
cat("\nA1 signature sets built successfully.\n")
cat("Signatures saved to:\n", signature_dir, "\n")
cat("Tables saved to:\n", table_dir, "\n")
cat("Plots saved to:\n", plot_dir, "\n")
cat("Intermediate RDS saved to:\n", intermediate_dir, "\n")