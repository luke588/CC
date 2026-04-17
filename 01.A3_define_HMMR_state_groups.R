rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
options(stringsAsFactors = FALSE)

# =========================================================
# 01.A3_define_HMMR_state_groups.R
# Define HMMR-high / middle / low groups
# in CA_HPV prioritized pooled epithelial cells (clusters 6 / 9 / 10)
# RNA-assay explicit version
# Final fixed version
# =========================================================

config_file <- "/Users/donglinlu/Desktop/CC1/9.scTenifoldKnk/05B.HMMR_KO_supporting_evidence/00.A3_config.R"

if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file)
}

source(config_file)

# -----------------------------
# Helper functions
# -----------------------------
extract_first_numeric <- function(x, candidate_names, default = NA_real_) {
  for (nm in candidate_names) {
    val <- x[[nm]]
    if (!is.null(val)) {
      val_num <- suppressWarnings(as.numeric(unlist(val)[1]))
      if (length(val_num) == 1 && is.finite(val_num)) {
        return(val_num)
      }
    }
  }
  default
}

extract_first_vector <- function(x, candidate_names) {
  for (nm in candidate_names) {
    val <- x[[nm]]
    if (!is.null(val)) {
      return(val)
    }
  }
  NULL
}

canonical_cluster_order <- c("9", "10", "6")

# -----------------------------
# 1. Check required input
# -----------------------------
check_required_files(c(epi_rdata_file))

# -----------------------------
# 2. Load Seurat object
# -----------------------------
seu <- load_first_seurat_object(epi_rdata_file)
meta <- seu@meta.data

if (nrow(meta) == 0) {
  stop("Metadata is empty in Seurat object.")
}

condition_col <- detect_condition_col(meta, target_condition)
cluster_col   <- detect_cluster_col(meta, target_clusters)

message("Detected condition column: ", condition_col)
message("Detected cluster column: ", cluster_col)

meta_work <- meta %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    .condition = as.character(.data[[condition_col]]),
    .cluster   = as.character(.data[[cluster_col]])
  )

# -----------------------------
# 3. Subset to target cells
# -----------------------------
target_meta <- meta_work %>%
  filter(
    .condition == target_condition,
    .cluster %in% target_clusters
  ) %>%
  distinct(cell, .keep_all = TRUE)

if (nrow(target_meta) == 0) {
  stop("No cells found after filtering to target condition and target clusters.")
}

target_meta$.cluster <- factor(as.character(target_meta$.cluster), levels = canonical_cluster_order)

target_cells <- target_meta$cell
seu_target <- subset(seu, cells = target_cells)

# -----------------------------
# 4. Extract HMMR expression
# -----------------------------
hmmr_expr <- extract_gene_expression(seu_target, hmmr_gene, assay_use = target_assay)

expr_df <- tibble(
  cell = names(hmmr_expr),
  HMMR_expr = as.numeric(hmmr_expr)
)

target_meta <- target_meta %>%
  left_join(expr_df, by = "cell")

if (!("HMMR_expr" %in% colnames(target_meta))) {
  stop("HMMR_expr column was not created successfully.")
}

if (all(is.na(target_meta$HMMR_expr))) {
  stop("HMMR expression is entirely NA after extraction.")
}

message("HMMR expression extracted successfully.")
message("Matched cells: ", sum(!is.na(target_meta$HMMR_expr)))
message("Median expression: ", round(median(target_meta$HMMR_expr, na.rm = TRUE), 4))
message("Max expression: ", round(max(target_meta$HMMR_expr, na.rm = TRUE), 4))

# -----------------------------
# 5. Define tertile-based state groups
# -----------------------------
state_res <- make_state_group(target_meta$HMMR_expr, method = state_method)

group_vec <- extract_first_vector(
  state_res,
  c("group", "groups", "state_group", "state", "label")
)

if (is.null(group_vec)) {
  stop("make_state_group() did not return a valid group vector.")
}

group_vec <- as.character(unlist(group_vec))

if (length(group_vec) != nrow(target_meta)) {
  stop("Length of state group vector does not match target_meta rows.")
}

group_vec[group_vec %in% c("low", "Low", "LOW")] <- "HMMR_low"
group_vec[group_vec %in% c("mid", "Mid", "MID", "middle", "Middle")] <- "middle"
group_vec[group_vec %in% c("high", "High", "HIGH")] <- "HMMR_high"

q_low_val <- extract_first_numeric(
  state_res,
  c("q_low", "low_cut", "cut_low", "low_cutoff", "lower_cutoff", "q1")
)

q_high_val <- extract_first_numeric(
  state_res,
  c("q_high", "high_cut", "cut_high", "high_cutoff", "upper_cutoff", "q2", "q3")
)

target_meta <- target_meta %>%
  mutate(
    HMMR_state = factor(group_vec, levels = c("HMMR_low", "middle", "HMMR_high"))
  )

cut_summary <- tibble(
  hmmr_gene = hmmr_gene,
  target_assay = target_assay,
  state_method = state_method,
  q_low = q_low_val,
  q_high = q_high_val,
  n_total = nrow(target_meta),
  n_low = sum(target_meta$HMMR_state == "HMMR_low", na.rm = TRUE),
  n_middle = sum(target_meta$HMMR_state == "middle", na.rm = TRUE),
  n_high = sum(target_meta$HMMR_state == "HMMR_high", na.rm = TRUE),
  condition_col = condition_col,
  cluster_col = cluster_col
)

# -----------------------------
# 6. Export tables
# -----------------------------
cell_manifest_df <- target_meta %>%
  transmute(
    cell = cell,
    condition = .condition,
    cluster = as.character(.cluster),
    HMMR_expr = HMMR_expr,
    HMMR_state = as.character(HMMR_state)
  ) %>%
  mutate(cluster = factor(cluster, levels = canonical_cluster_order)) %>%
  arrange(cluster, HMMR_state, desc(HMMR_expr)) %>%
  mutate(cluster = as.character(cluster))

group_summary_df <- target_meta %>%
  group_by(HMMR_state) %>%
  summarise(
    n_cells = n(),
    mean_expr = mean(HMMR_expr, na.rm = TRUE),
    median_expr = median(HMMR_expr, na.rm = TRUE),
    sd_expr = sd(HMMR_expr, na.rm = TRUE),
    min_expr = min(HMMR_expr, na.rm = TRUE),
    max_expr = max(HMMR_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    HMMR_state = factor(HMMR_state, levels = c("HMMR_low", "middle", "HMMR_high"))
  ) %>%
  arrange(HMMR_state)

cluster_state_summary_df <- target_meta %>%
  mutate(.cluster = factor(as.character(.cluster), levels = canonical_cluster_order)) %>%
  count(cluster = .cluster, HMMR_state) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(cluster, HMMR_state) %>%
  mutate(cluster = as.character(cluster))

write.csv(
  cell_manifest_df,
  file.path(table_dir, "01.A3_HMMR_state_cell_manifest.csv"),
  row.names = FALSE
)

write.csv(
  group_summary_df,
  file.path(table_dir, "02.A3_HMMR_state_group_summary.csv"),
  row.names = FALSE
)

write.csv(
  cut_summary,
  file.path(table_dir, "03.A3_HMMR_state_cut_summary.csv"),
  row.names = FALSE
)

write.csv(
  cluster_state_summary_df,
  file.path(table_dir, "04.A3_HMMR_state_cluster_breakdown.csv"),
  row.names = FALSE
)

# -----------------------------
# 7. Save intermediate objects
# -----------------------------
saveRDS(
  target_meta,
  file = file.path(intermediate_dir, "01.A3_target_meta_with_HMMR_states.rds")
)

saveRDS(
  seu_target,
  file = file.path(intermediate_dir, "02.A3_target_seurat_subset.rds")
)

state_vec <- target_meta$HMMR_state[match(colnames(seu_target), target_meta$cell)]
seu_target$A3_HMMR_state <- factor(
  state_vec,
  levels = c("HMMR_low", "middle", "HMMR_high")
)

saveRDS(
  seu_target,
  file = file.path(intermediate_dir, "03.A3_target_seurat_with_state_annotation.rds")
)

# -----------------------------
# 8. Plot 1: HMMR expression distribution
# -----------------------------
dist_df <- target_meta %>%
  mutate(
    HMMR_state = factor(HMMR_state, levels = c("HMMR_low", "middle", "HMMR_high"))
  )

p_hist <- ggplot(dist_df, aes(x = HMMR_expr)) +
  geom_histogram(
    bins = 45,
    fill = col_mid2,
    color = col_border,
    linewidth = 0.25,
    alpha = 0.90
  ) +
  labs(
    title = "HMMR expression distribution in prioritized CA_HPV epithelial cells",
    x = "Normalized HMMR expression",
    y = "Cell count"
  ) +
  theme_pub_cc()

if (is.finite(q_low_val)) {
  p_hist <- p_hist +
    geom_vline(xintercept = q_low_val, color = col_dark, linetype = "dashed", linewidth = 0.8) +
    annotate(
      "text",
      x = q_low_val,
      y = Inf,
      label = paste0("Low tertile cutoff = ", round(q_low_val, 3)),
      vjust = 1.5,
      hjust = 1.05,
      size = 3.8,
      color = col_dark
    )
}

if (is.finite(q_high_val)) {
  p_hist <- p_hist +
    geom_vline(xintercept = q_high_val, color = col_dark, linetype = "dashed", linewidth = 0.8) +
    annotate(
      "text",
      x = q_high_val,
      y = Inf,
      label = paste0("High tertile cutoff = ", round(q_high_val, 3)),
      vjust = 3.0,
      hjust = -0.05,
      size = 3.8,
      color = col_dark
    )
}

p_density <- ggplot(dist_df, aes(x = HMMR_expr, color = HMMR_state, fill = HMMR_state)) +
  geom_density(alpha = 0.18, linewidth = 1.0) +
  scale_color_manual(values = c("HMMR_low" = col_light2, "middle" = col_grey2, "HMMR_high" = col_dark)) +
  scale_fill_manual(values = c("HMMR_low" = col_light2, "middle" = col_grey2, "HMMR_high" = col_dark)) +
  labs(
    title = NULL,
    x = "Normalized HMMR expression",
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme_pub_cc() +
  theme(legend.position = "top")

if (is.finite(q_low_val)) {
  p_density <- p_density +
    geom_vline(xintercept = q_low_val, color = col_dark, linetype = "dashed", linewidth = 0.8)
}

if (is.finite(q_high_val)) {
  p_density <- p_density +
    geom_vline(xintercept = q_high_val, color = col_dark, linetype = "dashed", linewidth = 0.8)
}

p_dist <- p_hist / p_density + plot_layout(heights = c(2.2, 1.5))

save_plot_dual(
  p_dist,
  filename = "01.A3_HMMR_expression_distribution",
  width = 10.2,
  height = 8.6
)

# -----------------------------
# 9. Plot 2: HMMR state violin
# -----------------------------
p_violin <- ggplot(dist_df, aes(x = HMMR_state, y = HMMR_expr, fill = HMMR_state)) +
  geom_violin(trim = FALSE, color = col_border, linewidth = 0.35, alpha = 0.92) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", color = col_border, linewidth = 0.35) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.8, fill = "white", color = col_border, stroke = 0.45) +
  scale_fill_manual(values = c("HMMR_low" = col_light2, "middle" = col_grey, "HMMR_high" = col_dark)) +
  labs(
    title = "Definition of HMMR-high and HMMR-low state groups",
    x = NULL,
    y = "Normalized HMMR expression",
    fill = NULL
  ) +
  theme_pub_cc() +
  theme(legend.position = "none")

p_cluster_bar <- ggplot(
  cluster_state_summary_df %>%
    mutate(cluster = factor(cluster, levels = canonical_cluster_order)),
  aes(x = cluster, y = prop, fill = HMMR_state)
) +
  geom_col(position = "fill", width = 0.68, color = col_border, linewidth = 0.30) +
  scale_fill_manual(values = c("HMMR_low" = col_light2, "middle" = col_grey, "HMMR_high" = col_dark)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = NULL,
    x = "Cluster",
    y = "Proportion",
    fill = NULL
  ) +
  theme_pub_cc() +
  theme(legend.position = "top")

p_state <- p_violin / p_cluster_bar + plot_layout(heights = c(2.4, 1.3))

save_plot_dual(
  p_state,
  filename = "02.A3_HMMR_state_group_violin",
  width = 9.8,
  height = 8.4
)

# -----------------------------
# 10. Done
# -----------------------------
cat("\nA3 HMMR state groups defined successfully.\n")
cat("Detected condition column: ", condition_col, "\n")
cat("Detected cluster column: ", cluster_col, "\n")
cat("Target assay fixed to: ", target_assay, "\n")
cat("Target cell count: ", nrow(target_meta), "\n")

if (is.finite(q_low_val)) {
  cat("Low tertile cutoff: ", round(q_low_val, 4), "\n")
} else {
  cat("Low tertile cutoff: NA\n")
}

if (is.finite(q_high_val)) {
  cat("High tertile cutoff: ", round(q_high_val, 4), "\n")
} else {
  cat("High tertile cutoff: NA\n")
}

cat("Tables saved to:\n", table_dir, "\n")
cat("Plots saved to:\n", plot_dir, "\n")
cat("Intermediate objects saved to:\n", intermediate_dir, "\n")