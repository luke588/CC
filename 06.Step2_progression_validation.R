source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
# ============================================================
# Step 2: progression validation for high-risk epithelial states
# Author: OpenAI
# Version: single-script, run-ready
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
})

# =========================
# USER SETTINGS
# =========================
input_dir  <- "/Users/donglinlu/Desktop/CC1/4.epithelial subset"
output_dir <- "/Users/donglinlu/Desktop/CC1/5.progression_validation"

# fixed metadata columns based on your previous run
cluster_col <- "seurat_clusters"
stage_col   <- "group"
sample_col  <- "Sample"

# cluster-to-state mapping
cluster9_label <- "malignant_squamous_epithelial"
cluster10_label <- "EMT_like_epithelial"
cluster6_label <- "proliferation_core_epithelial"
other_label <- "other_epithelial_states"

# stage order
stage_levels <- c("NO_HPV", "N_HPV", "HSIL_HPV", "CA_HPV")

# marker files to export signatures
cluster6_marker_file  <- file.path(input_dir, "step1_5.cluster_6.markers_vs_rest.csv")
cluster9_marker_file  <- file.path(input_dir, "step1_5.cluster_9.markers_vs_rest.csv")
cluster10_marker_file <- file.path(input_dir, "step1_5.cluster_10.markers_vs_rest.csv")

# optional minimum filters for signature export
min_pct1 <- 0.20
max_pct2 <- 0.50
min_logfc <- 0.25
n_top_sig <- 50

# =========================
# PREPARE FOLDERS
# =========================
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "input"), showWarnings = FALSE)
dir.create(file.path(output_dir, "output"), showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "signatures"), showWarnings = FALSE)

# =========================
# HELPERS
# =========================
find_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

load_seurat_object <- function(rdata_path) {
  e <- new.env(parent = globalenv())
  loaded_names <- load(rdata_path, envir = e)
  if (length(loaded_names) == 0) stop("No objects found in RData: ", rdata_path)

  obj_names <- loaded_names[sapply(loaded_names, function(x) inherits(e[[x]], "Seurat"))]
  if (length(obj_names) == 0) {
    stop("No Seurat object found in: ", rdata_path)
  }

  obj <- e[[obj_names[1]]]
  attr(obj, "loaded_object_name") <- obj_names[1]
  obj
}

safe_write_csv <- function(x, path) {
  readr::write_csv(x, path)
}

make_signature <- function(marker_file, n_top = 50,
                           min_pct1 = 0.20,
                           max_pct2 = 0.50,
                           min_logfc = 0.25) {
  if (!file.exists(marker_file)) stop("Missing marker file: ", marker_file)

  df <- suppressMessages(readr::read_csv(marker_file, show_col_types = FALSE))

  need_cols <- c("gene", "avg_log2FC", "pct.1", "pct.2")
  if (!all(need_cols %in% colnames(df))) {
    stop("Marker file missing required columns: ", marker_file)
  }

  df2 <- df %>%
    filter(!is.na(gene)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    filter(avg_log2FC >= min_logfc,
           pct.1 >= min_pct1,
           pct.2 <= max_pct2) %>%
    arrange(desc(avg_log2FC), desc(pct.1), pct.2)

  if (nrow(df2) == 0) {
    df2 <- df %>%
      filter(!is.na(gene)) %>%
      distinct(gene, .keep_all = TRUE) %>%
      arrange(desc(avg_log2FC), desc(pct.1), pct.2)
  }

  sig <- head(df2$gene, n_top)
  sig <- sig[!is.na(sig) & sig != ""]
  unique(sig)
}

# =========================
# LOAD OBJECT
# =========================
obj_candidates <- c(
  file.path(input_dir, "05.highRisk_signature.Rdata"),
  file.path(input_dir, "02.Epithelial_subset.manualAnnotation.Rdata"),
  file.path(input_dir, "Epithelial_subset.Rdata"),
  file.path(input_dir, "04.cleanState.priority.Rdata"),
  file.path(input_dir, "Seurat.Rdata")
)

obj_file <- find_first_existing(obj_candidates)
if (is.na(obj_file)) {
  stop("No supported RData object found in input_dir: ", input_dir)
}

seu <- load_seurat_object(obj_file)
loaded_object_name <- attr(seu, "loaded_object_name")
meta <- seu@meta.data
meta$cell_id <- rownames(meta)

# save metadata columns for record
writeLines(colnames(meta), file.path(output_dir, "output", "06.0_meta_columns.txt"))

# check metadata columns
for (nm in c(cluster_col, stage_col, sample_col)) {
  if (!nm %in% colnames(meta)) stop("Missing metadata column: ", nm)
}

# =========================
# ADD HIGH-RISK STATE LABEL
# =========================
meta <- meta %>%
  mutate(
    cluster_chr = as.character(.data[[cluster_col]]),
    highRisk_state = case_when(
      cluster_chr == "9"  ~ cluster9_label,
      cluster_chr == "10" ~ cluster10_label,
      cluster_chr == "6"  ~ cluster6_label,
      TRUE ~ other_label
    )
  )

meta[[stage_col]] <- factor(meta[[stage_col]], levels = stage_levels)
meta$highRisk_state <- factor(
  meta$highRisk_state,
  levels = c(cluster9_label, cluster10_label, cluster6_label, other_label)
)

seu$highRisk_state <- meta$highRisk_state

# =========================
# 1. COUNT / PROPORTION BY STAGE (CELL LEVEL)
# =========================
state_count_by_stage <- meta %>%
  count(.data[[stage_col]], highRisk_state, name = "n_cells") %>%
  rename(stage = all_of(stage_col))

state_prop_by_stage <- state_count_by_stage %>%
  group_by(stage) %>%
  mutate(stage_total = sum(n_cells), proportion = n_cells / stage_total) %>%
  ungroup()

safe_write_csv(state_count_by_stage,
               file.path(output_dir, "output", "06.1_state_count_by_stage.csv"))
safe_write_csv(state_prop_by_stage,
               file.path(output_dir, "output", "06.2_state_proportion_by_stage.csv"))

# =========================
# 2. SAMPLE-LEVEL PROPORTION BY STAGE
# =========================
sample_state_count <- meta %>%
  count(.data[[sample_col]], .data[[stage_col]], highRisk_state, name = "n_cells") %>%
  rename(sample = all_of(sample_col), stage = all_of(stage_col))

sample_total <- meta %>%
  count(.data[[sample_col]], name = "sample_total") %>%
  rename(sample = all_of(sample_col))

sample_state_prop <- sample_state_count %>%
  left_join(sample_total, by = "sample") %>%
  mutate(proportion = n_cells / sample_total)

safe_write_csv(sample_state_prop,
               file.path(output_dir, "output", "06.3_sample_level_state_proportion_by_stage.csv"))

# =========================
# 3. FOCUS STATES ONLY
# =========================
focus_states <- c(cluster9_label, cluster10_label, cluster6_label)
focus_prop_by_stage <- state_prop_by_stage %>%
  filter(highRisk_state %in% focus_states)

safe_write_csv(focus_prop_by_stage,
               file.path(output_dir, "output", "06.4_focus_state_proportion_by_stage.csv"))

# =========================
# 4. CLUSTER-TO-STATE MAP
# =========================
cluster_to_state <- meta %>%
  count(.data[[cluster_col]], highRisk_state, name = "n_cells") %>%
  rename(cluster = all_of(cluster_col)) %>%
  arrange(as.numeric(as.character(cluster)), desc(n_cells))

safe_write_csv(cluster_to_state,
               file.path(output_dir, "output", "06.5_cluster_to_highRisk_state_map.csv"))

# =========================
# 5. EXPORT SIGNATURES
# =========================
cluster6_sig <- make_signature(cluster6_marker_file, n_top = n_top_sig,
                               min_pct1 = min_pct1, max_pct2 = max_pct2, min_logfc = min_logfc)
cluster9_sig <- make_signature(cluster9_marker_file, n_top = n_top_sig,
                               min_pct1 = min_pct1, max_pct2 = max_pct2, min_logfc = min_logfc)
cluster10_sig <- make_signature(cluster10_marker_file, n_top = n_top_sig,
                                min_pct1 = min_pct1, max_pct2 = max_pct2, min_logfc = min_logfc)

writeLines(cluster9_sig,
           file.path(output_dir, "signatures", "cluster9_malignant_squamous_epithelial_top50.txt"))
writeLines(cluster10_sig,
           file.path(output_dir, "signatures", "cluster10_EMT_like_epithelial_top50.txt"))
writeLines(cluster6_sig,
           file.path(output_dir, "signatures", "cluster6_proliferation_core_epithelial_top50.txt"))

sig_summary <- tibble(
  state = c(cluster9_label, cluster10_label, cluster6_label),
  n_genes = c(length(cluster9_sig), length(cluster10_sig), length(cluster6_sig)),
  file = c(
    "cluster9_malignant_squamous_epithelial_top50.txt",
    "cluster10_EMT_like_epithelial_top50.txt",
    "cluster6_proliferation_core_epithelial_top50.txt"
  )
)

safe_write_csv(sig_summary,
               file.path(output_dir, "output", "06.6_signature_export_summary.csv"))

# =========================
# 6. PLOTS
# =========================
p1 <- ggplot(state_prop_by_stage,
             aes(x = stage, y = proportion, fill = highRisk_state)) +
  geom_col(position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Proportion", fill = "highRisk_state",
       title = "High-risk state composition across stages") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "plots", "06.1_state_proportion_by_stage_barplot.pdf"),
       p1, width = 7, height = 5)

p2 <- ggplot(focus_prop_by_stage,
             aes(x = stage, y = proportion, fill = highRisk_state)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Proportion", fill = "highRisk_state",
       title = "Focus high-risk states across stages") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "plots", "06.2_focus_states_by_stage_barplot.pdf"),
       p2, width = 7, height = 5)

p3 <- ggplot(sample_state_prop %>% filter(highRisk_state %in% focus_states),
             aes(x = stage, y = proportion, fill = highRisk_state)) +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
  facet_wrap(~highRisk_state, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = "Sample-level proportion",
       title = "Sample-level proportions of focus states across stages") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(file.path(output_dir, "plots", "06.3_sample_level_focus_states_by_stage_boxplot.pdf"),
       p3, width = 7, height = 8)

# =========================
# 7. SAVE UPDATED OBJECT AND SESSION INFO
# =========================
saveRDS(seu, file.path(output_dir, "output", "06.Step2_progression_labeled_object.rds"))

session_txt <- c(
  paste0("Input object file: ", obj_file),
  paste0("Loaded object name: ", loaded_object_name),
  paste0("n cells: ", ncol(seu)),
  paste0("cluster_col: ", cluster_col),
  paste0("stage_col: ", stage_col),
  paste0("sample_col: ", sample_col),
  paste0("Output dir: ", output_dir),
  "",
  "State mapping:",
  paste0("cluster 9 -> ", cluster9_label),
  paste0("cluster 10 -> ", cluster10_label),
  paste0("cluster 6 -> ", cluster6_label),
  paste0("others -> ", other_label)
)
writeLines(session_txt, file.path(output_dir, "output", "06.7_run_info.txt"))

message("============================================================")
message("Step 2 finished successfully.")
message("Input object: ", obj_file)
message("Output dir: ", output_dir)
message("Key outputs:")
message(" - output/06.2_state_proportion_by_stage.csv")
message(" - output/06.3_sample_level_state_proportion_by_stage.csv")
message(" - signatures/cluster9_malignant_squamous_epithelial_top50.txt")
message(" - signatures/cluster10_EMT_like_epithelial_top50.txt")
message(" - signatures/cluster6_proliferation_core_epithelial_top50.txt")
message("============================================================")
