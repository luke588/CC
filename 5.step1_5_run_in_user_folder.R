rm(list = ls())
options(stringsAsFactors = FALSE)

source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(grid)
  library(stringr)
  library(tidyr)
  library(scales)
})

# =====================================
# 1. USER SETTINGS
# =====================================
base_dir <- "/Users/donglinlu/Desktop/CC1/4.epithelial subset"
outdir <- base_dir
object_name_in_RData <- NULL

seurat_candidates <- c(
  file.path(base_dir, "02.Epithelial_subset.manualAnnotation.Rdata"),
  file.path(base_dir, "04.cleanState.priority.Rdata"),
  file.path(base_dir, "05.highRisk_signature.Rdata"),
  file.path(base_dir, "Epithelial_subset.Rdata"),
  file.path(base_dir, "Seurat.Rdata")
)

manual_ann_candidates <- c(
  file.path(base_dir, "02.epiSubset.manualClusterAnn.txt"),
  file.path(base_dir, "01.epiSubset.clusterAnn.txt")
)

malig_marker_candidates <- c(
  file.path(base_dir, "05.cluster9_malignant_squamous_epithelial.markers.sig.txt"),
  file.path(base_dir, "05.cluster9_malignant_squamous_epithelial.markers.all.txt"),
  file.path(base_dir, "05.superficial_keratinizing_epithelial.markers.sig.txt"),
  file.path(base_dir, "05.superficial_keratinizing_epithelial.markers.all.txt"),
  file.path(base_dir, "04.superficial_keratinizing_epithelial.markers.txt")
)

emt_marker_candidates <- c(
  file.path(base_dir, "05.cluster10_EMT_like_epithelial.markers.sig.txt"),
  file.path(base_dir, "05.cluster10_EMT_like_epithelial.markers.all.txt"),
  file.path(base_dir, "05.EMT_like_epithelial.markers.sig.txt"),
  file.path(base_dir, "05.EMT_like_epithelial.markers.all.txt"),
  file.path(base_dir, "04.EMT_like_epithelial.markers.txt")
)

prolif_marker_candidates <- c(
  file.path(base_dir, "05.cluster6_proliferation_core_epithelial.markers.sig.txt"),
  file.path(base_dir, "05.cluster6_proliferation_core_epithelial.markers.all.txt"),
  file.path(base_dir, "05.proliferative_epithelial.markers.sig.txt"),
  file.path(base_dir, "05.proliferative_epithelial.markers.all.txt"),
  file.path(base_dir, "04.proliferative_epithelial.markers.txt")
)

sig_n <- 50
min_pct_diff <- 0.10
min_logfc <- 0.25

state_col_manual  <- NULL
sample_col_manual <- NULL
stage_col_manual  <- NULL

stage_levels <- c("NO_HPV", "N_HPV", "HSIL_HPV", "CA_HPV")

# =====================================
# 2. HELPERS
# =====================================
pick_first_existing <- function(x, label, allow_missing = FALSE) {
  hit <- x[file.exists(x)]
  if (length(hit) == 0) {
    if (allow_missing) return(NULL)
    stop("Cannot find ", label, ". Checked: ", paste(x, collapse = " | "))
  }
  hit[1]
}

load_seurat_object <- function(path, object_name_in_RData = NULL) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (!inherits(obj, "Seurat")) stop("RDS file does not contain a Seurat object.")
    return(obj)
  }

  if (grepl("\\.RData$|\\.rdata$", path, ignore.case = TRUE)) {
    e <- new.env(parent = emptyenv())
    load(path, envir = e)
    nms <- ls(e)

    if (!is.null(object_name_in_RData)) {
      if (!object_name_in_RData %in% nms) {
        stop("object_name_in_RData not found in RData: ", object_name_in_RData)
      }
      obj <- e[[object_name_in_RData]]
      if (!inherits(obj, "Seurat")) stop("Specified object is not a Seurat object.")
      return(obj)
    }

    seurat_nms <- nms[sapply(nms, function(x) inherits(e[[x]], "Seurat"))]
    if (length(seurat_nms) == 0) stop("No Seurat object found in: ", path)
    message("Auto-detected Seurat object: ", seurat_nms[1])
    return(e[[seurat_nms[1]]])
  }

  stop("Unsupported file format: ", path)
}

find_first_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

find_fc_col <- function(df) {
  hit <- intersect(c("avg_log2FC", "avg_logFC", "log2FC", "logFC"), colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_marker_table <- function(path) {
  x <- read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

  if (!"gene" %in% colnames(x)) {
    if (!is.null(rownames(x)) && all(rownames(x) != "")) {
      x$gene <- rownames(x)
    } else {
      stop("Missing gene column in marker file: ", path)
    }
  }

  fc_col <- find_fc_col(x)
  if (is.na(fc_col)) {
    stop("Cannot find logFC column in marker file: ", path)
  }

  if (!all(c("pct.1", "pct.2") %in% colnames(x))) {
    stop("Missing pct.1 / pct.2 in marker file: ", path)
  }

  x <- x %>%
    mutate(
      avg_log2FC_use = as.numeric(.data[[fc_col]]),
      pct_diff = as.numeric(pct.1) - as.numeric(pct.2)
    ) %>%
    arrange(desc(avg_log2FC_use))

  x
}

make_signature <- function(tab, sig_n = 50, min_pct_diff = 0.10, min_logfc = 0.25) {
  tab %>%
    filter(avg_log2FC_use >= min_logfc, pct_diff >= min_pct_diff) %>%
    arrange(desc(avg_log2FC_use), desc(pct_diff)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    slice_head(n = sig_n) %>%
    pull(gene)
}

auto_pick_col <- function(meta_cols, candidates, label) {
  hit <- intersect(candidates, meta_cols)
  if (length(hit) == 0) return(NULL)
  message("Auto-detected ", label, ": ", hit[1])
  hit[1]
}

read_manual_mapping <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ann <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("cluster", "label") %in% colnames(ann))) return(NULL)
  ann$cluster <- as.character(ann$cluster)
  ann$label <- as.character(ann$label)
  ann
}

ensure_state_labels <- function(obj, state_col, manual_map = NULL) {
  meta_df <- obj@meta.data

  if (!is.null(state_col) && state_col %in% colnames(meta_df)) {
    vals <- as.character(meta_df[[state_col]])
    if (length(vals) > 0 && any(grepl("[A-Za-z_]", vals))) {
      obj@meta.data[[state_col]] <- vals
      return(list(obj = obj, state_col = state_col))
    }
  }

  if ("manualCellType" %in% colnames(meta_df)) {
    obj@meta.data[["manualCellType"]] <- as.character(meta_df[["manualCellType"]])
    return(list(obj = obj, state_col = "manualCellType"))
  }

  if (!is.null(manual_map) && "seurat_clusters" %in% colnames(meta_df)) {
    cl <- as.character(meta_df$seurat_clusters)
    mapped <- manual_map$label[match(cl, manual_map$cluster)]
    if (sum(!is.na(mapped)) > 0) {
      obj@meta.data[["manualCellType"]] <- mapped
      return(list(obj = obj, state_col = "manualCellType"))
    }
  }

  ident_vals <- as.character(Idents(obj))
  if (length(ident_vals) > 0 && any(grepl("[A-Za-z_]", ident_vals))) {
    obj@meta.data[["active_ident_label"]] <- ident_vals
    return(list(obj = obj, state_col = "active_ident_label"))
  }

  stop("Cannot resolve state labels. Please provide valid state labels or manual annotation mapping.")
}

make_state_levels <- function(meta_df, state_col) {
  state_vals <- as.character(meta_df[[state_col]])

  if ("seurat_clusters" %in% colnames(meta_df)) {
    tmp <- data.frame(
      seurat_clusters = as.character(meta_df$seurat_clusters),
      state = state_vals,
      stringsAsFactors = FALSE
    ) %>%
      distinct() %>%
      mutate(seurat_clusters_num = suppressWarnings(as.numeric(seurat_clusters))) %>%
      arrange(seurat_clusters_num, state)

    lev <- unique(tmp$state)
    if (length(lev) > 0) return(lev)
  }

  unique(state_vals)
}

plot_violin <- function(obj, feature, group.by, file,
                        width = 13.2, height = 7.2,
                        title = NULL, ylab = NULL,
                        x_angle = 45) {
  if (is.null(title)) title <- feature
  if (is.null(ylab)) ylab <- feature

  state_levels <- make_state_levels(obj@meta.data, group.by)
  obj@meta.data[[group.by]] <- factor(as.character(obj@meta.data[[group.by]]), levels = state_levels)

  png_file <- sub("\\.pdf$", ".png", file)

  p <- VlnPlot(
    obj,
    features = feature,
    group.by = group.by,
    pt.size = 0
  ) +
    scale_x_discrete(limits = state_levels) +
    labs(
      title = title,
      x = NULL,
      y = ylab
    ) +
    theme_pub() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = x_angle,
        hjust = 1,
        vjust = 1,
        size = 12,
        colour = "black"
      )
    )

  save_pub(
    plot_obj = p,
    pdf_file = file,
    png_file = png_file,
    width = width,
    height = height,
    dpi = 600
  )
}

normalize_stage <- function(x) {
  x0 <- as.character(x)
  x1 <- toupper(trimws(x0))
  x1 <- gsub("[[:space:]]+", "", x1)
  x1 <- gsub("_EPITHELIALCELLS", "", x1)
  x1 <- gsub("EPITHELIALCELLS", "", x1)
  x1 <- gsub("[^A-Z0-9]+", "_", x1)

  out <- ifelse(
    grepl("NO_HPV|NOHPV", x1),
    "NO_HPV",
    ifelse(
      grepl("^N_HPV$|^NHPV$|N_HPV|NHPV", x1),
      "N_HPV",
      ifelse(
        grepl("HSIL_HPV|HSILHPV", x1),
        "HSIL_HPV",
        ifelse(
          grepl("CA_HPV|CAHPV", x1),
          "CA_HPV",
          NA_character_
        )
      )
    )
  )

  out
}

score_signature <- function(obj, genes, name) {
  genes_use <- intersect(unique(genes), rownames(obj))
  if (length(genes_use) < 5) {
    warning(name, ": matched genes < 5. Score set to NA.")
    obj[[paste0(name, "_score")]] <- NA_real_
    return(obj)
  }

  obj <- AddModuleScore(
    obj,
    features = list(genes_use),
    name = paste0(name, "_tmp_")
  )

  obj[[paste0(name, "_score")]] <- obj[[paste0(name, "_tmp_1")]]
  obj[[paste0(name, "_tmp_1")]] <- NULL
  obj
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

pick_available_state <- function(candidates, available_states, label) {
  hit <- intersect(candidates, available_states)
  if (length(hit) == 0) {
    stop("Cannot resolve state for ", label, ". Candidates: ", paste(candidates, collapse = " | "))
  }
  hit[1]
}

# =====================================
# 3. RESOLVE INPUTS
# =====================================
seurat_path <- pick_first_existing(seurat_candidates, "Seurat object file")
manual_ann_file <- pick_first_existing(manual_ann_candidates, "manual annotation mapping file", allow_missing = TRUE)

malig_marker_file  <- pick_first_existing(malig_marker_candidates, "malignant_squamous marker file")
emt_marker_file    <- pick_first_existing(emt_marker_candidates, "EMT_like marker file")
prolif_marker_file <- pick_first_existing(prolif_marker_candidates, "proliferation_core marker file")

message("Using Seurat file: ", seurat_path)
message("Using malignant_squamous markers: ", malig_marker_file)
message("Using EMT-like markers: ", emt_marker_file)
message("Using proliferation-core markers: ", prolif_marker_file)
if (!is.null(manual_ann_file)) {
  message("Using manual annotation mapping: ", manual_ann_file)
}

# =====================================
# 4. LOAD OBJECT AND RESOLVE METADATA
# =====================================
obj <- load_seurat_object(seurat_path, object_name_in_RData)
meta_cols <- colnames(obj@meta.data)
writeLines(meta_cols, file.path(outdir, "step1_5.meta_columns.txt"))

state_col <- state_col_manual
sample_col <- sample_col_manual
stage_col <- stage_col_manual

if (is.null(state_col)) {
  state_col <- auto_pick_col(
    meta_cols,
    c("manualCellType", "cleanCellType", "manualClusterAnn", "celltype", "CellType", "annotation", "cellAnn", "seurat_clusters"),
    "state_col"
  )
}
if (is.null(sample_col)) {
  sample_col <- auto_pick_col(
    meta_cols,
    c("sample", "orig.ident", "Sample", "sample_id", "patient", "Patient"),
    "sample_col"
  )
}
if (is.null(stage_col)) {
  stage_col <- auto_pick_col(
    meta_cols,
    c("group", "Group", "stage", "Stage", "disease_group", "condition", "Condition"),
    "stage_col"
  )
}

manual_map <- read_manual_mapping(manual_ann_file)
resolved <- ensure_state_labels(obj, state_col, manual_map = manual_map)
obj <- resolved$obj
state_col <- resolved$state_col

message("Final state column used for plotting: ", state_col)

available_states <- unique(as.character(obj@meta.data[[state_col]]))

program.info <- data.frame(
  cluster_id = c("9", "10", "6"),
  program_label = c(
    "malignant_squamous_epithelial",
    "EMT_like_epithelial",
    "proliferation_core_epithelial"
  ),
  output_prefix = c(
    "cluster9_malignant_squamous_epithelial",
    "cluster10_EMT_like_epithelial",
    "cluster6_proliferation_core_epithelial"
  ),
  score_prefix = c(
    "malignant_squamous_state",
    "emt_state",
    "proliferation_core_state"
  ),
  stringsAsFactors = FALSE
)

program.info$object_state <- c(
  pick_available_state(
    candidates = c("malignant_squamous_epithelial", "superficial_keratinizing_epithelial"),
    available_states = available_states,
    label = "cluster9 / malignant_squamous_epithelial"
  ),
  pick_available_state(
    candidates = c("EMT_like_epithelial"),
    available_states = available_states,
    label = "cluster10 / EMT_like_epithelial"
  ),
  pick_available_state(
    candidates = c("proliferation_core_epithelial", "proliferative_epithelial"),
    available_states = available_states,
    label = "cluster6 / proliferation_core_epithelial"
  )
)

program.info$marker_file <- c(
  malig_marker_file,
  emt_marker_file,
  prolif_marker_file
)

program.info$color <- c(
  "#D73027",
  "#4575B4",
  "#1A9850"
)

write_csv(
  program.info,
  file.path(outdir, "step1_5.program_mapping.csv")
)

# =====================================
# 5. BUILD SIGNATURES
# =====================================
marker_tables <- list()
signature_list <- list()

for (i in seq_len(nrow(program.info))) {
  marker_tab <- read_marker_table(program.info$marker_file[i])
  sig_genes <- make_signature(
    marker_tab,
    sig_n = sig_n,
    min_pct_diff = min_pct_diff,
    min_logfc = min_logfc
  )

  marker_tables[[program.info$program_label[i]]] <- marker_tab
  signature_list[[program.info$program_label[i]]] <- sig_genes

  writeLines(
    sig_genes,
    file.path(
      outdir,
      paste0("step1_5.", program.info$output_prefix[i], ".signature.top50.txt")
    )
  )
}

# =====================================
# 6. STEP 1.5-1 CELL CYCLE
# =====================================
cc_genes <- Seurat::cc.genes.updated.2019

obj <- CellCycleScoring(
  object = obj,
  s.features = cc_genes$s.genes,
  g2m.features = cc_genes$g2m.genes,
  set.ident = FALSE
)

plot_violin(
  obj,
  feature = "S.Score",
  group.by = state_col,
  file = file.path(outdir, "step1_5.01_S.Score_by_state.pdf"),
  title = "S phase score across epithelial states",
  ylab = "S.Score"
)

plot_violin(
  obj,
  feature = "G2M.Score",
  group.by = state_col,
  file = file.path(outdir, "step1_5.02_G2M.Score_by_state.pdf"),
  title = "G2/M phase score across epithelial states",
  ylab = "G2M.Score"
)

cc_summary <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(.data[[state_col]]) %>%
  summarise(
    n = n(),
    mean_S = safe_mean(S.Score),
    mean_G2M = safe_mean(G2M.Score),
    mean_cycle = safe_mean(S.Score + G2M.Score),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cycle))

write_csv(
  cc_summary,
  file.path(outdir, "step1_5.01_cell_cycle_summary_by_state.csv")
)

# =====================================
# 7. STEP 1.5-2 SIGNATURE SCORING
# =====================================
for (i in seq_len(nrow(program.info))) {
  obj <- score_signature(
    obj = obj,
    genes = signature_list[[program.info$program_label[i]]],
    name = program.info$score_prefix[i]
  )
}

plot_violin(
  obj,
  feature = "malignant_squamous_state_score",
  group.by = state_col,
  file = file.path(outdir, "step1_5.03_malignant_squamous_score_by_state.pdf"),
  title = "Malignant squamous signature score across epithelial states",
  ylab = "malignant_squamous_score"
)

plot_violin(
  obj,
  feature = "emt_state_score",
  group.by = state_col,
  file = file.path(outdir, "step1_5.04_EMT_like_score_by_state.pdf"),
  title = "EMT-like signature score across epithelial states",
  ylab = "EMT_like_score"
)

plot_violin(
  obj,
  feature = "proliferation_core_state_score",
  group.by = state_col,
  file = file.path(outdir, "step1_5.05_proliferation_core_score_by_state.pdf"),
  title = "Proliferation-core signature score across epithelial states",
  ylab = "proliferation_core_score"
)

sig_summary <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(.data[[state_col]]) %>%
  summarise(
    n = n(),
    mean_malignant_squamous = safe_mean(malignant_squamous_state_score),
    mean_EMT_like = safe_mean(emt_state_score),
    mean_proliferation_core = safe_mean(proliferation_core_state_score),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_proliferation_core))

write_csv(
  sig_summary,
  file.path(outdir, "step1_5.02_signature_score_summary_by_state.csv")
)

# =====================================
# 8. STEP 1.5-3 ALL-STATE STAGE PROPORTION
# =====================================
meta_df_plain <- obj@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    state_label = as.character(.data[[state_col]])
  )

if (!is.null(stage_col) && stage_col %in% colnames(meta_df_plain)) {
  meta_df_plain <- meta_df_plain %>%
    dplyr::mutate(
      stage_std = normalize_stage(.data[[stage_col]])
    )

  prop_df <- meta_df_plain %>%
    dplyr::filter(!is.na(stage_std), !is.na(state_label)) %>%
    dplyr::group_by(stage_std, state_label) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(stage_std) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      group = factor(stage_std, levels = stage_levels),
      state_label = state_label,
      n = n,
      prop = prop
    ) %>%
    dplyr::arrange(group, state_label)

  write_csv(
    prop_df,
    file.path(outdir, "step1_5.03_state_proportion_by_stage.csv")
  )
}

# =====================================
# 9. STEP 1.5-4 / 06 / 07 / 08 FOCUS PROGRAMS
# =====================================
obj@meta.data$focus_program_label <- NA_character_
obj@meta.data$focus_cluster_id <- NA_character_

for (i in seq_len(nrow(program.info))) {
  state_i <- program.info$object_state[i]
  program_i <- program.info$program_label[i]
  cluster_i <- program.info$cluster_id[i]

  hit <- as.character(obj@meta.data[[state_col]]) == state_i
  obj@meta.data$focus_program_label[hit] <- program_i
  obj@meta.data$focus_cluster_id[hit] <- cluster_i
}

focus_df <- obj@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    state_label = as.character(.data[[state_col]]),
    focus_program_label = as.character(focus_program_label),
    focus_cluster_id = as.character(focus_cluster_id)
  ) %>%
  dplyr::filter(!is.na(focus_program_label), !is.na(focus_cluster_id))

focus_summary <- focus_df %>%
  dplyr::group_by(focus_cluster_id, focus_program_label, state_label) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_S = safe_mean(S.Score),
    mean_G2M = safe_mean(G2M.Score),
    mean_cycle = safe_mean(S.Score + G2M.Score),
    mean_malignant_squamous = safe_mean(malignant_squamous_state_score),
    mean_EMT_like = safe_mean(emt_state_score),
    mean_proliferation_core = safe_mean(proliferation_core_state_score),
    .groups = "drop"
  ) %>%
  dplyr::transmute(
    seurat_clusters = focus_cluster_id,
    program_label = focus_program_label,
    state_label = state_label,
    n = n,
    mean_S = mean_S,
    mean_G2M = mean_G2M,
    mean_cycle = mean_cycle,
    mean_malignant_squamous = mean_malignant_squamous,
    mean_EMT_like = mean_EMT_like,
    mean_proliferation_core = mean_proliferation_core
  ) %>%
  dplyr::arrange(suppressWarnings(as.numeric(seurat_clusters)), program_label, state_label)

write_csv(
  focus_summary,
  file.path(outdir, "step1_5.06_focus_clusters_6_9_10_summary.csv")
)

if (!is.null(stage_col) && stage_col %in% colnames(focus_df)) {
  focus_df2 <- focus_df %>%
    dplyr::mutate(
      stage_std = normalize_stage(.data[[stage_col]])
    )

  focus_stage_prop <- focus_df2 %>%
    dplyr::filter(!is.na(stage_std)) %>%
    dplyr::group_by(stage_std, focus_cluster_id) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(stage_std) %>%
    dplyr::mutate(cluster_prop = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      group = factor(stage_std, levels = stage_levels),
      seurat_clusters = focus_cluster_id,
      n = n,
      cluster_prop = cluster_prop
    ) %>%
    dplyr::arrange(group, suppressWarnings(as.numeric(seurat_clusters)))

  write_csv(
    focus_stage_prop,
    file.path(outdir, "step1_5.07_focus_clusters_6_9_10_stage_proportion.csv")
  )

  focus_stage_prop_annot <- focus_df2 %>%
    dplyr::filter(!is.na(stage_std)) %>%
    dplyr::group_by(stage_std, focus_cluster_id, focus_program_label, state_label) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(stage_std) %>%
    dplyr::mutate(cluster_prop = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      group = factor(stage_std, levels = stage_levels),
      seurat_clusters = focus_cluster_id,
      program_label = focus_program_label,
      state_label = state_label,
      n = n,
      cluster_prop = cluster_prop
    ) %>%
    dplyr::arrange(group, suppressWarnings(as.numeric(seurat_clusters)), program_label, state_label)

  write_csv(
    focus_stage_prop_annot,
    file.path(outdir, "step1_5.07_focus_clusters_6_9_10_stage_proportion.annotated.csv")
  )
}

focus_composition <- focus_df %>%
  dplyr::group_by(focus_cluster_id, focus_program_label, state_label) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(focus_cluster_id) %>%
  dplyr::mutate(cellType_prop = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    seurat_clusters = focus_cluster_id,
    program_label = focus_program_label,
    state_label = state_label,
    n = n,
    cellType_prop = cellType_prop
  ) %>%
  dplyr::arrange(suppressWarnings(as.numeric(seurat_clusters)), program_label, state_label)

write_csv(
  focus_composition,
  file.path(outdir, "step1_5.08_focus_clusters_6_9_10_cellType_composition.csv")
)

# =====================================
# 10. SAVE OBJECT
# =====================================
saveRDS(
  obj,
  file.path(outdir, "step1_5.scored_object.rds")
)

# =====================================
# 11. DONE
# =====================================
message("Done. All outputs saved to: ", outdir)
message("Final state column used: ", state_col)
message("Current mainline used:")
message("cluster9 -> malignant_squamous_epithelial")
message("cluster10 -> EMT_like_epithelial")
message("cluster6 -> proliferation_core_epithelial")
message("")
message("Main outputs:")
message("step1_5.01_S.Score_by_state.pdf / .png")
message("step1_5.02_G2M.Score_by_state.pdf / .png")
message("step1_5.03_malignant_squamous_score_by_state.pdf / .png")
message("step1_5.04_EMT_like_score_by_state.pdf / .png")
message("step1_5.05_proliferation_core_score_by_state.pdf / .png")
message("step1_5.01_cell_cycle_summary_by_state.csv")
message("step1_5.02_signature_score_summary_by_state.csv")
message("step1_5.03_state_proportion_by_stage.csv")
message("step1_5.06_focus_clusters_6_9_10_summary.csv")
message("step1_5.07_focus_clusters_6_9_10_stage_proportion.csv")
message("step1_5.07_focus_clusters_6_9_10_stage_proportion.annotated.csv")
message("step1_5.08_focus_clusters_6_9_10_cellType_composition.csv")
message("step1_5.cluster9_malignant_squamous_epithelial.signature.top50.txt")
message("step1_5.cluster10_EMT_like_epithelial.signature.top50.txt")
message("step1_5.cluster6_proliferation_core_epithelial.signature.top50.txt")
message("step1_5.scored_object.rds")