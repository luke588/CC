# -----------------------------
# HMMR virtual knockout
# three-version framework
# cluster6_only / prioritized_6_9_10_pooled / all_epithelial_CA_HPV
# directly runnable
# -----------------------------

library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
set.seed(12345)

############################
# 0. User settings
############################
target_gene <- "HMMR"
target_stage <- "CA_HPV"
restrict_stage <- TRUE

work_dir <- "/Users/donglinlu/Desktop/CC1/9.scTenifoldKnk"
input_rdata <- file.path(work_dir, "Epithelial_subset.Rdata")
output_root <- file.path(work_dir, "HMMR_three_version_KO")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

# three KO modes
mode_table <- data.frame(
  mode = c(
    "cluster6_only",
    "prioritized_6_9_10_pooled",
    "all_epithelial_CA_HPV"
  ),
  nfeatures = c(3000, 3000, 3000),
  max_nc_nCells = c(300, 500, 800),
  stringsAsFactors = FALSE
)

############################
# 1. Helper functions
############################
theme_pub <- function(base_size = 18) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 5),
      plot.subtitle = element_text(hjust = 0.5, size = base_size + 1, margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = base_size + 2),
      axis.text = element_text(size = base_size, colour = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      plot.margin = margin(t = 14, r = 20, b = 12, l = 14)
    )
}

save_pdf <- function(plot_obj, path, width, height) {
  ggsave(
    filename = path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE
  )
}

pick_first_existing <- function(x, candidates) {
  hit <- candidates[candidates %in% x]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

find_col_case_insensitive <- function(df, candidates) {
  nms <- colnames(df)
  idx <- match(tolower(candidates), tolower(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NULL)
  nms[idx[1]]
}

get_counts_matrix <- function(seu, assay_use = "RNA") {
  mat <- tryCatch(
    GetAssayData(seu, assay = assay_use, layer = "counts"),
    error = function(e) NULL
  )

  if (is.null(mat)) {
    mat <- tryCatch(
      GetAssayData(seu, assay = assay_use, slot = "counts"),
      error = function(e) NULL
    )
  }

  if (is.null(mat)) {
    stop("Cannot extract counts matrix from Seurat object.")
  }

  return(mat)
}

safe_neglog10 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  valid_x <- x[is.finite(x) & !is.na(x) & x > 0]
  floor_x <- if (length(valid_x) > 0) min(valid_x) * 0.1 else 1e-300
  x[!is.finite(x) | is.na(x) | x <= 0] <- floor_x
  -log10(x)
}

write_gene_vector <- function(x, file) {
  x <- unique(as.character(x))
  x <- x[!is.na(x)]
  x <- x[x != ""]
  writeLines(x, con = file)
}

get_mode_cells <- function(meta_df, stage_col, cluster_col, mode_name, target_stage, restrict_stage = TRUE) {
  keep_idx <- rep(TRUE, nrow(meta_df))

  if (restrict_stage) {
    if (is.null(stage_col)) {
      stop("restrict_stage = TRUE, but no stage column found in metadata.")
    }
    keep_idx <- keep_idx & (as.character(meta_df[[stage_col]]) == target_stage)
  }

  if (mode_name == "cluster6_only") {
    keep_idx <- keep_idx & (as.character(meta_df[[cluster_col]]) %in% c("6"))
  }

  if (mode_name == "prioritized_6_9_10_pooled") {
    keep_idx <- keep_idx & (as.character(meta_df[[cluster_col]]) %in% c("6", "9", "10"))
  }

  if (mode_name == "all_epithelial_CA_HPV") {
    keep_idx <- keep_idx
  }

  rownames(meta_df)[keep_idx]
}

############################
# 2. Load Seurat object
############################
if (!file.exists(input_rdata)) {
  stop("Cannot find input file: ", input_rdata)
}

loaded_names <- load(input_rdata)
cat("Loaded objects:\n")
print(loaded_names)

seurat_obj_names <- loaded_names[sapply(loaded_names, function(x) inherits(get(x), "Seurat"))]

if (length(seurat_obj_names) == 0) {
  stop("No Seurat object found in Epithelial_subset.Rdata")
}

obj_name <- seurat_obj_names[1]
pbmc <- get(obj_name)

cat("\nUsing Seurat object:", obj_name, "\n")
cat("Object class:\n")
print(class(pbmc))
cat("Object dimension:\n")
print(dim(pbmc))

meta_cols <- colnames(pbmc@meta.data)
cat("\nMetadata columns:\n")
print(meta_cols)

stage_col <- pick_first_existing(meta_cols, c("stage", "Stage", "group", "Group"))
cluster_col <- pick_first_existing(
  meta_cols,
  c(
    "seurat_clusters",
    "cluster",
    "Cluster",
    "clean_cluster",
    "cleanCluster",
    "integrated_snn_res.0.8",
    "integrated_snn_res.0.6",
    "RNA_snn_res.0.8",
    "RNA_snn_res.0.6"
  )
)

if (is.null(cluster_col)) {
  stop("Cannot find cluster column in metadata.")
}

cat("\nDetected stage column:", ifelse(is.null(stage_col), "NULL", stage_col), "\n")
cat("Detected cluster column:", cluster_col, "\n")

if (!is.null(stage_col)) {
  cat("\nStage distribution:\n")
  print(table(pbmc@meta.data[[stage_col]]))
}

cat("\nCluster distribution:\n")
print(table(pbmc@meta.data[[cluster_col]]))

############################
# 3. Core runner
############################
run_knockout_mode <- function(mode_name, nfeatures_use, max_nc_nCells_use) {

  cat("\n==============================\n")
  cat("Running mode:", mode_name, "\n")
  cat("==============================\n")

  mode_outdir <- file.path(output_root, mode_name)
  dir.create(mode_outdir, recursive = TRUE, showWarnings = FALSE)

  meta_df <- pbmc@meta.data
  cells_use <- get_mode_cells(
    meta_df = meta_df,
    stage_col = stage_col,
    cluster_col = cluster_col,
    mode_name = mode_name,
    target_stage = target_stage,
    restrict_stage = restrict_stage
  )

  cat("Selected cells:", length(cells_use), "\n")

  if (length(cells_use) == 0) {
    stop("No cells found for mode: ", mode_name)
  }

  pbmc_sub <- subset(pbmc, cells = cells_use)

  if ("RNA" %in% Assays(pbmc_sub)) {
    DefaultAssay(pbmc_sub) <- "RNA"
  }

  assay_use <- DefaultAssay(pbmc_sub)
  countMat <- get_counts_matrix(pbmc_sub, assay_use = assay_use)

  if (!(target_gene %in% rownames(countMat))) {
    stop("Target gene not found in count matrix: ", target_gene)
  }

  pbmc_sub <- FindVariableFeatures(
    object = pbmc_sub,
    selection.method = "vst",
    nfeatures = nfeatures_use
  )

  hvgs <- VariableFeatures(pbmc_sub)
  feature_use <- unique(c(target_gene, hvgs))
  feature_use <- feature_use[feature_use %in% rownames(countMat)]

  data_use <- as.data.frame(as.matrix(countMat[feature_use, , drop = FALSE]))

  cat("Input matrix:", nrow(data_use), "genes x", ncol(data_use), "cells\n")

  use_nCells <- min(max_nc_nCells_use, ncol(data_use))
  cat("nc_nCells used:", use_nCells, "\n")

  result <- scTenifoldKnk(
    countMatrix = data_use,
    gKO = target_gene,
    qc_mtThreshold = 0.15,
    qc_minLSize = 1500,
    nc_nNet = 10,
    nc_nCells = use_nCells
  )

  save(
    result,
    file = file.path(mode_outdir, paste0(target_gene, "_", mode_name, "_scTenifoldKnk_result.Rdata"))
  )

  ############################
  # extract diffRegulation robustly
  ############################
  if (is.null(result$diffRegulation)) {
    stop("result$diffRegulation is NULL for mode: ", mode_name)
  }

  df <- as.data.frame(result$diffRegulation)

  if (nrow(df) == 0) {
    stop("diffRegulation result is empty for mode: ", mode_name)
  }

  if (!("gene" %in% colnames(df))) {
    df$gene <- rownames(df)
  }

  if (all(is.na(df$gene)) || all(df$gene == "")) {
    df$gene <- rownames(df)
  }

  z_col <- find_col_case_insensitive(df, c("Z", "z", "zscore", "z_score", "z.score"))
  padj_col <- find_col_case_insensitive(df, c("p.adj", "padj", "adj.p", "adj_p", "fdr", "FDR", "qvalue", "q.value"))
  p_col <- find_col_case_insensitive(df, c("p", "p.value", "p_value", "pvalue", "P.Value"))

  if (is.null(z_col)) {
    stop("Cannot find Z-score column in diffRegulation for mode: ", mode_name)
  }

  df$Z <- suppressWarnings(as.numeric(df[[z_col]]))

  if (!is.null(padj_col)) {
    df$p.adj <- suppressWarnings(as.numeric(df[[padj_col]]))
  } else {
    df$p.adj <- NA_real_
  }

  if (!is.null(p_col)) {
    df$p <- suppressWarnings(as.numeric(df[[p_col]]))
  } else {
    df$p <- NA_real_
  }

  if (all(is.na(df$p.adj)) && any(!is.na(df$p))) {
    df$p.adj <- p.adjust(df$p, method = "BH")
  }

  df <- df[!is.na(df$Z), , drop = FALSE]
  df <- df[df$gene != target_gene, , drop = FALSE]
  df <- df %>% arrange(desc(abs(Z)))

  write.table(
    df,
    file = file.path(mode_outdir, paste0(target_gene, "_", mode_name, "_all_diffRegulation.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  if (any(!is.na(df$p.adj))) {
    sig_df <- df %>% filter(!is.na(p.adj) & p.adj < 0.05)
  } else if (any(!is.na(df$p))) {
    sig_df <- df %>% filter(!is.na(p) & p < 0.05)
  } else {
    sig_df <- df[0, , drop = FALSE]
  }

  write.table(
    sig_df,
    file = file.path(mode_outdir, paste0(target_gene, "_", mode_name, "_sigDiff.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  if (any(!is.na(df$p))) {
    sig_raw_df <- df %>% filter(!is.na(p) & p < 0.05)
  } else {
    sig_raw_df <- df[0, , drop = FALSE]
  }

  write.table(
    sig_raw_df,
    file = file.path(mode_outdir, paste0(target_gene, "_", mode_name, "_rawP_sigDiff.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  top50_df <- df %>% slice_head(n = min(50, nrow(df)))
  write.table(
    top50_df,
    file = file.path(mode_outdir, paste0(target_gene, "_", mode_name, "_top50_perturbed_genes.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  ############################
  # top100 rank-based gene lists for downstream drug screening
  ############################
  top_up <- df %>%
    filter(Z > 0) %>%
    arrange(desc(Z), p.adj, p) %>%
    slice_head(n = min(100, sum(df$Z > 0, na.rm = TRUE)))

  top_down <- df %>%
    filter(Z < 0) %>%
    arrange(Z, p.adj, p) %>%
    slice_head(n = min(100, sum(df$Z < 0, na.rm = TRUE)))

  write_gene_vector(
    top_up$gene,
    file.path(mode_outdir, paste0(target_gene, "_", mode_name, "_rankTop100_up_genes.txt"))
  )

  write_gene_vector(
    top_down$gene,
    file.path(mode_outdir, paste0(target_gene, "_", mode_name, "_rankTop100_down_genes.txt"))
  )

  ############################
  # scatter plot
  ############################
  plot_p <- df$p.adj
  if (all(is.na(plot_p))) {
    plot_p <- df$p
  }

  df$log_p_plot <- safe_neglog10(plot_p)
  df$significant <- "Not significant"

  if (any(!is.na(df$p.adj))) {
    df$significant[!is.na(df$p.adj) & df$p.adj < 0.05] <- "Significant"
  } else if (any(!is.na(df$p))) {
    df$significant[!is.na(df$p) & df$p < 0.05] <- "Significant"
  }

  sig_df_plot <- df %>%
    filter(significant == "Significant") %>%
    arrange(desc(log_p_plot), desc(abs(Z)))

  x_min <- min(df$Z, na.rm = TRUE)
  x_max <- max(df$Z, na.rm = TRUE)
  y_max <- max(df$log_p_plot, na.rm = TRUE)

  x_pad <- (x_max - x_min) * 0.24
  y_pad <- max(0.35, y_max * 0.18)

  if (!is.finite(x_pad) || x_pad <= 0) x_pad <- 0.5
  if (!is.finite(y_pad) || y_pad <= 0) y_pad <- 0.5

  p1 <- ggplot(df, aes(x = Z, y = log_p_plot)) +
    geom_point(
      data = df %>% filter(significant == "Not significant"),
      color = "grey70",
      size = 2.3,
      alpha = 0.85
    ) +
    geom_point(
      data = df %>% filter(significant == "Significant"),
      color = "#D62728",
      size = 2.9,
      alpha = 0.95
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "#CC0000",
      linewidth = 0.9
    ) +
    labs(
      title = paste0(target_gene, " virtual knockout: ", mode_name),
      subtitle = paste0(
        "cells = ", ncol(pbmc_sub),
        " | genes = ", nrow(data_use),
        " | significant genes = ", nrow(sig_df)
      ),
      x = "Z-score",
      y = ifelse(any(!is.na(df$p.adj)), "-log10(adjusted p)", "-log10(raw p)")
    ) +
    theme_pub(base_size = 18) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 16, r = 70, b = 14, l = 16)
    ) +
    coord_cartesian(
      xlim = c(x_min - 0.2, x_max + x_pad),
      ylim = c(0, y_max + y_pad),
      clip = "off"
    )

  if (nrow(sig_df_plot) > 0) {
    p1 <- p1 +
      ggrepel::geom_text_repel(
        data = sig_df_plot,
        aes(label = gene),
        color = "#CC0000",
        size = 4.2,
        fontface = "bold",
        max.overlaps = Inf,
        box.padding = 0.55,
        point.padding = 0.25,
        segment.color = "grey60",
        segment.size = 0.45,
        min.segment.length = 0,
        force = 10,
        force_pull = 0.15,
        seed = 123
      )
  }

  save_pdf(
    p1,
    path = file.path(mode_outdir, paste0(target_gene, "_", mode_name, "_scatter.pdf")),
    width = 11.5,
    height = 8.8
  )

  ############################
  # top20 barplot
  ############################
  top_genes <- df %>%
    slice_head(n = min(20, nrow(df))) %>%
    mutate(gene = factor(gene, levels = rev(gene)))

  p2 <- ggplot(top_genes, aes(x = gene, y = Z, fill = Z)) +
    geom_col(width = 0.76) +
    coord_flip() +
    scale_fill_gradient2(
      low = "#5B8FD1",
      mid = "white",
      high = "#E10600",
      midpoint = 0
    ) +
    labs(
      title = paste0("Top 20 perturbed genes: ", mode_name),
      x = NULL,
      y = "Z-score"
    ) +
    theme_pub(base_size = 17) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 22)
    )

  save_pdf(
    p2,
    path = file.path(mode_outdir, paste0(target_gene, "_", mode_name, "_barplot.pdf")),
    width = 10.6,
    height = 8.4
  )

  ############################
  # summary
  ############################
  summary_tab <- data.frame(
    mode = mode_name,
    input_rdata = basename(input_rdata),
    seurat_object_used = obj_name,
    target_gene = target_gene,
    restrict_stage = restrict_stage,
    target_stage = ifelse(restrict_stage, target_stage, "ALL_STAGE"),
    selected_cells = ncol(pbmc_sub),
    input_genes = nrow(data_use),
    input_cells = ncol(data_use),
    nfeatures = nfeatures_use,
    nc_nNet = 10,
    nc_nCells = use_nCells,
    sig_gene_number = nrow(sig_df),
    raw_sig_gene_number = nrow(sig_raw_df),
    sig_gene_proportion = round(nrow(sig_df) / nrow(df), 4),
    stringsAsFactors = FALSE
  )

  write.table(
    summary_tab,
    file = file.path(mode_outdir, paste0(target_gene, "_", mode_name, "_run_summary.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  return(summary_tab)
}

############################
# 4. Run all modes
############################
all_summary <- list()

for (i in seq_len(nrow(mode_table))) {
  one_summary <- run_knockout_mode(
    mode_name = mode_table$mode[i],
    nfeatures_use = mode_table$nfeatures[i],
    max_nc_nCells_use = mode_table$max_nc_nCells[i]
  )
  all_summary[[length(all_summary) + 1]] <- one_summary
}

overall_summary <- do.call(rbind, all_summary)

write.table(
  overall_summary,
  file = file.path(output_root, paste0(target_gene, "_three_version_overall_summary.txt")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nAll three knockout modes finished successfully.\n")
cat("Output root:\n")
cat(output_root, "\n")