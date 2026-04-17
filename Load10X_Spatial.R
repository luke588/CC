options(stringsAsFactors = FALSE)

# =========================================================
# Final polished redraw for 10x cervical FFPE spatial module
# project: 15.main_spatial_support_10x_cervical_FFPE
# purpose:
# 1) formal cluster annotation mapping
# 2) redraw all main and supplementary plots in unified pubstyle
# 3) keep broad annotation for main text and annotated clusters for supplement
# final refined version: no subtitles + fixed title display
# =========================================================

# ----------------------------- packages -----------------------------
pkg_needed <- c(
  "Seurat", "SeuratObject", "Matrix", "ggplot2", "dplyr", "readr",
  "tibble", "tidyr", "patchwork", "stringr", "hdf5r"
)

pkg_to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if (length(pkg_to_install) > 0) {
  install.packages(pkg_to_install, repos = "https://cloud.r-project.org")
}
invisible(lapply(pkg_needed, library, character.only = TRUE))

set.seed(1234)

# ----------------------------- paths -----------------------------
root_dir <- "/Users/donglinlu/Desktop/CC1/15.main_spatial_support_10x_cervical_FFPE"

obj_candidates <- c(
  file.path(root_dir, "03.results", "04.redraw_pubstyle", "04.objects", "10x_cervical_FFPE_main_spatial_seurat_redraw_pubstyle.rds"),
  file.path(root_dir, "03.results", "01.objects", "10x_cervical_FFPE_main_spatial_seurat.rds")
)

obj_file <- obj_candidates[file.exists(obj_candidates)][1]
if (is.na(obj_file)) {
  stop("Cannot find Seurat object in expected locations.")
}

redraw_dir <- file.path(root_dir, "03.results", "05.final_pubstyle")
dir.create(redraw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(redraw_dir, "01.tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(redraw_dir, "02.plots_main"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(redraw_dir, "03.plots_supp"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(redraw_dir, "04.objects"), showWarnings = FALSE, recursive = TRUE)

# ----------------------------- fixed signature files -----------------------------
malignant_file <- "/Users/donglinlu/Desktop/CC1/4.epithelial subset/step1_5.cluster_9.markers_vs_rest.csv"
emt_file       <- "/Users/donglinlu/Desktop/CC1/4.epithelial subset/step1_5.cluster_10.markers_vs_rest.csv"
prolif_file    <- "/Users/donglinlu/Desktop/CC1/4.epithelial subset/step1_5.cluster_6.markers_vs_rest.csv"

if (!file.exists(malignant_file)) stop("Missing file: ", malignant_file)
if (!file.exists(emt_file)) stop("Missing file: ", emt_file)
if (!file.exists(prolif_file)) stop("Missing file: ", prolif_file)

# ----------------------------- helpers -----------------------------
theme_pub <- function(base_size = 17) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = "black",
        size = base_size + 3,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_blank(),
      axis.title = element_text(
        face = "bold",
        color = "black",
        size = base_size + 1
      ),
      axis.text = element_text(
        color = "black",
        size = base_size - 1
      ),
      axis.text.x = element_text(
        angle = 35,
        hjust = 1,
        size = base_size - 1,
        face = "bold"
      ),
      axis.text.y = element_text(
        size = base_size - 1,
        face = "bold"
      ),
      legend.title = element_text(
        face = "bold",
        color = "black",
        size = base_size
      ),
      legend.text = element_text(
        color = "black",
        size = base_size - 1
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        color = "black",
        size = base_size
      ),
      panel.border = element_blank(),
      plot.margin = margin(14, 18, 14, 18)
    )
}

theme_spatial_pub <- function(base_size = 17) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = "black",
        size = base_size + 3,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_blank(),
      legend.title = element_text(
        face = "bold",
        color = "black",
        size = base_size
      ),
      legend.text = element_text(
        color = "black",
        size = base_size - 1
      ),
      plot.margin = margin(18, 24, 16, 24)
    )
}

save_plot_pdf_png <- function(plot_obj, out_prefix, width = 8, height = 6, dpi = 320) {
  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE,
    useDingbats = FALSE
  )
  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    limitsize = FALSE,
    bg = "white"
  )
}

read_top_signature <- function(csv_file, top_n = 50) {
  df <- read.csv(csv_file, check.names = FALSE, stringsAsFactors = FALSE)

  gene_col_candidates <- c("gene", "Gene", "GENE", "symbol", "Symbol", "SYMBOL", "gene_symbol")
  fc_col_candidates   <- c("avg_log2FC", "avg_logFC", "avg_log2fc", "log2FC", "logFC", "avg_logfc")

  gene_col <- gene_col_candidates[gene_col_candidates %in% colnames(df)][1]
  fc_col   <- fc_col_candidates[fc_col_candidates %in% colnames(df)][1]

  if (is.na(gene_col)) {
    if (ncol(df) >= 1 && is.character(df[[1]])) {
      gene_col <- colnames(df)[1]
    } else {
      stop("Cannot identify gene column in:\n", csv_file)
    }
  }

  if (!is.na(fc_col)) {
    df <- df %>%
      dplyr::filter(!is.na(.data[[fc_col]])) %>%
      dplyr::filter(.data[[fc_col]] > 0) %>%
      dplyr::arrange(dplyr::desc(.data[[fc_col]]))
  }

  genes <- df[[gene_col]]
  genes <- genes[!is.na(genes)]
  genes <- genes[genes != ""]
  genes <- unique(genes)
  genes <- genes[seq_len(min(top_n, length(genes)))]

  return(genes)
}

get_expr_mat <- function(seu_obj, assay_name = "SCT") {
  assay_obj <- seu_obj[[assay_name]]

  if (inherits(assay_obj, "Assay5")) {
    layer_names <- SeuratObject::Layers(assay_obj)

    if ("scale.data" %in% layer_names) {
      expr_mat <- SeuratObject::LayerData(seu_obj, assay = assay_name, layer = "scale.data")
      used_layer <- "scale.data"
    } else if ("data" %in% layer_names) {
      expr_mat <- SeuratObject::LayerData(seu_obj, assay = assay_name, layer = "data")
      used_layer <- "data"
    } else if ("counts" %in% layer_names) {
      expr_mat <- SeuratObject::LayerData(seu_obj, assay = assay_name, layer = "counts")
      used_layer <- "counts"
    } else {
      stop("No usable layer found in assay: ", assay_name)
    }
  } else {
    if ("scale.data" %in% slotNames(assay_obj) && nrow(assay_obj@scale.data) > 0) {
      expr_mat <- assay_obj@scale.data
      used_layer <- "scale.data"
    } else if ("data" %in% slotNames(assay_obj) && nrow(assay_obj@data) > 0) {
      expr_mat <- assay_obj@data
      used_layer <- "data"
    } else if ("counts" %in% slotNames(assay_obj) && nrow(assay_obj@counts) > 0) {
      expr_mat <- assay_obj@counts
      used_layer <- "counts"
    } else {
      stop("No usable expression matrix found in assay: ", assay_name)
    }
  }

  return(list(mat = expr_mat, layer = used_layer))
}

compute_signature_score <- function(seu_obj, genes, score_name, assay_name = "SCT") {
  mat_info <- get_expr_mat(seu_obj, assay_name = assay_name)
  expr_mat <- mat_info$mat
  used_layer <- mat_info$layer

  use_genes <- intersect(genes, rownames(expr_mat))

  if (length(use_genes) < 5) {
    warning(score_name, ": fewer than 5 genes found in ", used_layer)
    seu_obj[[score_name]] <- NA_real_
    return(list(obj = seu_obj, used_genes = use_genes, used_layer = used_layer))
  }

  score_vec <- Matrix::colMeans(expr_mat[use_genes, , drop = FALSE])
  seu_obj[[score_name]] <- as.numeric(score_vec)

  return(list(obj = seu_obj, used_genes = use_genes, used_layer = used_layer))
}

# ----------------------------- load object -----------------------------
obj <- readRDS(obj_file)

if ("SCT" %in% Assays(obj)) {
  DefaultAssay(obj) <- "SCT"
} else {
  DefaultAssay(obj) <- DefaultAssay(obj)
}

# ----------------------------- formal mapping -----------------------------
cluster_map_df <- tibble::tribble(
  ~cluster, ~cluster_short, ~cluster_annotation, ~broad_annotation, ~annotation_confidence, ~annotation_basis,
  "0", "FS-interface",   "Fibro-stromal interface-like",   "Fibro-stromal",    "moderate", "mixed stromal / epithelial interface pattern",
  "1", "FS-matrix",      "Matrix fibro-stromal-like",      "Fibro-stromal",    "moderate", "COL1A1/COL1A2/DCN/LUM-rich stromal compartment",
  "2", "Basal-squamous", "Basal squamous epithelial-like", "Tumor epithelial", "high",     "KRT5/KRT14/TP63-enriched basal squamous program",
  "3", "Immune-rich",    "Immune-rich like",               "Immune",           "high",     "PTPRC/CD3D/CD3E/NKG7/LYZ immune marker enrichment",
  "4", "Vasc-peri",      "Vascular/perivascular-like",     "Fibro-stromal",    "moderate", "vascular/perivascular mixed marker pattern",
  "5", "Myofib-vasc",    "Myofibrovascular stromal-like",  "Fibro-stromal",    "moderate", "stromal + vascular + ACTA2/TAGLN-rich state",
  "6", "Sq-epi",         "Squamous tumor epithelial-like", "Tumor epithelial", "moderate", "squamous epithelial-related cluster",
  "7", "Sq-interface",   "Squamous-immune interface-like", "Tumor epithelial", "moderate", "interface state between squamous and immune compartments"
)

readr::write_csv(
  cluster_map_df,
  file.path(redraw_dir, "01.tables", "10x_cluster_number_to_annotation_formal_mapping.csv")
)

cluster_short_map  <- setNames(cluster_map_df$cluster_short, cluster_map_df$cluster)
cluster_detail_map <- setNames(cluster_map_df$cluster_annotation, cluster_map_df$cluster)
cluster_broad_map  <- setNames(cluster_map_df$broad_annotation, cluster_map_df$cluster)

obj$cluster_short <- unname(cluster_short_map[as.character(obj$seurat_clusters)])
obj$cluster_annotation <- unname(cluster_detail_map[as.character(obj$seurat_clusters)])
obj$broad_annotation_main <- unname(cluster_broad_map[as.character(obj$seurat_clusters)])

obj$cluster_short <- factor(obj$cluster_short, levels = cluster_map_df$cluster_short)
obj$cluster_annotation <- factor(obj$cluster_annotation, levels = cluster_map_df$cluster_annotation)
obj$broad_annotation_main <- factor(
  obj$broad_annotation_main,
  levels = c("Tumor epithelial", "Immune", "Fibro-stromal")
)

obj$cluster_display <- factor(
  paste0("C", as.character(obj$seurat_clusters), ": ", obj$cluster_short),
  levels = paste0("C", cluster_map_df$cluster, ": ", cluster_map_df$cluster_short)
)

# ----------------------------- colors -----------------------------
broad_colors <- c(
  "Tumor epithelial" = "#5B7BD5",
  "Immune" = "#5DAE61",
  "Fibro-stromal" = "#D97070"
)

cluster_colors <- c(
  "FS-interface" = "#C96E6E",
  "FS-matrix" = "#E48F8F",
  "Basal-squamous" = "#5877CE",
  "Immune-rich" = "#57A85D",
  "Vasc-peri" = "#53A5A9",
  "Myofib-vasc" = "#7AA59E",
  "Sq-epi" = "#8797E8",
  "Sq-interface" = "#B47EE5"
)

# ----------------------------- signatures -----------------------------
malignant_genes <- read_top_signature(malignant_file, top_n = 50)
emt_genes       <- read_top_signature(emt_file, top_n = 50)
prolif_genes    <- read_top_signature(prolif_file, top_n = 50)

if (!"malignant_squamous_score" %in% colnames(obj@meta.data)) {
  tmp1 <- compute_signature_score(obj, malignant_genes, "malignant_squamous_score", assay_name = DefaultAssay(obj))
  obj <- tmp1$obj
}
if (!"EMT_like_score" %in% colnames(obj@meta.data)) {
  tmp2 <- compute_signature_score(obj, emt_genes, "EMT_like_score", assay_name = DefaultAssay(obj))
  obj <- tmp2$obj
}
if (!"proliferation_core_score" %in% colnames(obj@meta.data)) {
  tmp3 <- compute_signature_score(obj, prolif_genes, "proliferation_core_score", assay_name = DefaultAssay(obj))
  obj <- tmp3$obj
}

# ----------------------------- MAIN FIGURE 1 -----------------------------
p_main_broad <- SpatialDimPlot(
  obj,
  group.by = "broad_annotation_main",
  cols = broad_colors,
  label = TRUE,
  label.size = 5.2,
  pt.size.factor = 1.8,
  image.alpha = 1
) &
  theme_spatial_pub(base_size = 17)

p_main_broad <- p_main_broad +
  labs(
    title = "Major spatial compartments"
  )

save_plot_pdf_png(
  p_main_broad,
  file.path(redraw_dir, "02.plots_main", "main_spatial_broad_annotation_pubstyle"),
  width = 11.2,
  height = 8.4
)

# ----------------------------- MAIN FIGURE 2 -----------------------------
DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else DefaultAssay(obj)

marker_panel_genes <- c(
  "KRT5", "KRT14", "TP63",
  "COL1A1", "DCN", "LUM",
  "PTPRC", "HLA-DRA", "PECAM1"
)
marker_panel_genes <- intersect(marker_panel_genes, rownames(obj))

if (length(marker_panel_genes) >= 6) {
  p_main_marker <- SpatialFeaturePlot(
    obj,
    features = marker_panel_genes,
    ncol = 3,
    pt.size.factor = 1.8,
    image.alpha = 1,
    min.cutoff = rep("q05", length(marker_panel_genes)),
    max.cutoff = rep("q95", length(marker_panel_genes))
  ) &
    theme_spatial_pub(base_size = 15)

  p_main_marker <- p_main_marker &
    theme(
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold")
    )

  save_plot_pdf_png(
    p_main_marker,
    file.path(redraw_dir, "02.plots_main", "main_selected_marker_spatial_panels_pubstyle"),
    width = 17.2,
    height = 12.2
  )
}

# ----------------------------- MAIN FIGURE 3 -----------------------------
sig_features <- c("malignant_squamous_score", "EMT_like_score", "proliferation_core_score")

p_main_sig_spatial <- SpatialFeaturePlot(
  obj,
  features = sig_features,
  ncol = 3,
  pt.size.factor = 1.8,
  image.alpha = 1,
  min.cutoff = c("q05", "q05", "q05"),
  max.cutoff = c("q95", "q95", "q95")
) &
  theme_spatial_pub(base_size = 15)

save_plot_pdf_png(
  p_main_sig_spatial,
  file.path(redraw_dir, "02.plots_main", "main_signature_score_spatial_panels_pubstyle"),
  width = 18.6,
  height = 6.9
)

# ----------------------------- MAIN FIGURE 4 -----------------------------
score_long <- obj@meta.data %>%
  tibble::rownames_to_column("spot_id") %>%
  dplyr::select(
    broad_annotation_main,
    malignant_squamous_score,
    EMT_like_score,
    proliferation_core_score
  ) %>%
  tidyr::pivot_longer(
    cols = c(malignant_squamous_score, EMT_like_score, proliferation_core_score),
    names_to = "signature",
    values_to = "score"
  ) %>%
  dplyr::mutate(
    broad_annotation_main = factor(
      broad_annotation_main,
      levels = c("Tumor epithelial", "Immune", "Fibro-stromal")
    ),
    signature = factor(
      signature,
      levels = c("malignant_squamous_score", "EMT_like_score", "proliferation_core_score"),
      labels = c("Malignant squamous", "EMT-like", "Proliferation core")
    )
  )

score_summary <- score_long %>%
  dplyr::group_by(broad_annotation_main, signature) %>%
  dplyr::summarise(
    n_spots = dplyr::n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  score_summary,
  file.path(redraw_dir, "01.tables", "10x_signature_score_by_broad_annotation_pubstyle.csv")
)

p_main_score <- ggplot(
  score_long,
  aes(x = broad_annotation_main, y = score, fill = broad_annotation_main)
) +
  geom_violin(scale = "width", color = "black", linewidth = 0.3, trim = TRUE) +
  geom_boxplot(width = 0.14, outlier.size = 0.2, fill = "white", color = "black") +
  facet_wrap(~signature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = broad_colors) +
  labs(
    title = "Program scores across spatial identities",
    x = NULL,
    y = "Average scaled expression score"
  ) +
  theme_pub(base_size = 17) +
  theme(legend.position = "none")

save_plot_pdf_png(
  p_main_score,
  file.path(redraw_dir, "02.plots_main", "main_signature_score_by_broad_identity_pubstyle"),
  width = 16.4,
  height = 6.3
)

# ----------------------------- SUPPLEMENT 1 -----------------------------
p_supp_cluster_spatial <- SpatialDimPlot(
  obj,
  group.by = "cluster_short",
  cols = cluster_colors,
  label = TRUE,
  label.size = 4.6,
  pt.size.factor = 1.8,
  image.alpha = 1
) &
  theme_spatial_pub(base_size = 17)

p_supp_cluster_spatial <- p_supp_cluster_spatial +
  labs(
    title = "Annotated spot clusters"
  )

save_plot_pdf_png(
  p_supp_cluster_spatial,
  file.path(redraw_dir, "03.plots_supp", "supp_spatial_cluster_annotation_pubstyle"),
  width = 11.0,
  height = 8.6
)

# ----------------------------- SUPPLEMENT 2 -----------------------------
p_supp_cluster_umap <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "cluster_short",
  cols = cluster_colors,
  label = TRUE,
  repel = TRUE,
  pt.size = 1.35
) +
  labs(
    title = "Annotated Visium spot clusters"
  ) +
  theme_pub(base_size = 17)

save_plot_pdf_png(
  p_supp_cluster_umap,
  file.path(redraw_dir, "03.plots_supp", "supp_UMAP_cluster_annotation_pubstyle"),
  width = 10.2,
  height = 7.7
)

# ----------------------------- SUPPLEMENT 3 -----------------------------
p_supp_broad_umap <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "broad_annotation_main",
  cols = broad_colors,
  label = TRUE,
  repel = TRUE,
  pt.size = 1.35
) +
  labs(
    title = "Broad spatial identities"
  ) +
  theme_pub(base_size = 17)

save_plot_pdf_png(
  p_supp_broad_umap,
  file.path(redraw_dir, "03.plots_supp", "supp_UMAP_broad_annotation_pubstyle"),
  width = 10.0,
  height = 7.5
)

# ----------------------------- SUPPLEMENT 4 -----------------------------
dotplot_features <- c(
  "EPCAM", "KRT5", "KRT14", "TP63", "KRT8", "KRT18",
  "COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "COL3A1",
  "PECAM1", "VWF", "KDR", "EMCN", "CLDN5",
  "PTPRC", "CD3D", "CD3E", "NKG7", "LYZ", "MS4A1", "HLA-DRA",
  "ACTA2", "TAGLN", "RGS5", "MCAM", "CSPG4"
)
dotplot_features <- intersect(dotplot_features, rownames(obj))

Idents(obj) <- obj$cluster_display

p_supp_dot <- DotPlot(
  obj,
  features = dotplot_features,
  assay = DefaultAssay(obj),
  dot.scale = 9
) +
  RotatedAxis() +
  labs(
    title = "Marker support for annotated spot clusters",
    x = NULL,
    y = NULL
  ) +
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      override.aes = list(alpha = 1)
    ),
    colour = guide_colorbar(
      title = "Average Expression",
      barheight = grid::unit(6.5, "cm"),
      barwidth = grid::unit(0.8, "cm")
    )
  ) +
  theme_pub(base_size = 18) +
  theme(
    plot.title = element_text(size = 21, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 14, angle = 35, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 13),
    legend.position = "right"
  )

save_plot_pdf_png(
  p_supp_dot,
  file.path(redraw_dir, "03.plots_supp", "supp_cluster_marker_dotplot_annotated_pubstyle"),
  width = 18.5,
  height = 10.8
)

# ----------------------------- SUPPLEMENT 5 -----------------------------
Idents(obj) <- obj$cluster_short

p_supp_sig_umap <- FeaturePlot(
  obj,
  features = sig_features,
  reduction = "umap",
  ncol = 3,
  pt.size = 1.35,
  min.cutoff = c("q05", "q05", "q05"),
  max.cutoff = c("q95", "q95", "q95")
) &
  theme_pub(base_size = 15)

save_plot_pdf_png(
  p_supp_sig_umap,
  file.path(redraw_dir, "03.plots_supp", "supp_signature_score_UMAP_panels_pubstyle"),
  width = 16.8,
  height = 6.4
)

# ----------------------------- SUPPLEMENT 6 -----------------------------
obj$sample_label <- "10x cervical FFPE SCC"

p_supp_qc_violin <- VlnPlot(
  obj,
  features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
  group.by = "sample_label",
  pt.size = 0,
  ncol = 3
) & theme_pub(base_size = 17)

save_plot_pdf_png(
  p_supp_qc_violin,
  file.path(redraw_dir, "03.plots_supp", "supp_QC_violin_pubstyle"),
  width = 14.8,
  height = 5.8
)

p_supp_qc_spatial <- SpatialFeaturePlot(
  obj,
  features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
  ncol = 3,
  pt.size.factor = 1.8,
  image.alpha = 1,
  min.cutoff = c("q05", "q05", "q05"),
  max.cutoff = c("q95", "q95", "q95")
) &
  theme_spatial_pub(base_size = 15)

save_plot_pdf_png(
  p_supp_qc_spatial,
  file.path(redraw_dir, "03.plots_supp", "supp_QC_spatial_panels_pubstyle"),
  width = 17.2,
  height = 6.2
)

# ----------------------------- save updated object -----------------------------
saveRDS(
  obj,
  file.path(redraw_dir, "04.objects", "10x_cervical_FFPE_main_spatial_seurat_final_pubstyle.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(redraw_dir, "01.tables", "sessionInfo_final_pubstyle.txt")
)

message("Final pubstyle redraw finished.")
message("Results saved to: ", redraw_dir)