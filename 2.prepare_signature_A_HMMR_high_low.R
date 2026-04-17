# ============================================================
# Prepare L1000CDS2 Signature A
# HMMR-high vs HMMR-low within CA_HPV prioritized epithelial cells
# prioritized clusters: 6 / 9 / 10
# input: Epithelial_subset.Rdata
# output: gene lists for L1000CDS2 reverse query
# ============================================================

library(Seurat)
library(dplyr)

rm(list = ls())
set.seed(12345)

# -----------------------------
# 0. Paths
# -----------------------------
input_rdata <- "/Users/donglinlu/Desktop/CC1/9.scTenifoldKnk/Epithelial_subset.Rdata"
outdir <- "/Users/donglinlu/Desktop/CC1/12.drugEnrich/L1000CDS2_main_screen"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Parameters
# -----------------------------
target_gene <- "HMMR"
target_stage <- "CA_HPV"
target_clusters <- c("6", "9", "10")
cluster_column_candidates <- c("seurat_clusters", "cluster", "Cluster", "cell_state", "state")
stage_column_candidates <- c("stage", "Stage", "group", "Group")

# use quartile split to generate a stronger transcriptional contrast
split_method <- "quartile"
upper_quantile <- 0.75
lower_quantile <- 0.25

# L1000 export size
n_top <- 100

# -----------------------------
# 2. Load Seurat object
# -----------------------------
obj_names <- load(input_rdata)

seurat_obj_name <- obj_names[sapply(obj_names, function(x) inherits(get(x), "Seurat"))][1]
if (is.na(seurat_obj_name)) {
  stop("No Seurat object found in Epithelial_subset.Rdata")
}

obj <- get(seurat_obj_name)
cat("Loaded Seurat object:", seurat_obj_name, "\n")
cat("Dimensions:", nrow(obj), "genes x", ncol(obj), "cells\n")

meta <- obj@meta.data

# -----------------------------
# 3. Detect metadata columns
# -----------------------------
cluster_col <- cluster_column_candidates[cluster_column_candidates %in% colnames(meta)][1]
stage_col <- stage_column_candidates[stage_column_candidates %in% colnames(meta)][1]

if (is.na(cluster_col)) stop("Cannot find cluster column in metadata")
if (is.na(stage_col)) stop("Cannot find stage column in metadata")
if (!(target_gene %in% rownames(obj))) stop("HMMR not found in object")

cat("Cluster column:", cluster_col, "\n")
cat("Stage column:", stage_col, "\n")

# -----------------------------
# 4. Subset prioritized CA_HPV cells
# -----------------------------
meta[[cluster_col]] <- as.character(meta[[cluster_col]])
meta[[stage_col]] <- as.character(meta[[stage_col]])

keep_cells <- rownames(meta)[
  meta[[stage_col]] == target_stage &
    meta[[cluster_col]] %in% target_clusters
]

obj_sub <- subset(obj, cells = keep_cells)

cat("Selected cells:", ncol(obj_sub), "\n")
print(table(obj_sub@meta.data[[stage_col]]))
print(table(as.character(obj_sub@meta.data[[cluster_col]])))

if (ncol(obj_sub) < 50) {
  stop("Too few cells after subsetting. Please re-check stage / cluster columns.")
}

# -----------------------------
# 5. Build HMMR-high and HMMR-low groups
# -----------------------------
DefaultAssay(obj_sub) <- "RNA"
expr_vec <- FetchData(obj_sub, vars = target_gene)[, 1]

if (split_method == "quartile") {
  q_low <- as.numeric(quantile(expr_vec, lower_quantile, na.rm = TRUE))
  q_high <- as.numeric(quantile(expr_vec, upper_quantile, na.rm = TRUE))

  hmmr_group <- ifelse(expr_vec <= q_low, "HMMR_low",
                       ifelse(expr_vec >= q_high, "HMMR_high", NA))
} else {
  med <- median(expr_vec, na.rm = TRUE)
  hmmr_group <- ifelse(expr_vec > med, "HMMR_high", "HMMR_low")
}

obj_sub$HMMR_group <- hmmr_group
obj_sigA <- subset(obj_sub, subset = !is.na(HMMR_group))

cat("Cells used for Signature A:\n")
print(table(obj_sigA$HMMR_group))

# -----------------------------
# 6. Differential expression
# -----------------------------
Idents(obj_sigA) <- "HMMR_group"

markers <- FindMarkers(
  object = obj_sigA,
  ident.1 = "HMMR_high",
  ident.2 = "HMMR_low",
  assay = "RNA",
  slot = "data",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

markers$gene <- rownames(markers)
markers <- markers %>%
  arrange(desc(avg_log2FC))

write.csv(
  markers,
  file = file.path(outdir, "SignatureA_HMMR_high_vs_low_all_DEGs.csv"),
  row.names = FALSE
)

# -----------------------------
# 7. Export L1000 gene lists
# -----------------------------
up_genes <- markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = n_top) %>%
  pull(gene)

# genes more highly expressed in HMMR_low (i.e. down in HMMR_high)
down_genes <- markers %>%
  filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
  arrange(avg_log2FC) %>%
  slice_head(n = n_top) %>%
  pull(gene)

write.table(
  up_genes,
  file = file.path(outdir, "SignatureA_HMMR_high_genes_top100.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  down_genes,
  file = file.path(outdir, "SignatureA_HMMR_low_genes_top100.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

max_len <- max(length(up_genes), length(down_genes))
submit_tab <- data.frame(
  up_genes = c(up_genes, rep(NA, max_len - length(up_genes))),
  down_genes = c(down_genes, rep(NA, max_len - length(down_genes)))
)

write.csv(
  submit_tab,
  file = file.path(outdir, "SignatureA_L1000_submission_table.csv"),
  row.names = FALSE
)

# -----------------------------
# 8. Summary
# -----------------------------
summary_tab <- data.frame(
  signature = "Signature_A",
  input_rdata = basename(input_rdata),
  seurat_object_used = seurat_obj_name,
  target_gene = target_gene,
  stage = target_stage,
  clusters = paste(target_clusters, collapse = ","),
  split_method = split_method,
  q_low = ifelse(split_method == "quartile", q_low, NA),
  q_high = ifelse(split_method == "quartile", q_high, NA),
  selected_cells_before_split = ncol(obj_sub),
  selected_cells_after_split = ncol(obj_sigA),
  n_HMMR_high = sum(obj_sigA$HMMR_group == "HMMR_high"),
  n_HMMR_low = sum(obj_sigA$HMMR_group == "HMMR_low"),
  n_sig_up = length(up_genes),
  n_sig_down = length(down_genes)
)

write.table(
  summary_tab,
  file = file.path(outdir, "SignatureA_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nSignature A preparation completed successfully.\n")
cat("Output directory:", outdir, "\n")
cat("Top up genes:", length(up_genes), "\n")
cat("Top down genes:", length(down_genes), "\n")
