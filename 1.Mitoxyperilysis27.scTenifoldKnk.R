# -----------------------------
# HMMR virtual knockout
# scTenifoldKnk formal script
# -----------------------------

library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
set.seed(12345)

target_gene <- "HMMR"
setwd("/Users/donglinlu/Desktop/CC1/9.scTenifoldKnk")

# 1. Load object
obj_names <- load("seurat.Rdata")
print(obj_names)

# 2. Basic checks
print(class(pbmc))
print(dim(pbmc))
print(colnames(pbmc@meta.data))

# 3. Select CA_HPV epithelial cells
pbmc_epi <- subset(
  pbmc,
  subset = stage == "CA_HPV" & cellType == "Epithelial cells"
)

cat("Selected CA_HPV epithelial cells:", ncol(pbmc_epi), "\n")
print(table(pbmc_epi@meta.data$stage))
print(table(pbmc_epi@meta.data$cellType))

# 4. Extract RNA counts
DefaultAssay(pbmc_epi) <- "RNA"
countMat <- GetAssayData(pbmc_epi, assay = "RNA", layer = "counts")

# 5. Check target gene
if (!(target_gene %in% rownames(countMat))) {
  stop(paste0("Target gene not found in count matrix: ", target_gene))
}
cat("Target gene found:", target_gene, "\n")

# 6. Select variable features
pbmc_epi <- FindVariableFeatures(
  object = pbmc_epi,
  selection.method = "vst",
  nfeatures = 3000
)

hvgs <- VariableFeatures(pbmc_epi)
feature_use <- unique(c(target_gene, hvgs))
feature_use <- feature_use[feature_use %in% rownames(countMat)]

data <- as.data.frame(countMat[feature_use, ])

cat("Input matrix dimension for scTenifoldKnk:",
    dim(data)[1], "genes x", dim(data)[2], "cells\n")

# 7. Choose nc_nCells adaptively
n_cells <- ncol(data)
use_nCells <- min(300, n_cells)
cat("nc_nCells used:", use_nCells, "\n")

# 8. Run virtual knockout
result <- scTenifoldKnk(
  countMatrix = data,
  gKO = target_gene,
  qc_mtThreshold = 0.15,
  qc_minLSize = 1500,
  nc_nNet = 10,
  nc_nCells = use_nCells
)

save(result, file = "HMMR_scTenifoldKnk_result.Rdata")

# 9. Extract differential regulation results
df <- result$diffRegulation
df <- df[df$gene != target_gene, ]
df <- df %>% arrange(desc(abs(Z)))

write.table(
  df,
  file = "HMMR_all_diffRegulation.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sig_df <- df %>% filter(p.adj < 0.05)

write.table(
  sig_df,
  file = "HMMR_sigDiff.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Significant perturbed genes (p.adj < 0.05):", nrow(sig_df), "\n")

# 10. Barplot: top 20 perturbed genes by |Z|
top_genes <- df %>% slice_head(n = 20)

p1 <- ggplot(top_genes, aes(x = reorder(gene, Z), y = Z, fill = Z)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_gradient2(
    low = "#3C78D8",
    mid = "white",
    high = "#CC0000"
  ) +
  labs(
    title = "Top 20 perturbed genes after HMMR virtual knockout",
    x = NULL,
    y = "Z-score"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.text.y = element_text(face = "bold")
  )

pdf(file = "HMMR_barplot.pdf", width = 7, height = 5.5)
print(p1)
dev.off()

# 11. Volcano plot
df$log_p.adj <- -log10(df$p.adj)
df$significant <- ifelse(df$p.adj < 0.05, "Significant", "Not significant")

label_genes <- df %>%
  filter(p.adj < 0.05) %>%
  slice_head(n = 15)

y_upper <- quantile(df$log_p.adj, 0.995, na.rm = TRUE)
if (!is.finite(y_upper)) y_upper <- max(df$log_p.adj, na.rm = TRUE)

p2 <- ggplot(df, aes(x = Z, y = log_p.adj, color = significant)) +
  geom_point(alpha = 0.75, size = 1.2) +
  scale_color_manual(values = c("Significant" = "#CC0000", "Not significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#CC0000") +
  geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 3.2,
    max.overlaps = 30,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey60"
  ) +
  labs(
    title = "HMMR virtual knockout",
    x = "Z-score",
    y = "-log10(adjusted p)"
  ) +
  coord_cartesian(ylim = c(0, y_upper)) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

pdf(file = "HMMR_volcano.pdf", width = 6.5, height = 5.5)
print(p2)
dev.off()

# 12. Summary file
summary_tab <- data.frame(
  target_gene = target_gene,
  selected_cells = ncol(pbmc_epi),
  input_genes = nrow(data),
  input_cells = ncol(data),
  nc_nNet = 10,
  nc_nCells = use_nCells,
  sig_gene_number = nrow(sig_df)
)

write.table(
  summary_tab,
  file = "HMMR_run_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("HMMR virtual knockout finished successfully.\n")