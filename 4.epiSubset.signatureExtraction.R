rm(list = ls())
options(stringsAsFactors = FALSE)

source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")

#######################00. libraries#######################
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(grid)
library(pheatmap)

#######################01. helper functions#######################
safe_scale_rows <- function(mat) {
  if (is.null(mat) || length(mat) == 0) return(NULL)
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (nrow(mat) == 0 || ncol(mat) == 0) return(NULL)

  keep_row <- apply(mat, 1, function(x) any(is.finite(x)))
  mat <- mat[keep_row, , drop = FALSE]
  if (nrow(mat) == 0 || ncol(mat) == 0) return(NULL)

  mat.scale <- t(scale(t(mat)))
  mat.scale[is.na(mat.scale)] <- 0
  mat.scale[is.infinite(mat.scale)] <- 0

  if (nrow(mat.scale) == 0 || ncol(mat.scale) == 0) return(NULL)
  if (all(mat.scale == 0)) return(NULL)

  mat.scale
}

#######################02. global settings#######################
logFCfilter <- 0.5
adjPvalFilter <- 0.05

#设置工作目录
workDir <- "/Users/donglinlu/Desktop/CC1/4.epithelial subset"
setwd(workDir)

#读取对象
load("02.Epithelial_subset.manualAnnotation.Rdata")

#检查对象
print(table(epi.clean$manualCellType))

#确保identity是clean状态
epi.clean$manualCellType <- as.character(epi.clean$manualCellType)

#状态顺序统一
state.order <- c(
  "EMT_like_epithelial",
  "secretory_like_epithelial",
  "keratinizing_squamous_epithelial",
  "aberrant_secretory_epithelial",
  "inflammatory_secretory_epithelial",
  "proliferative_epithelial",
  "basal_stromal_like_epithelial",
  "superficial_keratinizing_epithelial",
  "terminally_differentiated_squamous_epithelial",
  "basal_EMT_like_epithelial",
  "rare_epithelial_like_cells"
)

state.order <- state.order[state.order %in% unique(epi.clean$manualCellType)]
epi.clean$manualCellType <- factor(epi.clean$manualCellType, levels = state.order)
Idents(epi.clean) <- epi.clean$manualCellType

#当前重点候选状态
candidate.states <- c(
  "proliferative_epithelial",
  "EMT_like_epithelial",
  "inflammatory_secretory_epithelial"
)

#######################03.分别提取每个候选状态的marker#######################
all.marker.list <- list()
sig.marker.list <- list()

for (state in candidate.states) {
  cat("Running marker extraction for:", state, "\n")

  state.markers <- FindMarkers(
    object = epi.clean,
    ident.1 = state,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

  state.markers$gene <- rownames(state.markers)
  state.markers <- state.markers %>% arrange(desc(avg_log2FC))

  # 输出全部marker
  outFile1 <- paste0("05.", state, ".markers.all.txt")
  write.table(
    state.markers,
    file = outFile1,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # 筛显著marker
  sig.markers <- state.markers[
    (as.numeric(as.vector(state.markers$avg_log2FC)) > logFCfilter &
       as.numeric(as.vector(state.markers$p_val_adj)) < adjPvalFilter),
  ]

  outFile2 <- paste0("05.", state, ".markers.sig.txt")
  write.table(
    sig.markers,
    file = outFile2,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  all.marker.list[[state]] <- state.markers
  sig.marker.list[[state]] <- sig.markers
}

#######################04.构建signature gene list#######################
# 1. 各状态显著marker基因
prolif.genes <- sig.marker.list[["proliferative_epithelial"]]$gene
emt.genes    <- sig.marker.list[["EMT_like_epithelial"]]$gene
inflam.genes <- sig.marker.list[["inflammatory_secretory_epithelial"]]$gene

# 去掉NA
prolif.genes <- prolif.genes[!is.na(prolif.genes)]
emt.genes <- emt.genes[!is.na(emt.genes)]
inflam.genes <- inflam.genes[!is.na(inflam.genes)]

# 2. 联合集合：作为high-risk候选signature的大池子
highRisk.union <- unique(c(prolif.genes, emt.genes, inflam.genes))
write.table(
  data.frame(gene = highRisk.union),
  file = "05.highRisk_signature_union.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# 3. 交集集合：如果存在，代表多个危险状态共享程序
highRisk.intersect <- Reduce(intersect, list(prolif.genes, emt.genes, inflam.genes))
write.table(
  data.frame(gene = highRisk.intersect),
  file = "05.highRisk_signature_intersect.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# 4. 更适合后续bulk验证的topN signature
topN <- 20

prolif.top <- head(sig.marker.list[["proliferative_epithelial"]]$gene, topN)
emt.top <- head(sig.marker.list[["EMT_like_epithelial"]]$gene, topN)
inflam.top <- head(sig.marker.list[["inflammatory_secretory_epithelial"]]$gene, topN)

highRisk.topN <- unique(c(prolif.top, emt.top, inflam.top))
write.table(
  data.frame(gene = highRisk.topN),
  file = "05.highRisk_signature_topN.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

#######################05.高风险状态top gene可视化#######################
# 每个状态取前10个基因用于展示
top10.prolif <- head(all.marker.list[["proliferative_epithelial"]]$gene, 10)
top10.emt <- head(all.marker.list[["EMT_like_epithelial"]]$gene, 10)
top10.inflam <- head(all.marker.list[["inflammatory_secretory_epithelial"]]$gene, 10)

plot.genes <- unique(c(top10.prolif, top10.emt, top10.inflam))
plot.genes <- plot.genes[!is.na(plot.genes)]
plot.genes <- plot.genes[plot.genes %in% rownames(epi.clean)]

# 1. DotPlot（推荐主图之一）
if (length(plot.genes) > 0) {
  p_dot <- DotPlot(
    epi.clean,
    features = plot.genes,
    group.by = "manualCellType"
  ) +
    RotatedAxis() +
    labs(
      title = "Top marker expression across epithelial states",
      x = NULL,
      y = NULL
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = PUB_AXIS_TEXT_SIZE - 1,
        colour = "black"
      ),
      axis.text.y = element_text(
        size = PUB_AXIS_TEXT_SIZE - 1,
        colour = "black"
      ),
      legend.position = "right"
    )

  save_pub(
    plot_obj = p_dot,
    pdf_file = "05.highRisk.topGene.dotplot.vertical.pdf",
    png_file = "05.highRisk.topGene.dotplot.vertical.png",
    width = 13,
    height = 8.8,
    dpi = 600
  )
} else {
  cat("Skip DotPlot: no valid plot.genes\n")
}

# 2. 平均表达热图（推荐主图之一）
if (length(plot.genes) > 0) {
  avg.exp <- AverageExpression(
    epi.clean,
    features = plot.genes,
    group.by = "manualCellType",
    assays = "RNA",
    slot = "data"
  )

  if ("RNA" %in% names(avg.exp)) {
    heat.mat <- avg.exp$RNA
    state.order.heat <- state.order[state.order %in% colnames(heat.mat)]
    heat.mat <- heat.mat[, state.order.heat, drop = FALSE]

    heat.mat.scale <- safe_scale_rows(heat.mat)

    if (!is.null(heat.mat.scale)) {
      pdf(file = "05.highRisk.marker.avgExp.heatmap.pdf", width = 11.5, height = 8.8)
      pheatmap(
        heat.mat.scale,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        angle_col = 45,
        border_color = NA,
        fontsize_row = 11,
        fontsize_col = 11,
        main = "Average marker expression across epithelial states"
      )
      dev.off()
    } else {
      cat("Skip main heatmap: matrix is empty / non-finite / all zero after scaling\n")
    }
  } else {
    cat("Skip main heatmap: AverageExpression returned no RNA matrix\n")
  }
} else {
  cat("Skip main heatmap: no valid plot.genes\n")
}

# 3. 每个候选状态分别做平均表达热图（更清楚，带保护）
state.list <- c(
  "proliferative_epithelial",
  "EMT_like_epithelial",
  "inflammatory_secretory_epithelial"
)

for (state in state.list) {
  cat("Processing:", state, "\n")

  if (!(state %in% names(all.marker.list))) {
    cat("Skip:", state, "not found in all.marker.list\n")
    next
  }

  top10.state <- all.marker.list[[state]]$gene
  top10.state <- top10.state[!is.na(top10.state)]
  top10.state <- top10.state[top10.state != ""]
  top10.state <- unique(top10.state)
  top10.state <- head(top10.state, 10)
  top10.state <- top10.state[top10.state %in% rownames(epi.clean)]

  if (length(top10.state) == 0) {
    cat("Skip:", state, "has no valid genes\n")
    next
  }

  avg.state <- AverageExpression(
    epi.clean,
    features = top10.state,
    group.by = "manualCellType",
    assays = "RNA",
    slot = "data"
  )

  if (!"RNA" %in% names(avg.state)) {
    cat("Skip:", state, "AverageExpression result has no RNA matrix\n")
    next
  }

  mat.state <- avg.state$RNA
  valid.cols <- intersect(state.order, colnames(mat.state))
  mat.state <- mat.state[, valid.cols, drop = FALSE]

  if (nrow(mat.state) == 0 || ncol(mat.state) == 0) {
    cat("Skip:", state, "matrix is empty\n")
    next
  }

  mat.state.scale <- safe_scale_rows(mat.state)

  if (is.null(mat.state.scale)) {
    cat("Skip:", state, "scaled matrix is empty / non-finite / all zero\n")
    next
  }

  outFile <- paste0("05.", state, ".avgExp.heatmap.pdf")
  pdf(file = outFile, width = 11, height = 5.8)
  pheatmap(
    mat.state.scale,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    angle_col = 45,
    border_color = NA,
    fontsize_row = 11,
    fontsize_col = 11,
    main = paste0("Top markers of ", state)
  )
  dev.off()

  cat("Done:", state, "\n")
}

#######################06.候选状态在降维图上的可视化#######################
# 构建简单的group标签：主候选/其他
epi.clean$highRiskGroup <- "Other_states"
epi.clean$highRiskGroup[epi.clean$manualCellType == "proliferative_epithelial"] <- "proliferative_epithelial"
epi.clean$highRiskGroup[epi.clean$manualCellType == "EMT_like_epithelial"] <- "EMT_like_epithelial"
epi.clean$highRiskGroup[epi.clean$manualCellType == "inflammatory_secretory_epithelial"] <- "inflammatory_secretory_epithelial"

epi.clean$highRiskGroup <- factor(
  epi.clean$highRiskGroup,
  levels = c(
    "Other_states",
    "proliferative_epithelial",
    "EMT_like_epithelial",
    "inflammatory_secretory_epithelial"
  )
)

highRisk.cols <- c(
  "Other_states" = "grey80",
  "proliferative_epithelial" = "#1f78b4",
  "EMT_like_epithelial" = "#e31a1c",
  "inflammatory_secretory_epithelial" = "#33a02c"
)

# UMAP
p_umap_hr <- DimPlot(
  epi.clean,
  reduction = "umap",
  group.by = "highRiskGroup",
  pt.size = 1.0,
  cols = highRisk.cols
) +
  labs(title = "Candidate high-risk states on UMAP") +
  theme_pub()

save_pub(
  plot_obj = p_umap_hr,
  pdf_file = "05.highRisk.state.featurePlot.umap.pdf",
  png_file = "05.highRisk.state.featurePlot.umap.png",
  width = 11,
  height = 8.5,
  dpi = 600
)

# tSNE
p_tsne_hr <- TSNEPlot(
  epi.clean,
  group.by = "highRiskGroup",
  pt.size = 1.0,
  cols = highRisk.cols
) +
  labs(title = "Candidate high-risk states on tSNE") +
  theme_pub()

save_pub(
  plot_obj = p_tsne_hr,
  pdf_file = "05.highRisk.state.featurePlot.tsne.pdf",
  png_file = "05.highRisk.state.featurePlot.tsne.png",
  width = 11,
  height = 8.5,
  dpi = 600
)

#######################07.输出适合后续TCGA/GEO验证的基因表#######################
# 分别输出3个状态的top50基因
topN.bulk <- 50

bulk.gene.list <- list(
  proliferative_epithelial = head(all.marker.list[["proliferative_epithelial"]]$gene, topN.bulk),
  EMT_like_epithelial = head(all.marker.list[["EMT_like_epithelial"]]$gene, topN.bulk),
  inflammatory_secretory_epithelial = head(all.marker.list[["inflammatory_secretory_epithelial"]]$gene, topN.bulk)
)

for (nm in names(bulk.gene.list)) {
  outFile <- paste0("05.", nm, ".top", topN.bulk, ".geneList.txt")
  write.table(
    data.frame(gene = bulk.gene.list[[nm]]),
    file = outFile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

#######################08.保存对象#######################
save(
  epi.clean,
  all.marker.list,
  sig.marker.list,
  highRisk.union,
  highRisk.intersect,
  highRisk.topN,
  bulk.gene.list,
  file = "05.highRisk_signature.Rdata"
)

cat("\nDone.\n")
cat("Unified plotting style applied to 4.epiSubset.signatureExtraction.R\n")
cat("Main outputs saved:\n")
cat("05.highRisk.topGene.dotplot.vertical.pdf\n")
cat("05.highRisk.topGene.dotplot.vertical.png\n")
cat("05.highRisk.marker.avgExp.heatmap.pdf\n")
cat("05.highRisk.state.featurePlot.umap.pdf\n")
cat("05.highRisk.state.featurePlot.umap.png\n")
cat("05.highRisk.state.featurePlot.tsne.pdf\n")
cat("05.highRisk.state.featurePlot.tsne.png\n")