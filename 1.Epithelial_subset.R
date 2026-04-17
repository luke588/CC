rm(list = ls())
options(stringsAsFactors = FALSE)

source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")

#install.packages("clustree")
#install.packages("Seurat")
#install.packages("harmony")
#install.packages("gridExtra")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")

#install.packages("assertthat")
#install.packages("SCpubr")

logFCfilter <- 1
adjPvalFilter <- 0.05

#######################01.上皮细胞亚群提取和预处理#######################
#引用包
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(clustree)
library(harmony)
library(assertthat)
library(SCpubr)
library(gridExtra)
library(ggplot2)
library(grid)

#设置工作目录
workDir <- "/Users/donglinlu/Desktop/CC1/4.epithelial subset"
setwd(workDir)

#读取上一步全图谱分析保存的对象
load("/Users/donglinlu/Desktop/CC1/3.Seurat/Seurat.Rdata")

#查看当前cluster和细胞类型
print(table(pbmc$seurat_clusters))
print(table(pbmc$cellType))

#根据全图谱marker结果提取上皮相关cluster
#建议纳入的cluster：3,4,8,16,17,18,19,22,24,25
epi.cluster <- c("3", "4", "8", "16", "17", "18", "19", "22", "24", "25")

#提取上皮相关cluster
pbmc$seurat_clusters <- as.character(pbmc$seurat_clusters)
epi <- subset(pbmc, subset = seurat_clusters %in% epi.cluster)

#恢复身份为cluster
Idents(epi) <- epi$seurat_clusters

#重新计算线粒体比例
epi[["percent.mt"]] <- PercentageFeatureSet(object = epi, pattern = "^MT-")

#绘制QC图
p1 <- VlnPlot(epi, features = "nFeature_RNA", raster = FALSE) +
  labs(title = "nFeature_RNA") +
  theme_pub() +
  theme(legend.position = "none")

p2 <- VlnPlot(epi, features = "nCount_RNA", raster = FALSE) +
  labs(title = "nCount_RNA") +
  theme_pub() +
  theme(legend.position = "none")

p3 <- VlnPlot(epi, features = "percent.mt", raster = FALSE) +
  labs(title = "percent.mt") +
  theme_pub() +
  theme(legend.position = "none")

pdf(file = "02.epiSubset.featureViolin.pdf", width = 11.5, height = 7)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

plot1 <- FeatureScatter(
  object = epi,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  pt.size = 1.2
) +
  labs(title = "nCount_RNA vs nFeature_RNA") +
  theme_pub()

plot2 <- FeatureScatter(
  object = epi,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  pt.size = 1.2
) +
  labs(title = "nCount_RNA vs percent.mt") +
  theme_pub()

pdf(file = "02.epiSubset.featureCor.pdf", width = 14, height = 7.5)
print(plot1 + plot2)
dev.off()

#对上皮subset再次过滤
epi <- subset(
  x = epi,
  subset = nFeature_RNA > 300 &
    nFeature_RNA < 8000 &
    nCount_RNA < 100000 &
    percent.mt < 20
)

#标准化
epi <- NormalizeData(object = epi, normalization.method = "LogNormalize", scale.factor = 10000)

#高变基因
epi <- FindVariableFeatures(object = epi, selection.method = "vst", nfeatures = 1500)

#输出高变基因图
top10 <- head(VariableFeatures(object = epi), 10)

plot_var <- VariableFeaturePlot(object = epi) +
  labs(title = "Highly variable features") +
  theme_pub()

plot_var_label <- LabelPoints(plot = plot_var, points = top10, repel = TRUE) +
  theme_pub()

save_pub(
  plot_obj = plot_var_label,
  pdf_file = "02.epiSubset.featureVar.pdf",
  png_file = NULL,
  width = 11,
  height = 6.5
)

#######################02.PCA主成分分析#######################
epi <- ScaleData(epi)
epi <- RunPCA(epi, npcs = 30, features = VariableFeatures(object = epi))

#Harmony校正
epi <- RunHarmony(epi, "orig.ident")

#绘制每个PCA成分的特征基因
pdf(file = "03.epiSubset.pcaGene.pdf", width = 11, height = 8.5)
VizDimLoadings(object = epi, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

#绘制PCA图
p_pca <- DimPlot(object = epi, reduction = "pca", group.by = "stage") +
  labs(title = "PCA plot by stage") +
  theme_pub()

save_pub(
  plot_obj = p_pca,
  pdf_file = "03.epiSubset.PCA.pdf",
  png_file = NULL,
  width = 8.5,
  height = 6
)

#PCA热图
pdf(file = "03.epiSubset.pcaHeatmap.pdf", width = 11, height = 8.5)
DimHeatmap(object = epi, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

#JackStraw
epi <- JackStraw(object = epi, num.replicate = 100)
epi <- ScoreJackStraw(object = epi, dims = 1:20)

pdf(file = "03.epiSubset.pcaJackStraw.pdf", width = 8.5, height = 6.5)
JackStrawPlot(object = epi, dims = 1:20)
dev.off()

#ElbowPlot
p_elbow <- ElbowPlot(epi, ndims = 30) +
  labs(title = "Elbow plot") +
  theme_pub()

save_pub(
  plot_obj = p_elbow,
  pdf_file = "03.epiSubset.pcaElbow.pdf",
  png_file = NULL,
  width = 8.5,
  height = 6
)

#######################03.上皮细胞再聚类分析#######################
pcSelect <- 15

epi <- FindNeighbors(object = epi, reduction = "harmony", dims = 1:pcSelect)

#分辨率测试
epi <- FindClusters(epi, resolution = seq(0.2, 1.0, by = 0.1))
pdf(file = "03.epiSubset.clustree.pdf", width = 11, height = 8.5)
clustree(epi)
dev.off()

#根据clustree选择分辨率
epi <- FindClusters(object = epi, resolution = 0.4)

#降维可视化
epi <- RunUMAP(object = epi, reduction = "harmony", dims = 1:pcSelect)
epi <- RunTSNE(object = epi, reduction = "harmony", dims = 1:pcSelect)

p_umap_cluster <- DimPlot(
  epi,
  reduction = "umap",
  pt.size = 1.0,
  label = TRUE,
  label.size = 5.5,
  repel = TRUE
) +
  labs(title = "UMAP of epithelial subset") +
  theme_pub()

p_tsne_cluster <- TSNEPlot(
  object = epi,
  pt.size = 1.0,
  label = TRUE,
  label.size = 5.5,
  repel = TRUE
) +
  labs(title = "tSNE of epithelial subset") +
  theme_pub()

pdf(file = "01.epiSubset.cluster.pdf", width = 16, height = 7)
grid.arrange(p_umap_cluster, p_tsne_cluster, ncol = 2)
dev.off()

write.table(
  epi$seurat_clusters,
  file = "01.epiSubset.Cluster.txt",
  quote = FALSE,
  sep = "\t",
  col.names = FALSE
)

#统计各stage中的cluster数量和比例
tab.cluster.stage <- table(epi$seurat_clusters, epi$stage)
prop.cluster.stage <- prop.table(tab.cluster.stage, margin = 2)

write.csv(as.data.frame(tab.cluster.stage), file = "01.epiSubset.cluster_stage_count.csv", row.names = FALSE)
write.csv(as.data.frame(prop.cluster.stage), file = "01.epiSubset.cluster_stage_proportion.csv", row.names = FALSE)

#查找每个聚类的差异基因
epi.markers <- FindAllMarkers(
  object = epi,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.table(epi.markers, file = "01.epiSubset.clusterMarkers_all.txt", sep = "\t", row.names = FALSE, quote = FALSE)

sig.markers <- epi.markers[
  (abs(as.numeric(as.vector(epi.markers$avg_log2FC))) > logFCfilter &
     as.numeric(as.vector(epi.markers$p_val_adj)) < adjPvalFilter),
]

write.table(sig.markers, file = "01.epiSubset.clusterMarkers.txt", sep = "\t", row.names = FALSE, quote = FALSE)

top10 <- epi.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

pdf(file = "01.epiSubset.clusterHeatmap.pdf", width = 16, height = 15)
DoHeatmap(object = epi, features = unique(top10$gene)) + NoLegend()
dev.off()

#######################04.SingleR对上皮subset再注释#######################
#兼容不同Seurat版本提取表达矩阵
if ("data" %in% slotNames(epi[["RNA"]])) {
  epi_for_SingleR <- GetAssayData(epi, assay = "RNA", slot = "data")
} else {
  epi_for_SingleR <- GetAssayData(epi, assay = "RNA", layer = "data")
}

clusters <- epi@meta.data$seurat_clusters

#读取SingleR参考数据
ref1 <- get(load("/Users/donglinlu/Desktop/CC1/3.Seurat/ref_Human_all.RData"))
ref2 <- get(load("/Users/donglinlu/Desktop/CC1/3.Seurat/ref_Hematopoietic.RData"))
ref3 <- get(load("/Users/donglinlu/Desktop/CC1/3.Seurat/DatabaseImmuneCellExpressionData.Rdata"))
ref4 <- get(load("/Users/donglinlu/Desktop/CC1/3.Seurat/BlueprintEncode_bpe.se_human.RData"))
ref5 <- get(load("/Users/donglinlu/Desktop/CC1/3.Seurat/HumanPrimaryCellAtlas_hpca.se_human.RData"))
ref6 <- get(load("/Users/donglinlu/Desktop/CC1/3.Seurat/MonacoImmuneData.Rdata"))
ref7 <- get(load("/Users/donglinlu/Desktop/CC1/3.Seurat/NovershternHematopoieticData.Rdata"))

#运行SingleR
singler <- SingleR(
  test = epi_for_SingleR,
  ref = list(ref1, ref2, ref3, ref4, ref5, ref6, ref7),
  labels = list(
    ref1$label.main, ref2$label.main, ref3$label.main,
    ref4$label.main, ref5$label.main, ref6$label.main, ref7$label.main
  ),
  clusters = clusters
)

#统一命名
singler$labels <- gsub("_|-", " ", singler$labels)

clusterAnn <- as.data.frame(singler)
clusterAnn <- cbind(id = row.names(clusterAnn), clusterAnn)
clusterAnn <- clusterAnn[, c("id", "labels")]
write.table(clusterAnn, file = "01.epiSubset.clusterAnn.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#输出细胞注释结果
cellAnn <- clusterAnn[match(epi$seurat_clusters, clusterAnn[, 1]), 2]
cellAnnOut <- cbind(names(epi$seurat_clusters), cellAnn)
colnames(cellAnnOut) <- c("id", "labels")
write.table(cellAnnOut, file = "01.epiSubset.cellAnn.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#可视化注释结果
newLabels <- singler$labels
names(newLabels) <- levels(epi)
epi <- RenameIdents(epi, newLabels)
epi$cellType <- Idents(epi)

# 单独：UMAP
p_umap_ann <- DimPlot(
  epi,
  reduction = "umap",
  pt.size = 1.0,
  label = TRUE,
  label.size = 5.5,
  repel = TRUE
) +
  labs(title = "SingleR annotation: UMAP") +
  theme_pub()

save_pub(
  plot_obj = p_umap_ann,
  pdf_file = "01.epiSubset.clusterAnn.umap.pdf",
  png_file = NULL,
  width = 8.5,
  height = 7
)

# 单独：tSNE
p_tsne_ann <- TSNEPlot(
  object = epi,
  pt.size = 1.0,
  label = TRUE,
  label.size = 5.5,
  repel = TRUE
) +
  labs(title = "SingleR annotation: tSNE") +
  theme_pub()

save_pub(
  plot_obj = p_tsne_ann,
  pdf_file = "01.epiSubset.clusterAnn.tsne.pdf",
  png_file = NULL,
  width = 8.5,
  height = 7
)

# 单独：UMAP split by stage
p_umap_split <- DimPlot(
  epi,
  reduction = "umap",
  pt.size = 0.75,
  label = TRUE,
  label.size = 4.5,
  repel = TRUE,
  split.by = "stage"
) +
  labs(title = "SingleR annotation by stage: UMAP") +
  theme_pub()

save_pub(
  plot_obj = p_umap_split,
  pdf_file = "01.epiSubset.group.cellAnn.umap.pdf",
  png_file = NULL,
  width = 10.5,
  height = 7
)

# 单独：tSNE split by stage
p_tsne_split <- TSNEPlot(
  object = epi,
  pt.size = 0.75,
  label = TRUE,
  label.size = 4.5,
  repel = TRUE,
  split.by = "stage"
) +
  labs(title = "SingleR annotation by stage: tSNE") +
  theme_pub()

save_pub(
  plot_obj = p_tsne_split,
  pdf_file = "01.epiSubset.group.cellAnn.tsne.pdf",
  png_file = NULL,
  width = 10.5,
  height = 7
)

#######################05.阶段间差异分析#######################
groups <- paste0(epi$stage, "_", epi$cellType)
names(groups) <- colnames(epi)
epi <- AddMetaData(object = epi, metadata = groups, col.name = "group")

allGroupMarkers <- data.frame()

compare.list <- list(
  c("CA", "NO"),
  c("CA", "N"),
  c("CA", "HSIL"),
  c("HSIL", "N"),
  c("N", "NO")
)

for (cellName in unique(epi$cellType)) {
  for (comp in compare.list) {
    group1 <- paste0(comp[1], "_", cellName)
    group2 <- paste0(comp[2], "_", cellName)

    if ((length(groups[groups == group1]) > 5) & (length(groups[groups == group2]) > 5)) {
      tmp.markers <- FindMarkers(
        epi,
        ident.1 = group1,
        ident.2 = group2,
        group.by = "group",
        logfc.threshold = 0.1
      )
      tmp.markers$gene <- row.names(tmp.markers)
      tmp.markers$cellType <- cellName
      tmp.markers$compare <- paste0(comp[1], "_vs_", comp[2])
      allGroupMarkers <- rbind(allGroupMarkers, tmp.markers)

      sig.markersGroup <- tmp.markers[
        (abs(as.numeric(as.vector(tmp.markers$avg_log2FC))) > logFCfilter &
           as.numeric(as.vector(tmp.markers$p_val_adj)) < adjPvalFilter),
      ]
      sig.markersGroup <- cbind(Gene = row.names(sig.markersGroup), sig.markersGroup)

      outFile <- paste0("05.", cellName, ".", comp[1], "_vs_", comp[2], ".diffGene.txt")
      write.table(sig.markersGroup, file = outFile, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
}

write.table(allGroupMarkers, file = "05.allGroupMarkers.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#######################06.保存对象#######################
save(epi, cellAnn, clusterAnn, sig.markers, allGroupMarkers, file = "Epithelial_subset.Rdata")

cat("\nDone.\n")
cat("Unified plotting style applied to 1.Epithelial_subset.R\n")
cat("Separated files saved:\n")
cat("01.epiSubset.clusterAnn.umap.pdf\n")
cat("01.epiSubset.clusterAnn.tsne.pdf\n")
cat("01.epiSubset.group.cellAnn.umap.pdf\n")
cat("01.epiSubset.group.cellAnn.tsne.pdf\n")