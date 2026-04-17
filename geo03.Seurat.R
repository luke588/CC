#install.packages("clustree")
#install.packages("Seurat")
#install.packages("harmony")
#install.packages("gridExtra")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#BiocManager::install("scrapper")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("monocle")

#install.packages("devtools")
#devtools::install_github('immunogenomics/presto')

#install.packages("assertthat")
#install.packages("SCpubr")


#######################01. 数据前期处理和校正#######################
# 引用包
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(clustree)
library(harmony)
library(assertthat)
library(SCpubr)
library(gridExtra)

logFCfilter = 1          # 后续核心差异基因筛选阈值
adjPvalFilter = 0.05     # 校正后pvalue阈值

# 设置工作目录
workDir = "/Users/donglinlu/Desktop/CC1/3.Seurat"
setwd(workDir)

# 读取数据
# 只读取第一层样本文件夹，避免把多余子目录读进去
dirs_sample <- list.dirs(workDir, recursive = FALSE, full.names = TRUE)
names(dirs_sample) <- basename(dirs_sample)

counts <- Read10X(data.dir = dirs_sample)
pbmc <- CreateSeuratObject(counts = counts, min.cells = 5, min.features = 100)

# 增加分组信息（按列名提取）
stage <- gsub("(.*?)\\..*", "\\1", colnames(pbmc))
names(stage) <- colnames(pbmc)
pbmc <- AddMetaData(object = pbmc, metadata = stage, col.name = "stage")
pbmc$Sample <- pbmc$stage

# 计算线粒体基因比例
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")

# 绘制基因特征小提琴图
pdf(file = "01.featureViolin.pdf", width = 10, height = 6.5)
p1 <- VlnPlot(pbmc, features = "nFeature_RNA", raster = FALSE) + theme(legend.position = "none")
p2 <- VlnPlot(pbmc, features = "nCount_RNA", raster = FALSE) + theme(legend.position = "none")
p3 <- VlnPlot(pbmc, features = "percent.mt", raster = FALSE) + theme(legend.position = "none")
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

# 数据过滤
# 最小改动：增加上限过滤，避免极高UMI/极高基因数细胞影响聚类
pbmc <- subset(
  x = pbmc,
  subset = nFeature_RNA > 300 &
           nFeature_RNA < 8000 &
           nCount_RNA < 100000 &
           percent.mt < 20
)

# 绘制测序深度相关性图
pdf(file = "01.featureCor.pdf", width = 13, height = 7)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1.2)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 1.2)
print(plot1 + plot2)
dev.off()

# 对数据进行标准化
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 提取高变基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)

# 输出特征方差图
top10 <- head(VariableFeatures(object = pbmc), 10)
pdf(file = "01.featureVar.pdf", width = 10, height = 6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot1 + plot2)
dev.off()



#######################02.PCA主成分分析#######################
# PCA分析前标准化
pbmc <- ScaleData(pbmc)

# PCA分析
pbmc <- RunPCA(pbmc, npcs = 30, features = VariableFeatures(object = pbmc))

# Harmony校正
pbmc <- RunHarmony(pbmc, "orig.ident")

# 绘制每个PCA成分的特征基因
pdf(file = "02.pcaGene.pdf", width = 10, height = 8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

# 绘制PCA图
pdf(file = "02.PCA.pdf", width = 7.5, height = 5)
DimPlot(object = pbmc, reduction = "pca", group.by = "stage")
dev.off()

# PCA热图
pdf(file = "02.pcaHeatmap.pdf", width = 10, height = 8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

# JackStraw分析
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file = "02.pcaJackStraw.pdf", width = 8, height = 6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()

# ElbowPlot补充
pdf(file = "02.pcaElbow.pdf", width = 8, height = 6)
ElbowPlot(pbmc, ndims = 30)
dev.off()



#######################03.细胞聚类分析和marker基因#######################
# 聚类分析
pcSelect = 20

# 注意：这里改为真正使用Harmony结果
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:pcSelect)

# 分辨率测试（保留你原有思路）
pbmc <- FindClusters(pbmc, resolution = seq(0.5, 1.2, by = 0.1))

pdf(file = "03.clustree.pdf", width = 10, height = 8)
clustree(pbmc)
dev.off()

# 选择最终分辨率
pbmc <- FindClusters(pbmc, resolution = 0.6)

# 输出聚类图形
pdf(file = "03.cluster.pdf", width = 15, height = 6)

# UMAP
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:pcSelect)
p.umap <- DimPlot(pbmc, reduction = "umap", pt.size = 1.2, label = TRUE)

# tSNE
pbmc <- RunTSNE(pbmc, reduction = "harmony", dims = 1:pcSelect)
p.tsne <- TSNEPlot(object = pbmc, pt.size = 1.2, label = TRUE)

grid.arrange(p.umap, p.tsne, ncol = 2)
dev.off()

write.table(pbmc$seurat_clusters, file = "03.Cluster.txt", quote = FALSE, sep = "\t", col.names = FALSE)

# 统计各stage的cluster比例
tab.cluster.stage <- table(pbmc$seurat_clusters, pbmc$stage)
prop.cluster.stage <- prop.table(tab.cluster.stage, margin = 2)
write.csv(as.data.frame(tab.cluster.stage), file = "03.cluster_stage_count.csv", row.names = FALSE)
write.csv(as.data.frame(prop.cluster.stage), file = "03.cluster_stage_proportion.csv", row.names = FALSE)

## 查找每个聚类的差异基因
# 第一轮注释建议 only.pos=TRUE，更适合找代表性marker
pbmc.markers <- FindAllMarkers(
  object = pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.table(pbmc.markers, file = "03.clusterMarkers_all.txt", sep = "\t", row.names = FALSE, quote = FALSE)

sig.markers <- pbmc.markers[
  (abs(as.numeric(as.vector(pbmc.markers$avg_log2FC))) > logFCfilter &
     as.numeric(as.vector(pbmc.markers$p_val_adj)) < adjPvalFilter),
]
write.table(sig.markers, file = "03.clusterMarkers.txt", sep = "\t", row.names = FALSE, quote = FALSE)

top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.table(top10, file = "03.clusterTop10Markers.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 绘制marker在每个聚类的热图
pdf(file = "03.clusterHeatmap.pdf", width = 15, height = 15)
DoHeatmap(object = pbmc, features = unique(top10$gene)) + NoLegend()
dev.off()



#######################04.SingleR包注释细胞类型#######################
pbmc_for_SingleR <- GetAssayData(pbmc, slot = "data")
clusters <- pbmc@meta.data$seurat_clusters

ref1 <- get(load("ref_Human_all.RData"))
ref2 <- get(load("ref_Hematopoietic.RData"))
ref3 <- get(load("DatabaseImmuneCellExpressionData.Rdata"))
ref4 <- get(load("BlueprintEncode_bpe.se_human.RData"))
ref5 <- get(load("HumanPrimaryCellAtlas_hpca.se_human.RData"))
ref6 <- get(load("MonacoImmuneData.Rdata"))
ref7 <- get(load("NovershternHematopoieticData.Rdata"))

singler <- SingleR(
  test = pbmc_for_SingleR,
  ref = list(ref1, ref2, ref3, ref4, ref5, ref6, ref7),
  labels = list(
    ref1$label.main,
    ref2$label.main,
    ref3$label.main,
    ref4$label.main,
    ref5$label.main,
    ref6$label.main,
    ref7$label.main
  ),
  clusters = clusters
)

# 统一命名
singler$labels <- gsub("_|-", " ", singler$labels)
singler$labels[singler$labels == "T cells, CD4+"] <- "CD4+ T cells"
singler$labels[singler$labels == "T cells, CD8+"] <- "CD8+ T cells"
singler$labels[singler$labels == "Macrophage"] <- "Macrophages"
singler$labels[singler$labels == "B cell"] <- "B cells"

clusterAnn <- as.data.frame(singler)
clusterAnn <- cbind(id = row.names(clusterAnn), clusterAnn)
clusterAnn <- clusterAnn[, c("id", "labels")]
write.table(clusterAnn, file = "04.clusterAnn.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# 输出细胞注释结果
cellAnn <- clusterAnn[match(pbmc$seurat_clusters, clusterAnn[, 1]), 2]
cellAnnOut <- cbind(names(pbmc$seurat_clusters), cellAnn)
colnames(cellAnnOut) <- c("id", "labels")
write.table(cellAnnOut, file = "04.cellAnn.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# 细胞注释后的可视化
newLabels <- singler$labels
names(newLabels) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, newLabels)
pbmc$cellType <- Idents(pbmc)

pdf(file = "04.cellAnn.pdf", width = 15, height = 6)
p1 <- DimPlot(pbmc, reduction = "umap", pt.size = 1.2, label = TRUE)
p2 <- TSNEPlot(object = pbmc, pt.size = 1.2, label = TRUE)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# 分组可视化
pdf(file = "04.group.cellAnn.pdf", width = 24, height = 6)
p1 <- DimPlot(pbmc, reduction = "umap", pt.size = 0.8, label = TRUE, split.by = "stage")
p2 <- TSNEPlot(object = pbmc, pt.size = 0.8, label = TRUE, split.by = "stage")
grid.arrange(p1, p2, ncol = 2)
dev.off()

# 细胞类型marker分析
pbmc.markers.cell <- FindAllMarkers(
  object = pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.table(pbmc.markers.cell, file = "04.cellMarkers_all.txt", sep = "\t", row.names = FALSE, quote = FALSE)

sig.cellMarkers <- pbmc.markers.cell[
  (abs(as.numeric(as.vector(pbmc.markers.cell$avg_log2FC))) > logFCfilter &
     as.numeric(as.vector(pbmc.markers.cell$p_val_adj)) < adjPvalFilter),
]
write.table(sig.cellMarkers, file = "04.cellMarkers.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#######################05.不同阶段在同一种细胞类型中的差异分析#######################
# 原脚本这里是 Control vs Disease 模板
# 现在改成适合本项目的 stage 比较逻辑：NO / N / HSIL / CA

groups <- paste0(pbmc$stage, "_", pbmc$cellType)
names(groups) <- colnames(pbmc)
pbmc <- AddMetaData(object = pbmc, metadata = groups, col.name = "group")

allGroupMarkers <- data.frame()

# 需要比较的组合
compare.list <- list(
  c("CA", "NO"),
  c("CA", "N"),
  c("CA", "HSIL"),
  c("HSIL", "N"),
  c("N", "NO")
)

for(cellName in unique(pbmc$cellType)){
  for(comp in compare.list){
    group1 <- paste0(comp[1], "_", cellName)
    group2 <- paste0(comp[2], "_", cellName)

    if((length(groups[groups == group1]) > 5) & (length(groups[groups == group2]) > 5)){
      tmp.markers <- FindMarkers(
        pbmc,
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

# 保存单细胞对象
save(pbmc, cellAnn, clusterAnn, sig.markers, sig.cellMarkers, allGroupMarkers, file = "Seurat.Rdata")


