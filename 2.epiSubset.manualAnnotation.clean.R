rm(list = ls())
options(stringsAsFactors = FALSE)

source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")

#######################01.手动注释上皮细胞亚群#######################
#引用包
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(grid)

#设置工作目录
workDir <- "/Users/donglinlu/Desktop/CC1/4.epithelial subset"
setwd(workDir)

#读取上一步保存的上皮subset对象
load("Epithelial_subset.Rdata")

#确保当前身份为cluster
epi$seurat_clusters <- as.character(epi$seurat_clusters)
Idents(epi) <- epi$seurat_clusters

#查看cluster数量
print(table(epi$seurat_clusters))

#根据marker手动注释cluster
#注释依据见 01.epiSubset.clusterTop10Markers.txt
#0: APOB/MIA/FN1/CDH6，偏EMT-like epithelial
#1: TM4SF4/CTSE/LGALS4/CDH17，偏secretory-like epithelial
#2: SBSN/KRT1/DSG1/KRTDAP/KRT13，偏keratinizing squamous epithelial
#3: IGLC2/IGKC/IGHG3/JCHAIN，plasma contamination
#4: NTS/MMP13/ADH1C/CEL，aberrant secretory epithelial
#5: TCN1/FCGBP/DMBT1/CXCL5，inflammatory secretory epithelial
#6: KIF20A/GTSE1/HMMR/NEK2/DEPDC1/ASPM，proliferative epithelial
#7: IFNG/CCL4/GZMB/NKG7/PRF1，T/NK contamination
#8: CXCL14/PRAC1/DCN/ABI3BP，basal stromal-like epithelial
#9: KRT78/SPRR2D/SPRR2E/CEACAM7，superficial keratinizing epithelial
#10: CALML5/CASP14/UBD，terminally differentiated squamous epithelial
#11: IGFBP7/CHI3L1/NGFR/HEG1，basal EMT-like epithelial
#12: COL3A1/LUM/COL6A3/COL1A1/MGP，fibro-epithelial contamination
#13: TRPM5/GNG13/PROX1/LRMP，rare epithelial-like cells

new.cluster.ids <- c(
  "EMT_like_epithelial",                            # 0
  "secretory_like_epithelial",                     # 1
  "keratinizing_squamous_epithelial",              # 2
  "plasma_contamination",                          # 3
  "aberrant_secretory_epithelial",                 # 4
  "inflammatory_secretory_epithelial",             # 5
  "proliferative_epithelial",                      # 6
  "T_NK_contamination",                            # 7
  "basal_stromal_like_epithelial",                 # 8
  "superficial_keratinizing_epithelial",           # 9
  "terminally_differentiated_squamous_epithelial", # 10
  "basal_EMT_like_epithelial",                     # 11
  "fibro_epithelial_contamination",                # 12
  "rare_epithelial_like_cells"                     # 13
)

names(new.cluster.ids) <- levels(epi)
epi <- RenameIdents(epi, new.cluster.ids)
epi$manualCellType <- Idents(epi)

#输出cluster注释表
clusterAnn <- data.frame(
  cluster = names(new.cluster.ids),
  label = as.vector(new.cluster.ids),
  stringsAsFactors = FALSE
)

write.table(
  clusterAnn,
  file = "02.epiSubset.manualClusterAnn.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

#输出每个细胞的注释结果
cellAnnOut <- data.frame(
  id = colnames(epi),
  cluster = epi$seurat_clusters,
  label = epi$manualCellType,
  stringsAsFactors = FALSE
)

write.table(
  cellAnnOut,
  file = "02.epiSubset.manualCellAnn.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

#######################02.手动注释结果可视化#######################
# 1. UMAP单独图：保留标签，方便总览
p_umap_manual <- DimPlot(
  epi,
  reduction = "umap",
  pt.size = 1.0,
  label = TRUE,
  label.size = 5.5,
  repel = TRUE
) +
  NoLegend() +
  labs(title = "Manual annotation: UMAP") +
  theme_pub() +
  theme(legend.position = "none")

save_pub(
  plot_obj = p_umap_manual,
  pdf_file = "02.epiSubset.manualClusterAnn.umap.pdf",
  png_file = NULL,
  width = 11,
  height = 8.5
)

# 2. tSNE单独图：保留标签，方便总览
p_tsne_manual <- TSNEPlot(
  object = epi,
  pt.size = 1.0,
  label = TRUE,
  label.size = 5.5,
  repel = TRUE
) +
  NoLegend() +
  labs(title = "Manual annotation: tSNE") +
  theme_pub() +
  theme(legend.position = "none")

save_pub(
  plot_obj = p_tsne_manual,
  pdf_file = "02.epiSubset.manualClusterAnn.tsne.pdf",
  png_file = NULL,
  width = 11,
  height = 8.5
)

# 3. UMAP按stage分组单独图：不在图中写标签，只保留右侧图例
p_umap_manual_split <- DimPlot(
  epi,
  reduction = "umap",
  pt.size = 0.8,
  label = FALSE,
  split.by = "stage"
) +
  labs(title = "Manual annotation by stage: UMAP") +
  theme_pub() +
  theme(
    legend.position = "right",
    legend.key.size = unit(1.0, "lines")
  )

save_pub(
  plot_obj = p_umap_manual_split,
  pdf_file = "02.epiSubset.group.manualClusterAnn.umap.pdf",
  png_file = NULL,
  width = 20,
  height = 7
)

# 4. tSNE按stage分组单独图：不在图中写标签，只保留右侧图例
p_tsne_manual_split <- TSNEPlot(
  object = epi,
  pt.size = 0.8,
  label = FALSE,
  split.by = "stage"
) +
  labs(title = "Manual annotation by stage: tSNE") +
  theme_pub() +
  theme(
    legend.position = "right",
    legend.key.size = unit(1.0, "lines")
  )

save_pub(
  plot_obj = p_tsne_manual_split,
  pdf_file = "02.epiSubset.group.manualClusterAnn.tsne.pdf",
  png_file = NULL,
  width = 20,
  height = 7
)

#######################03.统计不同阶段各亚群比例#######################
#按手动注释后的亚群统计
tab.celltype.stage <- table(epi$manualCellType, epi$stage)
prop.celltype.stage <- prop.table(tab.celltype.stage, margin = 2)

write.csv(
  as.data.frame(tab.celltype.stage),
  file = "02.epiSubset.manualCellType_stage_count.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(prop.celltype.stage),
  file = "02.epiSubset.manualCellType_stage_proportion.csv",
  row.names = FALSE
)

#######################04.去除污染群，构建clean epithelial对象#######################
#明确去除污染群
remove.cluster <- c(
  "plasma_contamination",
  "T_NK_contamination",
  "fibro_epithelial_contamination"
)

epi.clean <- subset(epi, subset = !manualCellType %in% remove.cluster)

#重新设置identity
Idents(epi.clean) <- epi.clean$manualCellType

#输出clean对象的cluster表
write.table(
  data.frame(
    id = colnames(epi.clean),
    label = epi.clean$manualCellType,
    stringsAsFactors = FALSE
  ),
  file = "03.epiSubset.clean.manualCellAnn.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# 1. clean UMAP单独图：保留标签
p_umap_clean <- DimPlot(
  epi.clean,
  reduction = "umap",
  pt.size = 1.0,
  label = TRUE,
  label.size = 5.5,
  repel = TRUE
) +
  NoLegend() +
  labs(title = "Clean epithelial-state landscape: UMAP") +
  theme_pub() +
  theme(legend.position = "none")

save_pub(
  plot_obj = p_umap_clean,
  pdf_file = "03.epiSubset.clean.manualClusterAnn.umap.pdf",
  png_file = NULL,
  width = 11,
  height = 8.5
)

# 2. clean tSNE单独图：保留标签
p_tsne_clean <- TSNEPlot(
  object = epi.clean,
  pt.size = 1.0,
  label = TRUE,
  label.size = 5.5,
  repel = TRUE
) +
  NoLegend() +
  labs(title = "Clean epithelial-state landscape: tSNE") +
  theme_pub() +
  theme(legend.position = "none")

save_pub(
  plot_obj = p_tsne_clean,
  pdf_file = "03.epiSubset.clean.manualClusterAnn.tsne.pdf",
  png_file = NULL,
  width = 11,
  height = 8.5
)

# 3. clean UMAP按stage分组单独图：不在图中写标签，只保留右侧图例
p_umap_clean_split <- DimPlot(
  epi.clean,
  reduction = "umap",
  pt.size = 0.8,
  label = FALSE,
  split.by = "stage"
) +
  labs(title = "Clean epithelial states by stage: UMAP") +
  theme_pub() +
  theme(
    legend.position = "right",
    legend.key.size = unit(1.0, "lines")
  )

save_pub(
  plot_obj = p_umap_clean_split,
  pdf_file = "03.epiSubset.clean.group.manualClusterAnn.umap.pdf",
  png_file = NULL,
  width = 19,
  height = 7
)

# 4. clean tSNE按stage分组单独图：不在图中写标签，只保留右侧图例
p_tsne_clean_split <- TSNEPlot(
  object = epi.clean,
  pt.size = 0.8,
  label = FALSE,
  split.by = "stage"
) +
  labs(title = "Clean epithelial states by stage: tSNE") +
  theme_pub() +
  theme(
    legend.position = "right",
    legend.key.size = unit(1.0, "lines")
  )

save_pub(
  plot_obj = p_tsne_clean_split,
  pdf_file = "03.epiSubset.clean.group.manualClusterAnn.tsne.pdf",
  png_file = NULL,
  width = 19,
  height = 7
)

#统计clean对象不同阶段比例
tab.clean.stage <- table(epi.clean$manualCellType, epi.clean$stage)
prop.clean.stage <- prop.table(tab.clean.stage, margin = 2)

write.csv(
  as.data.frame(tab.clean.stage),
  file = "03.epiSubset.cleanCellType_stage_count.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(prop.clean.stage),
  file = "03.epiSubset.cleanCellType_stage_proportion.csv",
  row.names = FALSE
)

#######################05.保存对象#######################
save(
  epi, epi.clean, clusterAnn, cellAnnOut,
  file = "02.Epithelial_subset.manualAnnotation.Rdata"
)

cat("\nDone.\n")
cat("Unified plotting style applied to 2.epiSubset.manualAnnotation.clean.R\n")
cat("Main outputs saved:\n")
cat("02.epiSubset.manualClusterAnn.umap.pdf\n")
cat("02.epiSubset.manualClusterAnn.tsne.pdf\n")
cat("02.epiSubset.group.manualClusterAnn.umap.pdf\n")
cat("02.epiSubset.group.manualClusterAnn.tsne.pdf\n")
cat("03.epiSubset.clean.manualClusterAnn.umap.pdf\n")
cat("03.epiSubset.clean.manualClusterAnn.tsne.pdf\n")
cat("03.epiSubset.clean.group.manualClusterAnn.umap.pdf\n")
cat("03.epiSubset.clean.group.manualClusterAnn.tsne.pdf\n")