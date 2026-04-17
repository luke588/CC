rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
graphics.off()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# =========================================================
# Figure 1E | Functional summary of prioritized epithelial clusters
# Robust final version
# Fixed cluster order: 9 -> 10 -> 6
# =========================================================

# -----------------------------
# 0. paths
# -----------------------------
workDir <- "/Users/donglinlu/Desktop/CC1/4.epithelial subset"
setwd(workDir)

infile  <- file.path(workDir, "step1_5.06_focus_clusters_6_9_10_summary.csv")
out_pdf <- file.path(workDir, "Figure1E_functional_summary_final_unified.pdf")
out_png <- file.path(workDir, "Figure1E_functional_summary_final_unified.png")
out_csv <- file.path(workDir, "Figure1E_functional_summary_final_unified_plot_table.csv")

if (!file.exists(infile)) {
  stop("Cannot find input file: ", infile)
}

# -----------------------------
# 1. unified publication theme
# -----------------------------
theme_pub <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        size = base_size + 2,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10),
        colour = "black"
      ),
      plot.title.position = "plot",
      axis.title = element_blank(),
      axis.text.x = element_text(
        size = base_size - 2,
        colour = "black",
        face = "bold",
        angle = 30,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(
        size = base_size - 1,
        colour = "black",
        face = "bold"
      ),
      legend.title = element_text(
        size = base_size - 1,
        face = "bold",
        colour = "black"
      ),
      legend.text = element_text(
        size = base_size - 2,
        colour = "black"
      ),
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.8
      ),
      axis.line = element_blank(),
      plot.margin = margin(t = 14, r = 16, b = 14, l = 16)
    )
}

# -----------------------------
# 2. read input robustly
# -----------------------------
df <- read.csv(
  infile,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE
)

# 去掉列名两端空格和可能的 BOM
colnames(df) <- trimws(colnames(df))
colnames(df) <- sub("^\ufeff", "", colnames(df))

cat("Detected column names:\n")
print(colnames(df))
cat("\n")

# 如果没有标准 cluster 列，就强制把第一列当作 cluster
if (!("cluster" %in% colnames(df))) {
  colnames(df)[1] <- "cluster"
}

cat("Column names after fixing cluster column:\n")
print(colnames(df))
cat("\n")

required_cols <- c(
  "cluster",
  "mean_S",
  "mean_G2M",
  "mean_cycle",
  "mean_malignant_squamous",
  "mean_EMT_like",
  "mean_proliferation_core"
)

missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# -----------------------------
# 3. prepare plotting table
# -----------------------------
plot_df <- data.frame(
  cluster = as.character(df[["cluster"]]),
  mean_malignant_squamous = as.numeric(df[["mean_malignant_squamous"]]),
  mean_EMT_like = as.numeric(df[["mean_EMT_like"]]),
  mean_proliferation_core = as.numeric(df[["mean_proliferation_core"]]),
  mean_cycle = as.numeric(df[["mean_cycle"]]),
  mean_S = as.numeric(df[["mean_S"]]),
  mean_G2M = as.numeric(df[["mean_G2M"]]),
  stringsAsFactors = FALSE
)

plot_df <- plot_df[plot_df$cluster %in% c("9", "10", "6"), , drop = FALSE]

if (nrow(plot_df) == 0) {
  stop("No valid rows left after filtering clusters 9, 10, 6.")
}

plot_df$ClusterLabel <- ifelse(
  plot_df$cluster == "9", "Cluster 9\nMalignant-squamous",
  ifelse(
    plot_df$cluster == "10", "Cluster 10\nEMT-like",
    ifelse(
      plot_df$cluster == "6", "Cluster 6\nProliferation-core",
      NA
    )
  )
)

plot_df <- plot_df[!is.na(plot_df$ClusterLabel), , drop = FALSE]

plot_df$ClusterLabel <- factor(
  plot_df$ClusterLabel,
  levels = c(
    "Cluster 9\nMalignant-squamous",
    "Cluster 10\nEMT-like",
    "Cluster 6\nProliferation-core"
  )
)

# 按固定顺序排序
plot_df <- plot_df %>%
  arrange(ClusterLabel)

cat("plot_df:\n")
print(plot_df)
cat("\n")

# -----------------------------
# 4. reshape to long format
# -----------------------------
plot_long <- plot_df %>%
  pivot_longer(
    cols = c(
      mean_malignant_squamous,
      mean_EMT_like,
      mean_proliferation_core,
      mean_cycle,
      mean_S,
      mean_G2M
    ),
    names_to = "Metric",
    values_to = "Value"
  )

if (nrow(plot_long) == 0) {
  stop("pivot_longer produced an empty table.")
}

plot_long$Metric <- factor(
  plot_long$Metric,
  levels = c(
    "mean_malignant_squamous",
    "mean_EMT_like",
    "mean_proliferation_core",
    "mean_cycle",
    "mean_S",
    "mean_G2M"
  ),
  labels = c(
    "Malignant-squamous score",
    "EMT-like score",
    "Proliferation-core score",
    "Cell cycle score",
    "S phase score",
    "G2M score"
  )
)

cat("plot_long:\n")
print(plot_long)
cat("\n")

# -----------------------------
# 5. z-score within each metric
# -----------------------------
plot_long <- plot_long %>%
  group_by(Metric) %>%
  mutate(
    z = as.numeric(scale(Value))
  ) %>%
  ungroup()

plot_long$z[is.na(plot_long$z)] <- 0

plot_long <- plot_long %>%
  mutate(
    size_plot = rescale(abs(z), to = c(4.8, 13.5))
  )

write.csv(plot_long, out_csv, row.names = FALSE)

# -----------------------------
# 6. plot
# -----------------------------
p <- ggplot(
  plot_long,
  aes(x = ClusterLabel, y = Metric)
) +
  geom_point(
    aes(size = size_plot, fill = z),
    shape = 21,
    color = "black",
    stroke = 0.6
  ) +
  scale_size_identity() +
  scale_fill_gradient2(
    low = "#4C97D8",
    mid = "white",
    high = "#E76F6A",
    midpoint = 0,
    name = "Scaled score"
  ) +
  scale_y_discrete(
    limits = rev(levels(plot_long$Metric))
  ) +
  labs(
    title = "Functional summary of prioritized epithelial clusters"
  ) +
  theme_pub(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p)

# -----------------------------
# 7. save
# -----------------------------
ggsave(
  filename = out_pdf,
  plot = p,
  width = 10.2,
  height = 6.6,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  filename = out_png,
  plot = p,
  width = 10.2,
  height = 6.6,
  units = "in",
  dpi = 500
)

cat("Saved to:\n")
cat(out_pdf, "\n")
cat(out_png, "\n")
cat(out_csv, "\n")