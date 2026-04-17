rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
options(stringsAsFactors = FALSE)
# =====================================
# GLOBAL CANONICAL ORDER FOR PRIORITIZED CLUSTERS
# =====================================
CANONICAL_CLUSTER_IDS <- c("9", "10", "6")

CANONICAL_CLUSTER_LABELS <- c(
  "Cluster 9 | Malignant-squamous",
  "Cluster 10 | EMT-like",
  "Cluster 6 | Proliferation-core"
)

CANONICAL_PROGRAM_LABELS <- c(
  "malignant_squamous_epithelial",
  "EMT_like_epithelial",
  "proliferation_core_epithelial"
)

cluster_id_to_label <- c(
  "9"  = "Cluster 9 | Malignant-squamous",
  "10" = "Cluster 10 | EMT-like",
  "6"  = "Cluster 6 | Proliferation-core"
)
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(readr)
})

# =========================================================
# Figure 1E | Functional summary of prioritized epithelial clusters
# Final publication-style color version
# Input: step1_5.06_focus_clusters_6_9_10_summary.csv
# =========================================================

# -----------------------------
# 0. paths
# -----------------------------
workDir <- "/Users/donglinlu/Desktop/CC1/4.epithelial subset"
setwd(workDir)

infile  <- file.path(workDir, "step1_5.06_focus_clusters_6_9_10_summary.csv")
out_pdf <- file.path(workDir, "Figure1E_functional_summary_final.pdf")
out_png <- file.path(workDir, "Figure1E_functional_summary_final.png")

if (!file.exists(infile)) {
  stop("Cannot find input file: ", infile)
}

# -----------------------------
# 1. local theme refinement
# -----------------------------
theme_fig1e <- function(base_size = 16) {
  theme_pub(base_size = base_size) +
    theme(
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
      legend.position = "right",
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.8
      ),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(t = 14, r = 16, b = 14, l = 16)
    )
}

# -----------------------------
# 2. helpers
# -----------------------------
pick_first_col <- function(df, candidates, label) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    stop("Cannot find ", label, ". Checked: ", paste(candidates, collapse = " | "))
  }
  hit[1]
}

# -----------------------------
# 3. read input
# -----------------------------
df <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)

message("Input columns:")
print(colnames(df))
print(df)

# -----------------------------
# 4. detect columns
# -----------------------------
cluster_col <- colnames(df)[str_detect(tolower(colnames(df)), "cluster|seurat")][1]
if (is.na(cluster_col)) {
  cluster_col <- colnames(df)[1]
}

cycle_col <- pick_first_col(
  df,
  c("mean_cycle"),
  "mean cycle column"
)

s_col <- pick_first_col(
  df,
  c("mean_S", "mean_s"),
  "mean S-score column"
)

g2m_col <- pick_first_col(
  df,
  c("mean_G2M", "mean_g2m"),
  "mean G2M-score column"
)

prolif_col <- pick_first_col(
  df,
  c(
    "mean_prolif_state",
    "mean_proliferation_core",
    "mean_proliferation_core_score",
    "mean_proliferation",
    "mean_prolif"
  ),
  "proliferation metric column"
)

emt_col <- pick_first_col(
  df,
  c(
    "mean_emt_state",
    "mean_EMT_like",
    "mean_EMT_like_score",
    "mean_emt",
    "mean_EMT"
  ),
  "EMT-like metric column"
)

message("Detected columns:")
message("cluster_col = ", cluster_col)
message("cycle_col   = ", cycle_col)
message("s_col       = ", s_col)
message("g2m_col     = ", g2m_col)
message("prolif_col  = ", prolif_col)
message("emt_col     = ", emt_col)

# -----------------------------
# 5. build plotting dataframe
# -----------------------------
plot_df <- df %>%
  transmute(
    ClusterID = as.character(.data[[cluster_col]]),
    mean_cycle = as.numeric(.data[[cycle_col]]),
    mean_S = as.numeric(.data[[s_col]]),
    mean_G2M = as.numeric(.data[[g2m_col]]),
    mean_prolif_state = as.numeric(.data[[prolif_col]]),
    mean_emt_state = as.numeric(.data[[emt_col]])
  )

# -----------------------------
# 6. cluster mapping
# -----------------------------
plot_df <- plot_df %>%
  mutate(
    ClusterLabel = case_when(
      ClusterID == "10" ~ "Cluster 10 | EMT-like",
      ClusterID == "6"  ~ "Cluster 6 | Proliferation-core",
      ClusterID == "9"  ~ "Cluster 9 | Malignant-squamous",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(ClusterLabel))

plot_df$ClusterLabel <- factor(
  plot_df$ClusterLabel,
  levels = c(
    "Cluster 10 | EMT-like",
    "Cluster 6 | Proliferation-core",
    "Cluster 9 | Malignant-squamous"
  )
)

if (nrow(plot_df) == 0) {
  stop("No matched clusters found among 6 / 9 / 10 in input file.")
}

# -----------------------------
# 7. reshape
# -----------------------------
plot_long <- plot_df %>%
  pivot_longer(
    cols = c(mean_cycle, mean_S, mean_G2M, mean_prolif_state, mean_emt_state),
    names_to = "Metric",
    values_to = "Value"
  )

plot_long$Metric <- factor(
  plot_long$Metric,
  levels = c("mean_cycle", "mean_S", "mean_G2M", "mean_prolif_state", "mean_emt_state"),
  labels = c("Cell cycle", "S score", "G2M score", "Proliferation score", "EMT-like score")
)

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

# -----------------------------
# 8. plot
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
  labs(
    title = "Functional summary of prioritized epithelial clusters"
  ) +
  theme_fig1e(base_size = 16)

# -----------------------------
# 9. save
# -----------------------------
ggsave(
  filename = out_pdf,
  plot = p,
  width = 10.2,
  height = 6.6,
  units = "in",
  useDingbats = FALSE,
  bg = "white"
)

ggsave(
  filename = out_png,
  plot = p,
  width = 10.2,
  height = 6.6,
  units = "in",
  dpi = 500,
  bg = "white"
)

message("Saved to: ", out_pdf)
message("Saved to: ", out_png)