rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(grid)
  library(scales)
  library(tidyr)
})

# =========================
# 0. unified publication theme
# =========================
theme_pub <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = base_size + 4,
        margin = margin(b = 12)
      ),
      plot.title.position = "plot",
      axis.title = element_text(
        size = base_size + 1,
        face = "bold",
        colour = "black"
      ),
      axis.text = element_text(
        size = base_size - 2,
        colour = "black"
      ),
      axis.text.x = element_text(
        size = base_size - 2,
        colour = "black",
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      axis.text.y = element_text(
        size = base_size - 2,
        colour = "black"
      ),
      legend.title = element_text(
        size = base_size,
        face = "bold",
        colour = "black"
      ),
      legend.text = element_text(
        size = base_size - 1,
        colour = "black"
      ),
      legend.position = "right",
      legend.key = element_blank(),
      legend.key.size = unit(1.0, "lines"),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.9
      ),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        size = base_size,
        face = "bold",
        colour = "black"
      ),
      plot.margin = margin(t = 16, r = 20, b = 14, l = 18)
    )
}

save_plot_pdf_png <- function(plot_obj, pdf_file, width, height, dpi = 600) {
  png_file <- sub("\\.pdf$", ".png", pdf_file)

  ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE,
    bg = "white"
  )

  ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

normalize_stage <- function(x) {
  x0 <- as.character(x)
  x1 <- trimws(x0)
  x1 <- gsub("_Epithelial cells$", "", x1, ignore.case = TRUE)
  x1 <- gsub("Epithelial cells$", "", x1, ignore.case = TRUE)
  x1 <- gsub("[[:space:]]+", "", x1)
  x1 <- toupper(x1)

  out <- ifelse(
    grepl("NO_HPV|NOHPV", x1),
    "NO_HPV",
    ifelse(
      grepl("^N_HPV$|^NHPV$|N_HPV|NHPV", x1),
      "N_HPV",
      ifelse(
        grepl("HSIL_HPV|HSILHPV", x1),
        "HSIL_HPV",
        ifelse(
          grepl("CA_HPV|CAHPV", x1),
          "CA_HPV",
          NA_character_
        )
      )
    )
  )

  out
}

# =========================
# 1. path
# =========================
setwd("/Users/donglinlu/Desktop/CC1/4.epithelial subset")

infile <- "step1_5.07_focus_clusters_6_9_10_stage_proportion.csv"
outdir <- "prioritized_program_stage_progression"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 2. read data
# =========================
if (!file.exists(infile)) {
  stop("Cannot find input file: ", infile)
}

df <- fread(infile, data.table = FALSE)

req_cols <- c("group", "seurat_clusters", "cluster_prop")
miss_cols <- setdiff(req_cols, colnames(df))
if (length(miss_cols) > 0) {
  stop(
    "Missing required columns in input file: ",
    paste(miss_cols, collapse = ", ")
  )
}

# =========================
# 3. clean / reformat
# =========================
df2 <- df %>%
  dplyr::rename(
    cluster = seurat_clusters,
    proportion = cluster_prop
  ) %>%
  mutate(
    cluster = as.character(cluster),
    proportion = as.numeric(proportion),
    stage = normalize_stage(group),
    cluster_label = case_when(
      cluster == "9"  ~ "Malignant squamous",
      cluster == "10" ~ "EMT-like",
      cluster == "6"  ~ "Proliferation core",
      TRUE ~ NA_character_
    ),
    program_label = case_when(
      cluster == "9"  ~ "malignant_squamous_epithelial",
      cluster == "10" ~ "EMT_like_epithelial",
      cluster == "6"  ~ "proliferation_core_epithelial",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(stage),
    !is.na(cluster_label),
    !is.na(proportion)
  )

# =========================
# 4. set order
# =========================
stage_levels <- c("NO_HPV", "N_HPV", "HSIL_HPV", "CA_HPV")
cluster_levels <- c("Malignant squamous", "EMT-like", "Proliferation core")

df2$stage <- factor(df2$stage, levels = stage_levels)
df2$cluster_label <- factor(df2$cluster_label, levels = cluster_levels)

df2 <- df2 %>%
  arrange(cluster_label, stage)

# 补全缺失组合，避免后面画图丢组
df2 <- df2 %>%
  tidyr::complete(
    cluster_label = factor(cluster_levels, levels = cluster_levels),
    stage = factor(stage_levels, levels = stage_levels),
    fill = list(proportion = 0)
  ) %>%
  mutate(
    cluster = case_when(
      cluster_label == "Malignant squamous" ~ "9",
      cluster_label == "EMT-like" ~ "10",
      cluster_label == "Proliferation core" ~ "6",
      TRUE ~ NA_character_
    ),
    program_label = case_when(
      cluster_label == "Malignant squamous" ~ "malignant_squamous_epithelial",
      cluster_label == "EMT-like" ~ "EMT_like_epithelial",
      cluster_label == "Proliferation core" ~ "proliferation_core_epithelial",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(cluster_label, stage)

write.csv(
  df2,
  file = file.path(outdir, "prioritized_program_stage_progression.used_table.csv"),
  row.names = FALSE
)

write.csv(
  df2 %>%
    dplyr::select(stage, cluster, cluster_label, program_label, proportion),
  file = file.path(outdir, "prioritized_program_stage_progression.used_table.annotated.csv"),
  row.names = FALSE
)

df_wide <- df2 %>%
  dplyr::select(stage, program_label, proportion) %>%
  tidyr::pivot_wider(
    names_from = stage,
    values_from = proportion
  )

write.csv(
  df_wide,
  file = file.path(outdir, "prioritized_program_stage_progression.used_table.wide.csv"),
  row.names = FALSE
)

# =========================
# 5. colors
# =========================
my_cols <- c(
  "Malignant squamous" = "#D73027",
  "EMT-like" = "#4575B4",
  "Proliferation core" = "#1A9850"
)

# =========================
# 6. y-axis upper limit
# =========================
y_max <- max(df2$proportion, na.rm = TRUE)

if (y_max <= 1.05) {
  y_upper <- ceiling(y_max * 100 * 1.12) / 100
  y_lab_fun <- label_percent(accuracy = 1)
} else {
  y_upper <- y_max * 1.12
  y_lab_fun <- label_number(accuracy = 0.01)
}

if (!is.finite(y_upper) || y_upper <= 0) {
  y_upper <- 1
}

# =========================
# 7. main line plot
# =========================
p_line <- ggplot(
  df2,
  aes(
    x = stage,
    y = proportion,
    group = cluster_label,
    color = cluster_label
  )
) +
  geom_line(linewidth = 1.35) +
  geom_point(size = 3.6) +
  scale_color_manual(values = my_cols, breaks = cluster_levels) +
  scale_y_continuous(
    labels = y_lab_fun,
    expand = expansion(mult = c(0.01, 0.05)),
    limits = c(0, y_upper)
  ) +
  labs(
    title = "Progression-associated enrichment of\nprioritized epithelial programs",
    x = NULL,
    y = "Proportion within epithelial cells",
    color = NULL
  ) +
  theme_pub(base_size = 16)

save_plot_pdf_png(
  p_line,
  pdf_file = file.path(outdir, "prioritized_program_stage_progression.lineplot.main.pdf"),
  width = 12.2,
  height = 7.0
)

# =========================
# 8. grouped bar plot (backup)
# =========================
p_bar <- ggplot(
  df2,
  aes(
    x = stage,
    y = proportion,
    fill = cluster_label
  )
) +
  geom_col(
    position = position_dodge(width = 0.76),
    width = 0.64
  ) +
  scale_fill_manual(values = my_cols, breaks = cluster_levels) +
  scale_y_continuous(
    labels = y_lab_fun,
    expand = expansion(mult = c(0.01, 0.05)),
    limits = c(0, y_upper)
  ) +
  labs(
    title = "Stage-wise proportions of\nprioritized epithelial programs",
    x = NULL,
    y = "Proportion within epithelial cells",
    fill = NULL
  ) +
  theme_pub(base_size = 16)

save_plot_pdf_png(
  p_bar,
  pdf_file = file.path(outdir, "prioritized_program_stage_progression.barplot.backup.pdf"),
  width = 12.2,
  height = 7.0
)

# =========================
# 9. optional point-line backup
# =========================
p_point <- ggplot(
  df2,
  aes(
    x = stage,
    y = proportion,
    color = cluster_label
  )
) +
  geom_point(size = 4.0) +
  geom_path(
    aes(group = cluster_label),
    linewidth = 1.1
  ) +
  scale_color_manual(values = my_cols, breaks = cluster_levels) +
  scale_y_continuous(
    labels = y_lab_fun,
    expand = expansion(mult = c(0.01, 0.05)),
    limits = c(0, y_upper)
  ) +
  labs(
    title = "Stage-wise redistribution of\nprioritized epithelial programs",
    x = NULL,
    y = "Proportion within epithelial cells",
    color = NULL
  ) +
  theme_pub(base_size = 16)

save_plot_pdf_png(
  p_point,
  pdf_file = file.path(outdir, "prioritized_program_stage_progression.pointline.backup.pdf"),
  width = 12.2,
  height = 7.0
)

cat("\nDone.\n")
cat("Main figures saved as:\n")
cat(file.path(outdir, "prioritized_program_stage_progression.lineplot.main.pdf"), "\n")
cat(file.path(outdir, "prioritized_program_stage_progression.lineplot.main.png"), "\n")
cat(file.path(outdir, "prioritized_program_stage_progression.barplot.backup.pdf"), "\n")
cat(file.path(outdir, "prioritized_program_stage_progression.barplot.backup.png"), "\n")
cat(file.path(outdir, "prioritized_program_stage_progression.pointline.backup.pdf"), "\n")
cat(file.path(outdir, "prioritized_program_stage_progression.pointline.backup.png"), "\n")
cat("Supporting tables saved as:\n")
cat(file.path(outdir, "prioritized_program_stage_progression.used_table.csv"), "\n")
cat(file.path(outdir, "prioritized_program_stage_progression.used_table.annotated.csv"), "\n")
cat(file.path(outdir, "prioritized_program_stage_progression.used_table.wide.csv"), "\n")