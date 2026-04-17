rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(survival)
  library(stringr)
})

# ============================================================
# Figure 2E
# GSE44001 grouped Cox summary forest plot
# refined: fix right-side text overlap
# ============================================================

# -----------------------------
# 0. unified theme
# -----------------------------
theme_pub <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold", color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
      axis.line = element_blank(),
      plot.margin = margin(18, 110, 18, 18)
    )
}

# -----------------------------
# 1. helper
# -----------------------------
find_preferred_file <- function(project_root, candidate_paths = NULL, pattern = NULL,
                                prefer_keywords = c()) {
  if (!is.null(candidate_paths)) {
    for (p in candidate_paths) {
      if (file.exists(p)) return(normalizePath(p, winslash = "/"))
    }
  }

  if (!is.null(pattern)) {
    hits <- list.files(
      project_root,
      pattern = pattern,
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
    hits <- hits[!dir.exists(hits)]

    if (length(hits) == 0) return(NA_character_)

    score <- rep(0, length(hits))
    if (length(prefer_keywords) > 0) {
      for (kw in prefer_keywords) {
        score <- score + ifelse(str_detect(tolower(hits), tolower(kw)), 1, 0)
      }
    }

    hits <- hits[order(-score, nchar(hits))]
    return(normalizePath(hits[1], winslash = "/"))
  }

  return(NA_character_)
}

fmt_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    sprintf("%.3f", p)
  )
}

extract_grouped_cox <- function(df, group_col, display_name) {
  tmp <- df %>%
    filter(
      !is.na(.data[[group_col]]),
      !is.na(dfs_months),
      !is.na(dfs_status)
    ) %>%
    mutate(
      .group_tmp = factor(.data[[group_col]], levels = c("Low risk", "High risk"))
    )

  fit <- coxph(Surv(dfs_months, dfs_status) ~ .group_tmp, data = tmp)
  sm  <- summary(fit)

  coef_tab <- as.data.frame(sm$coefficients)
  ci_tab   <- as.data.frame(sm$conf.int)

  coef_tab$term <- rownames(coef_tab)
  ci_tab$term   <- rownames(ci_tab)

  idx <- grep("High risk", coef_tab$term)
  if (length(idx) == 0) idx <- 1 else idx <- idx[1]

  data.frame(
    signature = display_name,
    HR = as.numeric(ci_tab$`exp(coef)`[idx]),
    lower = as.numeric(ci_tab$`lower .95`[idx]),
    upper = as.numeric(ci_tab$`upper .95`[idx]),
    p = as.numeric(coef_tab$`Pr(>|z|)`[idx]),
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# 2. paths
# -----------------------------
project_root <- "/Users/donglinlu/Desktop/CC1"
figure_dir   <- file.path(project_root, "7.progression validation")
out_dir      <- file.path(figure_dir, "Figure2E_refined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidate_group_files <- c(
  file.path(project_root, "6.bulk_validation", "GSE44001", "02.signature_score_survival", "08.GSE44001_clinical_with_scores_and_groups.csv"),
  file.path(project_root, "6.bulk_validation", "GSE44001", "02.signature_score_survival", "05.GSE44001_clinical_with_scores_and_groups.csv"),
  file.path(project_root, "6.bulk_validation", "02.signature_score_survival", "08.GSE44001_clinical_with_scores_and_groups.csv"),
  file.path(project_root, "6.bulk_validation", "02.signature_score_survival", "05.GSE44001_clinical_with_scores_and_groups.csv")
)

group_file <- find_preferred_file(
  project_root = project_root,
  candidate_paths = candidate_group_files,
  pattern = "GSE44001_clinical_with_scores_and_groups\\.csv$",
  prefer_keywords = c("gse44001", "02.signature_score_survival")
)

if (is.na(group_file)) {
  stop("Cannot find GSE44001_clinical_with_scores_and_groups.csv automatically.")
}

message("Using input file: ", group_file)

# -----------------------------
# 3. read data
# -----------------------------
ana <- read.csv(group_file, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c(
  "dfs_months", "dfs_status",
  "malignant_squamous_group",
  "EMT_like_group",
  "proliferation_core_group"
)

missing_cols <- required_cols[!required_cols %in% colnames(ana)]
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

ana <- ana %>%
  mutate(
    dfs_months = as.numeric(dfs_months),
    dfs_status = as.numeric(dfs_status),
    malignant_squamous_group = factor(malignant_squamous_group, levels = c("Low risk", "High risk")),
    EMT_like_group = factor(EMT_like_group, levels = c("Low risk", "High risk")),
    proliferation_core_group = factor(proliferation_core_group, levels = c("Low risk", "High risk"))
  )

# -----------------------------
# 4. run grouped Cox
# fixed display order: top -> bottom
# -----------------------------
forest_df <- bind_rows(
  extract_grouped_cox(ana, "proliferation_core_group", "Proliferation-core"),
  extract_grouped_cox(ana, "EMT_like_group", "EMT-like"),
  extract_grouped_cox(ana, "malignant_squamous_group", "Malignant-squamous")
)

forest_df <- forest_df %>%
  mutate(
    y = c(3, 2, 1),
    hr_ci_label = sprintf("%.2f (%.2f\u2013%.2f)", HR, lower, upper),
    p_label = fmt_p(p)
  )

write.csv(
  forest_df,
  file.path(out_dir, "Figure2E_GSE44001_grouped_cox_summary_table.csv"),
  row.names = FALSE
)

# -----------------------------
# 5. layout positions
# -----------------------------
x_data_max <- max(forest_df$upper, na.rm = TRUE)
x_data_max <- max(8, ceiling(x_data_max * 1.10))

x_hr_col   <- x_data_max + 1.8
x_p_col    <- x_data_max + 5.0
x_plot_max <- x_data_max + 6.6

cap_half_height <- 0.08

# -----------------------------
# 6. plot
# -----------------------------
p <- ggplot(forest_df) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.9,
    color = "grey40"
  ) +
  geom_segment(
    aes(x = lower, xend = upper, y = y, yend = y, color = signature),
    linewidth = 1.2
  ) +
  geom_segment(
    aes(x = lower, xend = lower, y = y - cap_half_height, yend = y + cap_half_height, color = signature),
    linewidth = 1.2
  ) +
  geom_segment(
    aes(x = upper, xend = upper, y = y - cap_half_height, yend = y + cap_half_height, color = signature),
    linewidth = 1.2
  ) +
  geom_point(
    aes(x = HR, y = y, fill = signature),
    shape = 21,
    size = 5.2,
    color = "black",
    stroke = 0.8
  ) +
  geom_text(
    aes(x = x_hr_col, y = y, label = hr_ci_label),
    hjust = 0,
    size = 4.3,
    fontface = "bold"
  ) +
  geom_text(
    aes(x = x_p_col, y = y, label = p_label),
    hjust = 0,
    size = 4.3,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = x_hr_col,
    y = 3.60,
    label = "HR (95% CI)",
    hjust = 0,
    size = 4.8,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = x_p_col,
    y = 3.60,
    label = "P",
    hjust = 0,
    size = 4.8,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Proliferation-core" = "#55A868",
      "EMT-like" = "#4C97D8",
      "Malignant-squamous" = "#E76F6A"
    )
  ) +
  scale_color_manual(
    values = c(
      "Proliferation-core" = "#55A868",
      "EMT-like" = "#4C97D8",
      "Malignant-squamous" = "#E76F6A"
    )
  ) +
  scale_y_continuous(
    breaks = c(1, 2, 3),
    labels = c("Malignant-squamous", "EMT-like", "Proliferation-core"),
    limits = c(0.5, 3.7),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(0, x_plot_max),
    breaks = c(0, 1, 2, 4, 6, 8, 10)
  ) +
  labs(
    title = "GSE44001 grouped Cox summary",
    x = "Hazard ratio (High vs Low)",
    y = NULL
  ) +
  theme_pub(base_size = 15) +
  theme(
    legend.position = "none"
  )

print(p)

# -----------------------------
# 7. save
# -----------------------------
ggsave(
  filename = file.path(out_dir, "Figure2E_GSE44001_cox_summary_forest_refined.pdf"),
  plot = p,
  width = 11.6,
  height = 6.2,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  filename = file.path(out_dir, "Figure2E_GSE44001_cox_summary_forest_refined.png"),
  plot = p,
  width = 11.6,
  height = 6.2,
  units = "in",
  dpi = 500
)

cat("Done.\n")
cat("Input used:\n", group_file, "\n")
cat("Outputs saved in:\n", out_dir, "\n")