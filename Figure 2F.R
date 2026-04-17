rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(survival)
})

# ============================================================
# Figure 2F
# Refined cross-cohort role-assignment summary plot
# Cleaner / easier to read / no confusing bubble legends
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
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
      axis.line = element_blank(),
      plot.margin = margin(18, 18, 18, 18)
    )
}

# -----------------------------
# 1. helper functions
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

  NA_character_
}

extract_stage_p <- function(df, score_col) {
  tmp <- df %>%
    filter(!is.na(stage), !is.na(.data[[score_col]]))
  kruskal.test(tmp[[score_col]] ~ tmp$stage)$p.value
}

extract_grouped_cox <- function(df, group_col) {
  tmp <- df %>%
    filter(!is.na(.data[[group_col]]),
           !is.na(dfs_months),
           !is.na(dfs_status)) %>%
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
    HR = as.numeric(ci_tab$`exp(coef)`[idx]),
    lower = as.numeric(ci_tab$`lower .95`[idx]),
    upper = as.numeric(ci_tab$`upper .95`[idx]),
    p = as.numeric(coef_tab$`Pr(>|z|)`[idx]),
    stringsAsFactors = FALSE
  )
}

parse_gse63514_progression_p <- function(pdf_file) {
  # automatic parse from combined PDF if possible
  if (requireNamespace("pdftools", quietly = TRUE) && file.exists(pdf_file)) {
    txt <- paste(pdftools::pdf_text(pdf_file), collapse = " ")
    txt <- gsub("\n", " ", txt)
    txt <- gsub("−", "-", txt)

    p1 <- str_match(
      txt,
      "malignant_squamous_score\\s*KW p =\\s*([0-9eE\\.-]+)"
    )[, 2]
    p2 <- str_match(
      txt,
      "EMT_like_score\\s*KW p =\\s*([0-9eE\\.-]+)"
    )[, 2]
    p3 <- str_match(
      txt,
      "proliferation_core_score\\s*KW p =\\s*([0-9eE\\.-]+)"
    )[, 2]

    if (!any(is.na(c(p1, p2, p3)))) {
      return(data.frame(
        signature = c("Malignant-squamous", "EMT-like", "Proliferation-core"),
        progression_p = as.numeric(c(p1, p2, p3)),
        stringsAsFactors = FALSE
      ))
    }
  }

  # fallback to stable figure values if pdf parsing fails
  message("PDF parsing failed or pdftools unavailable. Using stable fallback values.")
  data.frame(
    signature = c("Malignant-squamous", "EMT-like", "Proliferation-core"),
    progression_p = c(5.1e-05, 0.0062, 1.4e-10),
    stringsAsFactors = FALSE
  )
}

status_progression <- function(p) {
  if (is.na(p)) return("limited")
  if (p < 1e-4) return("strong")
  if (p < 0.05) return("present")
  if (p < 0.1) return("trend")
  "limited"
}

status_stage <- function(p) {
  if (is.na(p)) return("limited")
  if (p < 0.05) return("present")
  if (p < 0.1) return("trend")
  "limited"
}

status_dfs <- function(p, hr) {
  if (is.na(p) || is.na(hr)) return("limited")
  if (p < 0.05 && hr > 1) return("risk")
  if (p < 0.1 && hr > 1) return("trend")
  "limited"
}

assign_role <- function(progression_status, dfs_status) {
  if (dfs_status == "risk") {
    return("Transferable\nrisk axis")
  }
  if (progression_status %in% c("strong", "present")) {
    return("Progression-\nlinked")
  }
  "Limited\nsupport"
}

# -----------------------------
# 2. paths
# -----------------------------
project_root <- "/Users/donglinlu/Desktop/CC1"
figure_dir   <- file.path(project_root, "7.progression validation")
out_dir      <- file.path(figure_dir, "Figure2F_refined")
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
  stop("Cannot find GSE44001 clinical_with_scores_and_groups.csv automatically.")
}

candidate_gse63514_pdf <- c(
  file.path(project_root, "7.progression validation", "GSE63514", "GSE63514_progression_validation_combined.pdf"),
  file.path(project_root, "7.progression validation", "GSE63514_progression_validation_combined.pdf")
)

gse63514_pdf <- find_preferred_file(
  project_root = project_root,
  candidate_paths = candidate_gse63514_pdf,
  pattern = "GSE63514_progression_validation_combined\\.pdf$",
  prefer_keywords = c("gse63514", "progression", "combined")
)

if (is.na(gse63514_pdf)) {
  message("Combined GSE63514 PDF not found. Fallback progression values will be used.")
}

message("Using GSE44001 input: ", group_file)
message("Using GSE63514 input: ", gse63514_pdf)

# -----------------------------
# 3. read GSE44001 data
# -----------------------------
ana <- read.csv(group_file, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c(
  "stage", "dfs_months", "dfs_status",
  "malignant_squamous_score", "EMT_like_score", "proliferation_core_score",
  "malignant_squamous_group", "EMT_like_group", "proliferation_core_group"
)

missing_cols <- required_cols[!required_cols %in% colnames(ana)]
if (length(missing_cols) > 0) {
  stop("Missing required columns in input file: ", paste(missing_cols, collapse = ", "))
}

ana <- ana %>%
  mutate(
    stage = factor(stage, levels = c("IA2", "IB1", "IB2", "IIA")),
    dfs_months = as.numeric(dfs_months),
    dfs_status = as.numeric(dfs_status),
    malignant_squamous_score = as.numeric(malignant_squamous_score),
    EMT_like_score = as.numeric(EMT_like_score),
    proliferation_core_score = as.numeric(proliferation_core_score),
    malignant_squamous_group = factor(malignant_squamous_group, levels = c("Low risk", "High risk")),
    EMT_like_group = factor(EMT_like_group, levels = c("Low risk", "High risk")),
    proliferation_core_group = factor(proliferation_core_group, levels = c("Low risk", "High risk"))
  )

# -----------------------------
# 4. collect cross-cohort evidence
# -----------------------------
progression_df <- parse_gse63514_progression_p(gse63514_pdf)

stage_df <- bind_rows(
  data.frame(signature = "Malignant-squamous",
             stage_p = extract_stage_p(ana, "malignant_squamous_score")),
  data.frame(signature = "EMT-like",
             stage_p = extract_stage_p(ana, "EMT_like_score")),
  data.frame(signature = "Proliferation-core",
             stage_p = extract_stage_p(ana, "proliferation_core_score"))
)

dfs_df <- bind_rows(
  cbind(
    data.frame(signature = "Malignant-squamous"),
    extract_grouped_cox(ana, "malignant_squamous_group")
  ),
  cbind(
    data.frame(signature = "EMT-like"),
    extract_grouped_cox(ana, "EMT_like_group")
  ),
  cbind(
    data.frame(signature = "Proliferation-core"),
    extract_grouped_cox(ana, "proliferation_core_group")
  )
)

summary_df <- progression_df %>%
  left_join(stage_df, by = "signature") %>%
  left_join(dfs_df, by = "signature") %>%
  mutate(
    progression_status = vapply(progression_p, status_progression, character(1)),
    stage_status = vapply(stage_p, status_stage, character(1)),
    dfs_status_label = mapply(status_dfs, p, HR),
    role_label = mapply(assign_role, progression_status, dfs_status_label)
  )

write.csv(
  summary_df,
  file.path(out_dir, "Figure2F_cross_cohort_role_assignment_table.csv"),
  row.names = FALSE
)

# -----------------------------
# 5. build plot table
# -----------------------------
plot_df <- bind_rows(
  summary_df %>%
    transmute(
      Program = signature,
      Column = "Progression\nGSE63514",
      Label = progression_status,
      FillClass = progression_status
    ),
  summary_df %>%
    transmute(
      Program = signature,
      Column = "Stage\nGSE44001",
      Label = stage_status,
      FillClass = stage_status
    ),
  summary_df %>%
    transmute(
      Program = signature,
      Column = "DFS risk\nGSE44001",
      Label = dfs_status_label,
      FillClass = dfs_status_label
    ),
  summary_df %>%
    transmute(
      Program = signature,
      Column = "Assigned role",
      Label = role_label,
      FillClass = case_when(
        role_label == "Transferable\nrisk axis" ~ "role_risk",
        role_label == "Progression-\nlinked" ~ "role_progression",
        TRUE ~ "role_limited"
      )
    )
)

plot_df$Program <- factor(
  plot_df$Program,
  levels = c("Malignant-squamous", "EMT-like", "Proliferation-core")
)

plot_df$Column <- factor(
  plot_df$Column,
  levels = c(
    "Progression\nGSE63514",
    "Stage\nGSE44001",
    "DFS risk\nGSE44001",
    "Assigned role"
  )
)

fill_colors <- c(
  "strong" = "#D73027",
  "present" = "#FC8D59",
  "trend" = "#FEE08B",
  "limited" = "#D9D9D9",
  "risk" = "#1A9850",
  "role_risk" = "#55A868",
  "role_progression" = "#A6CEE3",
  "role_limited" = "#E5E5E5"
)

# -----------------------------
# 6. plot
# -----------------------------
p <- ggplot(plot_df, aes(x = Column, y = Program, fill = FillClass)) +
  geom_tile(
    color = "black",
    linewidth = 0.8,
    width = 0.96,
    height = 0.88
  ) +
  geom_text(
    aes(label = Label),
    size = 5,
    fontface = "bold",
    lineheight = 0.95
  ) +
  scale_fill_manual(values = fill_colors, guide = "none") +
  labs(
    title = "Cross-cohort role assignment of prioritized epithelial programs",
    x = NULL,
    y = NULL
  ) +
  theme_pub(base_size = 15) +
  theme(
    axis.text.x = element_text(face = "bold", size = 13, lineheight = 0.95),
    axis.text.y = element_text(face = "bold", size = 14),
    panel.grid = element_blank()
  )

print(p)

# -----------------------------
# 7. save
# -----------------------------
ggsave(
  filename = file.path(out_dir, "Figure2F_cross_cohort_role_assignment_refined.pdf"),
  plot = p,
  width = 10.6,
  height = 5.8,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  filename = file.path(out_dir, "Figure2F_cross_cohort_role_assignment_refined.png"),
  plot = p,
  width = 10.6,
  height = 5.8,
  units = "in",
  dpi = 500
)

cat("Done.\n")
cat("GSE44001 input used:\n", group_file, "\n")
cat("GSE63514 input used:\n", gse63514_pdf, "\n")
cat("Outputs saved in:\n", out_dir, "\n")