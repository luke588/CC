source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# =========================================================
# Figure 2D
# Cross-cohort role-assignment bubble summary
# unified by canonical 00.manuscript_plot_theme.R
# =========================================================

# =========================
# USER SETTINGS
# =========================
base_dir <- "/Users/donglinlu/Desktop/CC1/6.bulk_validation"
out_dir  <- file.path(base_dir, "Figure2D_cross_cohort_role_assignment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(out_dir, "Figure2D_cross_cohort_role_assignment_bubble_clean.pdf")
out_png <- file.path(out_dir, "Figure2D_cross_cohort_role_assignment_bubble_clean.png")

# =========================
# DATA
# support_score:
# 1 = weak
# 2 = moderate
# 3 = strong
# canonical display order fixed as:
# Malignant squamous -> EMT-like -> Proliferation core
# =========================
plot_df <- tibble::tribble(
  ~program,                 ~domain,                       ~support_score,
  "Malignant squamous",     "Single-cell identity",        2.5,
  "Malignant squamous",     "Progression relevance",       2.2,
  "Malignant squamous",     "Clinical transferability",    1.2,

  "EMT-like",               "Single-cell identity",        2.4,
  "EMT-like",               "Progression relevance",       2.1,
  "EMT-like",               "Clinical transferability",    1.4,

  "Proliferation core",     "Single-cell identity",        2.6,
  "Proliferation core",     "Progression relevance",       2.9,
  "Proliferation core",     "Clinical transferability",    3.0
) %>%
  dplyr::mutate(
    domain = factor(
      domain,
      levels = c(
        "Single-cell identity",
        "Progression relevance",
        "Clinical transferability"
      )
    ),
    program = factor(
      program,
      levels = c(
        "Malignant squamous",
        "EMT-like",
        "Proliferation core"
      )
    )
  )

# =========================
# PLOT
# =========================
p <- plot_role_bubble_pub(
  plot_df = plot_df,
  title_txt = "Cross-cohort role assignment of prioritized epithelial programs",
  program_order = c(
    "Malignant squamous",
    "EMT-like",
    "Proliferation core"
  )
)

# =========================
# SAVE
# =========================
save_pub(
  plot_obj = p,
  pdf_file = out_pdf,
  png_file = out_png,
  width = 14.2,
  height = 8.0,
  dpi = 320
)

cat("Done.\n")
cat("Output files:\n")
cat(out_pdf, "\n")
cat(out_png, "\n")