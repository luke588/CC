rm(list = ls())
graphics.off()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(ggpubr)
  library(scales)
  library(broom)
})

# =========================================================
# TCGA-CESC validation of three prioritized epithelial programs
# Fixed-input, unified publication-style version
# Input/output root:
# /Users/donglinlu/Desktop/CC1/7.progression validation/TCGA
# =========================================================

# =========================
# USER SETTINGS
# =========================
base_dir <- "/Users/donglinlu/Desktop/CC1/7.progression validation/TCGA"

expr_file <- file.path(base_dir, "TCGA_CESC_geneMatrix_clean.txt")
clin_file <- file.path(base_dir, "TCGA_CESC_clinical_clean.txt")

sig9_file  <- file.path(base_dir, "cluster9_malignant_squamous_epithelial_top50.txt")
sig10_file <- file.path(base_dir, "cluster10_EMT_like_epithelial_top50.txt")
sig6_file  <- file.path(base_dir, "cluster6_proliferation_core_epithelial_top50.txt")

out_dir  <- file.path(base_dir, "results")
plot_dir <- file.path(base_dir, "plots")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# HELPERS
# =========================
standardize_tcga_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\.", "-", x)
  x <- gsub("\\s+", "", x)
  x
}

patient_id_12 <- function(x) {
  x <- standardize_tcga_id(x)
  substr(x, 1, 12)
}

read_signature <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[x != ""]
  unique(x)
}

make_unique_gene <- function(df, gene_col = "geneNames") {
  df <- df %>%
    dplyr::filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "")

  df %>%
    dplyr::group_by(.data[[gene_col]]) %>%
    dplyr::summarise(
      dplyr::across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

zscore_rows <- function(mat) {
  t(scale(t(mat)))
}

calc_signature_score_zmean <- function(expr_mat, genes) {
  genes_use <- intersect(genes, rownames(expr_mat))
  if (length(genes_use) == 0) {
    return(list(score = rep(NA_real_, ncol(expr_mat)), genes_used = character(0)))
  }

  submat <- expr_mat[genes_use, , drop = FALSE]
  zmat <- zscore_rows(submat)
  zmat[is.na(zmat)] <- 0
  score <- colMeans(zmat)

  list(score = score, genes_used = genes_use)
}

clean_stage_tcga <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NOT REPORTED", "NOT AVAILABLE", "[NOT AVAILABLE]")] <- NA
  x <- gsub("^STAGE\\s*", "", x)
  x <- gsub("\\s+", "", x)

  stage_levels <- c(
    "I", "IA", "IA1", "IA2", "IB", "IB1", "IB2",
    "II", "IIA", "IIB",
    "III", "IIIA", "IIIB", "IIIC",
    "IV", "IVA", "IVB"
  )

  factor(x, levels = stage_levels)
}

clean_grade_tcga <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NOT REPORTED", "NOT AVAILABLE", "[NOT AVAILABLE]")] <- NA
  x <- gsub("\\s+", "", x)
  factor(x, levels = c("G1", "G2", "G3", "G4", "GX"))
}

safe_stage_numeric <- function(x) {
  x <- as.character(x)
  x <- gsub("[ABC]$", "", x)
  map <- c(
    "I" = 1, "IA" = 1, "IA1" = 1, "IA2" = 1, "IB" = 1, "IB1" = 1, "IB2" = 1,
    "II" = 2, "IIA" = 2, "IIB" = 2,
    "III" = 3, "IIIA" = 3, "IIIB" = 3, "IIIC" = 3,
    "IV" = 4, "IVA" = 4, "IVB" = 4
  )
  unname(map[x])
}

cox_extract <- function(fit, var_name) {
  sm <- summary(fit)

  coef_tab <- as.data.frame(sm$coefficients)
  ci_tab   <- as.data.frame(sm$conf.int)

  coef_tab$term <- rownames(coef_tab)
  ci_tab$term   <- rownames(ci_tab)

  out <- dplyr::left_join(coef_tab, ci_tab, by = "term")
  colnames(out) <- make.names(colnames(out), unique = TRUE)
  out$variable <- var_name
  out
}

make_group <- function(x) {
  grp <- ifelse(x >= median(x, na.rm = TRUE), "High risk", "Low risk")
  factor(grp, levels = c("Low risk", "High risk"))
}

get_pairwise_comparisons <- function(groups) {
  lev <- unique(as.character(groups))
  lev <- lev[!is.na(lev)]
  if (length(lev) < 2) return(list())
  combn(lev, 2, simplify = FALSE)
}

pairwise_wilcox_one <- function(df, score_col, group_col) {
  df2 <- df %>%
    dplyr::filter(!is.na(.data[[group_col]]), !is.na(.data[[score_col]]))

  if (length(unique(df2[[group_col]])) < 2) {
    return(data.frame())
  }

  pw <- pairwise.wilcox.test(df2[[score_col]], df2[[group_col]], p.adjust.method = "BH")
  mat <- pw$p.value
  if (is.null(mat)) return(data.frame())

  out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(out) <- c("group1", "group2", "p.adj")
  out$signature <- score_col
  out$group_var <- group_col

  out %>% dplyr::filter(!is.na(p.adj))
}

format_p_value <- function(p, digits = 5) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) {
    return(format(p, scientific = TRUE, digits = 3))
  } else {
    out <- formatC(p, format = "f", digits = digits)
    out <- sub("0+$", "", out)
    out <- sub("\\.$", "", out)
    return(out)
  }
}

theme_pub <- function(base_size = 18) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = base_size + 4,
        margin = margin(b = 6),
        colour = "black"
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        face = "plain",
        size = base_size,
        margin = margin(b = 12),
        colour = "black"
      ),
      plot.title.position = "plot",
      axis.title = element_text(
        size = base_size + 1,
        face = "bold",
        colour = "black"
      ),
      axis.text = element_text(
        size = base_size - 1,
        colour = "black"
      ),
      axis.text.x = element_text(
        size = base_size - 1,
        face = "bold",
        colour = "black",
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(
        size = base_size - 1,
        face = "bold",
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
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.8
      ),
      axis.line = element_blank(),
      plot.margin = margin(t = 14, r = 16, b = 14, l = 16)
    )
}

theme_risk <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        hjust = 0,
        face = "bold",
        size = base_size + 1,
        colour = "black"
      ),
      axis.title = element_text(
        size = base_size,
        face = "bold",
        colour = "black"
      ),
      axis.text = element_text(
        size = base_size - 1,
        colour = "black"
      ),
      panel.border = element_blank(),
      axis.line = element_blank()
    )
}

save_plot_both <- function(plot_obj, pdf_file, png_file, width, height) {
  ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE
  )

  ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = 500
  )
}

make_assoc_boxplot <- function(df, score_col, group_col, title_txt, ylab_txt, out_prefix) {
  df2 <- df %>%
    dplyr::filter(!is.na(.data[[group_col]]), !is.na(.data[[score_col]]))

  if (nrow(df2) == 0 || length(unique(df2[[group_col]])) < 2) {
    return(NULL)
  }

  kw_p <- kruskal.test(df2[[score_col]] ~ df2[[group_col]])$p.value
  subtitle_txt <- paste0("Kruskal-Wallis p = ", format_p_value(kw_p, digits = 5))

  comps <- get_pairwise_comparisons(df2[[group_col]])

  levs <- unique(as.character(df2[[group_col]]))
  levs <- levs[!is.na(levs)]
  pal <- setNames(scales::hue_pal()(length(levs)), levs)

  ymax <- max(df2[[score_col]], na.rm = TRUE)
  ymin <- min(df2[[score_col]], na.rm = TRUE)
  yrng <- ymax - ymin
  if (!is.finite(yrng) || yrng == 0) yrng <- 1

  p <- ggplot(df2, aes(x = .data[[group_col]], y = .data[[score_col]], fill = .data[[group_col]])) +
    geom_boxplot(
      width = 0.68,
      color = "black",
      outlier.shape = 1,
      linewidth = 0.8,
      alpha = 0.9
    ) +
    geom_jitter(
      width = 0.12,
      size = 1.5,
      alpha = 0.55,
      color = "black"
    ) +
    scale_fill_manual(values = pal) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = group_col,
      y = ylab_txt
    ) +
    coord_cartesian(ylim = c(ymin, ymax + 0.35 * yrng)) +
    theme_pub(base_size = 18) +
    theme(
      legend.position = "none"
    )

  if (length(comps) > 0) {
    p <- p +
      stat_compare_means(
        comparisons = comps,
        method = "wilcox.test",
        label = "p.signif",
        step.increase = 0.08,
        size = 5
      )
  }

  save_plot_both(
    plot_obj = p,
    pdf_file = paste0(out_prefix, ".pdf"),
    png_file = paste0(out_prefix, ".png"),
    width = 10,
    height = 8
  )

  invisible(p)
}

plot_km_pub <- function(df, time_col, status_col, group_col, title_text, out_prefix) {
  surv_formula <- as.formula(
    paste0("Surv(", time_col, ", ", status_col, ") ~ ", group_col)
  )

  fit <- survfit(surv_formula, data = df)
  cox_fit <- coxph(surv_formula, data = df)

  sm <- summary(cox_fit)
  coef_df <- as.data.frame(sm$coefficients)
  ci_df   <- as.data.frame(sm$conf.int)

  coef_df$term <- rownames(coef_df)
  ci_df$term   <- rownames(ci_df)

  idx <- grep("High risk", coef_df$term)
  if (length(idx) == 0) idx <- 1 else idx <- idx[1]

  hr    <- ci_df$`exp(coef)`[idx]
  lower <- ci_df$`lower .95`[idx]
  upper <- ci_df$`upper .95`[idx]
  pval  <- coef_df$`Pr(>|z|)`[idx]

  ann_text <- paste0(
    "p = ", formatC(pval, format = "f", digits = 4), "\n",
    "HR = ", formatC(hr, format = "f", digits = 3), "\n",
    "95% CI: ", formatC(lower, format = "f", digits = 3),
    "–", formatC(upper, format = "f", digits = 3)
  )

  surv_df <- survminer::surv_summary(fit, data = df)
  surv_df$strata <- gsub(paste0(group_col, "="), "", surv_df$strata)
  surv_df$strata <- factor(surv_df$strata, levels = c("Low risk", "High risk"))
  censor_df <- surv_df %>% dplyr::filter(n.censor > 0)

  xmax <- max(df[[time_col]], na.rm = TRUE)
  xmax_break <- ceiling(xmax)

  curve_plot <- ggplot(surv_df, aes(x = time, y = surv, color = strata)) +
    geom_step(linewidth = 1.35) +
    geom_point(
      data = censor_df,
      aes(x = time, y = surv, color = strata),
      shape = 3,
      size = 3
    ) +
    scale_color_manual(
      values = c("Low risk" = "#2C7BE5", "High risk" = "#E74C3C"),
      name = "Risk"
    ) +
    coord_cartesian(xlim = c(0, xmax_break), ylim = c(0, 1.05)) +
    scale_x_continuous(breaks = pretty(c(0, xmax_break), n = 6)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = title_text,
      x = "OS time (years)",
      y = "Overall survival probability"
    ) +
    annotate(
      "text",
      x = xmax_break * 0.45, y = 0.18,
      label = ann_text,
      size = 6.2,
      fontface = "bold"
    ) +
    theme_pub(base_size = 18) +
    theme(
      legend.position = "top"
    )

  time_breaks <- pretty(c(0, xmax_break), n = 6)
  sf_sum <- summary(fit, times = time_breaks, extend = TRUE)

  risk_df <- data.frame(
    time = sf_sum$time,
    strata = gsub(".*=", "", sf_sum$strata),
    n.risk = sf_sum$n.risk,
    stringsAsFactors = FALSE
  )
  risk_df$strata <- factor(risk_df$strata, levels = c("Low risk", "High risk"))

  risk_plot <- ggplot(risk_df, aes(x = time, y = strata, label = n.risk, color = strata)) +
    geom_text(size = 6) +
    scale_color_manual(
      values = c("Low risk" = "#2C7BE5", "High risk" = "#E74C3C"),
      guide = "none"
    ) +
    scale_x_continuous(breaks = time_breaks, limits = c(min(time_breaks), max(time_breaks))) +
    labs(
      title = "Number at risk",
      x = "OS time (years)",
      y = "Risk"
    ) +
    theme_risk(base_size = 16)

  final_plot <- ggpubr::ggarrange(
    curve_plot,
    risk_plot,
    ncol = 1,
    heights = c(3.3, 1.2),
    align = "v"
  )

  save_plot_both(
    plot_obj = final_plot,
    pdf_file = paste0(out_prefix, ".pdf"),
    png_file = paste0(out_prefix, ".png"),
    width = 11,
    height = 9
  )

  invisible(list(fit = fit, cox_fit = cox_fit))
}

# =========================
# FILE CHECK
# =========================
cat("Expression file:\n", expr_file, "\n")
cat("Clinical file:\n", clin_file, "\n")
cat("Signature 9 file:\n", sig9_file, "\n")
cat("Signature 10 file:\n", sig10_file, "\n")
cat("Signature 6 file:\n", sig6_file, "\n\n")

need_files <- c(expr_file, clin_file, sig9_file, sig10_file, sig6_file)
missing_files <- need_files[!file.exists(need_files)]

if (length(missing_files) > 0) {
  stop("Missing file(s):\n", paste(missing_files, collapse = "\n"))
}

# =========================
# READ EXPRESSION
# =========================
expr <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!("geneNames" %in% colnames(expr))) {
  stop("Expression file must contain first column named 'geneNames'.")
}

colnames(expr)[1] <- "geneNames"
expr$geneNames <- trimws(expr$geneNames)

expr_num <- expr
for (j in 2:ncol(expr_num)) {
  expr_num[[j]] <- suppressWarnings(as.numeric(expr_num[[j]]))
}

expr_num <- make_unique_gene(expr_num, gene_col = "geneNames")

expr_mat <- as.matrix(expr_num[, -1])
rownames(expr_mat) <- expr_num$geneNames
mode(expr_mat) <- "numeric"

expr_ids <- standardize_tcga_id(colnames(expr_mat))
colnames(expr_mat) <- expr_ids

# =========================
# READ CLINICAL
# =========================
clin <- read.delim(
  clin_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

req_cols <- c("sample", "OS.time", "OS")
miss_cols <- req_cols[!req_cols %in% colnames(clin)]
if (length(miss_cols) > 0) {
  stop("Clinical file missing required columns: ", paste(miss_cols, collapse = ", "))
}

clin <- clin %>%
  dplyr::mutate(
    sample = standardize_tcga_id(sample),
    patient_id = patient_id_12(sample),
    OS.time = as.numeric(OS.time),
    OS = as.numeric(OS),
    OS.years = OS.time / 365,
    Age = if ("Age" %in% colnames(.)) as.numeric(Age) else NA_real_,
    Gender = if ("Gender" %in% colnames(.)) trimws(as.character(Gender)) else NA_character_,
    Grade = if ("Grade" %in% colnames(.)) clean_grade_tcga(Grade) else factor(NA),
    stage = if ("stage" %in% colnames(.)) clean_stage_tcga(stage) else factor(NA),
    T_stage = if ("T_stage" %in% colnames(.)) trimws(as.character(T_stage)) else NA_character_,
    M_stage = if ("M_stage" %in% colnames(.)) trimws(as.character(M_stage)) else NA_character_,
    N_stage = if ("N_stage" %in% colnames(.)) trimws(as.character(N_stage)) else NA_character_,
    stage_num = safe_stage_numeric(as.character(stage))
  )

# =========================
# MATCH SAMPLE ORDER
# =========================
common_samples <- intersect(colnames(expr_mat), clin$sample)

matching_mode <- "exact_sample_id"

if (length(common_samples) == 0) {
  expr_patient <- patient_id_12(colnames(expr_mat))
  keep_first <- !duplicated(expr_patient)
  expr_mat <- expr_mat[, keep_first, drop = FALSE]
  colnames(expr_mat) <- expr_patient[keep_first]

  clin <- clin %>%
    dplyr::arrange(sample) %>%
    dplyr::distinct(patient_id, .keep_all = TRUE)

  common_samples <- intersect(colnames(expr_mat), clin$patient_id)
  matching_mode <- "patient_id_12"
}

if (length(common_samples) == 0) {
  stop("No matched samples between expression and clinical data.")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]

if (matching_mode == "exact_sample_id") {
  clin <- clin %>% dplyr::filter(sample %in% common_samples)
  clin <- clin[match(common_samples, clin$sample), ]
  stopifnot(all(clin$sample == colnames(expr_mat)))
} else {
  clin <- clin %>% dplyr::filter(patient_id %in% common_samples)
  clin <- clin[match(common_samples, clin$patient_id), ]
  stopifnot(all(clin$patient_id == colnames(expr_mat)))
}

# =========================
# READ SIGNATURES
# =========================
sig9  <- read_signature(sig9_file)
sig10 <- read_signature(sig10_file)
sig6  <- read_signature(sig6_file)

sig9_res  <- calc_signature_score_zmean(expr_mat, sig9)
sig10_res <- calc_signature_score_zmean(expr_mat, sig10)
sig6_res  <- calc_signature_score_zmean(expr_mat, sig6)

score_df <- data.frame(
  sample = colnames(expr_mat),
  malignant_squamous_score = as.numeric(sig9_res$score),
  EMT_like_score = as.numeric(sig10_res$score),
  proliferation_core_score = as.numeric(sig6_res$score),
  stringsAsFactors = FALSE
)

write.csv(
  data.frame(
    signature = "cluster9_malignant_squamous_epithelial",
    genes_input = length(sig9),
    genes_used = length(sig9_res$genes_used),
    genes_used_list = paste(sig9_res$genes_used, collapse = ";"),
    signature_file = sig9_file
  ),
  file.path(out_dir, "01.signature_cluster9_usage.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    signature = "cluster10_EMT_like_epithelial",
    genes_input = length(sig10),
    genes_used = length(sig10_res$genes_used),
    genes_used_list = paste(sig10_res$genes_used, collapse = ";"),
    signature_file = sig10_file
  ),
  file.path(out_dir, "02.signature_cluster10_usage.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    signature = "cluster6_proliferation_core_epithelial",
    genes_input = length(sig6),
    genes_used = length(sig6_res$genes_used),
    genes_used_list = paste(sig6_res$genes_used, collapse = ";"),
    signature_file = sig6_file
  ),
  file.path(out_dir, "03.signature_cluster6_usage.csv"),
  row.names = FALSE
)

# =========================
# MERGE CLINICAL + SCORES
# =========================
if (matching_mode == "exact_sample_id") {
  ana <- dplyr::left_join(clin, score_df, by = "sample")
} else {
  score_df$patient_id <- score_df$sample
  ana <- dplyr::left_join(clin, score_df %>% dplyr::select(-sample), by = "patient_id")
}

write.csv(
  ana,
  file.path(out_dir, "04.TCGA_CESC_clinical_with_scores.csv"),
  row.names = FALSE
)

# =========================
# BASIC SUMMARIES
# =========================
score_summary <- ana %>%
  dplyr::summarise(
    n = dplyr::n(),
    malignant_squamous_mean = mean(malignant_squamous_score, na.rm = TRUE),
    EMT_like_mean = mean(EMT_like_score, na.rm = TRUE),
    proliferation_core_mean = mean(proliferation_core_score, na.rm = TRUE)
  )

write.csv(
  score_summary,
  file.path(out_dir, "05.score_summary.csv"),
  row.names = FALSE
)

# =========================
# GROUPING
# =========================
ana <- ana %>%
  dplyr::mutate(
    malignant_squamous_group = make_group(malignant_squamous_score),
    EMT_like_group = make_group(EMT_like_score),
    proliferation_core_group = make_group(proliferation_core_score)
  )

write.csv(
  ana,
  file.path(out_dir, "06.TCGA_CESC_clinical_with_scores_and_groups.csv"),
  row.names = FALSE
)

# =========================
# KM ANALYSIS
# =========================
plot_km_pub(
  df = ana,
  time_col = "OS.years",
  status_col = "OS",
  group_col = "malignant_squamous_group",
  title_text = "TCGA-CESC: OS by malignant-squamous score",
  out_prefix = file.path(plot_dir, "PUB_TCGA_KM_malignant_squamous_score")
)

plot_km_pub(
  df = ana,
  time_col = "OS.years",
  status_col = "OS",
  group_col = "EMT_like_group",
  title_text = "TCGA-CESC: OS by EMT-like score",
  out_prefix = file.path(plot_dir, "PUB_TCGA_KM_EMT_like_score")
)

plot_km_pub(
  df = ana,
  time_col = "OS.years",
  status_col = "OS",
  group_col = "proliferation_core_group",
  title_text = "TCGA-CESC: OS by proliferation-core score",
  out_prefix = file.path(plot_dir, "PUB_TCGA_KM_proliferation_core_score")
)

logrank_df <- dplyr::bind_rows(
  broom::tidy(survdiff(Surv(OS.years, OS) ~ malignant_squamous_group, data = ana)) %>% dplyr::mutate(signature = "malignant_squamous_score"),
  broom::tidy(survdiff(Surv(OS.years, OS) ~ EMT_like_group, data = ana)) %>% dplyr::mutate(signature = "EMT_like_score"),
  broom::tidy(survdiff(Surv(OS.years, OS) ~ proliferation_core_group, data = ana)) %>% dplyr::mutate(signature = "proliferation_core_score")
)

write.csv(
  logrank_df,
  file.path(out_dir, "07.logrank_results_OS.csv"),
  row.names = FALSE
)

# =========================
# COX ANALYSIS
# =========================
cox_uni_1 <- coxph(Surv(OS.time, OS) ~ malignant_squamous_score, data = ana)
cox_uni_2 <- coxph(Surv(OS.time, OS) ~ EMT_like_score, data = ana)
cox_uni_3 <- coxph(Surv(OS.time, OS) ~ proliferation_core_score, data = ana)

cox_uni_tab <- dplyr::bind_rows(
  cox_extract(cox_uni_1, "malignant_squamous_score"),
  cox_extract(cox_uni_2, "EMT_like_score"),
  cox_extract(cox_uni_3, "proliferation_core_score")
)

write.csv(
  cox_uni_tab,
  file.path(out_dir, "08.univariate_cox_results_OS.csv"),
  row.names = FALSE
)

cox_group_1 <- coxph(Surv(OS.time, OS) ~ malignant_squamous_group, data = ana)
cox_group_2 <- coxph(Surv(OS.time, OS) ~ EMT_like_group, data = ana)
cox_group_3 <- coxph(Surv(OS.time, OS) ~ proliferation_core_group, data = ana)

cox_group_tab <- dplyr::bind_rows(
  cox_extract(cox_group_1, "malignant_squamous_group"),
  cox_extract(cox_group_2, "EMT_like_group"),
  cox_extract(cox_group_3, "proliferation_core_group")
)

write.csv(
  cox_group_tab,
  file.path(out_dir, "09.grouped_cox_results_OS.csv"),
  row.names = FALSE
)

multivar_models <- list()

if (sum(!is.na(ana$Age)) > 10) {
  multivar_models[["model_malignant_squamous_age"]] <- coxph(Surv(OS.time, OS) ~ malignant_squamous_score + Age, data = ana)
  multivar_models[["model_EMT_like_age"]] <- coxph(Surv(OS.time, OS) ~ EMT_like_score + Age, data = ana)
  multivar_models[["model_proliferation_core_age"]] <- coxph(Surv(OS.time, OS) ~ proliferation_core_score + Age, data = ana)
}

if (sum(!is.na(ana$stage_num)) > 10) {
  multivar_models[["model_malignant_squamous_stage"]] <- coxph(Surv(OS.time, OS) ~ malignant_squamous_score + stage_num, data = ana)
  multivar_models[["model_EMT_like_stage"]] <- coxph(Surv(OS.time, OS) ~ EMT_like_score + stage_num, data = ana)
  multivar_models[["model_proliferation_core_stage"]] <- coxph(Surv(OS.time, OS) ~ proliferation_core_score + stage_num, data = ana)
}

if (length(multivar_models) > 0) {
  cox_multi_tab <- dplyr::bind_rows(
    lapply(names(multivar_models), function(nm) cox_extract(multivar_models[[nm]], nm))
  )

  write.csv(
    cox_multi_tab,
    file.path(out_dir, "10.multivariate_cox_results_OS.csv"),
    row.names = FALSE
  )
}

# =========================
# ASSOCIATION WITH STAGE / GRADE
# =========================
if (sum(!is.na(ana$stage)) > 0 && length(unique(na.omit(ana$stage))) >= 2) {
  stage_kruskal <- dplyr::bind_rows(
    broom::tidy(kruskal.test(malignant_squamous_score ~ stage, data = ana)) %>% dplyr::mutate(signature = "malignant_squamous_score"),
    broom::tidy(kruskal.test(EMT_like_score ~ stage, data = ana)) %>% dplyr::mutate(signature = "EMT_like_score"),
    broom::tidy(kruskal.test(proliferation_core_score ~ stage, data = ana)) %>% dplyr::mutate(signature = "proliferation_core_score")
  )

  write.csv(
    stage_kruskal,
    file.path(out_dir, "11.score_stage_kruskal.csv"),
    row.names = FALSE
  )

  stage_pairwise <- dplyr::bind_rows(
    pairwise_wilcox_one(ana, "malignant_squamous_score", "stage"),
    pairwise_wilcox_one(ana, "EMT_like_score", "stage"),
    pairwise_wilcox_one(ana, "proliferation_core_score", "stage")
  )

  write.csv(
    stage_pairwise,
    file.path(out_dir, "12.score_stage_pairwise_wilcox.csv"),
    row.names = FALSE
  )

  make_assoc_boxplot(
    ana,
    "malignant_squamous_score",
    "stage",
    "TCGA-CESC: malignant-squamous score by stage",
    "Signature score",
    file.path(plot_dir, "PUB_TCGA_stage_malignant_squamous_score")
  )

  make_assoc_boxplot(
    ana,
    "EMT_like_score",
    "stage",
    "TCGA-CESC: EMT-like score by stage",
    "Signature score",
    file.path(plot_dir, "PUB_TCGA_stage_EMT_like_score")
  )

  make_assoc_boxplot(
    ana,
    "proliferation_core_score",
    "stage",
    "TCGA-CESC: proliferation-core score by stage",
    "Signature score",
    file.path(plot_dir, "PUB_TCGA_stage_proliferation_core_score")
  )
}

if (sum(!is.na(ana$Grade)) > 0 && length(unique(na.omit(ana$Grade))) >= 2) {
  grade_kruskal <- dplyr::bind_rows(
    broom::tidy(kruskal.test(malignant_squamous_score ~ Grade, data = ana)) %>% dplyr::mutate(signature = "malignant_squamous_score"),
    broom::tidy(kruskal.test(EMT_like_score ~ Grade, data = ana)) %>% dplyr::mutate(signature = "EMT_like_score"),
    broom::tidy(kruskal.test(proliferation_core_score ~ Grade, data = ana)) %>% dplyr::mutate(signature = "proliferation_core_score")
  )

  write.csv(
    grade_kruskal,
    file.path(out_dir, "13.score_grade_kruskal.csv"),
    row.names = FALSE
  )

  grade_pairwise <- dplyr::bind_rows(
    pairwise_wilcox_one(ana, "malignant_squamous_score", "Grade"),
    pairwise_wilcox_one(ana, "EMT_like_score", "Grade"),
    pairwise_wilcox_one(ana, "proliferation_core_score", "Grade")
  )

  write.csv(
    grade_pairwise,
    file.path(out_dir, "14.score_grade_pairwise_wilcox.csv"),
    row.names = FALSE
  )

  make_assoc_boxplot(
    ana,
    "malignant_squamous_score",
    "Grade",
    "TCGA-CESC: malignant-squamous score by grade",
    "Signature score",
    file.path(plot_dir, "PUB_TCGA_grade_malignant_squamous_score")
  )

  make_assoc_boxplot(
    ana,
    "EMT_like_score",
    "Grade",
    "TCGA-CESC: EMT-like score by grade",
    "Signature score",
    file.path(plot_dir, "PUB_TCGA_grade_EMT_like_score")
  )

  make_assoc_boxplot(
    ana,
    "proliferation_core_score",
    "Grade",
    "TCGA-CESC: proliferation-core score by grade",
    "Signature score",
    file.path(plot_dir, "PUB_TCGA_grade_proliferation_core_score")
  )
}

# =========================
# SCORE SUMMARY BY STAGE / GRADE
# =========================
if (sum(!is.na(ana$stage)) > 0) {
  stage_score_table <- ana %>%
    dplyr::filter(!is.na(stage)) %>%
    dplyr::group_by(stage) %>%
    dplyr::summarise(
      n = dplyr::n(),
      malignant_squamous_mean = mean(malignant_squamous_score, na.rm = TRUE),
      malignant_squamous_sd = sd(malignant_squamous_score, na.rm = TRUE),
      EMT_like_mean = mean(EMT_like_score, na.rm = TRUE),
      EMT_like_sd = sd(EMT_like_score, na.rm = TRUE),
      proliferation_core_mean = mean(proliferation_core_score, na.rm = TRUE),
      proliferation_core_sd = sd(proliferation_core_score, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(
    stage_score_table,
    file.path(out_dir, "15.stage_wise_score_summary.csv"),
    row.names = FALSE
  )
}

if (sum(!is.na(ana$Grade)) > 0) {
  grade_score_table <- ana %>%
    dplyr::filter(!is.na(Grade)) %>%
    dplyr::group_by(Grade) %>%
    dplyr::summarise(
      n = dplyr::n(),
      malignant_squamous_mean = mean(malignant_squamous_score, na.rm = TRUE),
      malignant_squamous_sd = sd(malignant_squamous_score, na.rm = TRUE),
      EMT_like_mean = mean(EMT_like_score, na.rm = TRUE),
      EMT_like_sd = sd(EMT_like_score, na.rm = TRUE),
      proliferation_core_mean = mean(proliferation_core_score, na.rm = TRUE),
      proliferation_core_sd = sd(proliferation_core_score, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(
    grade_score_table,
    file.path(out_dir, "16.grade_wise_score_summary.csv"),
    row.names = FALSE
  )
}

# =========================
# RUN INFO
# =========================
report <- c(
  paste0("Base dir: ", base_dir),
  paste0("Expression file: ", expr_file),
  paste0("Clinical file: ", clin_file),
  paste0("Signature 9 file: ", sig9_file),
  paste0("Signature 10 file: ", sig10_file),
  paste0("Signature 6 file: ", sig6_file),
  paste0("Matching mode: ", matching_mode),
  paste0("Matched samples: ", nrow(ana)),
  paste0("cluster9 genes input/used: ", length(sig9), "/", length(sig9_res$genes_used)),
  paste0("cluster10 genes input/used: ", length(sig10), "/", length(sig10_res$genes_used)),
  paste0("cluster6 genes input/used: ", length(sig6), "/", length(sig6_res$genes_used)),
  paste0("Non-missing stage count: ", sum(!is.na(ana$stage))),
  paste0("Non-missing grade count: ", sum(!is.na(ana$Grade))),
  "Endpoint: TCGA-CESC overall survival",
  "KM style: custom ggplot + blue/red groups + HR/CI/p + risk table",
  "Grade/stage boxplots: Kruskal-Wallis p shown as subtitle below title",
  "IMPORTANT: use PUB_*.pdf files as final plots"
)

writeLines(
  report,
  file.path(out_dir, "17.run_info.txt")
)

cat("Done. Results saved in:\n", base_dir, "\n\n")
cat("Main outputs:\n")
cat("results/04.TCGA_CESC_clinical_with_scores.csv\n")
cat("results/08.univariate_cox_results_OS.csv\n")
cat("results/09.grouped_cox_results_OS.csv\n")
cat("results/17.run_info.txt\n\n")
cat("Plots:\n")
cat("plots/PUB_TCGA_KM_malignant_squamous_score.pdf\n")
cat("plots/PUB_TCGA_KM_EMT_like_score.pdf\n")
cat("plots/PUB_TCGA_KM_proliferation_core_score.pdf\n")
cat("plots/PUB_TCGA_stage_*.pdf (if stage available)\n")
cat("plots/PUB_TCGA_grade_*.pdf (if grade available)\n")