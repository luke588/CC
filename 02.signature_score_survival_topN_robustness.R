suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(ggpubr)
  library(broom)
})
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
# =========================================================
# 02.GSE44001_signature_score_survival_topN_robustness_refined.R
# top30 / top50 / top100 robustness validation
# refined publication-style script
# 1) stage boxplots: only global Kruskal-Wallis p, no pairwise ns/*
# 2) KM plots
# 3) grouped / continuous / multivariate Cox
# 4) cleaner grouped Cox forest plot with CI end caps
# =========================================================

# =========================
# USER SETTINGS
# =========================
base_dir <- "/Users/donglinlu/Desktop/CC1/6.bulk_validation"

expr_file <- file.path(base_dir, "geneMatrix.txt")
clin_file <- file.path(base_dir, "04.GSE44001_clinical_clean_matched.csv")

sig_dir <- "/Users/donglinlu/Desktop/CC1/5.progression_validation/signatures"
top_n_set <- c(30, 50, 100)

main_top_n_for_cleaner_forest <- 50

signature_info <- tibble::tribble(
  ~signature_prefix,                        ~score_col,                    ~group_col,
  "cluster9_malignant_squamous_epithelial", "malignant_squamous_score",   "malignant_squamous_group",
  "cluster10_EMT_like_epithelial",          "EMT_like_score",             "EMT_like_group",
  "cluster6_proliferation_core_epithelial", "proliferation_core_score",   "proliferation_core_group"
)

out_dir <- file.path(base_dir, "02.signature_score_survival_topN_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# HELPERS
# =========================
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
      dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
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

safe_stage_numeric <- function(x) {
  map <- c("IA2" = 1, "IB1" = 2, "IB2" = 3, "IIA" = 4)
  unname(map[x])
}

make_group <- function(x) {
  x_num <- as.numeric(x)

  if (all(is.na(x_num))) {
    return(factor(rep(NA_character_, length(x_num)), levels = c("Low risk", "High risk")))
  }

  med <- median(x_num, na.rm = TRUE)

  grp <- ifelse(
    is.na(x_num), NA_character_,
    ifelse(x_num >= med, "High risk", "Low risk")
  )

  if (length(unique(stats::na.omit(grp))) < 2) {
    ord <- order(x_num, na.last = TRUE)
    non_na_idx <- which(!is.na(x_num))
    n_non_na <- length(non_na_idx)
    low_n <- floor(n_non_na / 2)

    grp2 <- rep(NA_character_, length(x_num))
    grp2[ord[seq_len(low_n)]] <- "Low risk"
    grp2[ord[(low_n + 1):n_non_na]] <- "High risk"
    grp <- grp2
  }

  factor(grp, levels = c("Low risk", "High risk"))
}

safe_kruskal <- function(df, score_col) {
  tryCatch({
    fit <- kruskal.test(df[[score_col]] ~ df$stage)
    tibble(
      signature = score_col,
      statistic = as.numeric(fit$statistic),
      p.value = as.numeric(fit$p.value),
      method = fit$method
    )
  }, error = function(e) {
    tibble(
      signature = score_col,
      statistic = NA_real_,
      p.value = NA_real_,
      method = paste0("ERROR: ", e$message)
    )
  })
}

pairwise_wilcox_one <- function(df, score_col) {
  tryCatch({
    pw <- pairwise.wilcox.test(df[[score_col]], df$stage, p.adjust.method = "BH")
    mat <- pw$p.value

    if (is.null(mat)) return(tibble())

    res <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(res) <- c("group1", "group2", "p.adj")

    res %>%
      dplyr::filter(!is.na(p.adj)) %>%
      dplyr::mutate(signature = score_col) %>%
      dplyr::select(signature, group1, group2, p.adj)
  }, error = function(e) {
    tibble(
      signature = score_col,
      group1 = NA_character_,
      group2 = NA_character_,
      p.adj = NA_real_
    )
  })
}

safe_cox_fit <- function(formula_obj, data_df) {
  tryCatch(
    coxph(formula_obj, data = data_df),
    error = function(e) NULL,
    warning = function(w) suppressWarnings(coxph(formula_obj, data = data_df))
  )
}

extract_cox_table <- function(fit, model_name, top_n, signature_name) {
  if (is.null(fit)) {
    return(tibble(
      top_n = top_n,
      signature = signature_name,
      model = model_name,
      term = NA_character_,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      p.value = NA_real_
    ))
  }

  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::mutate(
      top_n = top_n,
      signature = signature_name,
      model = model_name
    ) %>%
    dplyr::select(top_n, signature, model, term, estimate, conf.low, conf.high, p.value)
}

safe_logrank_p <- function(df, time_col, status_col, group_col) {
  tryCatch({
    fml <- as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", group_col))
    sdiff <- survdiff(fml, data = df)
    1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1)
  }, error = function(e) {
    NA_real_
  })
}

format_p_text <- function(p, digits = 4) {
  if (is.na(p)) return("NA")
  if (p < 0.0001) return("< 0.0001")
  formatC(p, format = "f", digits = digits)
}

score_title_map <- c(
  "malignant_squamous_score" = "malignant squamous score",
  "EMT_like_score" = "EMT-like score",
  "proliferation_core_score" = "proliferation core score"
)

score_label_map <- c(
  "malignant_squamous_score" = "malignant squamous",
  "EMT_like_score" = "EMT-like",
  "proliferation_core_score" = "proliferation core"
)

score_color_map <- c(
  "proliferation core" = "#EE6352",
  "malignant squamous" = "#5DB7DE",
  "EMT-like" = "#43AA8B"
)

theme_pub <- function(base_size = 18) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = base_size + 5,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        face = "plain",
        size = base_size + 1,
        margin = margin(t = 0, b = 16)
      ),
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(t = 18, r = 20, b = 12, l = 12)
    )
}

make_stage_boxplot <- function(df, score_col, ylab_txt, title_txt, out_pdf, out_png = NULL) {
  stage_cols <- c(
    "IA2" = "#5DB7DE",
    "IB1" = "#43AA8B",
    "IB2" = "#EE6352",
    "IIA" = "#5B6C9E"
  )

  plot_df <- df %>%
    dplyr::filter(!is.na(stage), !is.na(.data[[score_col]])) %>%
    dplyr::mutate(
      stage = factor(as.character(stage), levels = c("IA2", "IB1", "IB2", "IIA"))
    )

  kw_res <- tryCatch(
    kruskal.test(plot_df[[score_col]] ~ plot_df$stage),
    error = function(e) NULL
  )

  kw_p_text <- if (is.null(kw_res)) {
    "Kruskal-Wallis, p = NA"
  } else {
    paste0("Kruskal-Wallis, p = ", formatC(kw_res$p.value, format = "f", digits = 4))
  }

  ymax <- max(plot_df[[score_col]], na.rm = TRUE)
  ymin <- min(plot_df[[score_col]], na.rm = TRUE)
  yrng <- ymax - ymin
  if (!is.finite(yrng) || yrng <= 0) yrng <- 1

  ylim_upper <- ymax + max(0.12 * yrng, 0.35)
  ylim_lower <- ymin - max(0.05 * yrng, 0.08)

  p <- ggplot(plot_df, aes(x = stage, y = .data[[score_col]], fill = stage)) +
    geom_boxplot(
      width = 0.68,
      color = "black",
      outlier.shape = 1,
      outlier.size = 2.2,
      size = 0.9,
      alpha = 0.95,
      na.rm = TRUE
    ) +
    geom_jitter(
      width = 0.12,
      size = 2.0,
      alpha = 0.55,
      color = "black",
      na.rm = TRUE
    ) +
    scale_fill_manual(values = stage_cols, drop = FALSE) +
    labs(
      title = title_txt,
      subtitle = kw_p_text,
      x = "Stage",
      y = ylab_txt
    ) +
    coord_cartesian(
      ylim = c(ylim_lower, ylim_upper),
      clip = "off"
    ) +
    theme_pub(base_size = 18) +
    theme(
      legend.position = "none"
    )

  pdf(out_pdf, width = 11, height = 9)
  print(p)
  dev.off()

  if (!is.null(out_png)) {
    ggsave(
      filename = out_png,
      plot = p,
      width = 11,
      height = 9,
      dpi = 300
    )
  }
}

plot_km_pub <- function(df, time_col, status_col, group_col, title_text, out_pdf, out_png = NULL) {
  surv_formula <- as.formula(
    paste0("Surv(", time_col, ", ", status_col, ") ~ ", group_col)
  )

  fit <- survfit(surv_formula, data = df)
  cox_fit <- safe_cox_fit(surv_formula, df)
  logrank_p <- safe_logrank_p(df, time_col, status_col, group_col)

  hr <- NA_real_
  lower <- NA_real_
  upper <- NA_real_
  cox_p <- NA_real_

  if (!is.null(cox_fit)) {
    td <- broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE)
    if (nrow(td) >= 1) {
      hr <- td$estimate[1]
      lower <- td$conf.low[1]
      upper <- td$conf.high[1]
      cox_p <- td$p.value[1]
    }
  }

  ann_text <- paste0(
    "log-rank p = ", ifelse(is.na(logrank_p), "NA", formatC(logrank_p, format = "f", digits = 4)), "\n",
    "HR = ", ifelse(is.na(hr), "NA", formatC(hr, format = "f", digits = 3)), "\n",
    "95% CI: ",
    ifelse(is.na(lower), "NA", formatC(lower, format = "f", digits = 3)),
    "–",
    ifelse(is.na(upper), "NA", formatC(upper, format = "f", digits = 3))
  )

  surv_df <- survminer::surv_summary(fit, data = df)
  surv_df$strata <- gsub(paste0(group_col, "="), "", surv_df$strata)
  surv_df$strata <- factor(surv_df$strata, levels = c("Low risk", "High risk"))
  censor_df <- surv_df %>% dplyr::filter(n.censor > 0)

  max_time <- max(df[[time_col]], na.rm = TRUE)
  max_time_axis <- max(5, ceiling(max_time))
  time_breaks <- 0:max_time_axis

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
      name = "Risk group"
    ) +
    coord_cartesian(xlim = c(0, max_time_axis + 0.2), ylim = c(0, 1.05)) +
    scale_x_continuous(breaks = time_breaks) +
    labs(
      title = title_text,
      x = "Time (years)",
      y = "Disease-free survival"
    ) +
    annotate(
      "text",
      x = max_time_axis * 0.45,
      y = 0.18,
      label = ann_text,
      size = 6.2,
      fontface = "bold"
    ) +
    theme_pub(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "top"
    )

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
    scale_x_continuous(breaks = time_breaks, limits = c(0, max_time_axis + 0.2)) +
    labs(
      title = "Number at risk",
      x = "Time (years)",
      y = "Risk group"
    ) +
    theme_pub(base_size = 17) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none"
    )

  final_plot <- ggpubr::ggarrange(
    curve_plot,
    risk_plot,
    ncol = 1,
    heights = c(3.3, 1.2),
    align = "v"
  )

  pdf(out_pdf, width = 11, height = 9)
  print(final_plot)
  dev.off()

  if (!is.null(out_png)) {
    ggsave(
      filename = out_png,
      plot = final_plot,
      width = 11,
      height = 9,
      dpi = 300
    )
  }

  invisible(
    list(
      fit = fit,
      cox_fit = cox_fit,
      logrank_p = logrank_p,
      cox_p = cox_p,
      hr = hr,
      conf.low = lower,
      conf.high = upper
    )
  )
}

create_cleaner_grouped_cox_forest <- function(cox_df, title_text, out_pdf, out_png = NULL) {
  plot_df <- cox_df %>%
    dplyr::filter(
      signature %in% c("proliferation_core_score", "malignant_squamous_score", "EMT_like_score"),
      !is.na(estimate),
      !is.na(conf.low),
      !is.na(conf.high)
    ) %>%
    dplyr::mutate(
      program = dplyr::case_when(
        signature == "proliferation_core_score" ~ "proliferation core",
        signature == "malignant_squamous_score" ~ "malignant squamous",
        signature == "EMT_like_score" ~ "EMT-like",
        TRUE ~ NA_character_
      )
    )

  order_df <- tibble::tibble(
    program = c("proliferation core", "malignant squamous", "EMT-like"),
    y = c(3, 2, 1)
  )

  plot_df <- plot_df %>%
    dplyr::left_join(order_df, by = "program") %>%
    dplyr::mutate(
      hr_ci = paste0(
        formatC(estimate, format = "f", digits = 3),
        " (",
        formatC(conf.low, format = "f", digits = 3),
        "–",
        formatC(conf.high, format = "f", digits = 3),
        ")"
      ),
      p_label = dplyr::case_when(
        is.na(p.value) ~ "NA",
        p.value < 0.0001 ~ "<0.0001",
        TRUE ~ formatC(p.value, format = "f", digits = 4)
      )
    )

  if (nrow(plot_df) == 0) {
    warning("No valid rows available for cleaner grouped Cox forest plot.")
    return(invisible(NULL))
  }

  color_map <- c(
    "proliferation core" = "#EE6352",
    "malignant squamous" = "#5DB7DE",
    "EMT-like" = "#43AA8B"
  )

  xmin_raw <- min(c(plot_df$conf.low, 1), na.rm = TRUE)
  xmax_raw <- max(c(plot_df$conf.high, 1), na.rm = TRUE)
  xrng <- xmax_raw - xmin_raw
  if (!is.finite(xrng) || xrng <= 0) xrng <- 1

  x_left   <- max(0, xmin_raw - 0.12 * xrng)
  x_hrcol  <- xmax_raw + 0.14 * xrng
  x_pcol   <- xmax_raw + 1.05 * xrng
  x_right  <- xmax_raw + 1.45 * xrng

  cap_h <- 0.10

  p <- ggplot(plot_df, aes(y = y)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.7, color = "grey40") +

    geom_segment(
      aes(x = conf.low, xend = conf.high, yend = y),
      linewidth = 1.15,
      color = "black"
    ) +
    geom_segment(
      aes(x = conf.low, xend = conf.low, y = y - cap_h / 2, yend = y + cap_h / 2),
      linewidth = 1.05,
      color = "black"
    ) +
    geom_segment(
      aes(x = conf.high, xend = conf.high, y = y - cap_h / 2, yend = y + cap_h / 2),
      linewidth = 1.05,
      color = "black"
    ) +

    geom_point(
      aes(x = estimate, fill = program),
      shape = 21,
      size = 4.8,
      stroke = 1,
      color = "black"
    ) +

    geom_text(
      aes(x = x_hrcol, label = hr_ci),
      hjust = 0,
      size = 5.4
    ) +
    geom_text(
      aes(x = x_pcol, label = p_label),
      hjust = 0,
      size = 5.4
    ) +

    annotate(
      "text",
      x = x_hrcol,
      y = 3.72,
      label = "HR (95% CI)",
      hjust = 0,
      fontface = "bold",
      size = 5.8
    ) +
    annotate(
      "text",
      x = x_pcol,
      y = 3.72,
      label = "P value",
      hjust = 0,
      fontface = "bold",
      size = 5.8
    ) +

    scale_fill_manual(values = color_map, guide = "none") +
    scale_y_continuous(
      breaks = c(3, 2, 1),
      labels = c("proliferation core", "malignant squamous", "EMT-like"),
      limits = c(0.5, 4.0)
    ) +
    coord_cartesian(
      xlim = c(x_left, x_right),
      clip = "off"
    ) +
    labs(
      title = title_text,
      x = "Hazard ratio",
      y = NULL
    ) +
    theme_pub(base_size = 18) +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(t = 18, r = 120, b = 12, l = 14)
    )

  pdf(out_pdf, width = 12.5, height = 6.4)
  print(p)
  dev.off()

  if (!is.null(out_png)) {
    ggsave(
      filename = out_png,
      plot = p,
      width = 12.5,
      height = 6.4,
      dpi = 300
    )
  }
}

# =========================
# FILE CHECK
# =========================
sig_files <- unlist(lapply(top_n_set, function(n_top) {
  c(
    file.path(sig_dir, paste0("cluster9_malignant_squamous_epithelial_top", n_top, ".txt")),
    file.path(sig_dir, paste0("cluster10_EMT_like_epithelial_top", n_top, ".txt")),
    file.path(sig_dir, paste0("cluster6_proliferation_core_epithelial_top", n_top, ".txt"))
  )
}))

need_files <- c(expr_file, clin_file, sig_files)
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

expr_num <- expr
for (j in 2:ncol(expr_num)) {
  expr_num[[j]] <- as.numeric(expr_num[[j]])
}

expr_num <- make_unique_gene(expr_num, gene_col = "geneNames")

expr_mat <- as.matrix(expr_num[, -1])
rownames(expr_mat) <- expr_num$geneNames
mode(expr_mat) <- "numeric"

# =========================
# READ CLINICAL
# =========================
clin <- read.csv(
  clin_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

req_cols <- c("sample", "stage", "largest_diameter", "dfs_months", "dfs_status")
miss_cols <- req_cols[!req_cols %in% colnames(clin)]

if (length(miss_cols) > 0) {
  stop("Clinical file missing required columns: ", paste(miss_cols, collapse = ", "))
}

clin <- clin %>%
  dplyr::mutate(
    sample = trimws(sample),
    stage = trimws(stage),
    stage = factor(stage, levels = c("IA2", "IB1", "IB2", "IIA")),
    largest_diameter = as.numeric(largest_diameter),
    dfs_months = as.numeric(dfs_months),
    dfs_status = as.numeric(dfs_status),
    dfs_years = dfs_months / 12
  )

# =========================
# MATCH SAMPLE ORDER
# =========================
common_samples <- intersect(colnames(expr_mat), clin$sample)

if (length(common_samples) == 0) {
  stop("No matched samples between expression and clinical data.")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
clin <- clin %>% dplyr::filter(sample %in% common_samples)
clin <- clin[match(common_samples, clin$sample), ]
stopifnot(all(clin$sample == colnames(expr_mat)))

# =========================
# STORAGE
# =========================
all_usage <- list()
all_score_summary <- list()
all_stage_kw <- list()
all_pairwise <- list()
all_continuous_univ <- list()
all_grouped_cox <- list()
all_multiv <- list()
all_stagewise <- list()
all_ana <- list()
all_km_summary <- list()

# =========================
# MAIN LOOP
# =========================
for (n_top in top_n_set) {
  cat("========================================================\n")
  cat("Running top", n_top, "...\n", sep = "")

  top_out_dir <- file.path(out_dir, paste0("top", n_top))
  dir.create(top_out_dir, recursive = TRUE, showWarnings = FALSE)

  top_plot_dir <- file.path(plot_dir, paste0("top", n_top))
  dir.create(top_plot_dir, recursive = TRUE, showWarnings = FALSE)

  sig9_file  <- file.path(sig_dir, paste0("cluster9_malignant_squamous_epithelial_top", n_top, ".txt"))
  sig10_file <- file.path(sig_dir, paste0("cluster10_EMT_like_epithelial_top", n_top, ".txt"))
  sig6_file  <- file.path(sig_dir, paste0("cluster6_proliferation_core_epithelial_top", n_top, ".txt"))

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

  usage_df <- dplyr::bind_rows(
    data.frame(
      top_n = n_top,
      signature = "cluster9_malignant_squamous_epithelial",
      score_col = "malignant_squamous_score",
      genes_input = length(sig9),
      genes_used = length(sig9_res$genes_used),
      genes_used_list = paste(sig9_res$genes_used, collapse = ";"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      top_n = n_top,
      signature = "cluster10_EMT_like_epithelial",
      score_col = "EMT_like_score",
      genes_input = length(sig10),
      genes_used = length(sig10_res$genes_used),
      genes_used_list = paste(sig10_res$genes_used, collapse = ";"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      top_n = n_top,
      signature = "cluster6_proliferation_core_epithelial",
      score_col = "proliferation_core_score",
      genes_input = length(sig6),
      genes_used = length(sig6_res$genes_used),
      genes_used_list = paste(sig6_res$genes_used, collapse = ";"),
      stringsAsFactors = FALSE
    )
  )

  write.csv(
    usage_df,
    file.path(top_out_dir, "01.signature_usage.csv"),
    row.names = FALSE
  )

  ana <- dplyr::left_join(clin, score_df, by = "sample") %>%
    dplyr::mutate(
      top_n = n_top,
      stage_num = safe_stage_numeric(as.character(stage)),
      malignant_squamous_group = make_group(malignant_squamous_score),
      EMT_like_group = make_group(EMT_like_score),
      proliferation_core_group = make_group(proliferation_core_score)
    )

  write.csv(
    ana,
    file.path(top_out_dir, "04.GSE44001_clinical_with_scores_and_groups.csv"),
    row.names = FALSE
  )

  score_summary <- ana %>%
    dplyr::summarise(
      top_n = n_top,
      n = dplyr::n(),
      malignant_squamous_mean = mean(malignant_squamous_score, na.rm = TRUE),
      EMT_like_mean = mean(EMT_like_score, na.rm = TRUE),
      proliferation_core_mean = mean(proliferation_core_score, na.rm = TRUE)
    )

  write.csv(
    score_summary,
    file.path(top_out_dir, "05.score_summary.csv"),
    row.names = FALSE
  )

  stage_assoc <- dplyr::bind_rows(
    safe_kruskal(ana, "malignant_squamous_score"),
    safe_kruskal(ana, "EMT_like_score"),
    safe_kruskal(ana, "proliferation_core_score")
  ) %>%
    dplyr::mutate(top_n = n_top) %>%
    dplyr::select(top_n, dplyr::everything())

  write.csv(
    stage_assoc,
    file.path(top_out_dir, "06.score_stage_kruskal.csv"),
    row.names = FALSE
  )

  pairwise_df <- dplyr::bind_rows(
    pairwise_wilcox_one(ana, "malignant_squamous_score"),
    pairwise_wilcox_one(ana, "EMT_like_score"),
    pairwise_wilcox_one(ana, "proliferation_core_score")
  ) %>%
    dplyr::mutate(top_n = n_top) %>%
    dplyr::select(top_n, dplyr::everything())

  write.csv(
    pairwise_df,
    file.path(top_out_dir, "07.score_stage_pairwise_wilcox.csv"),
    row.names = FALSE
  )

  stagewise_df <- ana %>%
    dplyr::select(sample, stage, malignant_squamous_score, EMT_like_score, proliferation_core_score) %>%
    tidyr::pivot_longer(
      cols = c(malignant_squamous_score, EMT_like_score, proliferation_core_score),
      names_to = "signature",
      values_to = "score"
    ) %>%
    dplyr::group_by(stage, signature) %>%
    dplyr::summarise(
      top_n = n_top,
      n = dplyr::n(),
      mean = mean(score, na.rm = TRUE),
      median = median(score, na.rm = TRUE),
      sd = sd(score, na.rm = TRUE),
      q1 = quantile(score, 0.25, na.rm = TRUE),
      q3 = quantile(score, 0.75, na.rm = TRUE),
      min = min(score, na.rm = TRUE),
      max = max(score, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(
    stagewise_df,
    file.path(top_out_dir, "13.stage_wise_score_summary.csv"),
    row.names = FALSE
  )

  # -------------------------
  # Stage boxplots
  # only global p, no pairwise ns/*
  # -------------------------
  make_stage_boxplot(
    df = ana,
    score_col = "malignant_squamous_score",
    ylab_txt = "Score",
    title_txt = paste0("malignant squamous score by stage (top", n_top, ")"),
    out_pdf = file.path(top_plot_dir, paste0("PUB_stage_malignant_squamous_score_top", n_top, ".pdf")),
    out_png = file.path(top_plot_dir, paste0("PUB_stage_malignant_squamous_score_top", n_top, ".png"))
  )

  make_stage_boxplot(
    df = ana,
    score_col = "EMT_like_score",
    ylab_txt = "Score",
    title_txt = paste0("EMT-like score by stage (top", n_top, ")"),
    out_pdf = file.path(top_plot_dir, paste0("PUB_stage_EMT_like_score_top", n_top, ".pdf")),
    out_png = file.path(top_plot_dir, paste0("PUB_stage_EMT_like_score_top", n_top, ".png"))
  )

  make_stage_boxplot(
    df = ana,
    score_col = "proliferation_core_score",
    ylab_txt = "Score",
    title_txt = paste0("proliferation core score by stage (top", n_top, ")"),
    out_pdf = file.path(top_plot_dir, paste0("PUB_stage_proliferation_core_score_top", n_top, ".pdf")),
    out_png = file.path(top_plot_dir, paste0("PUB_stage_proliferation_core_score_top", n_top, ".png"))
  )

  # -------------------------
  # KM + grouped Cox
  # -------------------------
  km_summary_list <- list()
  grouped_cox_list <- list()

  for (i in seq_len(nrow(signature_info))) {
    score_col <- signature_info$score_col[i]
    group_col <- signature_info$group_col[i]
    sig_name <- score_col

    km_title <- score_title_map[[score_col]]
    km_out_pdf <- file.path(top_plot_dir, paste0("PUB_KM_", score_col, "_top", n_top, ".pdf"))
    km_out_png <- file.path(top_plot_dir, paste0("PUB_KM_", score_col, "_top", n_top, ".png"))

    km_res <- plot_km_pub(
      df = ana,
      time_col = "dfs_years",
      status_col = "dfs_status",
      group_col = group_col,
      title_text = km_title,
      out_pdf = km_out_pdf,
      out_png = km_out_png
    )

    km_summary_list[[sig_name]] <- tibble(
      top_n = n_top,
      signature = sig_name,
      logrank_p = km_res$logrank_p,
      grouped_cox_p = km_res$cox_p,
      grouped_cox_hr = km_res$hr,
      grouped_cox_conf_low = km_res$conf.low,
      grouped_cox_conf_high = km_res$conf.high
    )

    grp_formula <- as.formula(paste0("Surv(dfs_years, dfs_status) ~ ", group_col))
    grp_fit <- safe_cox_fit(grp_formula, ana)

    grouped_cox_list[[sig_name]] <- extract_cox_table(
      fit = grp_fit,
      model_name = "grouped_univariate",
      top_n = n_top,
      signature_name = sig_name
    )
  }

  km_summary_df <- dplyr::bind_rows(km_summary_list)
  grp_df <- dplyr::bind_rows(grouped_cox_list)

  write.csv(
    km_summary_df,
    file.path(top_out_dir, "09.km_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    grp_df,
    file.path(top_out_dir, "11.grouped_cox_results.csv"),
    row.names = FALSE
  )

  # -------------------------
  # Cleaner grouped Cox forest plot for each topN
  # -------------------------
  create_cleaner_grouped_cox_forest(
    cox_df = grp_df,
    title_text = paste0("GSE44001 Cox summary (top", n_top, ")"),
    out_pdf = file.path(top_plot_dir, paste0("PUB_cleaner_grouped_cox_forest_top", n_top, ".pdf")),
    out_png = file.path(top_plot_dir, paste0("PUB_cleaner_grouped_cox_forest_top", n_top, ".png"))
  )

  # -------------------------
  # Continuous univariate Cox
  # -------------------------
  uni_list <- list()

  uni_formulas <- list(
    malignant_squamous_score = as.formula("Surv(dfs_years, dfs_status) ~ malignant_squamous_score"),
    EMT_like_score = as.formula("Surv(dfs_years, dfs_status) ~ EMT_like_score"),
    proliferation_core_score = as.formula("Surv(dfs_years, dfs_status) ~ proliferation_core_score")
  )

  for (nm in names(uni_formulas)) {
    fit_uni <- safe_cox_fit(uni_formulas[[nm]], ana)
    uni_list[[nm]] <- extract_cox_table(
      fit = fit_uni,
      model_name = "continuous_univariate",
      top_n = n_top,
      signature_name = nm
    )
  }

  uni_df <- dplyr::bind_rows(uni_list)

  write.csv(
    uni_df,
    file.path(top_out_dir, "10.univariate_cox_results.csv"),
    row.names = FALSE
  )

  # -------------------------
  # Continuous multivariate Cox
  # -------------------------
  mv_list <- list()

  mv_formulas <- list(
    malignant_squamous_score = as.formula("Surv(dfs_years, dfs_status) ~ malignant_squamous_score + stage_num + largest_diameter"),
    EMT_like_score = as.formula("Surv(dfs_years, dfs_status) ~ EMT_like_score + stage_num + largest_diameter"),
    proliferation_core_score = as.formula("Surv(dfs_years, dfs_status) ~ proliferation_core_score + stage_num + largest_diameter")
  )

  mv_df_input <- ana %>%
    dplyr::filter(
      !is.na(dfs_years),
      !is.na(dfs_status),
      !is.na(stage_num),
      !is.na(largest_diameter)
    )

  for (nm in names(mv_formulas)) {
    fit_mv <- safe_cox_fit(mv_formulas[[nm]], mv_df_input)
    mv_list[[nm]] <- extract_cox_table(
      fit = fit_mv,
      model_name = "continuous_multivariate",
      top_n = n_top,
      signature_name = nm
    )
  }

  mv_df <- dplyr::bind_rows(mv_list)

  write.csv(
    mv_df,
    file.path(top_out_dir, "12.multivariate_cox_results.csv"),
    row.names = FALSE
  )

  # -------------------------
  # Collect all
  # -------------------------
  all_usage[[paste0("top", n_top)]] <- usage_df
  all_score_summary[[paste0("top", n_top)]] <- score_summary
  all_stage_kw[[paste0("top", n_top)]] <- stage_assoc
  all_pairwise[[paste0("top", n_top)]] <- pairwise_df
  all_continuous_univ[[paste0("top", n_top)]] <- uni_df
  all_grouped_cox[[paste0("top", n_top)]] <- grp_df
  all_multiv[[paste0("top", n_top)]] <- mv_df
  all_stagewise[[paste0("top", n_top)]] <- stagewise_df
  all_ana[[paste0("top", n_top)]] <- ana
  all_km_summary[[paste0("top", n_top)]] <- km_summary_df
}

# =========================
# COMBINE ALL RESULTS
# =========================
usage_all <- dplyr::bind_rows(all_usage)
score_summary_all <- dplyr::bind_rows(all_score_summary)
stage_kw_all <- dplyr::bind_rows(all_stage_kw)
pairwise_all <- dplyr::bind_rows(all_pairwise)
uni_all <- dplyr::bind_rows(all_continuous_univ)
grp_all <- dplyr::bind_rows(all_grouped_cox)
mv_all <- dplyr::bind_rows(all_multiv)
stagewise_all <- dplyr::bind_rows(all_stagewise)
ana_all <- dplyr::bind_rows(all_ana)
km_summary_all <- dplyr::bind_rows(all_km_summary)

write.csv(
  usage_all,
  file.path(out_dir, "01.signature_usage_all_topN.csv"),
  row.names = FALSE
)

write.csv(
  score_summary_all,
  file.path(out_dir, "05.score_summary_all_topN.csv"),
  row.names = FALSE
)

write.csv(
  stage_kw_all,
  file.path(out_dir, "06.score_stage_kruskal_all_topN.csv"),
  row.names = FALSE
)

write.csv(
  pairwise_all,
  file.path(out_dir, "07.score_stage_pairwise_wilcox_all_topN.csv"),
  row.names = FALSE
)

write.csv(
  ana_all,
  file.path(out_dir, "08.GSE44001_clinical_with_scores_and_groups_all_topN.csv"),
  row.names = FALSE
)

write.csv(
  km_summary_all,
  file.path(out_dir, "09.km_summary_all_topN.csv"),
  row.names = FALSE
)

write.csv(
  uni_all,
  file.path(out_dir, "10.univariate_cox_results_all_topN.csv"),
  row.names = FALSE
)

write.csv(
  grp_all,
  file.path(out_dir, "11.grouped_cox_results_all_topN.csv"),
  row.names = FALSE
)

write.csv(
  mv_all,
  file.path(out_dir, "12.multivariate_cox_results_all_topN.csv"),
  row.names = FALSE
)

write.csv(
  stagewise_all,
  file.path(out_dir, "13.stage_wise_score_summary_all_topN.csv"),
  row.names = FALSE
)

# =========================
# ROBUSTNESS SUMMARY
# =========================
kw_small <- stage_kw_all %>%
  dplyr::select(top_n, signature, stage_kw_p = p.value)

usage_small <- usage_all %>%
  dplyr::select(top_n, score_col, genes_input, genes_used)

km_small <- km_summary_all %>%
  dplyr::select(
    top_n, signature, logrank_p,
    grouped_cox_p, grouped_cox_hr,
    grouped_cox_conf_low, grouped_cox_conf_high
  )

uni_small <- uni_all %>%
  dplyr::filter(term %in% c("malignant_squamous_score", "EMT_like_score", "proliferation_core_score")) %>%
  dplyr::select(
    top_n, signature,
    continuous_univ_hr = estimate,
    continuous_univ_conf_low = conf.low,
    continuous_univ_conf_high = conf.high,
    continuous_univ_p = p.value
  )

mv_small <- mv_all %>%
  dplyr::filter(term %in% c("malignant_squamous_score", "EMT_like_score", "proliferation_core_score")) %>%
  dplyr::select(
    top_n, signature,
    continuous_multiv_hr = estimate,
    continuous_multiv_conf_low = conf.low,
    continuous_multiv_conf_high = conf.high,
    continuous_multiv_p = p.value
  )

robustness_summary <- usage_small %>%
  dplyr::rename(signature = score_col) %>%
  dplyr::left_join(kw_small, by = c("top_n", "signature")) %>%
  dplyr::left_join(km_small, by = c("top_n", "signature")) %>%
  dplyr::left_join(uni_small, by = c("top_n", "signature")) %>%
  dplyr::left_join(mv_small, by = c("top_n", "signature")) %>%
  dplyr::arrange(
    match(signature, c("malignant_squamous_score", "EMT_like_score", "proliferation_core_score")),
    top_n
  )

write.csv(
  robustness_summary,
  file.path(out_dir, "14.robustness_summary.csv"),
  row.names = FALSE
)

# =========================
# MAIN MANUSCRIPT-STYLE TOP50 CLEANER FOREST
# =========================
grp_top50 <- grp_all %>%
  dplyr::filter(top_n == main_top_n_for_cleaner_forest)

write.csv(
  grp_top50,
  file.path(plot_dir, paste0("Figure2C_cleaner_grouped_cox_forest_top", main_top_n_for_cleaner_forest, "_source.csv")),
  row.names = FALSE
)

create_cleaner_grouped_cox_forest(
  cox_df = grp_top50,
  title_text = "GSE44001 Cox summary of prioritized epithelial programs",
  out_pdf = file.path(plot_dir, paste0("Figure2C_cleaner_grouped_cox_forest_top", main_top_n_for_cleaner_forest, ".pdf")),
  out_png = file.path(plot_dir, paste0("Figure2C_cleaner_grouped_cox_forest_top", main_top_n_for_cleaner_forest, ".png"))
)

# =========================
# RUN INFO
# =========================
run_info <- c(
  "GSE44001 topN robustness validation finished.",
  paste0("base_dir: ", base_dir),
  paste0("expr_file: ", expr_file),
  paste0("clin_file: ", clin_file),
  paste0("sig_dir: ", sig_dir),
  paste0("top_n_set: ", paste(top_n_set, collapse = ", ")),
  paste0("main_top_n_for_cleaner_forest: ", main_top_n_for_cleaner_forest),
  paste0("n matched samples: ", length(common_samples)),
  paste0("output dir: ", out_dir)
)

writeLines(
  run_info,
  file.path(out_dir, "00.run_info.txt")
)

cat("========================================================\n")
cat("Done.\n")
cat("Output directory:\n", out_dir, "\n\n")
cat("Stage plot files:\n")
cat("plots/top30/PUB_stage_*_top30.pdf/.png\n")
cat("plots/top50/PUB_stage_*_top50.pdf/.png\n")
cat("plots/top100/PUB_stage_*_top100.pdf/.png\n\n")
cat("Cleaner grouped Cox forest files:\n")
cat("plots/top30/PUB_cleaner_grouped_cox_forest_top30.pdf/.png\n")
cat("plots/top50/PUB_cleaner_grouped_cox_forest_top50.pdf/.png\n")
cat("plots/top100/PUB_cleaner_grouped_cox_forest_top100.pdf/.png\n\n")
cat("Main manuscript-style cleaner forest:\n")
cat(paste0("plots/Figure2C_cleaner_grouped_cox_forest_top", main_top_n_for_cleaner_forest, ".pdf/.png\n"))
cat("========================================================\n")