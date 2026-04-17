rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
options(stringsAsFactors = FALSE)

############################
# 0. Packages
############################
cran_pkgs <- c(
  "data.table", "dplyr", "tibble", "ggplot2", "tidyr",
  "survival", "ggpubr", "survminer"
)
bioc_pkgs <- c("GSVA")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(data.table)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(survival)
library(ggpubr)
library(survminer)
library(GSVA)

############################
# 1. Paths
############################
base_dir <- "/Users/donglinlu/Desktop/CC1/8.Hub Gene"
setwd(base_dir)

sig_file                 <- file.path(base_dir, "cluster6_proliferation_core_epithelial_top50.txt")
gse63514_expr_file       <- file.path(base_dir, "GSE63514_geneMatrix.txt")
gse63514_meta_file       <- file.path(base_dir, "GSE63514_metadata.txt")
gse44001_expr_file       <- file.path(base_dir, "GSE44001_geneMatrix.txt")
gse44001_score_meta_file <- file.path(base_dir, "GSE44001_clinical_score_metadata.csv")
gse44001_surv_file       <- file.path(base_dir, "gse44001_survival_metadata.txt")

out_dir <- file.path(base_dir, "HMMR_supporting_plots")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

lead_gene <- "HMMR"
marker_partners <- c("BIRC5", "MKI67", "TOP2A", "UBE2C")
progression_levels <- c("Normal", "LSIL", "HSIL", "Cancer")

############################
# 2. Helper functions
############################
fmt_p <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) {
      "NA"
    } else if (x < 1e-4) {
      formatC(x, format = "e", digits = 2)
    } else {
      sprintf("%.4f", x)
    }
  }, USE.NAMES = FALSE)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

read_signature <- function(file) {
  x <- fread(file, header = FALSE, data.table = FALSE)[, 1]
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & x != ""]
  x
}

prepare_expr_matrix <- function(expr_df, gene_col = "geneNames") {
  if (!gene_col %in% colnames(expr_df)) {
    stop(paste0("Expression matrix missing gene column: ", gene_col))
  }

  colnames(expr_df)[colnames(expr_df) == gene_col] <- "GeneSymbol"

  expr_df <- expr_df %>%
    filter(!is.na(GeneSymbol), GeneSymbol != "")

  sample_cols <- setdiff(colnames(expr_df), "GeneSymbol")
  expr_df[, sample_cols] <- lapply(expr_df[, sample_cols, drop = FALSE], as.numeric)

  expr_df <- expr_df[rowSums(is.na(expr_df[, sample_cols, drop = FALSE])) < length(sample_cols), ]

  expr_df <- expr_df %>%
    group_by(GeneSymbol) %>%
    summarise(across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  expr_mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
  rownames(expr_mat) <- expr_df$GeneSymbol
  expr_mat
}

find_first_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

find_gene_symbol <- function(target_gene, expr_mat) {
  rn <- rownames(expr_mat)

  if (target_gene %in% rn) return(target_gene)

  hit <- rn[toupper(rn) == toupper(target_gene)]
  if (length(hit) >= 1) return(hit[1])

  NA_character_
}

safe_spearman <- function(x, y) {
  tmp <- data.frame(x = x, y = y)
  tmp <- tmp[complete.cases(tmp), ]

  if (nrow(tmp) < 5) {
    return(list(rho = NA, p = NA, n = nrow(tmp)))
  }
  if (length(unique(tmp$x)) < 3 || length(unique(tmp$y)) < 2) {
    return(list(rho = NA, p = NA, n = nrow(tmp)))
  }

  res <- tryCatch(
    suppressWarnings(cor.test(tmp$x, tmp$y, method = "spearman")),
    error = function(e) NULL
  )

  if (is.null(res)) return(list(rho = NA, p = NA, n = nrow(tmp)))

  list(rho = unname(res$estimate), p = res$p.value, n = nrow(tmp))
}

safe_kw <- function(x, g) {
  tmp <- data.frame(x = x, g = g)
  tmp <- tmp[complete.cases(tmp), ]

  if (nrow(tmp) < 5) return(list(stat = NA, p = NA))
  if (length(unique(tmp$g)) < 2) return(list(stat = NA, p = NA))

  res <- tryCatch(
    kruskal.test(x ~ g, data = tmp),
    error = function(e) NULL
  )

  if (is.null(res)) return(list(stat = NA, p = NA))
  list(stat = unname(res$statistic), p = res$p.value)
}

safe_cox <- function(time, status, gene_exp) {
  tmp <- data.frame(time = time, status = status, gene = gene_exp)
  tmp <- tmp[complete.cases(tmp), ]

  if (nrow(tmp) < 10) {
    return(list(hr = NA, p = NA, low = NA, high = NA, n = nrow(tmp)))
  }
  if (length(unique(tmp$gene)) < 3) {
    return(list(hr = NA, p = NA, low = NA, high = NA, n = nrow(tmp)))
  }
  if (length(unique(tmp$status)) < 2) {
    return(list(hr = NA, p = NA, low = NA, high = NA, n = nrow(tmp)))
  }

  fit <- tryCatch(
    coxph(Surv(time, status) ~ gene, data = tmp),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(hr = NA, p = NA, low = NA, high = NA, n = nrow(tmp)))
  }

  s <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(s)) {
    return(list(hr = NA, p = NA, low = NA, high = NA, n = nrow(tmp)))
  }

  hr <- unname(s$coefficients[1, "exp(coef)"])
  p  <- unname(s$coefficients[1, "Pr(>|z|)"])
  ci <- tryCatch(s$conf.int[1, c("lower .95", "upper .95")], error = function(e) c(NA, NA))

  list(
    hr = hr,
    p = p,
    low = unname(ci[1]),
    high = unname(ci[2]),
    n = nrow(tmp)
  )
}

safe_grouped_cox <- function(time, status, group_factor) {
  tmp <- data.frame(time = time, status = status, group = group_factor)
  tmp <- tmp[complete.cases(tmp), ]

  if (nrow(tmp) < 10) return(list(hr = NA, p = NA, low = NA, high = NA))
  if (length(unique(tmp$group)) < 2) return(list(hr = NA, p = NA, low = NA, high = NA))
  if (length(unique(tmp$status)) < 2) return(list(hr = NA, p = NA, low = NA, high = NA))

  fit <- tryCatch(
    coxph(Surv(time, status) ~ group, data = tmp),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(hr = NA, p = NA, low = NA, high = NA))

  s <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(s)) return(list(hr = NA, p = NA, low = NA, high = NA))

  hr <- unname(s$coefficients[1, "exp(coef)"])
  p  <- unname(s$coefficients[1, "Pr(>|z|)"])
  ci <- tryCatch(s$conf.int[1, c("lower .95", "upper .95")], error = function(e) c(NA, NA))

  list(
    hr = hr,
    p = p,
    low = unname(ci[1]),
    high = unname(ci[2])
  )
}

make_binary_group <- function(x, low_label = "Low", high_label = "High") {
  x_num <- as.numeric(x)
  if (all(is.na(x_num))) {
    return(factor(rep(NA_character_, length(x_num)), levels = c(low_label, high_label)))
  }

  med <- median(x_num, na.rm = TRUE)
  grp <- ifelse(is.na(x_num), NA_character_,
                ifelse(x_num >= med, high_label, low_label))

  if (length(unique(stats::na.omit(grp))) < 2) {
    ord <- order(x_num, na.last = TRUE)
    non_na_idx <- which(!is.na(x_num))
    n_non_na <- length(non_na_idx)
    low_n <- floor(n_non_na / 2)

    grp2 <- rep(NA_character_, length(x_num))
    grp2[ord[seq_len(low_n)]] <- low_label
    grp2[ord[(low_n + 1):n_non_na]] <- high_label
    grp <- grp2
  }

  factor(grp, levels = c(low_label, high_label))
}

theme_pub <- function(base_size = 18) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 5),
      plot.subtitle = element_text(hjust = 0.5, size = base_size + 1, margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = base_size + 2),
      axis.text = element_text(size = base_size, colour = "black"),
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      plot.margin = margin(t = 10, r = 18, b = 10, l = 10)
    )
}

save_pdf <- function(plot_obj, path, width, height) {
  ggsave(
    filename = path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE
  )
}

make_pairwise_order <- function(levels_vec) {
  cmb <- combn(levels_vec, 2, simplify = FALSE)
  bind_rows(lapply(cmb, function(z) {
    data.frame(group1 = z[1], group2 = z[2], stringsAsFactors = FALSE)
  }))
}

plot_corr_pub <- function(df, xvar, yvar, title_text, xlab, ylab, out_pdf,
                          point_color = "#4E79A7",
                          add_lm = TRUE,
                          x_breaks = NULL,
                          x_labels = NULL) {
  tmp <- df[, c(xvar, yvar)]
  tmp <- tmp[complete.cases(tmp), , drop = FALSE]

  cor_res <- safe_spearman(tmp[[xvar]], tmp[[yvar]])
  sub_txt <- paste0(
    "Spearman rho = ", fmt_num(cor_res$rho, 3),
    "; p = ", fmt_p(cor_res$p)
  )

  p <- ggplot(df, aes_string(x = xvar, y = yvar)) +
    geom_point(size = 4.0, alpha = 0.78, color = point_color)

  if (add_lm) {
    p <- p +
      geom_smooth(
        method = "lm",
        se = TRUE,
        color = "#E41A1C",
        fill = "grey75",
        linewidth = 1.25
      )
  }

  if (!is.null(x_breaks) && !is.null(x_labels)) {
    p <- p + scale_x_continuous(breaks = x_breaks, labels = x_labels)
  }

  p <- p +
    labs(
      title = title_text,
      subtitle = sub_txt,
      x = xlab,
      y = ylab
    ) +
    theme_pub(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )

  save_pdf(p, out_pdf, width = 10.4, height = 8.0)
}

plot_km_pub <- function(df, time_col, status_col, group_col, title_text, out_pdf) {
  surv_formula <- as.formula(
    paste0("Surv(", time_col, ", ", status_col, ") ~ ", group_col)
  )

  fit <- survfit(surv_formula, data = df)

  sdiff <- tryCatch(
    survdiff(surv_formula, data = df),
    error = function(e) NULL
  )

  logrank_p <- NA_real_
  if (!is.null(sdiff)) {
    logrank_p <- 1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1)
  }

  grp_cox <- safe_grouped_cox(
    time = df[[time_col]],
    status = df[[status_col]],
    group_factor = df[[group_col]]
  )

  ann_text <- paste0(
    "log-rank p = ", fmt_p(logrank_p), "\n",
    "HR = ", fmt_num(grp_cox$hr, 3), "\n",
    "95% CI: ", fmt_num(grp_cox$low, 3), "–", fmt_num(grp_cox$high, 3)
  )

  surv_df <- survminer::surv_summary(fit, data = df)
  surv_df$strata <- gsub(paste0(group_col, "="), "", surv_df$strata)

  group_levels <- levels(df[[group_col]])
  if (is.null(group_levels)) group_levels <- unique(surv_df$strata)
  surv_df$strata <- factor(surv_df$strata, levels = group_levels)

  censor_df <- surv_df %>% dplyr::filter(n.censor > 0)

  max_time <- max(df[[time_col]], na.rm = TRUE)
  max_time_axis <- max(5, ceiling(max_time))
  time_breaks <- 0:max_time_axis

  curve_cols <- c(
    "Low HMMR" = "#2C7BE5",
    "High HMMR" = "#E74C3C"
  )

  curve_plot <- ggplot(surv_df, aes(x = time, y = surv, color = strata)) +
    geom_step(linewidth = 1.35) +
    geom_point(
      data = censor_df,
      aes(x = time, y = surv, color = strata),
      shape = 3,
      size = 3
    ) +
    scale_color_manual(values = curve_cols, name = "Group") +
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
  risk_df$strata <- factor(risk_df$strata, levels = group_levels)

  risk_plot <- ggplot(risk_df, aes(x = time, y = strata, label = n.risk, color = strata)) +
    geom_text(size = 6) +
    scale_color_manual(values = curve_cols, guide = "none") +
    scale_x_continuous(breaks = time_breaks, limits = c(0, max_time_axis + 0.2)) +
    labs(
      title = "Number at risk",
      x = "Time (years)",
      y = "Group"
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

  save_pdf(final_plot, out_pdf, width = 11.0, height = 9.2)

  invisible(
    list(
      logrank_p = logrank_p,
      hr = grp_cox$hr,
      low = grp_cox$low,
      high = grp_cox$high,
      cox_p = grp_cox$p
    )
  )
}

############################
# 3. Read files
############################
top50_genes <- read_signature(sig_file)

gse63514_expr_df  <- fread(gse63514_expr_file, data.table = FALSE)
gse63514_meta     <- fread(gse63514_meta_file, data.table = FALSE)
gse44001_expr_df  <- fread(gse44001_expr_file, data.table = FALSE)
gse44001_score_df <- fread(gse44001_score_meta_file, data.table = FALSE)
gse44001_surv_df  <- fread(gse44001_surv_file, data.table = FALSE)

############################
# 4. Check headers
############################
required_63514_meta_cols <- c("sample_id", "group", "group_order")
if (!all(required_63514_meta_cols %in% colnames(gse63514_meta))) {
  stop("GSE63514_metadata.txt is missing required columns: sample_id, group, group_order")
}

score_sample_col <- find_first_col(gse44001_score_df, c("sample", "Sample", "sample_id", "SampleID", "ID"))
if (is.na(score_sample_col)) {
  stop("GSE44001_clinical_score_metadata.csv is missing sample/ID column")
}
if (!("proliferation_core_score" %in% colnames(gse44001_score_df))) {
  stop("GSE44001_clinical_score_metadata.csv is missing proliferation_core_score column")
}

surv_sample_col <- find_first_col(gse44001_surv_df, c("ID", "sample", "Sample", "sample_id", "SampleID"))
surv_time_col <- find_first_col(
  gse44001_surv_df,
  c("disease_free_survival_(dfs)_(months)", "dfs_months", "DFS_months", "dfs_time", "time_months", "time")
)
surv_status_col <- find_first_col(
  gse44001_surv_df,
  c("status_of_dfs", "dfs_status", "DFS_status", "status", "event")
)
surv_stage_col <- find_first_col(gse44001_surv_df, c("Stage", "stage"))
surv_diameter_col <- find_first_col(gse44001_surv_df, c("largest diameter", "largest_diameter"))

if (is.na(surv_sample_col) || is.na(surv_time_col) || is.na(surv_status_col)) {
  stop("gse44001_survival_metadata.txt is missing ID / disease_free_survival_(dfs)_(months) / status_of_dfs columns")
}

############################
# 5. Standardize metadata names
############################
gse44001_score_df <- gse44001_score_df %>%
  rename(sample = all_of(score_sample_col)) %>%
  mutate(
    sample = as.character(sample),
    proliferation_core_score = as.numeric(proliferation_core_score)
  ) %>%
  select(sample, proliferation_core_score)

gse44001_surv_df <- gse44001_surv_df %>%
  rename(
    sample = all_of(surv_sample_col),
    dfs_months = all_of(surv_time_col),
    dfs_status = all_of(surv_status_col)
  ) %>%
  mutate(
    sample = as.character(sample),
    dfs_months = as.numeric(dfs_months),
    dfs_status = as.numeric(dfs_status),
    dfs_years = dfs_months / 12
  )

if (!is.na(surv_stage_col)) {
  gse44001_surv_df <- gse44001_surv_df %>%
    rename(stage = all_of(surv_stage_col))
}
if (!is.na(surv_diameter_col)) {
  gse44001_surv_df <- gse44001_surv_df %>%
    rename(largest_diameter = all_of(surv_diameter_col))
}

############################
# 6. Prepare expression matrices
############################
expr63514 <- prepare_expr_matrix(gse63514_expr_df, gene_col = "geneNames")
expr44001 <- prepare_expr_matrix(gse44001_expr_df, gene_col = "geneNames")

############################
# 7. Prepare metadata
############################
gse63514_meta <- gse63514_meta %>%
  mutate(
    group = factor(group, levels = progression_levels),
    group_order = as.numeric(group_order)
  )

common_63514_samples <- intersect(colnames(expr63514), gse63514_meta$sample_id)
expr63514 <- expr63514[, common_63514_samples, drop = FALSE]
gse63514_meta <- gse63514_meta %>%
  filter(sample_id %in% common_63514_samples) %>%
  arrange(group_order, group, sample_id)
expr63514 <- expr63514[, gse63514_meta$sample_id, drop = FALSE]

common_44001_score_samples <- intersect(colnames(expr44001), gse44001_score_df$sample)
score_meta_44001 <- gse44001_score_df %>%
  filter(sample %in% common_44001_score_samples)
expr44001_score <- expr44001[, score_meta_44001$sample, drop = FALSE]

common_44001_surv_samples <- intersect(colnames(expr44001), gse44001_surv_df$sample)
surv_meta_44001 <- gse44001_surv_df %>%
  filter(sample %in% common_44001_surv_samples)
expr44001_surv <- expr44001[, surv_meta_44001$sample, drop = FALSE]

############################
# 8. Compute GSE63514 proliferation_core_score by ssGSEA
############################
gene_sets <- list(proliferation_core_score = intersect(top50_genes, rownames(expr63514)))

if (length(gene_sets[[1]]) < 5) {
  stop("Too few signature genes matched in GSE63514 for stable ssGSEA calculation")
}

ssgsea_par_63514 <- ssgseaParam(
  exprData = expr63514,
  geneSets = gene_sets,
  alpha = 0.25,
  normalize = TRUE,
  minSize = 5,
  maxSize = 5000
)

score63514_mat <- gsva(ssgsea_par_63514, verbose = TRUE)

score63514_df <- data.frame(
  sample_id = colnames(score63514_mat),
  proliferation_core_score = as.numeric(score63514_mat[1, ]),
  stringsAsFactors = FALSE
)

score63514_df <- gse63514_meta %>%
  left_join(score63514_df, by = "sample_id")

############################
# 9. Build merged HMMR tables
############################
lead_gene_symbol_63514 <- find_gene_symbol(lead_gene, expr63514)
lead_gene_symbol_44001 <- find_gene_symbol(lead_gene, expr44001)

if (is.na(lead_gene_symbol_63514)) stop(paste0("Cannot find ", lead_gene, " in GSE63514 expression matrix"))
if (is.na(lead_gene_symbol_44001)) stop(paste0("Cannot find ", lead_gene, " in GSE44001 expression matrix"))

partner_63514 <- setNames(
  lapply(marker_partners, function(g) find_gene_symbol(g, expr63514)),
  marker_partners
)
partner_44001 <- setNames(
  lapply(marker_partners, function(g) find_gene_symbol(g, expr44001)),
  marker_partners
)

gse63514_hmmr <- gse63514_meta %>%
  left_join(score63514_df %>% select(sample_id, proliferation_core_score), by = "sample_id") %>%
  mutate(
    HMMR = as.numeric(expr63514[lead_gene_symbol_63514, sample_id]),
    group = factor(group, levels = progression_levels),
    group_order_plot = c(1, 2, 3, 4)[match(as.character(group), progression_levels)]
  )

for (g in marker_partners) {
  if (!is.na(partner_63514[[g]])) {
    gse63514_hmmr[[g]] <- as.numeric(expr63514[partner_63514[[g]], gse63514_hmmr$sample_id])
  } else {
    gse63514_hmmr[[g]] <- NA_real_
  }
}

write.csv(
  gse63514_hmmr,
  file.path(out_dir, "GSE63514_HMMR_merged.csv"),
  row.names = FALSE
)

merged_44001_meta <- surv_meta_44001 %>%
  left_join(score_meta_44001, by = "sample")

common_44001_merge_samples <- intersect(colnames(expr44001), merged_44001_meta$sample)

gse44001_hmmr <- merged_44001_meta %>%
  filter(sample %in% common_44001_merge_samples) %>%
  mutate(
    HMMR = as.numeric(expr44001[lead_gene_symbol_44001, sample]),
    HMMR_group = make_binary_group(HMMR, low_label = "Low HMMR", high_label = "High HMMR")
  )

for (g in marker_partners) {
  if (!is.na(partner_44001[[g]])) {
    gse44001_hmmr[[g]] <- as.numeric(expr44001[partner_44001[[g]], gse44001_hmmr$sample])
  } else {
    gse44001_hmmr[[g]] <- NA_real_
  }
}

write.csv(
  gse44001_hmmr,
  file.path(out_dir, "GSE44001_HMMR_merged.csv"),
  row.names = FALSE
)

############################
# 10. Save input logic
############################
writeLines(
  c(
    "Part 2 script inputs:",
    paste0("sig_file: ", sig_file),
    paste0("gse63514_expr_file: ", gse63514_expr_file),
    paste0("gse63514_meta_file: ", gse63514_meta_file),
    paste0("gse44001_expr_file: ", gse44001_expr_file),
    paste0("gse44001_score_meta_file: ", gse44001_score_meta_file),
    paste0("gse44001_surv_file: ", gse44001_surv_file),
    "",
    "Logic used in Part 2:",
    "GSE63514 uses top50 signature to compute proliferation_core_score by ssGSEA.",
    "GSE44001_clinical_score_metadata.csv is used only for proliferation_core_score correlation.",
    "gse44001_survival_metadata.txt is used for KM and Cox.",
    "",
    "Detected survival columns:",
    paste0("sample = ", surv_sample_col),
    paste0("dfs_months = ", surv_time_col),
    paste0("dfs_status = ", surv_status_col),
    paste0("stage = ", ifelse(is.na(surv_stage_col), "NA", surv_stage_col)),
    paste0("largest_diameter = ", ifelse(is.na(surv_diameter_col), "NA", surv_diameter_col))
  ),
  file.path(out_dir, "00_input_files_used_part2.txt")
)

############################
# 11. GSE63514 statistics
############################
hmmr_kw <- safe_kw(gse63514_hmmr$HMMR, gse63514_hmmr$group)

writeLines(
  c(
    paste0("Gene: ", lead_gene),
    paste0("Kruskal-Wallis statistic: ", fmt_num(hmmr_kw$stat, 3)),
    paste0("Kruskal-Wallis p: ", fmt_p(hmmr_kw$p))
  ),
  file.path(out_dir, "GSE63514_HMMR_kruskal.txt")
)

pairwise_order_df <- make_pairwise_order(progression_levels)

pairwise_hmmr <- ggpubr::compare_means(
  HMMR ~ group,
  data = gse63514_hmmr,
  method = "wilcox.test",
  p.adjust.method = "BH"
) %>%
  dplyr::select(group1, group2, p, p.adj, p.format, p.signif, method)

pairwise_hmmr <- pairwise_order_df %>%
  left_join(pairwise_hmmr, by = c("group1", "group2"))

write.table(
  pairwise_hmmr,
  file.path(out_dir, "GSE63514_HMMR_pairwise_wilcox.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################
# 12. GSE63514 plots
############################
y_max <- max(gse63514_hmmr$HMMR, na.rm = TRUE)
y_min <- min(gse63514_hmmr$HMMR, na.rm = TRUE)
y_range <- y_max - y_min
if (!is.finite(y_range) || y_range <= 0) y_range <- 1

pairwise_hmmr$y.position <- y_max + seq(0.10, by = 0.10, length.out = nrow(pairwise_hmmr)) * y_range

plot_corr_pub(
  df = gse63514_hmmr,
  xvar = "group_order_plot",
  yvar = "HMMR",
  title_text = "GSE63514: HMMR and progression order",
  xlab = "Progression group",
  ylab = "HMMR expression",
  out_pdf = file.path(out_dir, "GSE63514_HMMR_vs_group_order.pdf"),
  point_color = "#4E79A7",
  add_lm = TRUE,
  x_breaks = 1:4,
  x_labels = progression_levels
)

kw_subtitle <- paste0("Kruskal-Wallis, p = ", fmt_p(hmmr_kw$p))

p_hmmr_group <- ggplot(
  gse63514_hmmr,
  aes(x = group, y = HMMR, fill = group)
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.72,
    color = NA,
    width = 0.92
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.82,
    color = "black",
    linewidth = 0.9
  ) +
  geom_jitter(
    width = 0.10,
    size = 3.2,
    alpha = 0.68,
    color = "black"
  ) +
  stat_pvalue_manual(
    pairwise_hmmr,
    label = "p.signif",
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    tip.length = 0.01,
    step.increase = 0,
    size = 6
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#F3B5AF",
      "LSIL"   = "#B8CF88",
      "HSIL"   = "#8FD3D8",
      "Cancer" = "#D4B1EE"
    ),
    drop = FALSE
  ) +
  coord_cartesian(ylim = c(y_min, y_max + 0.78 * y_range)) +
  labs(
    title = "GSE63514: HMMR across progression groups",
    subtitle = kw_subtitle,
    x = NULL,
    y = "HMMR expression"
  ) +
  theme_pub(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

save_pdf(
  p_hmmr_group,
  file.path(out_dir, "GSE63514_HMMR_progression_violin.pdf"),
  width = 10.8,
  height = 8.8
)

plot_corr_pub(
  df = gse63514_hmmr,
  xvar = "HMMR",
  yvar = "proliferation_core_score",
  title_text = "GSE63514: HMMR and proliferation_core_score",
  xlab = "HMMR expression",
  ylab = "proliferation_core_score",
  out_pdf = file.path(out_dir, "GSE63514_HMMR_vs_proliferation_core_score.pdf"),
  point_color = "#4E79A7",
  add_lm = TRUE
)

for (g in marker_partners) {
  if (!all(is.na(gse63514_hmmr[[g]]))) {
    plot_corr_pub(
      df = gse63514_hmmr,
      xvar = "HMMR",
      yvar = g,
      title_text = paste0("GSE63514: HMMR and ", g),
      xlab = "HMMR expression",
      ylab = paste0(g, " expression"),
      out_pdf = file.path(out_dir, paste0("GSE63514_HMMR_vs_", g, ".pdf")),
      point_color = "#4E79A7",
      add_lm = TRUE
    )
  }
}

############################
# 13. GSE44001 plots
############################
plot_corr_pub(
  df = gse44001_hmmr,
  xvar = "HMMR",
  yvar = "proliferation_core_score",
  title_text = "GSE44001: HMMR and proliferation_core_score",
  xlab = "HMMR expression",
  ylab = "proliferation_core_score",
  out_pdf = file.path(out_dir, "GSE44001_HMMR_vs_proliferation_core_score.pdf"),
  point_color = "#4E79A7",
  add_lm = TRUE
)

for (g in marker_partners) {
  if (!all(is.na(gse44001_hmmr[[g]]))) {
    plot_corr_pub(
      df = gse44001_hmmr,
      xvar = "HMMR",
      yvar = g,
      title_text = paste0("GSE44001: HMMR and ", g),
      xlab = "HMMR expression",
      ylab = paste0(g, " expression"),
      out_pdf = file.path(out_dir, paste0("GSE44001_HMMR_vs_", g, ".pdf")),
      point_color = "#4E79A7",
      add_lm = TRUE
    )
  }
}

km_hmmr <- plot_km_pub(
  df = gse44001_hmmr,
  time_col = "dfs_years",
  status_col = "dfs_status",
  group_col = "HMMR_group",
  title_text = "GSE44001: HMMR disease-free survival",
  out_pdf = file.path(out_dir, "GSE44001_HMMR_KM.pdf")
)

cox_hmmr <- safe_cox(
  time = gse44001_hmmr$dfs_months,
  status = gse44001_hmmr$dfs_status,
  gene_exp = gse44001_hmmr$HMMR
)

cox_hmmr_df <- data.frame(
  gene = lead_gene,
  n = cox_hmmr$n,
  HR = cox_hmmr$hr,
  conf.low = cox_hmmr$low,
  conf.high = cox_hmmr$high,
  p.value = cox_hmmr$p,
  HR_CI = ifelse(
    is.na(cox_hmmr$hr),
    NA_character_,
    paste0(
      fmt_num(cox_hmmr$hr, 3), " (",
      fmt_num(cox_hmmr$low, 3), "–",
      fmt_num(cox_hmmr$high, 3), ")"
    )
  ),
  p.value.formatted = fmt_p(cox_hmmr$p),
  stringsAsFactors = FALSE
)

write.csv(
  cox_hmmr_df,
  file.path(out_dir, "GSE44001_HMMR_univariate_cox.csv"),
  row.names = FALSE
)

############################
# 14. Correlation summary table
############################
summary_list <- list()

res_63514_order <- safe_spearman(gse63514_hmmr$HMMR, gse63514_hmmr$group_order)
summary_list[[length(summary_list) + 1]] <- data.frame(
  dataset = "GSE63514",
  comparison = "HMMR_vs_group_order",
  rho = res_63514_order$rho,
  p.value = res_63514_order$p,
  n = res_63514_order$n,
  stringsAsFactors = FALSE
)

res_63514_core <- safe_spearman(gse63514_hmmr$HMMR, gse63514_hmmr$proliferation_core_score)
summary_list[[length(summary_list) + 1]] <- data.frame(
  dataset = "GSE63514",
  comparison = "HMMR_vs_proliferation_core_score",
  rho = res_63514_core$rho,
  p.value = res_63514_core$p,
  n = res_63514_core$n,
  stringsAsFactors = FALSE
)

res_44001_core <- safe_spearman(gse44001_hmmr$HMMR, gse44001_hmmr$proliferation_core_score)
summary_list[[length(summary_list) + 1]] <- data.frame(
  dataset = "GSE44001",
  comparison = "HMMR_vs_proliferation_core_score",
  rho = res_44001_core$rho,
  p.value = res_44001_core$p,
  n = res_44001_core$n,
  stringsAsFactors = FALSE
)

for (g in marker_partners) {
  if (!all(is.na(gse63514_hmmr[[g]]))) {
    res <- safe_spearman(gse63514_hmmr$HMMR, gse63514_hmmr[[g]])
    summary_list[[length(summary_list) + 1]] <- data.frame(
      dataset = "GSE63514",
      comparison = paste0("HMMR_vs_", g),
      rho = res$rho,
      p.value = res$p,
      n = res$n,
      stringsAsFactors = FALSE
    )
  }

  if (!all(is.na(gse44001_hmmr[[g]]))) {
    res <- safe_spearman(gse44001_hmmr$HMMR, gse44001_hmmr[[g]])
    summary_list[[length(summary_list) + 1]] <- data.frame(
      dataset = "GSE44001",
      comparison = paste0("HMMR_vs_", g),
      rho = res$rho,
      p.value = res$p,
      n = res$n,
      stringsAsFactors = FALSE
    )
  }
}

hmmr_corr_summary <- bind_rows(summary_list) %>%
  mutate(
    rho_fmt = fmt_num(rho, 3),
    p_fmt = fmt_p(p.value)
  )

write.csv(
  hmmr_corr_summary,
  file.path(out_dir, "HMMR_correlation_summary.csv"),
  row.names = FALSE
)

cat("Part 2 finished successfully.\n")
cat("Outputs saved in: ", out_dir, "\n")