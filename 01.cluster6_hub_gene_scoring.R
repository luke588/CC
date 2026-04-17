rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
options(stringsAsFactors = FALSE)

############################
# 0. Packages
############################
cran_pkgs <- c(
  "data.table", "dplyr", "tibble", "ggplot2", "tidyr",
  "survival", "ggrepel"
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
library(ggrepel)
library(GSVA)

############################
# 1. Paths
############################
base_dir <- "/Users/donglinlu/Desktop/CC1/8.Hub Gene"
setwd(base_dir)

marker_file              <- file.path(base_dir, "cluster6_marker_table.csv")
sig_file                 <- file.path(base_dir, "cluster6_proliferation_core_epithelial_top50.txt")
gse63514_expr_file       <- file.path(base_dir, "GSE63514_geneMatrix.txt")
gse63514_meta_file       <- file.path(base_dir, "GSE63514_metadata.txt")
gse44001_expr_file       <- file.path(base_dir, "GSE44001_geneMatrix.txt")
gse44001_score_meta_file <- file.path(base_dir, "GSE44001_clinical_score_metadata.csv")
gse44001_surv_file       <- file.path(base_dir, "gse44001_survival_metadata.txt")

out_dir <- file.path(base_dir, "output_hub_gene")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

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
  return(x)
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
  return(expr_mat)
}

find_first_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
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
      plot.margin = margin(t = 16, r = 40, b = 14, l = 16)
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

############################
# 3. Read files
############################
marker_df <- fread(marker_file, data.table = FALSE)
top50_genes <- read_signature(sig_file)

gse63514_expr_df  <- fread(gse63514_expr_file, data.table = FALSE)
gse63514_meta     <- fread(gse63514_meta_file, data.table = FALSE)
gse44001_expr_df  <- fread(gse44001_expr_file, data.table = FALSE)
gse44001_score_df <- fread(gse44001_score_meta_file, data.table = FALSE)
gse44001_surv_df  <- fread(gse44001_surv_file, data.table = FALSE)

############################
# 4. Check headers
############################
required_marker_cols <- c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
if (!all(required_marker_cols %in% colnames(marker_df))) {
  stop("cluster6_marker_table.csv is missing required columns: gene, avg_log2FC, pct.1, pct.2, p_val_adj")
}

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
    group = factor(group, levels = c("Normal", "LSIL", "HSIL", "Cancer")),
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

write.csv(score63514_df, file.path(out_dir, "GSE63514_proliferation_core_scores.csv"), row.names = FALSE)

############################
# 9. Build candidate gene table
############################
marker_df <- marker_df %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(desc(avg_log2FC)) %>%
  mutate(
    marker_rank = row_number(),
    pct_delta = pct.1 - pct.2,
    in_top50 = ifelse(gene %in% top50_genes, 1, 0),
    in_top10 = ifelse(gene %in% top50_genes[1:min(10, length(top50_genes))], 1, 0)
  )

candidate_genes <- marker_df$gene
candidate_genes <- candidate_genes[
  candidate_genes %in% rownames(expr63514) &
    candidate_genes %in% rownames(expr44001)
]

candidate_info <- marker_df %>%
  filter(gene %in% candidate_genes)

write.csv(candidate_info, file.path(out_dir, "cluster6_candidate_genes_intersected.csv"), row.names = FALSE)

############################
# 10. Gene-level scoring
############################
res_list <- lapply(candidate_genes, function(g) {

  sc_row <- candidate_info[candidate_info$gene == g, , drop = FALSE]

  gene_63514 <- as.numeric(expr63514[g, gse63514_meta$sample_id])
  kw_63514   <- safe_kw(gene_63514, gse63514_meta$group)
  tr_63514   <- safe_spearman(gene_63514, gse63514_meta$group_order)
  cor_63514  <- safe_spearman(gene_63514, score63514_df$proliferation_core_score)

  gene_44001_cor <- as.numeric(expr44001_score[g, score_meta_44001$sample])
  cor_44001      <- safe_spearman(gene_44001_cor, score_meta_44001$proliferation_core_score)

  gene_44001_surv_gene <- as.numeric(expr44001_surv[g, surv_meta_44001$sample])
  cox_44001 <- safe_cox(
    time = surv_meta_44001$dfs_months,
    status = surv_meta_44001$dfs_status,
    gene_exp = gene_44001_surv_gene
  )

  data.frame(
    gene = g,
    marker_rank = sc_row$marker_rank,
    p_val = sc_row$p_val,
    avg_log2FC = sc_row$avg_log2FC,
    pct.1 = sc_row$pct.1,
    pct.2 = sc_row$pct.2,
    pct_delta = sc_row$pct_delta,
    p_val_adj = sc_row$p_val_adj,
    in_top50 = sc_row$in_top50,
    in_top10 = sc_row$in_top10,
    GSE63514_kw_p = kw_63514$p,
    GSE63514_trend_rho = tr_63514$rho,
    GSE63514_trend_p = tr_63514$p,
    GSE63514_cor_score_rho = cor_63514$rho,
    GSE63514_cor_score_p = cor_63514$p,
    GSE44001_cor_score_rho = cor_44001$rho,
    GSE44001_cor_score_p = cor_44001$p,
    GSE44001_cox_HR = cox_44001$hr,
    GSE44001_cox_p = cox_44001$p,
    stringsAsFactors = FALSE
  )
})

hub_score_df <- bind_rows(res_list)

############################
# 11. Rule-based scoring system
############################
hub_score_df <- hub_score_df %>%
  mutate(
    score_sc_logFC = ifelse(!is.na(avg_log2FC) & avg_log2FC >= 2, 1, 0),
    score_sc_pct_delta = ifelse(!is.na(pct_delta) & pct_delta >= 0.30, 1, 0),
    score_top50 = ifelse(in_top50 == 1, 1, 0),
    score_top10 = ifelse(in_top10 == 1, 1, 0),
    score_63514_kw = ifelse(!is.na(GSE63514_kw_p) & GSE63514_kw_p < 0.05, 1, 0),
    score_63514_trend = ifelse(!is.na(GSE63514_trend_rho) & !is.na(GSE63514_trend_p) & GSE63514_trend_rho > 0 & GSE63514_trend_p < 0.05, 1, 0),
    score_63514_cor = ifelse(!is.na(GSE63514_cor_score_rho) & !is.na(GSE63514_cor_score_p) & GSE63514_cor_score_rho > 0.30 & GSE63514_cor_score_p < 0.05, 1, 0),
    score_44001_cor = ifelse(!is.na(GSE44001_cor_score_rho) & !is.na(GSE44001_cor_score_p) & GSE44001_cor_score_rho > 0.30 & GSE44001_cor_score_p < 0.05, 1, 0),
    score_44001_cox = ifelse(!is.na(GSE44001_cox_HR) & !is.na(GSE44001_cox_p) & GSE44001_cox_HR > 1 & GSE44001_cox_p < 0.05, 1, 0)
  ) %>%
  mutate(
    total_score = score_sc_logFC +
      score_sc_pct_delta +
      score_top50 +
      score_top10 +
      score_63514_kw +
      score_63514_trend +
      score_63514_cor +
      score_44001_cor +
      score_44001_cox
  ) %>%
  arrange(
    desc(total_score),
    desc(score_44001_cox),
    desc(score_63514_cor),
    desc(score_44001_cor),
    desc(avg_log2FC)
  )

write.csv(hub_score_df, file.path(out_dir, "cluster6_hub_gene_scoring_table.csv"), row.names = FALSE)

############################
# 12. Export top tables
############################
top20_df <- hub_score_df %>% slice_head(n = 20)
top10_df <- hub_score_df %>% slice_head(n = 10)
top5_df  <- hub_score_df %>% slice_head(n = 5)

write.csv(top20_df, file.path(out_dir, "cluster6_top20_hub_genes.csv"), row.names = FALSE)
write.csv(top10_df, file.path(out_dir, "cluster6_top10_hub_genes.csv"), row.names = FALSE)
write.csv(top5_df,  file.path(out_dir, "cluster6_top5_hub_genes.csv"), row.names = FALSE)

summary_df <- data.frame(
  n_marker_input = nrow(marker_df),
  n_top50_signature = length(top50_genes),
  n_candidate_after_intersection = length(candidate_genes),
  n_top20 = nrow(top20_df),
  n_top10 = nrow(top10_df),
  n_top5 = nrow(top5_df)
)
write.csv(summary_df, file.path(out_dir, "cluster6_hub_gene_summary.csv"), row.names = FALSE)

############################
# 13. Plot 1: hub gene ranking barplot
############################
plot_bar_df <- top20_df %>%
  mutate(gene = factor(gene, levels = rev(gene)))

p_bar <- ggplot(plot_bar_df, aes(x = gene, y = total_score)) +
  geom_col(fill = "#5B8FD1", width = 0.72) +
  coord_flip() +
  labs(
    title = "Cluster 6 hub gene prioritization",
    subtitle = "Integrated evidence from single-cell and external cohorts",
    x = NULL,
    y = "Total score"
  ) +
  theme_pub(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(face = "bold")
  )

save_pdf(
  p_bar,
  file.path(out_dir, "cluster6_top20_hub_gene_barplot.pdf"),
  width = 8.2,
  height = 8.6
)

############################
# 14. Plot 2: bubble plot
# top20 circles + all 20 gene labels
############################
bubble_df <- top20_df %>%
  mutate(
    neglog10_cox_p = -log10(GSE44001_cox_p),
    neglog10_cox_p = ifelse(is.infinite(neglog10_cox_p), NA, neglog10_cox_p)
  )

x_min <- min(bubble_df$GSE63514_cor_score_rho, na.rm = TRUE)
x_max <- max(bubble_df$GSE63514_cor_score_rho, na.rm = TRUE)
y_min <- min(bubble_df$GSE44001_cor_score_rho, na.rm = TRUE)
y_max <- max(bubble_df$GSE44001_cor_score_rho, na.rm = TRUE)

x_pad <- (x_max - x_min) * 0.35
y_pad <- (y_max - y_min) * 0.30

if (!is.finite(x_pad) || x_pad == 0) x_pad <- 0.08
if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.08

p_bubble <- ggplot(
  bubble_df,
  aes(
    x = GSE63514_cor_score_rho,
    y = GSE44001_cor_score_rho,
    size = total_score,
    color = neglog10_cox_p
  )
) +
  geom_point(alpha = 0.92) +
  ggrepel::geom_text_repel(
    aes(label = gene),
    size = 4.7,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.65,
    point.padding = 0.35,
    force = 12,
    force_pull = 0.2,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.size = 0.45,
    seed = 123,
    show.legend = FALSE
  ) +
  scale_color_gradient(low = "#9FD3C7", high = "#D95F5F", na.value = "grey70") +
  scale_size_continuous(range = c(2.6, 8.2)) +
  scale_x_continuous(
    limits = c(x_min - 0.03, x_max + x_pad),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(y_min - 0.03, y_max + y_pad),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Top20 hub genes: cross-cohort concordance",
    x = "GSE63514",
    y = "GSE44001",
    color = expression(-log[10]("Cox p")),
    size = "Total score"
  ) +
  theme_pub(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.margin = margin(t = 18, r = 85, b = 14, l = 18)
  )

save_pdf(
  p_bubble,
  file.path(out_dir, "cluster6_top20_hub_gene_bubbleplot.pdf"),
  width = 12.2,
  height = 8.6
)

############################
# 15. Plot 3: top20 metric heatmap
############################
heat_df <- top20_df %>%
  select(
    gene,
    avg_log2FC,
    pct_delta,
    GSE63514_trend_rho,
    GSE63514_cor_score_rho,
    GSE44001_cor_score_rho,
    GSE44001_cox_HR,
    total_score
  ) %>%
  mutate(
    gene = factor(gene, levels = rev(gene))
  )

heat_mat <- heat_df %>%
  column_to_rownames("gene") %>%
  as.matrix()

heat_mat_scaled <- scale(heat_mat)

heat_long <- as.data.frame(heat_mat_scaled) %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "metric",
    values_to = "zscore"
  ) %>%
  mutate(
    gene = factor(gene, levels = rev(rownames(heat_mat_scaled))),
    metric = factor(
      metric,
      levels = colnames(heat_mat_scaled),
      labels = c(
        "SC_logFC",
        "SC_pctΔ",
        "63514_trend",
        "63514_scoreCor",
        "44001_scoreCor",
        "44001_CoxHR",
        "Total"
      )
    )
  )

p_heat <- ggplot(heat_long, aes(x = metric, y = gene, fill = zscore)) +
  geom_tile(color = "grey70", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#5B8FD1",
    mid = "#F3F2D0",
    high = "#E53935",
    midpoint = 0,
    name = "Z score"
  ) +
  labs(
    title = "Top20 hub gene metric heatmap",
    subtitle = "Scaled evidence across prioritization metrics",
    x = NULL,
    y = NULL
  ) +
  theme_pub(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )

save_pdf(
  p_heat,
  file.path(out_dir, "cluster6_top20_metric_heatmap.pdf"),
  width = 9.2,
  height = 9.4
)

############################
# 16. Save manuscript-ready table
############################
manuscript_ready_df <- hub_score_df %>%
  mutate(
    GSE63514_kw_p_fmt = fmt_p(GSE63514_kw_p),
    GSE63514_trend_p_fmt = fmt_p(GSE63514_trend_p),
    GSE63514_cor_score_p_fmt = fmt_p(GSE63514_cor_score_p),
    GSE44001_cor_score_p_fmt = fmt_p(GSE44001_cor_score_p),
    GSE44001_cox_p_fmt = fmt_p(GSE44001_cox_p)
  )

write.csv(
  manuscript_ready_df,
  file.path(out_dir, "cluster6_hub_gene_scoring_table_manuscript_ready.csv"),
  row.names = FALSE
)

############################
# 17. Save input file record
############################
writeLines(
  c(
    "Part 1 script inputs:",
    paste0("marker_file: ", marker_file),
    paste0("sig_file: ", sig_file),
    paste0("gse63514_expr_file: ", gse63514_expr_file),
    paste0("gse63514_meta_file: ", gse63514_meta_file),
    paste0("gse44001_expr_file: ", gse44001_expr_file),
    paste0("gse44001_score_meta_file: ", gse44001_score_meta_file),
    paste0("gse44001_surv_file: ", gse44001_surv_file),
    "",
    "Logic used in Part 1:",
    "GSE63514 is used for progression trend and score correlation.",
    "GSE44001_clinical_score_metadata.csv is used for proliferation_core_score correlation.",
    "gse44001_survival_metadata.txt is used for Cox contribution in the ranking."
  ),
  file.path(out_dir, "00_input_files_used_part1.txt")
)

cat("Part 1 finished successfully.\n")
cat("Outputs saved in: ", out_dir, "\n")