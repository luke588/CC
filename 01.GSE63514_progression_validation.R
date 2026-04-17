rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
options(stringsAsFactors = FALSE)

############################
# 0. Packages
############################
cran_pkgs <- c(
  "data.table", "dplyr", "tibble", "ggplot2", "ggpubr",
  "reshape2", "tidyr"
)

bioc_pkgs <- c("GSVA", "limma")

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
library(ggpubr)
library(GSVA)
library(limma)
library(reshape2)
library(tidyr)

############################
# 1. Paths
############################
base_dir <- "/Users/donglinlu/Desktop/CC1/7.progression validation/GSE63514"
setwd(base_dir)

expr_file  <- file.path(base_dir, "geneMatrix.txt")
meta_file  <- file.path(base_dir, "GSE63514_metadata.txt")

sig6_file  <- file.path(base_dir, "cluster6_proliferation_core_epithelial_top50.txt")
sig9_file  <- file.path(base_dir, "cluster9_malignant_squamous_epithelial_top50.txt")
sig10_file <- file.path(base_dir, "cluster10_EMT_like_epithelial_top50.txt")

out_dir <- base_dir

############################
# 2. Read expression matrix
############################
expr_df <- fread(expr_file, data.table = FALSE)

if (!"geneNames" %in% colnames(expr_df)) {
  stop("表达矩阵第一列列名必须是 geneNames")
}

colnames(expr_df)[colnames(expr_df) == "geneNames"] <- "GeneSymbol"

expr_df <- expr_df %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "")

sample_cols <- setdiff(colnames(expr_df), "GeneSymbol")
expr_df[, sample_cols] <- lapply(expr_df[, sample_cols, drop = FALSE], as.numeric)

expr_df <- expr_df[
  rowSums(is.na(expr_df[, sample_cols, drop = FALSE])) < length(sample_cols),
]

expr_df <- expr_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

expr_mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
rownames(expr_mat) <- expr_df$GeneSymbol

cat("Expression matrix dimension:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

############################
# 3. Read metadata
############################
meta <- fread(meta_file, data.table = FALSE)

required_meta_cols <- c("sample_id", "group", "group_order")
if (!all(required_meta_cols %in% colnames(meta))) {
  stop("metadata 文件必须包含列: sample_id, group, group_order")
}

meta$group <- factor(meta$group, levels = c("Normal", "LSIL", "HSIL", "Cancer"))
meta$group_order <- as.numeric(meta$group_order)

common_samples <- intersect(colnames(expr_mat), meta$sample_id)

if (length(common_samples) == 0) {
  stop("表达矩阵列名与 metadata$sample_id 没有匹配到任何样本，请检查样本名")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
meta <- meta %>% filter(sample_id %in% common_samples)

meta <- meta %>% arrange(group_order, group, sample_id)
expr_mat <- expr_mat[, meta$sample_id, drop = FALSE]

cat("Matched samples:", ncol(expr_mat), "\n")
cat("Group counts:\n")
print(table(meta$group))

############################
# 4. Read signatures
############################
read_signature <- function(file) {
  x <- fread(file, header = FALSE, data.table = FALSE)[, 1]
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & x != ""]
  x
}

sig6  <- read_signature(sig6_file)
sig9  <- read_signature(sig9_file)
sig10 <- read_signature(sig10_file)

# Canonical order fixed as 9 -> 10 -> 6
score_names <- c(
  "malignant_squamous_score",
  "EMT_like_score",
  "proliferation_core_score"
)

gene_sets <- list(
  malignant_squamous_score = intersect(sig9, rownames(expr_mat)),
  EMT_like_score = intersect(sig10, rownames(expr_mat)),
  proliferation_core_score = intersect(sig6, rownames(expr_mat))
)

sig_overlap <- data.frame(
  signature = names(gene_sets),
  n_genes_in_signature = c(length(sig9), length(sig10), length(sig6)),
  n_genes_matched = sapply(gene_sets, length)
)

write.csv(
  sig_overlap,
  file = file.path(out_dir, "signature_overlap_summary.csv"),
  row.names = FALSE
)

print(sig_overlap)

if (any(sig_overlap$n_genes_matched < 5)) {
  warning("某些 signature 命中的基因少于5个，请检查基因名格式")
}

############################
# 5. ssGSEA scoring
############################
ssgsea_par <- ssgseaParam(
  exprData = expr_mat,
  geneSets = gene_sets,
  alpha = 0.25,
  normalize = TRUE,
  minSize = 5,
  maxSize = 5000
)

score_mat <- tryCatch(
  {
    gsva(ssgsea_par, verbose = TRUE)
  },
  error = function(e) {
    stop("ssGSEA failed: ", e$message)
  }
)

score_df <- as.data.frame(t(score_mat))
score_df$sample_id <- rownames(score_df)

score_df <- meta %>%
  left_join(score_df, by = "sample_id")

write.csv(
  score_df,
  file.path(out_dir, "GSE63514_ssGSEA_scores.csv"),
  row.names = FALSE
)

############################
# 6. Statistics
############################
kw_res_list <- list()
trend_res_list <- list()
pairwise_res_list <- list()

for (sc in score_names) {
  kw <- kruskal.test(as.formula(paste(sc, "~ group")), data = score_df)
  kw_res_list[[sc]] <- data.frame(
    signature = sc,
    method = "Kruskal-Wallis",
    statistic = as.numeric(kw$statistic),
    p_value = kw$p.value
  )

  sp <- suppressWarnings(cor.test(score_df[[sc]], score_df$group_order, method = "spearman"))
  trend_res_list[[sc]] <- data.frame(
    signature = sc,
    method = "Spearman trend",
    rho = as.numeric(sp$estimate),
    p_value = sp$p.value
  )

  pw <- pairwise.wilcox.test(
    x = score_df[[sc]],
    g = score_df$group,
    p.adjust.method = "BH"
  )

  pw_df <- as.data.frame(as.table(pw$p.value))
  colnames(pw_df) <- c("group1", "group2", "p_adj")
  pw_df$signature <- sc
  pairwise_res_list[[sc]] <- pw_df
}

kw_res <- bind_rows(kw_res_list)
trend_res <- bind_rows(trend_res_list)
pairwise_res <- bind_rows(pairwise_res_list)

write.csv(
  kw_res,
  file.path(out_dir, "kruskal_wallis_results.csv"),
  row.names = FALSE
)

write.csv(
  trend_res,
  file.path(out_dir, "trend_test_results.csv"),
  row.names = FALSE
)

write.csv(
  pairwise_res,
  file.path(out_dir, "pairwise_wilcox_results.csv"),
  row.names = FALSE
)

############################
# 7. Display helpers
############################
display_title <- function(score_col) {
  dplyr::case_when(
    score_col == "malignant_squamous_score" ~ "Malignant squamous",
    score_col == "EMT_like_score" ~ "EMT-like",
    score_col == "proliferation_core_score" ~ "Proliferation core",
    TRUE ~ score_col
  )
}

display_ylabel <- function(score_col) {
  dplyr::case_when(
    score_col == "malignant_squamous_score" ~ "Malignant squamous score",
    score_col == "EMT_like_score" ~ "EMT-like score",
    score_col == "proliferation_core_score" ~ "Proliferation core score",
    TRUE ~ score_col
  )
}

get_kw_p <- function(score_col, kw_table) {
  kw_table$p_value[match(score_col, kw_table$signature)]
}

progression_fill <- c(
  "Normal" = "#8EC5E8",
  "LSIL"   = "#9FD3A8",
  "HSIL"   = "#F2B36D",
  "Cancer" = "#E97A6E"
)

############################
# 8. Plot function
############################
plot_signature_box <- function(df, score_col, out_file, kw_table) {
  comp_list <- list(
    c("Normal", "LSIL"),
    c("LSIL", "HSIL"),
    c("HSIL", "Cancer"),
    c("Normal", "Cancer")
  )

  y_max <- max(df[[score_col]], na.rm = TRUE)
  y_min <- min(df[[score_col]], na.rm = TRUE)
  y_range <- y_max - y_min
  if (!is.finite(y_range) || y_range <= 0) y_range <- 1

  kw_p <- get_kw_p(score_col, kw_table)
  subtitle_text <- paste0("Kruskal–Wallis p = ", fmt_p(kw_p))

  p <- ggplot(df, aes(x = group, y = .data[[score_col]], fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.38, color = NA) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      alpha = 0.90,
      color = "black",
      linewidth = 0.75
    ) +
    geom_jitter(
      width = 0.12,
      size = 1.7,
      alpha = 0.65,
      color = "black"
    ) +
    stat_compare_means(
      comparisons = comp_list,
      method = "wilcox.test",
      label = "p.signif",
      step.increase = 0.11,
      size = 5.6
    ) +
    scale_fill_manual(values = progression_fill, drop = FALSE) +
    coord_cartesian(
      ylim = c(y_min - 0.03 * y_range, y_max + 0.34 * y_range),
      clip = "off"
    ) +
    labs(
      x = NULL,
      y = display_ylabel(score_col),
      title = display_title(score_col),
      subtitle = subtitle_text
    ) +
    theme_pub() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = PUB_AXIS_TEXT_SIZE),
      axis.text.y = element_text(size = PUB_AXIS_TEXT_SIZE),
      plot.margin = margin(t = 18, r = 16, b = 12, l = 14)
    )

  ggsave(out_file, plot = p, width = 7, height = 6)
  return(p)
}

############################
# 9. Individual plots
############################
plot_signature_box(
  score_df,
  "malignant_squamous_score",
  file.path(out_dir, "GSE63514_malignant_squamous_score_by_group.pdf"),
  kw_res
)

plot_signature_box(
  score_df,
  "EMT_like_score",
  file.path(out_dir, "GSE63514_EMT_like_score_by_group.pdf"),
  kw_res
)

plot_signature_box(
  score_df,
  "proliferation_core_score",
  file.path(out_dir, "GSE63514_proliferation_core_score_by_group.pdf"),
  kw_res
)

############################
# 10. Combined plot
############################
score_long <- score_df %>%
  select(sample_id, group, group_order, all_of(score_names)) %>%
  pivot_longer(
    cols = all_of(score_names),
    names_to = "signature",
    values_to = "score"
  ) %>%
  mutate(
    signature = factor(signature, levels = score_names),
    program = dplyr::case_when(
      signature == "malignant_squamous_score" ~ "Malignant squamous",
      signature == "EMT_like_score" ~ "EMT-like",
      signature == "proliferation_core_score" ~ "Proliferation core",
      TRUE ~ as.character(signature)
    )
  )

facet_lab_df <- kw_res %>%
  mutate(
    signature = factor(signature, levels = score_names),
    program = dplyr::case_when(
      signature == "malignant_squamous_score" ~ "Malignant squamous",
      signature == "EMT_like_score" ~ "EMT-like",
      signature == "proliferation_core_score" ~ "Proliferation core",
      TRUE ~ as.character(signature)
    ),
    facet_label = paste0(program, "\nKruskal–Wallis p = ", fmt_p(p_value))
  ) %>%
  arrange(signature)

score_long <- score_long %>%
  left_join(facet_lab_df[, c("signature", "facet_label")], by = "signature")

score_long$facet_label <- factor(
  score_long$facet_label,
  levels = facet_lab_df$facet_label
)

p_all <- ggplot(score_long, aes(x = group, y = score, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.38, color = NA) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.90,
    color = "black",
    linewidth = 0.75
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.0,
    alpha = 0.60,
    color = "black"
  ) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = progression_fill, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
  labs(
    x = NULL,
    y = "ssGSEA score",
    title = NULL
  ) +
  theme_pub() +
  theme(
    legend.position = "none",
    strip.text = element_text(
      face = "bold",
      size = PUB_STRIP_TEXT,
      lineheight = 1.08
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = PUB_AXIS_TEXT_SIZE
    ),
    axis.text.y = element_text(size = PUB_AXIS_TEXT_SIZE),
    plot.margin = margin(t = 10, r = 16, b = 12, l = 14)
  )

ggsave(
  file.path(out_dir, "GSE63514_progression_validation_combined.pdf"),
  plot = p_all,
  width = 14,
  height = 5.8
)

############################
# 11. Group summary
############################
group_summary <- score_df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    malignant_squamous_mean = mean(malignant_squamous_score, na.rm = TRUE),
    malignant_squamous_sd   = sd(malignant_squamous_score, na.rm = TRUE),
    EMT_like_mean           = mean(EMT_like_score, na.rm = TRUE),
    EMT_like_sd             = sd(EMT_like_score, na.rm = TRUE),
    proliferation_core_mean = mean(proliferation_core_score, na.rm = TRUE),
    proliferation_core_sd   = sd(proliferation_core_score, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  group_summary,
  file.path(out_dir, "group_score_summary.csv"),
  row.names = FALSE
)

############################
# 12. Session info
############################
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo_GSE63514.txt")
)

cat("All analysis finished successfully.\n")
cat("Results saved in:", out_dir, "\n")