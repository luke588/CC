# ==========================================
# TCGA-CESC minimal validation (clean version)
# proliferation_core_score + HMMR + 3 hub genes
# ==========================================

rm(list = ls())
source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")
options(stringsAsFactors = FALSE)

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(GSVA)

base_dir <- "/Users/donglinlu/Desktop/CC1/10.TCGA_CESC"
setwd(base_dir)

out_dir <- file.path(base_dir, "TCGA_results")
dir.create(out_dir, showWarnings = FALSE)

# -----------------------------
# 0. helper functions
# -----------------------------
safe_extract_gene <- function(expr_df, gene_symbol, gene_col = "geneNames") {
  sub <- expr_df[expr_df[[gene_col]] == gene_symbol, , drop = FALSE]
  if (nrow(sub) == 0) stop(paste0("Gene not found: ", gene_symbol))
  if (nrow(sub) > 1) {
    message("Multiple rows found for ", gene_symbol, ", keeping the first one.")
    sub <- sub[1, , drop = FALSE]
  }
  return(sub)
}

expr_row_to_long <- function(expr_row, gene_symbol, gene_col = "geneNames") {
  x <- expr_row[, setdiff(colnames(expr_row), gene_col), drop = FALSE]
  data.frame(
    sample = colnames(x),
    expr = as.numeric(x[1, ]),
    gene = gene_symbol,
    stringsAsFactors = FALSE
  )
}

make_cor_plot <- function(df, xvar, yvar, title_txt, xlab_txt, ylab_txt, out_file) {
  df2 <- df %>%
    dplyr::filter(!is.na(.data[[xvar]]), !is.na(.data[[yvar]]))

  cor_res <- suppressWarnings(cor.test(df2[[xvar]], df2[[yvar]], method = "spearman"))
  rho_txt <- sprintf("Spearman rho = %.3f\np = %.3g", cor_res$estimate, cor_res$p.value)

  p <- ggplot(df2, aes_string(x = xvar, y = yvar)) +
    geom_point(size = 2.3, alpha = 0.8, color = "#4C78A8") +
    geom_smooth(method = "lm", se = TRUE, color = "#D62728", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, label = rho_txt,
             hjust = 1.1, vjust = 1.5, size = 4.3) +
    labs(title = title_txt, x = xlab_txt, y = ylab_txt) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(out_file, p, width = 6.2, height = 5.2)
}

make_group_plot <- function(df, group_var, yvar, title_txt, ylab_txt, out_file, group_levels = NULL) {
  df2 <- df %>%
    dplyr::filter(!is.na(.data[[group_var]]), !is.na(.data[[yvar]]))

  if (nrow(df2) == 0) {
    warning(paste("No data left for", group_var, "vs", yvar))
    return(NULL)
  }

  if (!is.null(group_levels)) {
    df2[[group_var]] <- factor(df2[[group_var]], levels = group_levels)
    df2 <- df2 %>% dplyr::filter(!is.na(.data[[group_var]]))
  }

  if (length(unique(df2[[group_var]])) < 2) {
    warning(paste("Not enough groups for", group_var))
    return(NULL)
  }

  kw_res <- kruskal.test(df2[[yvar]] ~ df2[[group_var]])

  p <- ggplot(df2, aes_string(x = group_var, y = yvar, fill = group_var)) +
    geom_violin(trim = FALSE, alpha = 0.45, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 1.5, alpha = 0.6) +
    annotate("text",
             x = 1.1,
             y = max(df2[[yvar]], na.rm = TRUE) * 1.03,
             label = paste0("Kruskal-Wallis p = ", format(kw_res$p.value, digits = 3)),
             hjust = 0, size = 4.5) +
    labs(title = title_txt, x = NULL, y = ylab_txt) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  ggsave(out_file, p, width = 7.0, height = 5.6)
}

# -----------------------------
# 1. read files
# -----------------------------
expr_df <- fread("TCGA_CESC_geneMatrix_clean.txt", data.table = FALSE)
cli_df  <- fread("TCGA_CESC_clinical_clean.txt", data.table = FALSE)
sig_df  <- fread("cluster6_proliferation_core_epithelial_top50.txt", header = FALSE, data.table = FALSE)

colnames(expr_df)[1] <- "geneNames"
colnames(sig_df)[1] <- "gene"

required_cli <- c("sample", "OS.time", "OS", "stage", "T_stage", "N_stage")
missing_cli <- setdiff(required_cli, colnames(cli_df))
if (length(missing_cli) > 0) {
  stop(paste("Missing clinical columns:", paste(missing_cli, collapse = ", ")))
}

# -----------------------------
# 2. calculate proliferation_core_score by ssGSEA
# -----------------------------
sig_genes <- unique(sig_df$gene)
sig_genes <- sig_genes[sig_genes %in% expr_df$geneNames]

cat("Signature genes found in TCGA matrix:", length(sig_genes), "\n")

expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
rownames(expr_mat) <- expr_df$geneNames
mode(expr_mat) <- "numeric"

gene_sets <- list(proliferation_core = sig_genes)

# GSVA新版写法
par_obj <- ssgseaParam(
  exprData = expr_mat,
  geneSets = gene_sets,
  alpha = 0.25,
  normalize = TRUE
)

score_mat <- gsva(par_obj, verbose = FALSE)

score_df <- data.frame(
  sample = colnames(score_mat),
  proliferation_core_score = as.numeric(score_mat["proliferation_core", ]),
  stringsAsFactors = FALSE
)

tcga_df <- cli_df %>%
  left_join(score_df, by = "sample")

write.csv(tcga_df, file.path(out_dir, "TCGA_CESC_ssGSEA_scores.csv"), row.names = FALSE)

# -----------------------------
# 3. extract HMMR and classic proliferation genes
# -----------------------------
genes_interest <- c("HMMR", "NUF2", "CDKN3", "MKI67", "TOP2A", "UBE2C", "BIRC5")

for (g in genes_interest) {
  g_row <- safe_extract_gene(expr_df, g, gene_col = "geneNames")
  g_long <- expr_row_to_long(g_row, g, gene_col = "geneNames")
  colnames(g_long)[2] <- paste0(g, "_expr")
  tcga_df <- tcga_df %>% left_join(g_long[, c("sample", paste0(g, "_expr"))], by = "sample")
}

write.csv(tcga_df, file.path(out_dir, "TCGA_CESC_clinical_with_scores_and_genes.csv"), row.names = FALSE)

# -----------------------------
# 4. Task 1: proliferation_core_score relations
# -----------------------------
make_cor_plot(
  tcga_df, "HMMR_expr", "proliferation_core_score",
  "TCGA-CESC: HMMR vs proliferation_core_score",
  "HMMR expression", "proliferation_core_score",
  file.path(out_dir, "TCGA_HMMR_vs_proliferation_core_score.pdf")
)

for (g in c("MKI67", "TOP2A", "UBE2C", "BIRC5")) {
  make_cor_plot(
    tcga_df, paste0(g, "_expr"), "proliferation_core_score",
    paste0("TCGA-CESC: ", g, " vs proliferation_core_score"),
    paste0(g, " expression"), "proliferation_core_score",
    file.path(out_dir, paste0("TCGA_", g, "_vs_proliferation_core_score.pdf"))
  )
}

make_group_plot(
  tcga_df, "stage", "proliferation_core_score",
  "TCGA-CESC: proliferation_core_score by stage",
  "proliferation_core_score",
  file.path(out_dir, "TCGA_proliferation_core_score_by_stage.pdf")
)

make_group_plot(
  tcga_df, "T_stage", "proliferation_core_score",
  "TCGA-CESC: proliferation_core_score by T stage",
  "proliferation_core_score",
  file.path(out_dir, "TCGA_proliferation_core_score_by_T_stage.pdf")
)

make_group_plot(
  tcga_df, "N_stage", "proliferation_core_score",
  "TCGA-CESC: proliferation_core_score by N stage",
  "proliferation_core_score",
  file.path(out_dir, "TCGA_proliferation_core_score_by_N_stage.pdf")
)

# -----------------------------
# 5. Task 2: HMMR validation in TCGA
# -----------------------------
for (g in c("MKI67", "TOP2A", "UBE2C", "BIRC5")) {
  make_cor_plot(
    tcga_df, "HMMR_expr", paste0(g, "_expr"),
    paste0("TCGA-CESC: HMMR vs ", g),
    "HMMR expression", paste0(g, " expression"),
    file.path(out_dir, paste0("TCGA_HMMR_vs_", g, ".pdf"))
  )
}

# HMMR KM
tcga_df$HMMR_group <- ifelse(
  tcga_df$HMMR_expr >= median(tcga_df$HMMR_expr, na.rm = TRUE),
  "High", "Low"
)

fit_km <- survfit(Surv(OS.time, OS) ~ HMMR_group, data = tcga_df)

p_km <- ggsurvplot(
  fit_km,
  data = tcga_df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("#4C78A8", "#D62728"),
  title = "TCGA-CESC: OS by HMMR expression",
  xlab = "OS time",
  ylab = "Overall survival probability",
  ggtheme = theme_classic(base_size = 14)
)

pdf(file.path(out_dir, "TCGA_HMMR_KM.pdf"), width = 7.2, height = 6.8)
print(p_km)
dev.off()

# HMMR univariate Cox
cox_uni <- coxph(Surv(OS.time, OS) ~ HMMR_expr, data = tcga_df)
cox_sum <- summary(cox_uni)

cox_out <- data.frame(
  gene = "HMMR",
  HR = cox_sum$coefficients[1, "exp(coef)"],
  lower95 = cox_sum$conf.int[1, "lower .95"],
  upper95 = cox_sum$conf.int[1, "upper .95"],
  p_value = cox_sum$coefficients[1, "Pr(>|z|)"]
)

write.csv(cox_out, file.path(out_dir, "TCGA_HMMR_univariate_cox.csv"), row.names = FALSE)

# -----------------------------
# 6. Task 3: 3 hub genes expression and correlations
# -----------------------------
hub_cor_summary <- data.frame()

hub_genes <- c("HMMR", "NUF2", "CDKN3")

# pairwise hub-gene correlations
pairs_list <- list(
  c("HMMR", "NUF2"),
  c("HMMR", "CDKN3"),
  c("NUF2", "CDKN3")
)

for (pp in pairs_list) {
  g1 <- pp[1]
  g2 <- pp[2]

  cor_res <- suppressWarnings(cor.test(tcga_df[[paste0(g1, "_expr")]],
                                       tcga_df[[paste0(g2, "_expr")]],
                                       method = "spearman"))

  hub_cor_summary <- rbind(
    hub_cor_summary,
    data.frame(
      comparison = paste0(g1, "_vs_", g2),
      rho = as.numeric(cor_res$estimate),
      p_value = cor_res$p.value
    )
  )

  make_cor_plot(
    tcga_df,
    paste0(g1, "_expr"),
    paste0(g2, "_expr"),
    paste0("TCGA-CESC: ", g1, " vs ", g2),
    paste0(g1, " expression"),
    paste0(g2, " expression"),
    file.path(out_dir, paste0("TCGA_", g1, "_vs_", g2, ".pdf"))
  )
}

# each hub gene vs proliferation_core_score
for (g in hub_genes) {
  cor_res <- suppressWarnings(cor.test(tcga_df[[paste0(g, "_expr")]],
                                       tcga_df$proliferation_core_score,
                                       method = "spearman"))

  hub_cor_summary <- rbind(
    hub_cor_summary,
    data.frame(
      comparison = paste0(g, "_vs_proliferation_core_score"),
      rho = as.numeric(cor_res$estimate),
      p_value = cor_res$p.value
    )
  )

  make_cor_plot(
    tcga_df,
    paste0(g, "_expr"),
    "proliferation_core_score",
    paste0("TCGA-CESC: ", g, " vs proliferation_core_score"),
    paste0(g, " expression"),
    "proliferation_core_score",
    file.path(out_dir, paste0("TCGA_", g, "_vs_proliferation_core_score.pdf"))
  )
}

write.csv(hub_cor_summary, file.path(out_dir, "TCGA_hub_gene_correlation_summary.csv"), row.names = FALSE)

# optional: hub genes by stage
for (g in hub_genes) {
  make_group_plot(
    tcga_df, "stage", paste0(g, "_expr"),
    paste0("TCGA-CESC: ", g, " expression by stage"),
    paste0(g, " expression"),
    file.path(out_dir, paste0("TCGA_", g, "_by_stage.pdf"))
  )
}

cat("TCGA-CESC validation finished.\n")