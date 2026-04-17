# ============================================================
# Prepare L1000CDS2 Signature B
# HMMR virtual knockout signature from prioritized_6_9_10_pooled
# input: prioritized_6_9_10_pooled scTenifoldKnk result
# output: gene lists for L1000CDS2 mimic query
# ============================================================

library(dplyr)

rm(list = ls())
set.seed(12345)

# -----------------------------
# 0. Paths
# -----------------------------
# If your actual file name differs, only edit this one line.
ko_rdata <- "/Users/donglinlu/Desktop/CC1/12.drugEnrich/L1000CDS2_main_screen/HMMR_prioritized_6_9_10_pooled_scTenifoldKnk_result.Rdata"
outdir <- "/Users/donglinlu/Desktop/CC1/12.drugEnrich/L1000CDS2_main_screen"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Parameters
# -----------------------------
target_gene <- "HMMR"
n_top <- 100
adj_cutoff <- 0.05

# -----------------------------
# 2. Load result object
# -----------------------------
obj_names <- load(ko_rdata)

result_obj_name <- obj_names[sapply(obj_names, function(x) {
  obj <- get(x)
  is.list(obj) && "diffRegulation" %in% names(obj)
})][1]

if (is.na(result_obj_name)) {
  stop("No scTenifoldKnk result object with diffRegulation found in KO Rdata")
}

res <- get(result_obj_name)
cat("Loaded KO result object:", result_obj_name, "\n")

# -----------------------------
# 3. Extract and standardize table
# -----------------------------
df <- res$diffRegulation

if (!("gene" %in% colnames(df))) {
  stop("diffRegulation table does not contain a gene column")
}
if (!("Z" %in% colnames(df))) {
  stop("diffRegulation table does not contain a Z column")
}

if (!("p.adj" %in% colnames(df))) {
  df$p.adj <- NA_real_
}
if (!("p" %in% colnames(df))) {
  df$p <- NA_real_
}

df <- df %>%
  filter(gene != target_gene) %>%
  mutate(
    absZ = abs(Z),
    significance_class = case_when(
      !is.na(p.adj) & p.adj < adj_cutoff & Z > 0 ~ "positive_Z_adjSig",
      !is.na(p.adj) & p.adj < adj_cutoff & Z < 0 ~ "negative_Z_adjSig",
      TRUE ~ "not_adjSig"
    )
  ) %>%
  arrange(desc(absZ))

write.csv(
  df,
  file = file.path(outdir, "SignatureB_HMMR_virtual_KO_all_diffRegulation.csv"),
  row.names = FALSE
)

# -----------------------------
# 4. Export L1000 gene lists
# -----------------------------
positive_z_genes <- df %>%
  filter(!is.na(p.adj), p.adj < adj_cutoff, Z > 0) %>%
  arrange(desc(Z)) %>%
  slice_head(n = n_top) %>%
  pull(gene)

negative_z_genes <- df %>%
  filter(!is.na(p.adj), p.adj < adj_cutoff, Z < 0) %>%
  arrange(Z) %>%
  slice_head(n = n_top) %>%
  pull(gene)

write.table(
  positive_z_genes,
  file = file.path(outdir, "SignatureB_positiveZ_top100.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  negative_z_genes,
  file = file.path(outdir, "SignatureB_negativeZ_top100.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

max_len <- max(length(positive_z_genes), length(negative_z_genes))
submit_tab <- data.frame(
  up_genes = c(positive_z_genes, rep(NA, max_len - length(positive_z_genes))),
  down_genes = c(negative_z_genes, rep(NA, max_len - length(negative_z_genes)))
)

write.csv(
  submit_tab,
  file = file.path(outdir, "SignatureB_L1000_submission_table.csv"),
  row.names = FALSE
)

# -----------------------------
# 5. Summary
# -----------------------------
summary_tab <- data.frame(
  signature = "Signature_B",
  ko_rdata = basename(ko_rdata),
  result_object_used = result_obj_name,
  target_gene = target_gene,
  adj_cutoff = adj_cutoff,
  total_genes_in_diffRegulation = nrow(df),
  adj_sig_positive_Z = sum(!is.na(df$p.adj) & df$p.adj < adj_cutoff & df$Z > 0),
  adj_sig_negative_Z = sum(!is.na(df$p.adj) & df$p.adj < adj_cutoff & df$Z < 0),
  exported_positive_Z = length(positive_z_genes),
  exported_negative_Z = length(negative_z_genes)
)

write.table(
  summary_tab,
  file = file.path(outdir, "SignatureB_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nSignature B preparation completed successfully.\n")
cat("Output directory:", outdir, "\n")
cat("Positive Z genes exported:", length(positive_z_genes), "\n")
cat("Negative Z genes exported:", length(negative_z_genes), "\n")
