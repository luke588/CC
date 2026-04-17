rm(list = ls())
options(stringsAsFactors = FALSE)

# =========================================================
# Figure 6A
# Cleaned Venn diagram of candidate drugs:
# left circle = HMMR-high
# right circle = HMMR-KO
# =========================================================

# -----------------------------#
# paths
# -----------------------------#
root_dir  <- "/Users/donglinlu/Desktop/CC1/L1000CDS2_main_screen"
input_dir <- file.path(root_dir, "input")
out_dir   <- file.path(input_dir, "Venn_compare_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------#
# display settings
# -----------------------------#
left_label  <- "HMMR-high"
right_label <- "HMMR-KO"
main_title  <- "Convergent candidate drugs"

# -----------------------------#
# packages
# -----------------------------#
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
library(VennDiagram)
library(grid)

# -----------------------------#
# helper functions
# -----------------------------#
find_first_match <- function(search_root, patterns) {
  all_files <- list.files(search_root, recursive = TRUE, full.names = TRUE)
  for (pt in patterns) {
    hit <- all_files[grepl(pt, basename(all_files), ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  return(NA_character_)
}

read_drug_set_raw <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)

  df <- tryCatch(
    read.csv(file, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  )

  preferred_cols <- grep(
    "drug|compound|pert|name|candidate",
    colnames(df),
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(preferred_cols) > 0) {
    vals <- unique(unlist(df[, preferred_cols, drop = FALSE], use.names = FALSE))
  } else {
    non_num <- names(df)[!sapply(df, is.numeric)]
    use_col <- if (length(non_num) > 0) non_num[1] else names(df)[1]
    vals <- df[[use_col]]
  }

  vals <- as.character(vals)
  vals <- trimws(vals)
  vals <- vals[!is.na(vals) & vals != ""]
  vals
}

is_valid_drug_entry <- function(x) {
  x2 <- trimws(as.character(x))
  x2u <- toupper(x2)

  if (x2 == "" || is.na(x2)) return(FALSE)

  # remove obvious non-drug entries
  if (x2u %in% c("NONE", "NA", "N/A", "NULL")) return(FALSE)
  if (grepl("^HTTPS?://", x2u)) return(FALSE)
  if (grepl("^WWW\\.", x2u)) return(FALSE)
  if (grepl("PUBCHEM|LIFE\\.CCS|SUMMARY\\.CGI|SOURCE=|INPUT=|CID=", x2u)) return(FALSE)

  # very long strings are usually URLs / annotations, not clean drug names
  if (nchar(x2u) > 60) return(FALSE)

  return(TRUE)
}

clean_drug_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  x <- gsub("\\s+", " ", x)
  x <- gsub("^['\"]|['\"]$", "", x)
  x <- toupper(x)
  unique(x)
}

# -----------------------------#
# locate Signature A and B files
# -----------------------------#
sigA_file <- find_first_match(
  root_dir,
  patterns = c(
    "^SignatureA_single_drug_deduplicated\\.csv$",
    "^SignatureA_drugs_cleaned\\.csv$",
    "^SignatureA.*deduplicated.*\\.csv$",
    "^SignatureA.*cleaned.*\\.csv$"
  )
)

sigB_file <- find_first_match(
  root_dir,
  patterns = c(
    "^SignatureB_single_drug_deduplicated\\.csv$",
    "^SignatureB_drugs_cleaned\\.csv$",
    "^SignatureB.*deduplicated.*\\.csv$",
    "^SignatureB.*cleaned.*\\.csv$"
  )
)

if (is.na(sigA_file)) stop("Cannot find Signature A drug file under: ", root_dir)
if (is.na(sigB_file)) stop("Cannot find Signature B drug file under: ", root_dir)

# -----------------------------#
# read raw sets
# -----------------------------#
set_A_raw <- read_drug_set_raw(sigA_file)
set_B_raw <- read_drug_set_raw(sigB_file)

# -----------------------------#
# clean sets
# -----------------------------#
set_A_removed <- set_A_raw[!vapply(set_A_raw, is_valid_drug_entry, logical(1))]
set_B_removed <- set_B_raw[!vapply(set_B_raw, is_valid_drug_entry, logical(1))]

set_A <- clean_drug_name(set_A_raw[vapply(set_A_raw, is_valid_drug_entry, logical(1))])
set_B <- clean_drug_name(set_B_raw[vapply(set_B_raw, is_valid_drug_entry, logical(1))])

if (length(set_A) == 0) stop("No valid drug names left in Signature A after cleaning.")
if (length(set_B) == 0) stop("No valid drug names left in Signature B after cleaning.")

overlap_drugs <- intersect(set_A, set_B)

# -----------------------------#
# save cleaned outputs
# -----------------------------#
write.csv(
  data.frame(drug = sort(set_A), stringsAsFactors = FALSE),
  file.path(out_dir, "01.SignatureA_HMMR_high_drugs_cleaned.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(drug = sort(set_B), stringsAsFactors = FALSE),
  file.path(out_dir, "02.SignatureB_HMMR_KO_drugs_cleaned.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(drug = sort(overlap_drugs), stringsAsFactors = FALSE),
  file.path(out_dir, "03.HMMR_high_vs_KO_overlap_drugs_cleaned.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(removed_entry = unique(set_A_removed), stringsAsFactors = FALSE),
  file.path(out_dir, "04.SignatureA_removed_non_drug_entries.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(removed_entry = unique(set_B_removed), stringsAsFactors = FALSE),
  file.path(out_dir, "05.SignatureB_removed_non_drug_entries.csv"),
  row.names = FALSE,
  quote = FALSE
)

summary_lines <- c(
  paste0("Signature A file: ", sigA_file),
  paste0("Signature B file: ", sigB_file),
  paste0("Left circle label: ", left_label),
  paste0("Right circle label: ", right_label),
  paste0("Raw entries in Signature A: ", length(set_A_raw)),
  paste0("Raw entries in Signature B: ", length(set_B_raw)),
  paste0("Valid drugs in Signature A after cleaning: ", length(set_A)),
  paste0("Valid drugs in Signature B after cleaning: ", length(set_B)),
  paste0("Overlap drugs after cleaning: ", length(overlap_drugs))
)

writeLines(
  summary_lines,
  con = file.path(out_dir, "06.venn_summary.txt")
)

# -----------------------------#
# draw venn
# left = HMMR-high, right = HMMR-KO
# -----------------------------#
venn.plot <- venn.diagram(
  x = list(
    "LEFT"  = set_A,
    "RIGHT" = set_B
  ),
  filename = NULL,
  scaled = FALSE,

  # user-preferred colors
  fill = c("#79B6D8", "#F0B65A"),
  alpha = 0.90,
  col = "black",
  lwd = 2,
  lty = "solid",

  # numbers
  cex = 2.2,
  fontface = "bold",
  fontfamily = "sans",

  # do NOT use default category names, we add them manually
  category.names = c("", ""),
  cat.cex = 0.01,
  cat.col = c("white", "white"),

  margin = 0.08
)

# -----------------------------#
# save pdf
# -----------------------------#
pdf(
  file = file.path(out_dir, "07.HMMR_high_vs_KO_drug_overlap_venn.pdf"),
  width = 7.8,
  height = 6.2
)
grid.newpage()
grid.draw(venn.plot)

# main title
grid.text(
  label = main_title,
  x = 0.5, y = 0.965,
  gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "sans")
)

# left/right labels written directly over each circle
grid.text(
  label = left_label,
  x = 0.30, y = 0.80,
  gp = gpar(fontsize = 15, fontface = "bold", fontfamily = "sans")
)
grid.text(
  label = right_label,
  x = 0.70, y = 0.80,
  gp = gpar(fontsize = 15, fontface = "bold", fontfamily = "sans")
)

dev.off()

# -----------------------------#
# save png
# -----------------------------#
png(
  filename = file.path(out_dir, "07.HMMR_high_vs_KO_drug_overlap_venn.png"),
  width = 2400,
  height = 1900,
  res = 300
)
grid.newpage()
grid.draw(venn.plot)

grid.text(
  label = main_title,
  x = 0.5, y = 0.965,
  gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "sans")
)
grid.text(
  label = left_label,
  x = 0.30, y = 0.80,
  gp = gpar(fontsize = 15, fontface = "bold", fontfamily = "sans")
)
grid.text(
  label = right_label,
  x = 0.70, y = 0.80,
  gp = gpar(fontsize = 15, fontface = "bold", fontfamily = "sans")
)

dev.off()

# -----------------------------#
# also show in current plotting window
# -----------------------------#
grid.newpage()
grid.draw(venn.plot)
grid.text(
  label = main_title,
  x = 0.5, y = 0.965,
  gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "sans")
)
grid.text(
  label = left_label,
  x = 0.30, y = 0.80,
  gp = gpar(fontsize = 15, fontface = "bold", fontfamily = "sans")
)
grid.text(
  label = right_label,
  x = 0.70, y = 0.80,
  gp = gpar(fontsize = 15, fontface = "bold", fontfamily = "sans")
)

cat("Finished.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Valid overlap drugs after cleaning: ", length(overlap_drugs), "\n", sep = "")
print(sort(overlap_drugs))