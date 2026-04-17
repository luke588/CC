rm(list = ls())
options(stringsAsFactors = FALSE)

source("/Users/donglinlu/Desktop/CC1/00.manuscript_plot_theme.R")

#######################00. libraries#######################
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

#######################01. helper functions#######################
safe_numeric <- function(x) {
  x <- gsub(",", "", as.character(x))
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  x
}

sanitize_names <- function(x) {
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[[:space:][:punct:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

find_first_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

normalize_stage <- function(x) {
  x0 <- as.character(x)
  x1 <- toupper(trimws(x0))
  x1 <- gsub("[[:space:]]+", "", x1)
  x1 <- gsub("_EPITHELIALCELLS", "", x1)
  x1 <- gsub("EPITHELIALCELLS", "", x1)
  x1 <- gsub("[^A-Z0-9]+", "_", x1)

  out <- ifelse(
    grepl("NO_HPV|NOHPV", x1),
    "NO_HPV",
    ifelse(
      grepl("^N_HPV$|^NHPV$|N_HPV|NHPV", x1),
      "N_HPV",
      ifelse(
        grepl("HSIL_HPV|HSILHPV", x1),
        "HSIL_HPV",
        ifelse(
          grepl("CA_HPV|CAHPV", x1),
          "CA_HPV",
          NA_character_
        )
      )
    )
  )

  out
}

#######################02. global settings#######################
workDir <- "/Users/donglinlu/Desktop/CC1/4.epithelial subset"
setwd(workDir)

stage_levels <- c("NO_HPV", "N_HPV", "HSIL_HPV", "CA_HPV")
stage_labels_heatmap <- c(
  "NO_HPV" = "NO\nHPV",
  "N_HPV" = "N\nHPV",
  "HSIL_HPV" = "HSIL\nHPV",
  "CA_HPV" = "CA\nHPV"
)

state_order <- c(
  "EMT_like_epithelial",
  "secretory_like_epithelial",
  "keratinizing_squamous_epithelial",
  "aberrant_secretory_epithelial",
  "inflammatory_secretory_epithelial",
  "proliferative_epithelial",
  "basal_stromal_like_epithelial",
  "superficial_keratinizing_epithelial",
  "terminally_differentiated_squamous_epithelial",
  "basal_EMT_like_epithelial",
  "rare_epithelial_like_cells"
)

state_colors <- c(
  "EMT_like_epithelial" = "#E76F61",
  "secretory_like_epithelial" = "#D98C00",
  "keratinizing_squamous_epithelial" = "#B3A100",
  "aberrant_secretory_epithelial" = "#67B300",
  "inflammatory_secretory_epithelial" = "#12B85B",
  "proliferative_epithelial" = "#19B7A5",
  "basal_stromal_like_epithelial" = "#1CA9C9",
  "superficial_keratinizing_epithelial" = "#2296E3",
  "terminally_differentiated_squamous_epithelial" = "#9D78E3",
  "basal_EMT_like_epithelial" = "#D85BE1",
  "rare_epithelial_like_cells" = "#EE5AA8"
)

# canonical priority order fixed as 9 -> 10 -> 6
priority_levels <- c(
  "malignant_squamous_epithelial",
  "EMT_like_epithelial",
  "proliferation_core_epithelial"
)

priority_colors <- c(
  "malignant_squamous_epithelial" = "#D73027",
  "EMT_like_epithelial" = "#4575B4",
  "proliferation_core_epithelial" = "#1A9850"
)

#######################03. load object#######################
load("02.Epithelial_subset.manualAnnotation.Rdata")

#######################04. read clean-state proportion table#######################
clean_file <- "03.epiSubset.cleanCellType_stage_proportion.csv"
clean.raw <- read.csv(clean_file, stringsAsFactors = FALSE)
colnames(clean.raw) <- sanitize_names(colnames(clean.raw))

cell_col <- find_first_col(
  clean.raw,
  c("celltype", "cell_type", "state", "cell_state", "cellstate")
)

stage_col <- find_first_col(
  clean.raw,
  c("stage", "group", "progression_stage")
)

prop_col <- find_first_col(
  clean.raw,
  c(
    "proportion", "prop", "frequency", "freq", "ratio",
    "fraction", "percentage", "percent", "value"
  )
)

if (is.na(cell_col) | is.na(stage_col) | is.na(prop_col)) {
  if (ncol(clean.raw) >= 3) {
    colnames(clean.raw)[1:3] <- c("celltype", "stage", "proportion")
    cell_col <- "celltype"
    stage_col <- "stage"
    prop_col <- "proportion"
  } else {
    stop("Cannot identify CellType / Stage / Proportion columns in 03.epiSubset.cleanCellType_stage_proportion.csv")
  }
}

prop.df <- data.frame(
  CellType = as.character(clean.raw[[cell_col]]),
  Stage = normalize_stage(clean.raw[[stage_col]]),
  Proportion = safe_numeric(clean.raw[[prop_col]]),
  stringsAsFactors = FALSE
)

prop.df.sum <- prop.df %>%
  filter(
    !is.na(CellType),
    CellType != "",
    CellType != "NA",
    !is.na(Stage)
  ) %>%
  mutate(
    CellType = as.character(CellType),
    Stage = as.character(Stage)
  ) %>%
  group_by(CellType, Stage) %>%
  summarise(Proportion = sum(Proportion, na.rm = TRUE), .groups = "drop") %>%
  filter(
    !is.na(CellType),
    CellType != "",
    CellType != "NA"
  )

state_breaks <- state_order[state_order %in% unique(prop.df.sum$CellType)]
state_breaks <- state_breaks[!is.na(state_breaks)]

prop.df.sum <- prop.df.sum %>%
  filter(CellType %in% state_breaks) %>%
  complete(
    CellType = state_breaks,
    Stage = stage_levels,
    fill = list(Proportion = 0)
  )

prop.df.sum$Stage <- factor(prop.df.sum$Stage, levels = stage_levels)
prop.df.sum$CellType <- factor(prop.df.sum$CellType, levels = state_breaks)

write.table(
  prop.df.sum,
  file = "04.cleanState.proportionTable.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#######################05. Figure: stage-wise composition barplot#######################
p_bar <- ggplot(prop.df.sum, aes(x = Stage, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", width = 0.72, colour = NA) +
  scale_fill_manual(values = state_colors[state_breaks], breaks = state_breaks) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Stage-wise composition of clean epithelial states",
    x = NULL,
    y = "Relative proportion",
    fill = "Cell state"
  ) +
  theme_pub()

save_pub(
  plot_obj = p_bar,
  pdf_file = "04.cleanState.proportion.barplot.pdf",
  png_file = "04.cleanState.proportion.barplot.png",
  width = 11.5,
  height = 7.8,
  dpi = 600
)

#######################06. Figure: all clean states lineplot#######################
p_line_all <- ggplot(
  prop.df.sum,
  aes(x = Stage, y = Proportion, group = CellType, color = CellType)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.8) +
  scale_color_manual(values = state_colors[state_breaks], breaks = state_breaks) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Progression-associated enrichment of clean epithelial states",
    x = NULL,
    y = "Proportion",
    color = "Cell state"
  ) +
  theme_pub()

save_pub(
  plot_obj = p_line_all,
  pdf_file = "04.cleanState.allStates.lineplot.pdf",
  png_file = "04.cleanState.allStates.lineplot.png",
  width = 13.0,
  height = 8.0,
  dpi = 600
)

#######################07. Figure: heatmap (ggplot version)#######################
heat.df <- prop.df.sum %>%
  filter(
    !is.na(CellType),
    as.character(CellType) != "NA"
  ) %>%
  mutate(
    Stage = factor(as.character(Stage), levels = stage_levels),
    CellType = factor(as.character(CellType), levels = rev(state_breaks))
  ) %>%
  droplevels()

print(unique(as.character(heat.df$CellType)))

p_heat <- ggplot(heat.df, aes(x = Stage, y = CellType, fill = Proportion)) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_gradientn(
    colours = c("#3B73B9", "#DCE9F6", "#F4E7A1", "#E64B35"),
    name = "Proportion"
  ) +
  scale_x_discrete(labels = stage_labels_heatmap) +
  labs(
    title = "Stage-wise proportion heatmap of clean epithelial states",
    x = NULL,
    y = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.9,
      size = PUB_AXIS_TEXT_SIZE
    )
  )

save_pub(
  plot_obj = p_heat,
  pdf_file = "04.cleanState.proportion.heatmap.pdf",
  png_file = "04.cleanState.proportion.heatmap.png",
  width = 11.8,
  height = 8.2,
  dpi = 600
)

#######################08. read priority file for Figure 1C#######################
priority.file <- "step1_5.07_focus_clusters_6_9_10_stage_proportion.csv"

if (file.exists(priority.file)) {
  focus.raw <- read.csv(priority.file, stringsAsFactors = FALSE)
  colnames(focus.raw) <- sanitize_names(colnames(focus.raw))

  stage_col2 <- find_first_col(
    focus.raw,
    c("group", "stage", "progression_stage")
  )

  cluster_col2 <- find_first_col(
    focus.raw,
    c("seurat_clusters", "cluster", "focus_cluster")
  )

  prop_col2 <- find_first_col(
    focus.raw,
    c(
      "cluster_prop", "proportion", "prop", "frequency", "freq",
      "ratio", "fraction", "percentage", "percent", "value"
    )
  )

  if (is.na(stage_col2)) stop("Cannot find stage column in step1_5.07_focus_clusters_6_9_10_stage_proportion.csv")
  if (is.na(cluster_col2)) stop("Cannot find cluster column in step1_5.07_focus_clusters_6_9_10_stage_proportion.csv")
  if (is.na(prop_col2)) stop("Cannot find proportion column in step1_5.07_focus_clusters_6_9_10_stage_proportion.csv")

  priority.df <- data.frame(
    Stage = normalize_stage(focus.raw[[stage_col2]]),
    cluster_id = as.character(focus.raw[[cluster_col2]]),
    Proportion = safe_numeric(focus.raw[[prop_col2]]),
    stringsAsFactors = FALSE
  )

  priority.df$PriorityState <- case_when(
    priority.df$cluster_id == "9"  ~ "malignant_squamous_epithelial",
    priority.df$cluster_id == "10" ~ "EMT_like_epithelial",
    priority.df$cluster_id == "6"  ~ "proliferation_core_epithelial",
    TRUE ~ NA_character_
  )

  priority.df <- priority.df %>%
    filter(!is.na(Stage), !is.na(PriorityState)) %>%
    group_by(PriorityState, Stage) %>%
    summarise(Proportion = sum(Proportion, na.rm = TRUE), .groups = "drop")

} else {
  warning(
    "step1_5.07_focus_clusters_6_9_10_stage_proportion.csv not found. Using clean-state fallback mapping."
  )

  priority.df <- prop.df.sum %>%
    mutate(
      PriorityState = case_when(
        as.character(CellType) == "superficial_keratinizing_epithelial" ~ "malignant_squamous_epithelial",
        as.character(CellType) == "EMT_like_epithelial" ~ "EMT_like_epithelial",
        as.character(CellType) == "proliferative_epithelial" ~ "proliferation_core_epithelial",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(PriorityState)) %>%
    group_by(PriorityState, Stage) %>%
    summarise(Proportion = sum(Proportion, na.rm = TRUE), .groups = "drop")
}

priority.df <- priority.df %>%
  complete(
    PriorityState = priority_levels,
    Stage = stage_levels,
    fill = list(Proportion = 0)
  )

priority.df$Stage <- factor(priority.df$Stage, levels = stage_levels)
priority.df$PriorityState <- factor(priority.df$PriorityState, levels = priority_levels)

priority.df <- priority.df %>%
  arrange(PriorityState, Stage)

write.table(
  priority.df,
  file = "04.cleanState.priority3.proportionTable.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#######################09. Figure 1C: prioritized programs lineplot#######################
p_priority <- ggplot(
  priority.df,
  aes(x = Stage, y = Proportion, group = PriorityState, color = PriorityState)
) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3.4) +
  scale_color_manual(values = priority_colors, breaks = priority_levels) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Progression-associated enrichment of prioritized epithelial programs",
    x = NULL,
    y = "Proportion",
    color = "Prioritized program"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = PUB_AXIS_TEXT_SIZE
    )
  )

save_pub(
  plot_obj = p_priority,
  pdf_file = "04.cleanState.progression.lineplot.pdf",
  png_file = "04.cleanState.progression.lineplot.png",
  width = 12.8,
  height = 7.4,
  dpi = 600
)

#######################10. summary table of prioritized programs#######################
priority.summary <- priority.df %>%
  mutate(
    Stage = as.character(Stage),
    PriorityState = as.character(PriorityState)
  ) %>%
  pivot_wider(
    names_from = Stage,
    values_from = Proportion,
    values_fill = 0
  )

priority.summary$NO_HPV <- safe_numeric(priority.summary$NO_HPV)
priority.summary$N_HPV <- safe_numeric(priority.summary$N_HPV)
priority.summary$HSIL_HPV <- safe_numeric(priority.summary$HSIL_HPV)
priority.summary$CA_HPV <- safe_numeric(priority.summary$CA_HPV)

priority.summary <- priority.summary %>%
  mutate(
    delta_CA_vs_NO = CA_HPV - NO_HPV,
    delta_CA_vs_HSIL = CA_HPV - HSIL_HPV,
    mean_preCA = rowMeans(cbind(NO_HPV, N_HPV, HSIL_HPV), na.rm = TRUE),
    CA_enrichment = CA_HPV / (mean_preCA + 1e-6)
  )

write.table(
  priority.summary,
  file = "04.cleanState.priority3.summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#######################11. save objects#######################
save_objs <- c("prop.df.sum", "priority.df", "priority.summary")
if (exists("epi")) save_objs <- c(save_objs, "epi")
if (exists("epi.clean")) save_objs <- c(save_objs, "epi.clean")

save(
  list = save_objs,
  file = "04.cleanState.priority.Rdata"
)

#######################12. done#######################
cat("\nDone.\n")
cat("Main outputs saved:\n")
cat("04.cleanState.proportionTable.txt\n")
cat("04.cleanState.proportion.barplot.pdf\n")
cat("04.cleanState.proportion.barplot.png\n")
cat("04.cleanState.proportion.heatmap.pdf\n")
cat("04.cleanState.proportion.heatmap.png\n")
cat("04.cleanState.allStates.lineplot.pdf\n")
cat("04.cleanState.allStates.lineplot.png\n")
cat("04.cleanState.priority3.proportionTable.txt\n")
cat("04.cleanState.progression.lineplot.pdf\n")
cat("04.cleanState.progression.lineplot.png\n")
cat("04.cleanState.priority3.summary.txt\n")
cat("04.cleanState.priority.Rdata\n")