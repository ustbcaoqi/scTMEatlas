# ==============================================================================
# Script Name: Figure1.R
# Description: TME-related analyses for Figure 1 and ED Figures
# Author: Zhang Lab
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(tidytext)
  library(scales)
  library(ggradar)
  library(patchwork)
  library(geomtextpath)
  library(forestplot)
  library(grid)
})

args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript Figure2.R <output_prefix>")
}
pat <- args[1]

check_file <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
}

check_columns <- function(df, required_cols, object_name = "data frame") {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

# ==============================================================================
# Fig. 1c. Identification of cellular modules
# ==============================================================================

# Input:
#   seu.rds: integrated Seurat object with metadata columns:
#     CancerType, SampleID, FineCellType, MajorLineage
# Output:
#   <pat>_TME_lineage_normalized_abundance.rds
#   <pat>_Cellular_Module_Assignments.csv
#   <pat>_Figure2_TME_Subtyping_Heatmap.pdf

scObject <- readRDS("allcell_seu.rds")

required_meta_cols <- c("CancerType", "SampleID", "FineCellType", "MajorLineage")
check_columns(scObject@meta.data, required_meta_cols, "Seurat metadata")

# Step 1. Sample-level quality control
# Exclude blood-related tumor samples due to their distinct systemic microenvironment.
sc_filtered <- subset(scObject, subset = CancerType != "Blood_related_tumors")

# Retain samples containing at least 4,000 cells.
cell_counts_per_sample <- table(sc_filtered$SampleID)
valid_samples <- names(cell_counts_per_sample)[cell_counts_per_sample >= 4000]
sc_tme <- subset(sc_filtered, subset = SampleID %in% valid_samples)

if (ncol(sc_tme) == 0) {
  stop("No cells remain after sample-level filtering.")
}

# Step 2. Lineage-normalized relative abundance calculation
cell_matrix <- as.matrix(table(sc_tme$FineCellType, sc_tme$SampleID))

cell_metadata <- sc_tme@meta.data %>%
  select(FineCellType, MajorLineage) %>%
  distinct() %>%
  remove_rownames()

abundance_list <- list()

for (lineage in unique(cell_metadata$MajorLineage)) {
  target_subtypes <- cell_metadata$FineCellType[cell_metadata$MajorLineage == lineage]
  target_subtypes <- intersect(target_subtypes, rownames(cell_matrix))
  
  if (length(target_subtypes) == 0) next
  
  sub_matrix <- cell_matrix[target_subtypes, , drop = FALSE]
  lineage_totals <- colSums(sub_matrix)
  
  sub_norm <- sweep(
    sub_matrix,
    2,
    ifelse(lineage_totals == 0, NA, lineage_totals),
    FUN = "/"
  )
  
  sub_norm[is.na(sub_norm)] <- 0
  abundance_list[[lineage]] <- sub_norm
}

if (length(abundance_list) == 0) {
  stop("No lineage-specific abundance matrices were generated.")
}

tme_abundance <- do.call(rbind, abundance_list)
saveRDS(tme_abundance, paste0(pat, "_TME_lineage_normalized_abundance.rds"))

# Step 3. Pairwise Pearson correlation
cor_matrix <- cor(t(tme_abundance), method = "pearson")
cor_matrix[is.na(cor_matrix)] <- 0
dist_matrix <- as.dist(1 - cor_matrix)

# Step 4. Hierarchical clustering
# Note: correlation distance is used here to group cell subpopulations with similar
# sample-level abundance patterns.
hc_modules <- hclust(dist_matrix, method = "ward.D2")
module_assignments <- cutree(hc_modules, k = 5)

module_df <- data.frame(
  CellSubpopulation = names(module_assignments),
  CellularModule = paste0("CM", module_assignments),
  stringsAsFactors = FALSE
) %>%
  mutate(
    CM_Label = factor(CellularModule, levels = paste0("CM", 1:5)),
    DendrogramOrder = match(CellSubpopulation, hc_modules$labels[hc_modules$order])
  ) %>%
  arrange(CM_Label, DendrogramOrder)

write.csv(
  module_df %>% select(CellSubpopulation, CellularModule),
  paste0(pat, "_Cellular_Module_Assignments.csv"),
  row.names = FALSE
)

# Step 5. Visualization using ComplexHeatmap
cm_colors <- c(
  CM1 = "#33b39f",
  CM2 = "#f5af98",
  CM3 = "#eb6f5d",
  CM4 = "#6376a0",
  CM5 = "#71c9dd",
  Tumor = "#bebebe"
)

plot_mat <- cor_matrix[module_df$CellSubpopulation, module_df$CellSubpopulation]
split_factor <- module_df$CM_Label

ha_top <- HeatmapAnnotation(
  Module = split_factor,
  col = list(Module = cm_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE
)

ha_left <- rowAnnotation(
  Module = split_factor,
  col = list(Module = cm_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE
)

col_fun <- colorRamp2(c(-2,-1, 0, 1,2), c("#4d9ac6","#d4e6f0","white",  "#cb4e40", "#b5222e"))

final_heatmap <- Heatmap(
  plot_mat,
  name = "Correlation",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  row_split = split_factor,
  column_split = split_factor,
  row_gap = unit(1.2, "mm"),
  column_gap = unit(1.2, "mm"),
  top_annotation = ha_top,
  left_annotation = ha_left,
  show_row_names = FALSE,
  show_column_names = FALSE,
  border = TRUE
)

output_pdf <- paste0(pat, "_Figure2_TME_Subtyping_Heatmap.pdf")
pdf(output_pdf, width = 7, height = 6.8)
draw(final_heatmap, merge_legends = TRUE)
dev.off()

# ==============================================================================
# Fig. 1h. Radar plot
# ==============================================================================

# Inputs:
#   TME_marker_ring_data.csv
#     Required columns: x, y, group
#   TME_marker_genes.csv
#     Required column: gene
#   TME_radar_scores.csv
#     First column should be the group/TME subtype column required by ggradar;
#     remaining columns should be scaled numeric features from 0 to 1.
# Output:
#   <pat>_TME_subtypes_marker_radar.pdf

check_file("TME_marker_ring_data.csv")
check_file("TME_marker_genes.csv")
check_file("TME_radar_scores.csv")

tmp <- read.csv("TME_marker_ring_data.csv", check.names = FALSE)
gene_df <- read.csv("TME_marker_genes.csv", check.names = FALSE)
df_radar <- read.csv("TME_radar_scores.csv", check.names = FALSE)

check_columns(tmp, c("x", "y", "group"), "TME_marker_ring_data.csv")
check_columns(gene_df, c("gene"), "TME_marker_genes.csv")

genes <- gene_df$gene

tme_colors <- c(
  TLS = "#33b39f",
  Lymphoid = "#f5af98",
  Myeloid = "#eb6f5d",
  Steroid = "#6376a0",
  Stromal = "#71c9dd",
  Desert = "#bebebe"
)

missing_tme_colors <- setdiff(unique(tmp$group), names(tme_colors))
if (length(missing_tme_colors) > 0) {
  stop("Missing TME colors for: ", paste(missing_tme_colors, collapse = ", "))
}

p1 <- ggplot() +
  geom_bar(
    data = tmp,
    aes(x = x, y = y, fill = group),
    stat = "identity",
    position = "dodge"
  ) +
  geom_textpath(
    data = tmp,
    aes(x = 1, y = 1.1, label = group, group = group),
    position = position_dodge(width = 0.9),
    size = 5,
    color = "white",
    vjust = 0,
    upright = TRUE
  ) +
  geom_text(
    aes(
      x = rep(1, length(genes)),
      y = rep(3.5, length(genes)),
      label = genes,
      group = seq_along(genes)
    ),
    color = "black",
    size = 5,
    fontface = "italic",
    position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = tme_colors) +
  ylim(-7, 4.5) +
  coord_polar(start = -0.5) +
  theme_void() +
  theme(legend.position = "none")

p2 <- ggradar(
  df_radar,
  axis.labels = rep(NA, ncol(df_radar) - 1),
  grid.min = 0,
  grid.mid = 0.5,
  grid.max = 1,
  group.line.width = 1,
  group.point.size = 2,
  group.colours = tme_colors,
  background.circle.colour = "white",
  gridline.mid.colour = "#2b8c96",
  legend.position = "none",
  label.gridline.min = FALSE,
  label.gridline.mid = FALSE,
  label.gridline.max = FALSE
) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank()
  )

final_radar_plot <- p1 + inset_element(p2, left = 0, bottom = 0, right = 1, top = 1)

ggsave(
  filename = paste0(pat, "_TME_subtypes_marker_radar.pdf"),
  plot = final_radar_plot,
  width = 5,
  height = 5,
  dpi = 300
)





# ==============================================================================
# Figure 1i. Circular barplot
# ==============================================================================

# Input:
#   CancerType_6TME_sample_matrix.csv
#     Rows: cancer types
#     Columns: TME subtypes, expected to include:
#       TLS, Lymphoid, Myeloid, Steroid, Stromal, Desert
# Output:
#   <pat>_6TME_circular_plot.pdf

check_file("CancerType_6TME_sample_matrix.csv")
mat <- read.csv("CancerType_6TME_sample_matrix.csv", row.names = 1, check.names = FALSE)

mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
mat_norm <- sweep(mat, 1, rowSums(mat, na.rm = TRUE), FUN = "/")

df <- as.data.frame(mat_norm)
df$CancerType <- rownames(df)

df_long <- df %>%
  pivot_longer(
    cols = -CancerType,
    names_to = "TME_subtype",
    values_to = "value"
  ) %>%
  filter(value > 0) %>%
  group_by(TME_subtype) %>%
  mutate(
    rank = rank(-value, ties.method = "first"),
    label_text = ifelse(rank <= 3, as.character(CancerType), "")
  ) %>%
  ungroup() %>%
  mutate(
    CancerType_raw = CancerType,
    CancerType = reorder_within(CancerType, -value, TME_subtype)
  )

colors_custom <- c(
  AC = "#A65628", LEU = "#13ae68", BRCA = "#f8b62d", BTC = "#74a272",
  CC = "#afd0ee", CHC = "#9c5193", CRC = "#bcbd95", CTCL = "#6b8ec1",
  ESCC = "#efcfd3", FTC = "#a48b78", GBM = "#c89696", STAD = "#c7848f",
  HC = "#13393e", HCC = "#dea34d", MELA = "#585656", MM = "#ba5a67",
  NB = "#4994c4", NPC = "#96bd5b", NSCLC = "#7d4460", OS = "#e4cd54",
  OSCC = "#dfb5cf", OV = "#f3c2b2", PC = "#768876", PDAC = "#eeac85",
  RCC = "#5c9d9a", THCA = "#ea5515", UC = "#f39800", UCEC = "#b79ec6"
)

missing_colors <- setdiff(unique(df_long$CancerType_raw), names(colors_custom))
if (length(missing_colors) > 0) {
  stop("Missing colors for cancer types: ", paste(missing_colors, collapse = ", "))
}

df_plot <- df_long %>%
  group_by(TME_subtype) %>%
  arrange(TME_subtype, desc(value)) %>%
  mutate(
    n_current = n(),
    id = row_number(),
    angle = 90 - 360 * (id - 0.5) / n_current,
    hjust = ifelse(angle < -90, 1, 0),
    angle = ifelse(angle < -90, angle + 180, angle)
  ) %>%
  ungroup()

df_plot$TME_subtype <- factor(
  df_plot$TME_subtype,
  levels = c("TLS", "Lymphoid", "Myeloid", "Steroid", "Stromal", "Desert")
)

max_val <- max(df_plot$value, na.rm = TRUE)
center_y <- -0.8 * max_val

p <- ggplot(df_plot, aes(x = CancerType, y = value, fill = CancerType_raw)) +
  geom_bar(stat = "identity", width = 0.85, color = NA) +
  geom_text(
    aes(
      y = value + (max_val * 0.05),
      label = label_text,
      angle = angle,
      hjust = hjust
    ),
    size = 4,
    show.legend = FALSE
  ) +
  geom_text(
    data = df_plot %>% group_by(TME_subtype) %>% slice(1),
    aes(x = 1, y = center_y, label = TME_subtype),
    size = 5.5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  coord_polar(start = 0) +
  ylim(center_y, max_val * 1.5) +
  facet_wrap(~TME_subtype, scales = "free", nrow = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = colors_custom) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5, "pt"),
    panel.spacing = unit(0.2, "lines")
  )

ggsave(
  filename = paste0(pat, "_6TME_circular_plot.pdf"),
  plot = p,
  width = 20,
  height = 5,
  dpi = 300
)








# ==============================================================================
# ED Fig. 6c-d, Gini importance Bar plots
# ==============================================================================


library(ggplot2)
library(dplyr)
library(tidytext)

#Dataset Initialization and Logic Transformation
input_file  <- "data/processed/TME_importance_analysis_top20.csv"
output_pdf  <- "results/figures/TME_Subcluster_Contribution_Barplot.pdf"

if (!file.exists(input_file)) {
  stop(paste("Input data file missing:", input_file))
}

plot_data <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(Subtype = factor(Subtype, levels = c(
    "TLS", "Lymphoid", "Myeloid", "Steroid", "Stromal", "Desert"
  ))) %>%
  mutate(Direction = if_else(Contribution_Score > 0, "Up", "Down")) %>%
  mutate(Stars = case_when(
    Pvalue < 0.001 ~ "***",
    Pvalue < 0.01  ~ "**",
    Pvalue < 0.05  ~ "*",
    TRUE           ~ ""
  )) %>%
  mutate(CellType = reorder_within(CellType, Contribution_Score, Subtype))

#Multi-panel Visualization Configuration
p <- ggplot(plot_data, aes(x = Contribution_Score, y = CellType, fill = Direction)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.1) +
  
  # Fixed Alignment Bug: Embedded hjust expression within aes() to prevent facet re-indexing mismatch
  geom_text(
    aes(
      label = Stars, 
      hjust = if_else(Contribution_Score >= 0, -0.2, 1.2)
    ), 
    size = 3.5, 
    family = "sans"
  ) +
  
  facet_wrap(~Subtype, scales = "free_y", ncol = 6) +
  scale_y_reordered() + 
  scale_fill_manual(values = c("Up" = "#da6763", "Down" = "#4991c5")) +
  
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),     
    panel.grid.minor = element_blank(),     
    panel.border     = element_rect(linewidth = 0.8, color = "black", fill = NA),
    axis.text        = element_text(color = "black", size = 9),
    axis.title       = element_text(size = 11),
    strip.text       = element_text(size = 13, face = "bold"),
    strip.background = element_blank(),     
    legend.position  = "none"        
  ) +
  
  # Reserve symmetric margins for flanking significance labels
  scale_x_continuous(expand = expansion(mult = c(0.25, 0.25))) +
  labs(x = "Contribution Score", y = NULL)

#Save Figure with Publication-grade Parameters
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(output_pdf, p, width = 18, height = 3.5, device = "pdf")









# ==============================================================================
# ED Fig. 8 Forest plots
# ==============================================================================

# Input:
#   TME_survival_results.csv
#     Required columns: MP, cancer, HR, lower, upper, p
# Output:
#   TME_survival_forest_plots/Forest_Plot_<pat>_<MP>.pdf

check_file("TME_survival_results.csv")
all_res <- read.csv("TME_survival_results.csv", check.names = FALSE)

check_columns(
  all_res,
  c("MP", "cancer", "HR", "lower", "upper", "p"),
  "TME_survival_results.csv"
)

out_dir <- "TME_survival_forest_plots"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

mp_list_names <- unique(all_res$MP)

for (current_mp in mp_list_names) {
  df_mp <- all_res %>% filter(MP == current_mp)
  
  df_single <- df_mp %>%
    filter(cancer != "PanCancer") %>%
    arrange(HR)
  
  df_pan <- df_mp %>% filter(cancer == "PanCancer")
  
  if (nrow(df_pan) != 1) {
    stop("Expected exactly one PanCancer row for MP = ", current_mp)
  }
  
  plot_data <- bind_rows(df_single, df_pan)
  
  table_text <- matrix(
    c("Cancer type", "HR", "HR range", "p value"),
    nrow = 1
  )
  
  for (i in seq_len(nrow(plot_data))) {
    p_str <- formatC(plot_data$p[i], format = "e", digits = 2)
    
    row_text <- c(
      plot_data$cancer[i],
      sprintf("%.2f", plot_data$HR[i]),
      sprintf("(%.2f-%.2f)", plot_data$lower[i], plot_data$upper[i]),
      p_str
    )
    
    table_text <- rbind(table_text, row_text)
  }
  
  mean_val <- c(NA, plot_data$HR)
  lower_val <- c(NA, plot_data$lower)
  upper_val <- c(NA, plot_data$upper)
  
  is_summary <- c(TRUE, rep(FALSE, nrow(plot_data) - 1), TRUE)
  
  pdf_path <- file.path(out_dir, paste0("Forest_Plot_", pat, "_", current_mp, ".pdf"))
  pdf(pdf_path, width = 6, height = 6.5, onefile = FALSE)
  
  p <- forestplot(
    labeltext = table_text,
    mean = mean_val,
    lower = lower_val,
    upper = upper_val,
    is.summary = is_summary,
    lwd.ci = 1.6,
    lwd.xaxis = 1.5,
    lwd.zero = 1.5,
    hrzl_lines = list("2" = gpar(lwd = 1.5, col = "#222222")),
    clip = c(0.4, 2.0),
    xlog = FALSE,
    grid = structure(c(1), gp = gpar(lty = 2, col = "red")),
    boxsize = 0.3,
    xticks = c(0.5, 1.0, 1.5, 2.0),
    col = fpColors(
      box = "black",
      lines = "black",
      summary = "#3174af",
      zero = "black"
    ),
    txt_gp = fpTxtGp(
      label = gpar(fontsize = 12),
      ticks = gpar(fontsize = 12),
      xlab = gpar(fontsize = 18),
      title = gpar(fontsize = 12)
    ),
    xlab = "",
    title = paste("TME -", current_mp)
  )
  
  p <- fp_set_zebra_style(p, "#F5F5F5", "#FFFFFF")
  
  print(p)
  dev.off()
}



