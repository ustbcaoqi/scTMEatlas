# ==============================================================================
# Script Name: Figure1.R
# Description: Tumor-related analyses 
# Author: Zhang Lab
# ==============================================================================



##--- Load Dependencies
library(Seurat)
library(infercnv)
library(Matrix)
library(NMF)
library(future)
library(ggradar)
library(patchwork)
library(tidyverse)
library(vegan)
library(ggnewscale)
library(scales)
library(tidytext)

##--- Parallel Computing Configuration
plan("multicore", workers = 4)
options(future.globals.maxSize = 50 * 1024^3)
options(future.rng.onMisuse = "ignore")

##--- Global Input Arguments
args <- commandArgs(TRUE)
pat  <- args[1] # Patient ID
num  <- as.numeric(args[2]) # NMF Rank K

# Load Core Seurat Object
scObject <- readRDS("seu.rds")


# ==============================================================================
#  Identification of malignant epithelial cells (inferCNV)
# ==============================================================================
message(">>> Start inferCNV analysis for ", pat)

out_dir <- paste0(pat, "_inferCNV")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## Step 1. Restrict to epithelial population
sc_epi <- subset(scObject, subset = FinalCellType == "Epithelial")

## Step 2. Define reference vs observation
sc_epi$anno <- "Observation"
sc_epi$anno[sc_epi$TissueType == "Adjacent_NonTumor"] <- "Reference"

## Step 3. Write annotation file
annotation_file <- file.path(out_dir, paste0(pat, "_infercnv_anno.txt"))
write.table(
  data.frame(cell = colnames(sc_epi), group = sc_epi$anno),
  file = annotation_file,
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

## Step 4. Define reference group names
ref_group_names <- "Reference"

## Step 5. Create inferCNV object
if (Version(sc_epi)$major >= 5) {
  sc_epi <- JoinLayers(sc_epi, layers = "counts")
  counts_matrix <- LayerData(sc_epi, layer = "counts")
} else {
  counts_matrix <- GetAssayData(sc_epi, assay = "RNA", slot = "counts")
}

gene_order_file <- system.file("extdata", "hg38_gencode_v27.txt", package = "infercnv")
gene_pos <- read.table(gene_order_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
counts_matrix <- counts_matrix[rownames(counts_matrix) %in% gene_pos$V1, ]

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix       = counts_matrix,
  annotations_file        = annotation_file,
  delim                   = "\t",
  gene_order_file         = gene_order_file,
  ref_group_names         = ref_group_names,
  min_max_counts_per_cell = c(100, Inf),
  chr_exclude             = c("chrX", "chrY", "chrM")
)

## Step 6. Run inferCNV
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff                            = 0.1,
  out_dir                           = out_dir,
  cluster_by_groups                 = FALSE,
  denoise                           = TRUE,
  HMM                               = TRUE,
  analysis_mode                     = "subclusters",
  tumor_subcluster_partition_method = "leiden",
  leiden_resolution                 = 0.01,
  output_format                     = "pdf",
  num_threads                       = 16
)


# ==============================================================================
# Identification of tumor programs using NMF
# ==============================================================================
message(">>> Start NMF analysis for ", pat, " (K=", num, ")")

## Step 1. Load malignant cells & Extract Counts
metadata <- scObject@meta.data
cells    <- rownames(metadata)[metadata$malignant_overlap == "Both_malignant"]

if (Version(scObject)$major >= 5) {
  scObject <- JoinLayers(scObject, layers = "counts")
  expr_tumor_pure <- LayerData(scObject, layer = "counts")[, cells]
} else {
  expr_tumor_pure <- GetAssayData(scObject, slot = "counts")[, cells]
}

## Step 2. Filter malignant cells (< 500 detected genes excluded)
detected_genes  <- colSums(expr_tumor_pure > 0)
expr_tumor_pure <- expr_tumor_pure[, detected_genes >= 500]

## Step 3. Remove mitochondrial and ribosomal genes
keep_genes      <- !grepl("^MT-|^RPL|^RPS", rownames(expr_tumor_pure), ignore.case = FALSE)
expr_tumor_pure <- expr_tumor_pure[keep_genes, ]

## Step 4. Select top 5,000 genes by mean expression
gene_mean       <- Matrix::rowMeans(expr_tumor_pure)
top_genes       <- names(sort(gene_mean, decreasing = TRUE))[1:min(5000, length(gene_mean))]
expr_tumor_pure <- expr_tumor_pure[top_genes, ]

if (ncol(expr_tumor_pure) < 50) {
  stop(paste0("Sample ", pat, " contains fewer than 50 valid malignant cells after filtering. Exiting..."))
}

## Step 5. CPM normalization (10^6) and log2 transformation
cpm_matrix <- t(t(expr_tumor_pure) / colSums(expr_tumor_pure)) * 1000000
log_matrix <- log2(cpm_matrix + 1)

## Step 6. Mean centering & Set negative values to zero
row_means       <- rowMeans(log_matrix)
centered_matrix <- t(t(log_matrix) - row_means)
centered_matrix[centered_matrix < 0] <- 0

## Step 7. Run NMF (snmf/r algorithm)
nmf_matrix <- as.matrix(centered_matrix)
nmf_obj <- nmf(nmf_matrix, rank = num, method = "snmf/r", seed = 1)

##--- Save NMF Outputs
saveRDS(nmf_obj, paste0(pat, ".rank", num, ".NMF.rds"))
saveRDS(nmf_matrix, paste0(pat, ".rank", num, ".Matrix.rds"))
message(">>> Finish NMF for ", pat, " (K=", num, ")")


# ==============================================================================
# 10MPs-6TME Radar plots visualization
# ==============================================================================
message(">>> Start Figure 3E Radar Plot visualization...")

## Step 1. Define Color Palettes and Factor Levels
tme_colors <- c(
  "TLS" = "#33b39f", "Lymphoid" = "#f5af98", "Myeloid" = "#eb6f5d",
  "Steroid" = "#6376a0", "Stromal" = "#71c9dd", "Desert" = "#bebebe"
)
tme_levels <- names(tme_colors)

mp_colors <- c(
  "EMT" = "#e05555", "TGFB" = "#f68d36", "Hormone_response" = "#dc87ba",
  "IFN&MHC" = "#3186bd", "Metabolism" = "#fcbb78", "Angiogenesis" = "#f59796",
  "WNT&Stem" = "#2ca249", "Hypoxia" = "#65c6cc", "Proliferation" = "#bfc03b",
  "NFKB" = "#f3bd29"
)
mp_names <- names(mp_colors)

## Step 2. Data Preparation
df_radar <- as.data.frame(df_scaled)
colnames(df_radar)[1] <- "group"
df_radar$group        <- factor(df_radar$group, levels = tme_levels)

tmp <- data.frame(
  sector = factor(mp_names, levels = mp_names), y = 1,
  group  = factor(mp_names, levels = mp_names)
)

## Step 3. Define Outer Color Ring Plot
p_ring <- ggplot() +
  geom_bar(data = tmp, aes(x = sector, y = y, fill = group), stat = "identity", width = 1) +
  scale_fill_manual(values = mp_colors) +
  ylim(-11, 2.3) +
  coord_polar(start = -0.3) +
  theme_void() +
  theme(legend.position = "none")

## Step 4. Define Single Radar Plot Function (With strict bound alignment)
make_single_tme_radar <- function(df, tme_name) {
  df_one <- df %>% filter(group == tme_name)
  
  ggradar(
    df_one,
    axis.labels               = rep(NA, 10),
    grid.min                  = 0, 
    grid.mid                  = 0.5, 
    grid.max                  = 1,
    group.line.width          = 0.6,
    group.point.size          = 1,
    group.colours             = tme_colors[tme_name],
    fill                      = TRUE,
    fill.alpha                = 0.25,
    background.circle.colour  = "white",
    gridline.mid.colour       = "#2b8c96",
    legend.position           = "none",
    label.gridline.min        = FALSE,
    label.gridline.mid        = FALSE,
    label.gridline.max        = FALSE
  ) +
    theme(
      plot.background  = element_blank(),
      panel.background = element_blank(),
      plot.title       = element_text(hjust = 0.5, size = 14)
    ) +
    ggtitle(tme_name)
}

## Step 5. Combine Radar Plots and Outer Rings via Patchwork
radar_list <- lapply(tme_levels, function(tme) {
  p_radar <- make_single_tme_radar(df_radar, tme)
  p_ring + inset_element(p_radar, left = 0, bottom = 0, right = 1, top = 1.108)
})
names(radar_list) <- tme_levels
final_radar_plot   <- wrap_plots(radar_list, ncol = 3)

## Step 6. Save PDF Output
output_pdf <- paste0(pat, "_6TME_cluster_radar_chart.pdf")
pdf(output_pdf, width = 10, height = 6)
print(final_radar_plot)
dev.off()
message(">>> Figure 3E Radar Plot successfully saved to: ", output_pdf)


# ==============================================================================
#  Shannon’s index visualization (Lollipop Plot)
# ==============================================================================
message(">>> Start Figure S7G Shannon's index analysis...")

## Step 1. Load Functional Count Matrix
input_csv <- paste0(pat, "_cancertype_by_10Function_count_matrix.csv")
if (!file.exists(input_csv)) {
  input_csv <- "Tu_cancertype_by_10Function_count_matrix.csv"
}

df      <- read.csv(input_csv, header = TRUE, row.names = 1, check.names = FALSE)
df_norm <- df / rowSums(df)

## Step 2. Calculate Shannon Diversity Index for each MP
shannon_index <- apply(df_norm, 2, function(x) {
  vegan::diversity(x, index = "shannon")
})

plot_df <- data.frame(
  Function = names(shannon_index),
  Shannon  = as.numeric(shannon_index),
  stringsAsFactors = FALSE
)
plot_df          <- plot_df[order(plot_df$Shannon), ]
plot_df$Function <- factor(plot_df$Function, levels = plot_df$Function)

## Step 3. Generate Lollipop Plot using ggplot2
p_shannon <- ggplot(plot_df, aes(x = Shannon, y = Function, color = Function)) + 
  geom_segment(aes(x = 0, xend = Shannon, y = Function, yend = Function), linewidth = 2.5) +
  geom_point(size = 7) +
  scale_color_manual(values = mp_colors) + 
  labs(x = "Shannon Index", y = "") +
  theme_bw() +
  theme(
    legend.position   = "none", panel.border = element_blank(),
    axis.line         = element_line(color = "black"), panel.spacing = unit(0, "mm"),
    axis.text.x       = element_text(size = 15, color = "black"),
    axis.text.y       = element_text(size = 16, color = "black"),
    axis.title.x      = element_text(size = 18, color = "black"),
    panel.grid.major  = element_blank(), panel.grid.minor = element_blank(),
    plot.margin       = margin(t = 5, r = 5, b = 5, l = 5)
  )

## Step 4. Save Exportation to PDF
output_shannon_pdf <- paste0(pat, "_Shannon_plot.pdf")
ggsave(output_shannon_pdf, plot = p_shannon, width = 4.5, height = 5, useDingbats = FALSE)
message(">>> Figure S7G Shannon's index plot successfully saved to: ", output_shannon_pdf)


# ==============================================================================
#  Spearman’s correlation between MPs and TME Cell Types
# ==============================================================================
message(">>> Start Figure S8I Spearman's correlation analysis...")

## Step 1. Load Proportion/Expression Matrix
input_matrix_csv <- paste0(pat, "_Sample_clusters_Tu_normalized_matrix.csv")
if (!file.exists(input_matrix_csv)) {
  input_matrix_csv <- "New_Sample_clusters_Tu_normalized_matrix.csv"
}

mat <- read.csv(input_matrix_csv, row.names = 1, check.names = FALSE)
mat <- mat[complete.cases(mat), ]
mat <- as.matrix(mat)
storage.mode(mat) <- "numeric"

target_cells <- mp_names 

## Step 2. Compute Spearman Correlation
cor_result <- expand.grid(CellType1 = colnames(mat), CellType2 = colnames(mat), stringsAsFactors = FALSE) %>%
  filter(CellType1 != CellType2) %>%
  filter(!(CellType1 %in% target_cells & CellType2 %in% target_cells)) %>% 
  rowwise() %>%
  mutate(
    test = list(cor.test(mat[, CellType1], mat[, CellType2], method = "spearman")),
    R_value = as.numeric(test$estimate),
    P_value = test$p.value
  ) %>%
  ungroup() %>%
  mutate(
    P_adj      = p.adjust(P_value, method = "BH"),
    log10P     = -log10(P_value),
    log10P_adj = -log10(P_adj)
  ) %>%
  select(CellType1, CellType2, R_value, P_value, P_adj, log10P, log10P_adj)

## Step 3. Extract Top 20 Correlated Features for Each MP
top20_list <- lapply(target_cells, function(ct) {
  cor_result %>%
    filter(CellType1 == ct, !(CellType2 %in% target_cells)) %>%
    arrange(desc(abs(R_value))) %>%
    slice_head(n = 20) %>%
    mutate(Target = ct)
})
top20_df <- bind_rows(top20_list) %>%
  mutate(CellType1 = factor(CellType1, levels = target_cells))

## Step 4. Build Layered Multicolored Bar Plot using ggplot2
p_cor <- ggplot()

for (sub in names(mp_colors)) {
  sub_df   <- top20_df %>% filter(CellType1 == sub)
  base_col <- mp_colors[sub]
  
  p_cor <- p_cor +
    geom_bar(
      data      = sub_df,
      aes(x = R_value, y = tidytext::reorder_within(CellType2, R_value, CellType1), fill = log10P_adj),
      stat      = "identity", width = 0.7, color = "black", linewidth = 0.3
    ) +
    scale_fill_gradient(low = alpha(base_col, 0.15), high = base_col, guide = "none") +
    ggnewscale::new_scale_fill()
}

# Apply layouts
p_cor <- p_cor +
  facet_wrap(~ CellType1, nrow = 2, scales = "free") +
  tidytext::scale_y_reordered() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 12) +
  theme(
    panel.border       = element_blank(), axis.line = element_line(color = "black"),
    panel.spacing      = unit(2, "mm"),
    axis.text.x        = element_text(size = 10, color = "black"),
    axis.text.y        = element_text(size = 12, color = "black"),
    strip.text         = element_text(size = 14, face = "bold", color = "black"),
    strip.background   = element_blank(), legend.position = "none",
    panel.grid.major   = element_blank(), panel.grid.minor = element_blank(),
    plot.margin        = margin(t = 5, r = 5, b = 5, l = 5)
  ) +
  labs(x = "Spearman's Rho (R value)", y = "", title = "")

## Step 5. Save Exportation to PDF
output_cor_pdf <- paste0(pat, "_MPs_TME_correlation_barplot.pdf")
ggsave(filename = output_cor_pdf, plot = p_cor, width = 22, height = 9, useDingbats = FALSE)
message(">>> Figure S8I Correlation plot successfully saved to: ", output_cor_pdf)
