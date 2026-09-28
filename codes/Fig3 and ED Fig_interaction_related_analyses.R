# ==============================================================================
# Script Name: Figure4.R
# Description: Ligand-receptor interactions across tumor microenvironment
# Author: Zhang Lab
# ==============================================================================

# ==============================================================================
# ED Fig. 14b: Chord plots illustrating cell-cell interactions across TME subtypes
# ==============================================================================

library(CellChat)

# 1. Path Definitions (Using Relative Paths for GitHub Reproducibility)
# Note: Ensure the files follow consistent case-sensitivity for Linux compatibility
file_list <- c(
  TLS      = "data/processed/rank_TLS_interaction_tier3.csv",
  Lymphoid = "data/processed/rank_Lymphoid_interaction_tier3.csv",
  Myeloid  = "data/processed/rank_Myeloid_interaction_tier3.csv",
  Steroid  = "data/processed/rank_Steroid_interaction_tier3.csv",
  Stromal  = "data/processed/rank_Stromal_interaction_tier3.csv", 
  Desert   = "data/processed/rank_Desert_interaction_tier3.csv"
)

output_pdf <- "results/figures/New_Rank_6TME_network_plot.pdf"
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

# 2. Canonical Cell Type Ordering
celltype_order <- c(
  "B", "CD4", "CD8", "NK", "DC", "Macrophage", 
  "Neutrophil", "Mast", "EC", "Fb", "Tumor"
)

# 3. Defined Color Palette for Cell Lineages
cell_colors <- c(
  "B"          = "#f08940",
  "CD4"        = "#64ba9f",
  "CD8"        = "#c9e3c1",
  "NK"         = "#f7e673",
  "DC"         = "#f6cbde",
  "Macrophage" = "#b37cb3",
  "Neutrophil" = "#d6251f",
  "Mast"       = "#f4af63",
  "EC"         = "#bcb7d6",
  "Fb"         = "#653b90",
  "Tumor"      = "#2572a9"
)

# 4. Multi-panel Visualization
pdf(output_pdf, width = 14, height = 3.5)

# Configure a 1-row, 6-column layout matrix
par(
  mfrow = c(1, 6),
  mar = c(1, 2.5, 4, 1) # Bottom, Left, Top, Right
)

for (nm in names(file_list)) {
  
  if (!file.exists(file_list[nm])) {
    stop(paste("Required input file not found for subtype:", nm, "-", file_list[nm]))
  }
  
  # Load interaction strength matrix
  interaction_raw <- read.csv(file_list[nm], row.names = 1, check.names = FALSE)
  interaction_matrix <- as.matrix(interaction_raw)
  
  # Robust filtering: Keep only cell types present in both the matrix and the canonical list
  valid_cells <- celltype_order[celltype_order %in% rownames(interaction_matrix) & 
                                  celltype_order %in% colnames(interaction_matrix)]
  
  if (length(valid_cells) == 0) {
    stop(paste("No matching canonical cell types found in dataset for:", nm))
  }
  
  # Align and subset matrix to the specified order safely
  sub_matrix <- interaction_matrix[valid_cells, valid_cells, drop = FALSE]
  color_use  <- cell_colors[rownames(sub_matrix)]
  
  # Render Circle Plot via CellChat
  netVisual_circle(
    sub_matrix,
    weight.scale      = TRUE,
    edge.width.max    = 4,
    arrow.size        = 0.35,
    label.edge        = FALSE,
    vertex.label.cex  = 1.3,
    color.use         = color_use,
    title.name        = paste0("TME-", nm)
  )
}

dev.off()







# ==============================================================================
# Figure S9B-E and H. Half-circle dot plots
# ==============================================================================


library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)

#' 批量计算多个基因在各细胞类型中的表达分布与贡献度
#' @param obj Seurat对象
#' @param genes 基因名称向量 (如 c("CCL19", "CXCL13", "CD8A"))
#' @param celltype_col 细胞类型的元数据列名
#' @param exclude_types 需要剔除的细胞类型，默认为 NULL
#' @param assay 使用的检测类型，默认为 "RNA"
#' @param slot 使用的数据槽，默认为 "counts"

calculate_gene_contribution <- function(obj, 
                                        genes, 
                                        celltype_col = "major_celltypes_tier3", 
                                        exclude_types = NULL, 
                                        assay = "RNA", 
                                        slot = "count") {
  
  # 1. 检查基因是否存在，并筛选出存在的基因
  existing_genes <- genes[genes %in% rownames(GetAssayData(obj, assay = assay, slot = slot))]
  missing_genes <- setdiff(genes, existing_genes)
  
  if (length(missing_genes) > 0) {
    warning(paste("以下基因未在对象中找到:", paste(missing_genes, collapse = ", ")))
  }
  
  if (length(existing_genes) == 0) {
    stop("提供的基因列表均不在对象中，请检查输入。")
  }
  
  # 2. 提取所有基因的表达数据
  # FetchData 返回一个 Data.Frame，行是 Cell ID，列是 Gene
  exp_matrix <- FetchData(obj, vars = existing_genes, slot = slot)
  
  # 3. 合并元数据并转换为长矩阵格式进行统一处理
  plot_df_long <- obj@meta.data %>%
    select(all_of(celltype_col)) %>%
    bind_cols(exp_matrix) %>%
    # 转换为长表：Gene | Expression
    pivot_longer(cols = all_of(existing_genes), names_to = "gene", values_to = "expression")
  
  # 4. 核心逻辑计算
  result <- plot_df_long %>%
    # 过滤掉不需要的细胞类型
    filter(!(.data[[celltype_col]] %in% exclude_types)) %>%
    # 判定阳性细胞 (Expression > 0)
    mutate(is_positive = ifelse(expression > 0, 1, 0)) %>%
    # 按照 基因 + 细胞类型 + 是否阳性 分组统计
    group_by(gene, .data[[celltype_col]], is_positive) %>%
    summarise(n = n(), .groups = 'drop') %>%
    # --- 计算指标 1: 组内阳性率 (Within-type Percentage) ---
    group_by(gene, .data[[celltype_col]]) %>%
    mutate(
      type_total = sum(n),
      pct_within_type = n / type_total * 100
    ) %>%
    # 只分析阳性细胞部分
    filter(is_positive == 1) %>%
    # --- 计算指标 2: 全局贡献度 (Across-type Contribution) ---
    group_by(gene) %>%
    mutate(
      gene_positive_total = sum(n),
      pct_contribution = n / gene_positive_total * 100
    ) %>%
    ungroup() %>%
    # 整理输出列
    select(
      gene,
      cell_type = !!sym(celltype_col),
      positive_cells = n,
      type_total_cells = type_total,
      pct_within_type,      # 该基因在该类群中的阳性率
      pct_contribution      # 该类群在所有阳性细胞中的占比
    ) %>%
    arrange(gene, desc(pct_contribution))
  
  return(result)
}

# --- 使用方式 ---
my_genes <- c("CCL19", "CCL21", "CXCL13")
res_multi <- calculate_gene_contribution(
  obj = seu_pan_0.5M,
  genes = my_genes,
  celltype_col = "major_celltypes_tier3",
  exclude_types = "Epithelial"
)







# ==============================================================================
# Figure S10. dot plots
# ==============================================================================



calculate_gene_TPK <- function(obj, 
                               genes, 
                               celltype_col = "major_celltypes_tier3", 
                               exclude_types = NULL, 
                               assay = "RNA", 
                               slot = "count") {
  
  # 加载必要的包
  library(dplyr)
  library(tidyr)
  library(Seurat)
  
  # 1. 检查基因是否存在，并筛选出存在的基因
  existing_genes <- genes[genes %in% rownames(GetAssayData(obj, assay = assay, slot = slot))] %>% unique()
  missing_genes <- setdiff(genes, existing_genes)
  
  if (length(missing_genes) > 0) {
    warning(paste("以下基因未在对象中找到:", paste(missing_genes, collapse = ", ")))
  }
  
  if (length(existing_genes) == 0) {
    stop("提供的基因列表均不在对象中，请检查输入。")
  }
  
  # 2. 提取所有基因的表达数据
  exp_matrix <- FetchData(obj, vars = existing_genes, slot = slot)
  
  # 3. 合并元数据
  meta_data <- obj@meta.data %>% 
    select(all_of(celltype_col))
  
  plot_df_long <- meta_data %>%
    bind_cols(exp_matrix) %>%
    # 转换为长表格式
    pivot_longer(cols = all_of(existing_genes), 
                 names_to = "gene", 
                 values_to = "count") %>% 
    mutate(gene = factor(gene, levels = existing_genes))
  
  # 4. 计算TPK
  result <- plot_df_long %>%
    # 过滤不需要的细胞类型
    {if(!is.null(exclude_types)) 
      filter(., !(.data[[celltype_col]] %in% exclude_types)) 
      else .} %>%
    # 按照基因和细胞类型分组
    group_by(gene, .data[[celltype_col]]) %>%
    summarise(
      # 计算该基因在该细胞类型中的总表达量
      n_cells_celltype = n(),
      count_celltype = sum(count, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    # 可选：添加每细胞平均表达量
    group_by(gene) %>%
    mutate(
      # 计算细胞数量
      n_cells = sum(n_cells_celltype),
      total_count = sum(count_celltype),
      total_count_per_1000cells = (total_count / n_cells) * 1000,
      count_celltype_per_100cells = count_celltype / n_cells * 100,
      # 计算TPK：每1000个细胞中的表达量
      TPK = (count_celltype / n_cells) * 1000,
      TPK_pct = TPK / total_count_per_1000cells * 100
    ) %>%
    # 重新排序列
    select(
      gene,
      cell_type = !!sym(celltype_col),
      n_cells_celltype,
      n_cells,
      total_count,
      total_count_per_1000cells,
      count_celltype,
      count_celltype_per_100cells,
      TPK,
      TPK_pct
    ) %>%
    # 按照基因和TPK降序排列
    arrange(gene, desc(TPK))
  
  return(result)
}

# --- 使用示例 ---
my_genes <- c("CCL19", "CCL21", "CXCL13")
tpk_result <- calculate_gene_TPK(
  obj = seu_pan_0.5M,
  genes = my_genes,
  celltype_col = "major_celltypes_tier3",
  exclude_types = "Epithelial"
)
tpk_result
# 结果解读：
# - n_cells_celltype: 该细胞类型的细胞总数
# - n_cells: 所有细胞类型的细胞总数
# - total_count_celltype: 该基因在该细胞类型中的总表达分子数
# - avg_count_per_cell: 平均每个细胞的表达量
# - TPK: 每1000个细胞中,特定细胞类型的表达分子数；所有加起来等于count_per_1000cells










