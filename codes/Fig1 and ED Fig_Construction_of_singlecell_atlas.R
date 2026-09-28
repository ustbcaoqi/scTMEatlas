# ==============================================================================
# Script Name: Figure1.R
# Description: Construction of single-cell atlas
# Author: Zhang Lab
# ==============================================================================




# ==============================================================================
# ED Fig. 5a. Circos plot
# ==============================================================================



### circos 

cancer_order <- c("MM","NPC","THCA","LEU","BRCA","CC","UCEC","NSCLC","CRC","GC","BTC","OSCC","PDAC",
                  "ESCC","HCC","RCC","FTC", "OV","OS","UC","NB","GBM","PC","AC","MELA","CTCL")


cancer_circos_r <- data.frame(cancer = cancer_order) %>% 
  mutate( n = 1:n(),
          r1 = 0.86 - 0.02*n,
          r2 = 0.84 - 0.02*n,
          # n_r = ifelse(cancer %in% c("LEU",""),1,0),
          r2 = case_when( n != 25 & n %% 5 == 0 ~ r2 + 0.006, .default = r2),
          r1 = case_when(n != 1 & n !=26 & n %% 5 == 1 ~ r1 - 0.006, .default = r1)) %>% 
  arrange(cancer)

options(digits = 3)
cancer_circos_r



#### circos --------------


# 1.准备数据

## 所有细胞类型，cluster含量，然后cancer_type富集比例，tissue Ro/e热图

seu_pan_meta_atlas <- seu_pan_meta %>% 
  filter(low_quality == "N") %>% 
  filter(major_celltypes_tier1 != "Epi_Tumor") %>% 
  filter(major_celltypes_tier3 != "Other") %>% 
  filter(major_celltypes_tier3 != "Epithelial") 
nrow(seu_pan_meta_atlas) 
seu_pan_meta_atlas[1:5,]



sam_clus <- seu_pan_meta_atlas %>% 
  group_by(Cohort,Sample,major_celltypes_tier3,clusters,CancerType_final,Tissue_final) %>% 
  summarise(n=n()) %>% 
  group_by(Sample) %>%
  mutate(all=sum(n),prop=n/all*100)

table(sam_clus$Sample) %>% names() %>% is.na() %>% sum
table(sam_clus$clusters) 
table(sam_clus$Cohort) 
table(sam_clus$clusters) %>% length()
sam_clus[1:5,]

table(sam_clus$Tissue_final)
tumor_clus_cancer <- sam_clus %>% 
  filter(Tissue_final %in% c("Tu")) %>% 
  select(Sample, clusters, prop) %>% 
  pivot_wider(names_from = Sample, values_from = prop,values_fill = 0) %>% 
  pivot_longer(cols = -1, names_to = "Sample",values_to = "prop") %>% 
  left_join(sam_info[,c("Sample","CancerType_final")], by = "Sample") %>% 
  group_by(CancerType_final,clusters) %>% 
  summarise(prop_median = median(prop) + 0.1,
            prop_mean = mean(prop)) %>% 
  group_by(clusters) %>% 
  mutate(mean=mean(prop_median),sd=sd(prop_median),
         prop_median_scale = (prop_median - mean)/sd,
         prop_median_scale = case_when(is.na(prop_median_scale) ~ 0.01,.default = prop_median_scale),
         prop_median_scale = case_when(prop_median_scale > 1.5 ~ 1.5,
                                       prop_median_scale < (-1.5) ~ -1.5,
                                       .default = prop_median_scale))

tumor_clus_cancer
write.csv(tumor_clus_cancer,"2.seu_summary/res.tumor_clus_cancer.prop.csv")


tumor_clus_cancer2 <- tumor_clus_cancer %>% 
  select(clusters,CancerType_final,prop_median_scale) %>% 
  pivot_wider(names_from = CancerType_final, values_from = prop_median_scale, values_fill = 0)


ggplot(tumor_clus_cancer, aes(x=clusters,y=CancerType_final)) + 
  geom_point(aes(size=abs(prop_median_scale),color=prop_median_scale)) +
  cowplot::theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
        axis.text = element_text(color='black')) + 
  scale_color_distiller(palette = "RdBu")




### karyocyte.txt 

"chr - hs1 1 0 248956422 chr1
chr - hs2 2 0 242193529 chr2
# 格式：chr - 染色体ID 编号 起始 终止 标签"

karyotype_atlas <- seu_pan_meta_atlas %>% 
  mutate(major_celltypes_tier3 = factor(major_celltypes_tier3, levels = names(tier3_color_panel))) %>% 
  group_by(major_celltypes_tier3) %>% 
  summarise(n=n()) %>% 
  mutate(chr="chr -", 
         id = 1:n(),
         start = 0,
         color=tier3_color_panel[major_celltypes_tier3] %>% gsub("#","",.) %>% tolower()) %>% 
  mutate(final = paste0(chr," ",major_celltypes_tier3," ",major_celltypes_tier3," ",start," ", n," ", 
                        major_celltypes_tier3," ","color=",color))
karyotype_atlas

write.table(karyotype_atlas$final,"4.TME/circos_atlas/karyotype.txt",quote = F,row.names = F, col.names = F)





## data_atlas_clus

data_atlas_clus <- seu_pan_meta_atlas %>% 
  mutate(major_celltypes_tier3 = factor(major_celltypes_tier3, levels = names(tier3_color_panel)),
         clusters = factor(clusters, levels = names(cluster_color_panel))) %>% 
  group_by(major_celltypes_tier3, clusters) %>% 
  summarise(cell_count=n()) %>% 
  # filter(cell_count > 20) %>%
  left_join(cluster_color_df[,c("clusters","cluster_colors","clus_id")], by="clusters") %>% 
  mutate(start = 0,
         end=cumsum(cell_count),
         start=end-cell_count,
         label=ifelse(cell_count> 1000, clusters,"" ),
         color=cluster_colors %>% gsub("#","",.) %>% tolower()) %>% 
  mutate(final = paste0(major_celltypes_tier3,"\t",start,"\t", end,"\t",clus_id,"\tcolor=",color),
         final2 = paste0(major_celltypes_tier3,"\t",start,"\t", end,"\t",clus_id,"\tcolor=black"))
data_atlas_clus


write.table(data_atlas_clus$final,"4.TME/circos_atlas/data_clus.txt",quote = F,row.names = F, col.names = F)
write.table(data_atlas_clus$final2,"4.TME/circos_atlas/data_clus2.txt",quote = F,row.names = F, col.names = F)





## cancer scaled median proportion 
data_atlas_clus_cancer <- data_atlas_clus %>% 
  left_join(tumor_clus_cancer2, by = "clusters")

cancers <- colnames(data_atlas_clus_cancer)[11:36]
for(cancer_coi in cancers){
  data_atlas_clus_cancer_coi <- data_atlas_clus_cancer %>% 
    mutate(final = paste0(major_celltypes_tier3,"\t",start,"\t", end,"\t",.data[[cancer_coi]]))
  data_atlas_clus_cancer_coi$final
  write.table(data_atlas_clus_cancer_coi$final,paste0("4.TME/circos_atlas/data_line_",cancer_coi,".txt"),quote = F,row.names = F, col.names = F)
}

data_atlas_clus_cancer[data_atlas_clus_cancer$clusters == "Tumor",]



## tissue Ro/e
Roe_df <- read.csv("2.seu_summary/res.Roe_clus_tissue.csv",row.names = "X")

tissues <- c("BM","PB","AN","Tu")
for(coi in tissues){
  data_atlas_clus_tissue <- data_atlas_clus %>% 
    left_join(Roe_df[Roe_df$Tissue_final == coi,c("clusters","Roe")], by = "clusters") %>% 
    mutate(final = paste0(major_celltypes_tier3,"\t",start,"\t", end,"\t",Roe))
  write.table(data_atlas_clus_tissue$final, paste0("4.TME/circos_atlas/data_heatmap_",coi,".txt"),quote = F,row.names = F, col.names = F)
}


cancers <- colnames(data_atlas_clus_cancer)[11:36]
for(cancer_coi in cancers){
  data_atlas_clus_cancer_coi <- data_atlas_clus_cancer %>% 
    mutate(final = paste0(major_celltypes_tier3,"\t",start,"\t", end,"\t",.data[[cancer_coi]]))
  data_atlas_clus_cancer_coi$final
  write.table(data_atlas_clus_cancer_coi$final,paste0("4.TME/circos_atlas/data_line_",cancer_coi,".txt"),quote = F,row.names = F, col.names = F)
}













# ==============================================================================
# ED Fig. 3e. Scatter plot for tissue distribution of 119 cell subpopulations
# ==============================================================================



### tumor enrichment and infiltrating level

cluster_info <- read.csv("2.seu_summary/cluster_info_sum_v2.csv",row.names = "X")

cluster_info <- cluster_info %>% 
  filter(! major_celltypes_tier3 %in% c("Other","Epithelial"))
cluster_info[1:5,]

table(seu_pan_meta$Tissue_final)
table(seu_pan_meta$major_celltypes_tier3)
table(seu_pan_meta$Cohort)
table(sam_info$Tissue_final)

seu_pan_meta[1:5,]

clus_sample_prop <- seu_pan_meta %>% 
  filter(Tissue_final %in% c("Tu","Me","TT")) %>%
  mutate(Tissue_final = case_when( Tissue_final %in% c("Me","TT") ~ "Tu", .default = Tissue_final)) %>% 
  filter(! major_celltypes_tier3 %in% c("Other","Epithelial")) %>%
  group_by(Cohort, Sample, CancerType_final,clusters) %>% 
  summarise(Cell_count = n() ) %>% 
  group_by(Cohort, Sample, CancerType_final) %>% 
  mutate(cell_sum=sum(Cell_count)) %>% 
  filter(cell_sum > 1000) %>% 
  mutate(prop=Cell_count / sum(Cell_count)*100)  %>% 
  left_join(cluster_info[,c("clusters","major_celltypes_tier3")],by="clusters")
clus_sample_prop[1:5,]
table(clus_sample_prop$Sample) %>% length()
table(clus_sample_prop$clusters) %>% length()
table(clus_sample_prop$CancerType_final) %>% length() # 21
write.csv(clus_sample_prop, "2.seu_summary/res.clus_sample_prop.csv")


clus_sample_prop <- read.csv("2.seu_summary/res.clus_sample_prop.csv",row.names = "X")

clus_cancer_infilt <- clus_sample_prop %>% 
  group_by(CancerType_final, clusters) %>% 
  summarise(prop_median = median(prop)) %>% 
  group_by(clusters) %>% 
  mutate(prop_scale = scale(prop_median)) %>% 
  left_join(cluster_info[,c("clusters","major_celltypes_tier3")],by="clusters") %>% 
  mutate(major_celltypes_tier3 = factor(major_celltypes_tier3, levels = names(tier3_color_panel))) %>% 
  filter(major_celltypes_tier3 != "Tumor")
clus_cancer_infilt[1:5,]
table(clus_cancer_infilt$CancerType_final) %>% length()
table(clus_cancer_infilt$clusters) %>% length()


clus_tu_infilt <- clus_cancer_infilt %>% 
  group_by(clusters) %>% 
  summarise(prop_median_all_avg = median(prop_median)) %>% 
  ungroup()


Roe_df <- read.csv("2.seu_summary/res.Roe_clus119_tissue.csv",row.names = "X")
Roe_df[1:5,]
Roe_df_tu <- Roe_df %>% filter(Tissue_final == "Tu")
Roe_df_tu


clus_tu_infilt_roe <- clus_tu_infilt %>% 
  left_join(Roe_df_tu[,c("clusters","Roe","major_celltypes_tier3")],by = "clusters") %>% 
  mutate(major_celltypes_tier3 = factor(major_celltypes_tier3, levels = names(tier3_color_panel)),
         show = case_when( clusters %in% c( 
           "CD4_12_Treg_TNFRSF9",
           "CD8_03_Teff_GZMK",
           "B_04_Bm_GPR183",
           "DC_02_cDC2_CD1C",
           "DC_04_mreg_LAMP3",
           "Peri_01_ADGRF5",
           "Fb_01_myCAF_FAP",
           "Mph_08_TREM2",
           "Mph_09_SPP1",  
           "Mph_07_CCL4",
           "Neu_16_FKBP5",
           "Neu_17_CCL4",
           "Neu_19_VEGFA" ) ~ "Y",
           .default = "N"))
table(clus_tu_infilt_roe$show)
length(clus_tu_infilt_roe$clusters)

library(ggrepel)
ggplot(clus_tu_infilt_roe, aes(y=log10(prop_median_all_avg), x=(Roe))) + 
  geom_point(aes(color=major_celltypes_tier3),size=3, show.legend = F) + 
  geom_text_repel(data = clus_tu_infilt_roe   %>% filter(show == "Y"), 
                  aes(label=clusters),
                  force = 100,          
                  force_pull = 1,       
                  max.overlaps = Inf,   
                  min.segment.length = 0, 
                  seed = 123,           
                  max.iter = 10000,     
                  box.padding = 1,      
                  point.padding = 1,   
                  size = 3) + 
  scale_color_manual(values = tier3_color_panel) + 
  labs(y="Median proportion of clusters\n in tumor (%)", x="Roe in tumor") + 
  cowplot::theme_cowplot() 

ggsave("0.figure/Fig1_dot_clus_tu_infilt_roe.pdf", width = 6,height = 4)





# ==============================================================================
# ED Fig. 4. Characteristics of tumor-associated cell subpopulations
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggh4x)
})

# 1. Path Definition and Safety Checks
input_main  <- "data/processed/figure_s3_main_summary.csv"
input_roe   <- "data/processed/figure_s3_roe_tissue.csv"
input_prop  <- "data/processed/figure_s3_cancer_proportion.csv"
output_pdf  <- "results/figures/FigureS3_tumor_associated_cell_subpopulations.pdf"

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input_main), file.exists(input_roe), file.exists(input_prop))

panel_label <- c(
  tissue     = "Tissue preference",
  cancer     = "Distribution across cancer types",
  spec       = "Cancer type specificity", 
  tcga_prog  = "Prognosis in TCGA",
  tcga_stage = "Progression in TCGA",
  study_prog = "Prognosis in this study",
  study_stage= "Progression in this study"
)

to_numeric <- function(x) {
  x_clean <- if_else(as.character(x) %in% c("#N/A", "NA", ""), NA_character_, as.character(x))
  as.numeric(x_clean)
}

star_label <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# 2. Main Dataset Processing
main_raw <- read.csv(input_main, stringsAsFactors = FALSE, check.names = FALSE)
if (ncol(main_raw) < 12) stop("Main input table must contain at least 12 columns.")

main_df <- tibble(
  cell          = main_raw[[1]],
  facet         = main_raw[[3]],
  retain        = main_raw[[2]],
  spec          = to_numeric(main_raw[[4]]),
  tcga_prog     = to_numeric(main_raw[[5]]),
  tcga_prog_p   = to_numeric(main_raw[[6]]),
  tcga_stage    = to_numeric(main_raw[[7]]),
  tcga_stage_p  = to_numeric(main_raw[[8]]),
  study_prog    = to_numeric(main_raw[[9]]),
  study_prog_p  = to_numeric(main_raw[[10]]),
  study_stage   = to_numeric(main_raw[[11]]),
  study_stage_p = to_numeric(main_raw[[12]]),
  order         = row_number()
) %>%
  filter(retain == "Y") %>%
  arrange(order) %>%
  mutate(
    cell   = factor(cell, levels = rev(unique(cell))),
    facet  = factor(facet, levels = unique(facet)),
    s_tcga_prog   = star_label(tcga_prog_p),
    s_tcga_stage  = star_label(tcga_stage_p),
    s_study_prog  = star_label(study_prog_p),
    s_study_stage = star_label(study_stage_p),
    c_tcga_prog   = if_else(tcga_prog > 0, log2(tcga_prog), NA_real_),
    c_tcga_stage  = tcga_stage,
    c_study_prog  = if_else(study_prog > 0, log2(study_prog), NA_real_),
    c_study_stage = study_stage
  )

cell_levels  <- levels(main_df$cell)
facet_levels <- levels(main_df$facet)

theme_box <- theme_bw() +
  theme(
    panel.grid       = element_blank(),
    strip.text.y     = element_blank(),
    strip.background = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.y     = element_blank(),
    axis.text.x      = element_text(size = 7, color = "black"),
    axis.text.x.top  = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_text(size = 10, margin = margin(b = 10)),
    panel.spacing.y  = unit(0.25, "lines"),
    plot.margin      = margin(2, 2, 2, 2)
  )

# 3. Core Plotting Function
plot_metric <- function(data, x, color = NULL, star = NULL, label, center = 0, bar = FALSE) {
  p <- ggplot(data, aes(y = cell))
  
  if (bar) {
    p <- p +
      geom_col(aes(x = .data[[x]], fill = .data[[x]])) +
      scale_fill_gradient(
        low = "#f0f6fb", high = "#2C7CBB",
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",
                               barwidth = unit(0.7, "cm"), barheight = unit(3.2, "cm"), title.position = "top")
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.05)), sec.axis = dup_axis(name = label))
  } else {
    valid_vals <- data[[x]][!is.na(data[[x]])]
    max_d <- if (length(valid_vals) > 0) max(abs(valid_vals - center), na.rm = TRUE) else 1
    if (!is.finite(max_d) || max_d == 0) max_d <- 1
    
    p <- p +
      geom_vline(xintercept = center, linetype = "dashed", color = "grey50") +
      geom_segment(aes(x = center, xend = .data[[x]], yend = cell), color = "grey80") +
      geom_point(aes(x = .data[[x]], color = .data[[color]]), size = 3.5) +
      geom_text(aes(x = .data[[x]], label = .data[[star]]), vjust = 0.8, size = 3) +
      scale_color_gradient2(
        low = "#3B6FB6", mid = "white", high = "#B63B3B", midpoint = 0,
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",
                               barwidth = unit(0.7, "cm"), barheight = unit(3.2, "cm"), title.position = "top")
      ) +
      scale_x_continuous(limits = c(center - max_d, center + max_d),
                         expand = expansion(mult = c(0.05, 0.05)), sec.axis = dup_axis(name = label))
  }
  
  p + labs(x = NULL, y = NULL) +
    facet_grid(facet ~ ., scales = "free_y", space = "free_y") +
    theme_box
}

p_bar <- plot_metric(main_df, "spec", label = panel_label["spec"], bar = TRUE)
p1    <- plot_metric(main_df, "tcga_prog", "c_tcga_prog", "s_tcga_prog", panel_label["tcga_prog"], center = 1)
p2    <- plot_metric(main_df, "tcga_stage", "c_tcga_stage", "s_tcga_stage", panel_label["tcga_stage"], center = 0)
p3    <- plot_metric(main_df, "study_prog", "c_study_prog", "s_study_prog", panel_label["study_prog"], center = 1)
p4    <- plot_metric(main_df, "study_stage", "c_study_stage", "s_study_stage", panel_label["study_stage"], center = 0)

# 4. Tissue Preference (RO/E) Plot
strip_colors <- c(
  CD8T = "#CCEBC5", CD4T = "#66C2A5", NK = "#FFED6F", ILC = "#B3DE69", B = "#FD8D3C",
  Monocyte = "#4DAF4A", Neutrophil = "#E41A1C", Macrophage = "#BC80BD", DC = "#FCCDE5",
  Mast = "#FDB462", EC = "#BEBADA", Fb = "#6A3D9A"
)

if (any(!facet_levels %in% names(strip_colors))) {
  stop("Missing strip colors for: ", paste(setdiff(facet_levels, names(strip_colors)), collapse = ", "))
}

roe_raw <- read.csv(input_roe, stringsAsFactors = FALSE, check.names = FALSE)
if (ncol(roe_raw) < 5) stop("Tissue preference input table must contain at least 5 columns.")

roe_df <- tibble(
  cell   = roe_raw[[1]],
  tissue = roe_raw[[2]],
  facet  = roe_raw[[3]],
  roe    = roe_raw[[4]],
  label  = roe_raw[[5]]
) %>%
  filter(cell %in% cell_levels) %>%
  mutate(
    cell   = factor(cell, levels = cell_levels),
    tissue = factor(tissue, levels = c("BM", "PB", "AN", "Tu")),
    facet  = factor(facet, levels = facet_levels)
  )

p_roe <- ggplot(roe_df, aes(x = tissue, y = cell)) +
  geom_tile(aes(fill = roe)) +
  geom_text(aes(label = label), size = 2.5) +
  scale_x_discrete(expand = c(0, 0), sec.axis = dup_axis(name = panel_label["tissue"])) +
  scale_fill_stepsn(
    name = "RO/E",
    colors = c("#f0f6fb", "#85B7D9", "#2C7CBB", "#08306B"),
    breaks = c(0, 1, 2, 3, 5),
    labels = c("", "+/-", "+", "++", "+++"),
    guide = guide_colorsteps(frame.colour = "black", ticks.colour = "black",
                             barheight = unit(3.2, "cm"), barwidth = unit(0.7, "cm"),
                             ticks = FALSE, even.steps = TRUE, show.limits = FALSE)
  ) +
  facet_grid2(
    facet ~ ., scales = "free_y", space = "free_y", switch = "y",
    strip = strip_themed(background_y = elem_list_rect(fill = unname(strip_colors[facet_levels]), color = "black", linewidth = 0.6))
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    panel.grid          = element_blank(),
    panel.border        = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text.y         = element_text(size = 9, color = "black"),
    axis.text.x         = element_text(size = 9, color = "black"),
    axis.text.x.top     = element_blank(),
    axis.ticks.x.top    = element_blank(),
    axis.title.x.top    = element_text(size = 10, margin = margin(b = 10)),
    strip.placement     = "inside",
    strip.background    = element_blank(),
    strip.text.y.left   = element_text(angle = 90, size = 8, margin = margin(0, 0, 0, 0), color = "black"),
    panel.spacing.y     = unit(0.25, "lines"),
    legend.text         = element_text(size = 10, face = "bold"),
    plot.margin         = margin(2, 2, 2, 2)
  )

# 5. Cancer Proportion Distribution Plot
prop_raw <- read.csv(input_prop, stringsAsFactors = FALSE, check.names = FALSE)
if (ncol(prop_raw) < 5) stop("Cancer distribution input table must contain at least 5 columns.")

cancer_order <- c(
  "MM", "NPC", "THCA", "LEU", "BRCA", "CC", "UCEC", "NSCLC", "CRC", "GC",
  "BTC", "OSCC", "PDAC", "ESCC", "HCC", "RCC", "FTC", "OV", "OS", "UC",
  "NB", "GBM", "PC", "AC", "MELA", "CTCL"
)

prop_df <- tibble(
  cell        = prop_raw[[1]],
  cancer      = prop_raw[[2]],
  facet       = prop_raw[[3]],
  prop        = prop_raw[[4]],
  scaled_prop = prop_raw[[5]]
) %>%
  filter(cell %in% cell_levels) %>%
  mutate(
    cell        = factor(cell, levels = cell_levels),
    cancer      = factor(cancer, levels = cancer_order),
    facet       = factor(facet, levels = facet_levels)
  )

p_dot <- ggplot(prop_df, aes(x = cancer, y = cell)) +
  geom_point(aes(size = prop, color = scaled_prop)) +
  scale_color_gradient(
    low = "#dae8f5", high = "#2a79b9", name = "Scaled\nProportion",
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barwidth = unit(0.7, "cm"), barheight = unit(3.2, "cm"), title.position = "top")
  ) +
  scale_size(
    range = c(0.2, 5), name = "Median\nProportion", breaks = c(1, 10, 20, 25),
    guide = guide_legend(title.position = "top", override.aes = list(shape = 21, fill = "transparent", color = "black", stroke = 0.8), order = 1)
  ) +
  scale_x_discrete(sec.axis = dup_axis(name = panel_label["cancer"])) +
  facet_grid(facet ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_box +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, color = "black"),
    legend.key       = element_blank(),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    legend.spacing.y = unit(0.2, "cm")
  )

# 6. Figure Assembly and Export
final_plot <- (p_roe + p_dot + p_bar + p1 + p2 + p3 + p4) +
  plot_layout(ncol = 7, widths = c(0.8, 3, 0.6, 1, 1, 1, 1), guides = "collect")

ggsave(output_pdf, final_plot, width = 12, height = 15.5)

