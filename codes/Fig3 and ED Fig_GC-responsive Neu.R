# ==============================================================================
# Script Name: Figure3.R
# Description: GC-responsive neutrophils
# Author: Zhang Lab
# ==============================================================================




# ==============================================================================
# Fig. 3k. Volcanol plot
# ==============================================================================



degs <- read.csv("2.seu_neu/res_degs_Neu_GCneuvsOthers.csv")
degs

# 添加显著性标签列
degs$significance <- ifelse(
  degs$pvals_adj < 0.05 & abs(degs$logfoldchanges) > 1,
  ifelse(degs$logfoldchanges > 0.5, "Up", "Down"),
  "Not Sig"
)

# 绘制火山图
# 选择要标记的基因 - 这里选择最显著的10个基因
genes_to_label <- degs %>%
  arrange(-abs(scores)) %>%
  filter( abs(logfoldchanges) > 1) %>% 
  head(20) 


# 绘制火山图
library(ggrepel)
ggplot(degs, aes(x = logfoldchanges, y = abs(scores), color = significance)) +
  geom_point(alpha = 1, size = 1) +
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "Not Sig" = "gray")) +
  geom_text_repel(
    data = genes_to_label,
    aes(label = gene),
    size = 2.5,
    box.padding = 0.35,  # 控制标签周围的空白区域
    point.padding = 0.3, # 控制标签与点之间的空白区域
    segment.color = 'grey50',  # 连接线的颜色
    segment.size = 0.2,        # 连接线的粗细
    max.overlaps = Inf,        # 允许无限重叠尝试
    min.segment.length = 0,    # 总是绘制连接线
    force = 1,                # 调整标签间的排斥力
    nudge_x = ifelse(genes_to_label$logfoldchanges > 0, 0.5, -0.5), # 根据方向调整初始位置
    direction = "both"           # 主要沿y轴方向调整
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "Abs scores",
    color = "Significant",
    # title = "Volcano Plot of Differential Expression"
  ) + 
  # 添加阈值线
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(-6,6)) + 
  # 美化主题
  theme_bw() +
  theme(
    # 添加边框
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    aspect.ratio = 1,
    legend.position = "top",
    # panel.grid.major = element_line(color = "grey90"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

ggsave(paste0(fdir,"Fig_neu_DEG_GCneuvsOther.pdf"),
       width = 4,
       height = 4)









# ==============================================================================
# ED Fig. 12d. TF plot
# ==============================================================================




seu_Pan_neu_coi <- seu_Pan_neu@meta.data %>% 
  filter(Tissue_final %in% c("Tu","Me")) %>% 
  filter( MP_neu %in% names(MP_neu_color_panel)[3:6] ) %>% 
  mutate(MP_GC = case_when(MP_neu != "MP5_GC" ~ "Other", .default = "MP5_GC"),
         MP_neu = factor(MP_neu, levels = names(MP_neu_color_panel)[3:6]))

seu_Pan_neu_MP <- seu_Pan_neu[,rownames(seu_Pan_neu_coi)]
seu_Pan_neu_MP$MP_GC <- seu_Pan_neu_coi$MP_GC

Idents(seu_Pan_neu_MP) <- seu_Pan_neu_MP$MP_GC
markers <- FindAllMarkers(seu_Pan_neu_MP, only.pos = F)
neu_GC_DEG <- markers
neu_GC_DEG[1:5,]


TF_hs <- read.table("cisTarget/hs_hgn_curated_tfs.txt") 

neu_GC_DEG_TF  <- neu_GC_DEG %>% 
  filter(cluster == "MP5_GC") %>% 
  filter(gene %in% TF_hs$V1) %>% 
  arrange(-avg_log2FC) %>% 
  mutate( id=1:n() ) %>% 
  mutate( label = case_when((avg_log2FC) > 1.4 ~ gene, 
                            (avg_log2FC) < -1.7 ~ gene, 
                            .default = NA) )

neu_GC_DEG_TF[1:20,]
nrow(neu_GC_DEG_TF)

ggplot(neu_GC_DEG_TF, aes(x=id, y=avg_log2FC)) + 
  geom_point() + 
  geom_label_repel(aes(label = label),
                   box.padding = 0.5,
                   max.overlaps = 20,
                   na.rm = T,
                   size = 3,
                   direction = "both",         # 允许双向移动
                   segment.size = 0.2,         # 连接线粗细
                   show.legend = FALSE) + 
  labs(x="Transcriptional Factor",y="Log2 Fold Change") + 
  cowplot::theme_cowplot()

ggsave("0.figure/Fig_neu_GC_TF_log2FC.pdf", width = 3, height = 4)


