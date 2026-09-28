# ==============================================================================
# Script Name: Figure4.R
# Description: mouse neutrophils
# Author: Zhang Lab
# ==============================================================================


# ==============================================================================
# Figure 4. tissue distribution
# ==============================================================================


seu_mm_neu <- rbind(seu_mm_neu_umap@meta.data[,c("cluster_final_neu2","tissue_combined_final2")],
                    seu_mm_neu_hly@meta.data[,c("cluster_final_neu2","tissue_combined_final2")])


plot_df <- seu_mm_neu%>% 
  group_by(cluster_final_neu2, tissue_combined_final2) %>% 
  summarise(n=n()) %>% 
  group_by(cluster_final_neu2) %>% 
  mutate(prop= n/sum(n) *100 ) %>% 
  mutate(tissue_combined_final2 = factor(tissue_combined_final2, levels = names(mm_tissue_color_panel) %>% rev ))


ggplot(plot_df, aes( y=cluster_final_neu2, x=prop,fill=tissue_combined_final2)) + 
  geom_bar(stat = "identity") + 
  labs(x="Proportion",y=NULL, fill ="Tissue") + 
  scale_y_discrete(limits= rev) + 
  scale_fill_manual(values = mm_tissue_color_panel) + 
  cowplot::theme_cowplot() + 
  theme(
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "plain"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_blank(),
    # aspect.ratio = 1,  # 设置纵横比,
    legend.position = "top",
    # legend.position = c(1, 0), legend.justification = c(1, 0)
  )


ggsave("0.figure/Fig_mm_neu_bar_tissue.pdf",width = 3.2,height = 4.2)






