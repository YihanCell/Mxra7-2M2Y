library(RColorBrewer)
library(ggplot2)
setwd(dir = '/data/yihan/Mxra7_2m/')
combinedh_afterHMCluster.2y = readRDS("/data/yihan/Mxra7_2Y/02.Merge.Result/RDS/combinedh_afterHMCluster.2y.rds")
meta.2y = combinedh_afterHMCluster.2y@meta.data
head(meta.2y)

new.cluster.ids = c("Granulocytes","Granulocytes","Granulocytes","Monocytes","Granulocytes",
                     "Monocytes","Monocytes","Granulocytes","T cells","Granulocytes",
                     "B cells","Granulocytes","B cells","Erythrocytes","Granulocytes",
                     "Monocytes","Monocytes","B cells","B cells","Granulocytes",
                     "Granulocytes","Monocytes","Granulocytes","Granulocytes","Granulocytes",
                     "T cells","Monocytes","Granulocytes","Monocytes","Monocytes",
                     "Monocytes","B cells","Granulocytes","Monocytes","Granulocytes")
names(new.cluster.ids) = levels(combinedh_afterHMCluster.2y)
combinedh_afterHMCluster.2y = RenameIdents(combinedh_afterHMCluster.2y, new.cluster.ids)
DimPlot(combinedh_afterHMCluster.2y, reduction = "umap") + scale_color_manual(values = brewer.pal(6, "Paired"))
ggsave('./Result_for_2Y/2Y_Annotation.pdf', plot = DimPlot(combinedh_afterHMCluster.2y, reduction = "umap") + scale_color_manual(values = brewer.pal(6, "Paired")), height = 5, width = 5, dpi = 300)
DimPlot(combinedh_afterHMCluster.2y, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Set2")) + NoLegend()
ggsave('./Result_for_2Y/2Y_Annotation_split_nolegend.pdf', plot = DimPlot(combinedh_afterHMCluster.2y, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Paired")) + NoLegend(), height = 5, width = 10, dpi = 300)
DimPlot(combinedh_afterHMCluster.2y, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Set2"))
ggsave('./Result_for_2Y/2Y_Annotation_split.pdf', plot = DimPlot(combinedh_afterHMCluster.2y, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Paired")) , height = 5, width = 5, dpi = 300)
combinedh_afterHMCluster.2y@meta.data$predicted_celltype = combinedh_afterHMCluster.2y@active.ident
ggsave('./Result_for_2Y/2Y_Annotation_nolegend.pdf', plot = DimPlot(combinedh_afterHMCluster.2y, reduction = "umap") + scale_color_manual(values = brewer.pal(6, "Paired")) + NoLegend(), height = 5, width = 5, dpi = 300)


sub = table(combinedh_afterHMCluster.2y@meta.data$split, combinedh_afterHMCluster.2y@meta.data$predicted_celltype)
sub
prop = prop.table(sub, margin = 1)
prop
cell_count = as.data.frame(table(combinedh_afterHMCluster.2y@meta.data$split, combinedh_afterHMCluster.2y@meta.data$predicted_celltype))
cell_count = as.data.frame(prop)
head(cell_count)
names(cell_count) = c("Group", "CellType", "Percentage")

Per1 = ggplot(cell_count, aes(x=Group, y=Percentage, fill=CellType)) +
  geom_col(width = 0.6) +
  scale_fill_brewer(palette = "Paired") +
  xlab(NULL) + ylab(NULL) +
  theme_minimal() +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
ggsave('./Result_for_2Y/2Y_CelltypePercentage_nolegend.pdf', plot = Per1, height = 5, width = 5, dpi = 300)

Per2 = ggplot(cell_count, aes(x=Group, y=Percentage, fill=CellType)) +
  geom_col() +
  scale_fill_brewer(palette = "Paired") +
  ggtitle("Cell type distribution in two groups") +
  xlab("Group") + ylab("Cell Percentage") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        legend.title = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
ggsave('./Result_for_2Y/2Y_CelltypePercentage.pdf', plot = Per2, height = 5, width = 5, dpi = 300)


#my_col = brewer.pal(9, "Oranges")[3:6]
Per3 = ggplot(cell_count, aes(x=CellType, y=Group, fill=Percentage)) +
  geom_tile() +
  #scale_fill_gradientn(colours = my_col) +
  ggtitle("Percentage of Cells by Type") +
  xlab("Cell Type") + ylab("Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        legend.position = "right", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
Per3
ggsave('./Result_for_2Y/2Y_CelltypePercentage_Compare.pdf', plot = Per3, height = 5, width = 5, dpi = 300)

saveRDS(combinedh_afterHMCluster.2y,file = './03.Annotation.Result/RDS/combinedh_afterAnnotation.rds')
