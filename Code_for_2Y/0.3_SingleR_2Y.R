library(SingleR)
library(RColorBrewer)
library(ggplot2)
Mouse.se = MouseRNAseqData()
setwd('/data/yihan/Mxra7/')
meta.2y = combinedh.2y@meta.data
head(meta.2y)
combinedh.2y_for_SingleR = GetAssayData(combinedh.2y, slot="data")

combined.2y.hesc.mm = SingleR(test = combinedh.2y_for_SingleR, ref = Mouse.se, labels = Mouse.se$label.main)
combined.2y.hesc.mm
combined.2y.hesc.main = SingleR(test = combinedh.2y_for_SingleR, ref = Mouse.se, labels = Mouse.se$label.main)
table(combined.2y.hesc.mm$labels, meta.2y$seurat_clusters)
combinedh.2y@meta.data$labels = combined.2y.hesc.main$labels
UMAP1 = DimPlot(combinedh.2y, group.by = "labels",reduction = "umap") + scale_color_manual(values = brewer.pal(19, "Paired")) + ggtitle("")
UMAP2 = DimPlot(combined.2y, group.by = "labels",reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(12, "Paired")) + NoLegend() + ggtitle("")
UMAP = grid.arrange(UMAP1, UMAP2, ncol = 1)
ggsave("./PLOT/UMAP2Y.svg", plot = UMAP, width = 10, height = 10, units = "in", dpi = 600)

new.cluster.ids <- c("Granulocytes","Granulocytes","Granulocytes","Monocytes","Granulocytes",
                     "Monocytes","Monocytes","Granulocytes","T cells","Granulocytes",
                     "B cells","Granulocytes","B cells","Erythrocytes","Granulocytes",
                     "Monocytes","Monocytes","B cells","B cells","Granulocytes",
                     "Granulocytes","Monocytes","Granulocytes","Granulocytes","Granulocytes",
                     "NK cells","Monocytes","Granulocytes","Monocytes","Monocytes",
                     "Monocytes","B cells","Granulocytes","Monocytes","Granulocytes")
names(new.cluster.ids) <- levels(combinedh.2y)
combinedh.2y <- RenameIdents(combinedh.2y, new.cluster.ids)
DimPlot(combinedh.2y, reduction = "umap") + scale_color_manual(values = brewer.pal(6, "Paired"))
ggsave('./03.Annotation.Result/Plot/2Y_Annotation.pdf', plot = DimPlot(combinedh.2y, reduction = "umap") + scale_color_manual(values = brewer.pal(6, "Paired")), height = 10, width = 10, dpi = 500)
DimPlot(combinedh.2y, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Set2")) + NoLegend()
ggsave('./03.Annotation.Result/Plot/2Y_Annotation_split_nolegend.pdf', plot = DimPlot(combinedh.2y, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Paired")) + NoLegend(), height = 10, width = 10, dpi = 500)
DimPlot(combinedh.2y, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Set2"))
ggsave('./03.Annotation.Result/Plot/2Y_Annotation_split.pdf', plot = DimPlot(combinedh.2y, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Paired")) , height = 10, width = 10, dpi = 500)
combinedh.2y@meta.data$predicted_celltype <- combinedh.2y@active.ident

sub = table(combinedh.2y@meta.data$split, combinedh.2y@meta.data$predicted_celltype)
sub
prop <- prop.table(sub, margin = 1)
prop
cell_count = as.data.frame(table(combinedh.2y@meta.data$split, combinedh.2y@meta.data$predicted_celltype))
cell_count = as.data.frame(prop)
head(cell_count)
names(cell_count) <- c("Group", "CellType", "Percentage")

Per1 = ggplot(cell_count, aes(x=Group, y=Percentage, fill=CellType)) +
  geom_col(width = 0.6) +
  scale_fill_brewer(palette = "Paired") +
  xlab(NULL) + ylab(NULL) +
  theme_minimal() +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
ggsave('./03.Annotation.Result/Plot/2Y_CelltypePercentage_nolegend.pdf', plot = Per1, height = 10, width = 10, dpi = 500)

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
ggsave('./03.Annotation.Result/Plot/2Y_CelltypePercentage.pdf', plot = Per2, height = 10, width = 10, dpi = 500)


#my_col <- brewer.pal(9, "Oranges")[3:6]

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
ggsave('./03.Annotation.Result/Plot/2Y_CelltypePercentage_Compare.pdf', plot = Per3, height = 10, width = 10, dpi = 500)

saveRDS(combinedh.2y,file = './03.Annotation.Result/RDS/combinedh_afterAnnotation.svg')
