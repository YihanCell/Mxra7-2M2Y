library(SingleR)
library(RColorBrewer)
library(ggplot2)
Mouse.se = MouseRNAseqData()

setwd('/data/yihan/Mxra7_2m/')
combinedh.2m = readRDS("./02.Merge.Result/RDS/2M_combinedh_afterHMCluster.rds")
meta.2y = combinedh.2m@meta.data
head(meta.2y)

combinedh.2m_for_SingleR = GetAssayData(combinedh.2m, slot="data")

combined.2y.hesc.mm = SingleR(test = combinedh.2m_for_SingleR, ref = Mouse.se, labels = Mouse.se$label.main)
combined.2y.hesc.mm
table(combined.2y.hesc.mm$labels, meta.2y$seurat_clusters)
combinedh.2m@meta.data$labels = combined.2y.hesc.mm$labels
DimPlot(combinedh.2m, group.by = "labels",reduction = "umap") + scale_color_manual(values = brewer.pal(19, "Paired")) + ggtitle("")
DimPlot(combinedh.2m, group.by = "labels",reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(12, "Paired")) + NoLegend() + ggtitle("")
ggsave("./PLOT/UMAP2Y.svg", plot = UMAP, width = 10, height = 10, units = "in", dpi = 600)
new.cluster.ids <- c("Granulocytes","Granulocytes","Granulocytes","T cells","Granulocytes",
                     "Granulocytes","B cells", "Monocytes","Monocytes","B cells",
                     "B cells","Granulocytes","B cells","Granulocytes","Monocytes",
                     "Monocytes","Granulocytes","Erythrocytes","Granulocytes","Granulocytes",
                     "Granulocytes","Granulocytes","Monocytes","B cells","Monocytes",
                     "B cells","B cells","Granulocytes","Monocytes","Granulocytes",
                     "Monocytes","Granulocytes")
names(new.cluster.ids) <- levels(combinedh.2m)
combinedh.2m <- RenameIdents(combinedh.2m, new.cluster.ids)
DimPlot(combinedh.2m, reduction = "umap") + scale_color_manual(values = brewer.pal(5, "Paired"))
ggsave('./03.Annotation.Result/Plot/2M_Annotation.pdf', plot = DimPlot(combinedh.2m, reduction = "umap") + scale_color_manual(values = brewer.pal(5, "Paired")), height = 5, width = 7)
DimPlot(combinedh.2m, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(5, "Set2")) + NoLegend()
ggsave('./03.Annotation.Result/Plot/2M_Annotation_split_nolegend.pdf', plot = DimPlot(combinedh.2m, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(5, "Paired")) + NoLegend(), height = 5, width = 5)
DimPlot(combinedh.2m, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(5, "Set2"))
ggsave('./03.Annotation.Result/Plot/2M_Annotation_split.pdf', plot = DimPlot(combinedh.2m, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Paired")) , height = 5, width = 7)
combinedh.2m@meta.data$predicted_celltype <- combinedh.2m@active.ident

sub = table(combinedh.2m@meta.data$split, combinedh.2m@meta.data$predicted_celltype)
sub
prop <- prop.table(sub, margin = 1)
prop
cell_count = as.data.frame(table(combinedh.2m@meta.data$split, combinedh.2m@meta.data$predicted_celltype))
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
Per1
ggsave('./03.Annotation.Result/Plot/2M_CelltypePercentage_nolegend.pdf', plot = Per1, height = 5, width = 5)

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
Per2
ggsave('./03.Annotation.Result/Plot/2M_CelltypePercentage.pdf', plot = Per2, height = 5, width = 8)


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
ggsave('./03.Annotation.Result/Plot/2M_CelltypePercentage_Compare.pdf', plot = Per3, height = 5, width = 6)

saveRDS(combinedh.2m,file = './03.Annotation.Result/RDS/2M_combinedh_afterAnnotation.RDS')
