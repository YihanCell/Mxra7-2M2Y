library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
setwd('/data/yihan/Mxra7/')

All.Markers = FindAllMarkers(combinedh.2y, min.pct = 0.5, logfc.threshold = 0.5)
head(All.Markers)
write.csv(All.Markers, file = "./03.Annotation.Result/RDS/significant.markers.csv")
write.csv(All.Markers, file = "./03.Annotation.Result/RDS/significant.markers.txt")

All.Markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(combinedh.2y, features = top4$gene) + NoLegend()
ggsave('./03.Annotation.Result/Plot/MarkersHeatmap.svg', plot = DoHeatmap(combinedh.2y, features = top10$gene) + NoLegend() , height = 20, width = 20, dpi = 500)
ggsave('./03.Annotation.Result/Plot/MarkersHeatmap.jpg', plot = DoHeatmap(combinedh.2y, features = top10$gene) + NoLegend(), height = 20, width = 20, dpi = 500)
ggsave('./03.Annotation.Result/Plot/MarkersHeatmap.pdf', plot = DoHeatmap(combinedh.2y, features = top10$gene) + NoLegend(), height = 20, width = 20, dpi = 500)

plot12=DoHeatmap(combinedh.2y, features = top4$gene, group.colors = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C"),angle = 90) + NoLegend() +
  scale_fill_gradientn(colors = c("white","#F0F0F0","#FFD92F"), guide = FALSE) + theme(axis.text.y = element_text(size = 18))
ggsave('./03.Annotation.Result/Plot/MarkersHeatmapNew4cell.jpg', plot = plot12, height = 20, width = 20, dpi = 500)
ggsave('./03.Annotation.Result/Plot/MarkersHeatmapNew4cell.svg', plot = plot12, height = 20, width = 20, dpi = 500)
ggsave('./03.Annotation.Result/Plot/MarkersHeatmapNew5cell.pdf', plot = plot12, height = 15, width = 15, dpi = 500)

saveRDS(All.Markers,file = './03.Annotation.Result/RDS/Markers.rds')

