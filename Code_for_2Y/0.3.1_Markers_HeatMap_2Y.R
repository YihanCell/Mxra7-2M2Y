library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(MAST)

setwd('/data/yihan/Mxra7_2m')

All.Markers = FindAllMarkers(combinedh_afterHMCluster.2y, min.pct = 0.5, logfc.threshold = 0.5)
head(All.Markers)
write.csv(All.Markers, file = "./Result_for_2Y/significant.markers.csv")

All.Markers %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC) -> top7
DoHeatmap(combinedh_afterHMCluster.2y, features = top4$gene) + NoLegend()
ggsave('./Result_for_2Y/MarkersHeatmap.pdf', plot = DoHeatmap(combinedh_afterHMCluster.2y, features = top10$gene) + NoLegend(), height = 8, width = 8, dpi = 300)

plot17=DoHeatmap(combinedh_afterHMCluster.2y, features = top7$gene, group.colors = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99")) + NoLegend() +
  scale_fill_gradientn(colors = c("white","#F0F0F0","#FFD92F"), guide = FALSE) + theme(axis.text.y = element_text(size = 13))
ggsave('./Result_for_2Y/MarkersHeatmapNew7cell.pdf', plot = plot17, height = 8, width = 8, dpi = 300)

saveRDS(All.Markers,file = './Result_for_2Y/Markers.rds')

