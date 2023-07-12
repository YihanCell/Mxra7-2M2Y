library(SingleR)
library(RColorBrewer)
library(ggplot2)
library(Seurat)
Mouse.se = MouseRNAseqData()
Imm.se = celldex::ImmGenData()

setwd('/data/yihan/Mxra7_2m/')
combinedh.X = readRDS('/data/yihan/Mxra7_2m/Result_for_X/02.Merge.Result/RDS/X_combinedh_afterHMCluster.rds')
meta.X = combinedh.X@meta.data
head(meta.X)
combinedh.X_for_SingleR = GetAssayData(combinedh.X, slot="data")

combined.X.hesc.mm = SingleR(test = combinedh.X_for_SingleR, ref = Mouse.se, labels = Mouse.se$label.fine)
combined.X.hesc.mm
saveRDS(object = combined.X.hesc.mm, file = "./Result_for_X/04.Cell_Subtype.Result/SingleR_result/RDS/combined.X.hesc.mm.RDS")
table(combined.X.hesc.mm$labels, meta.X$seurat_clusters)
combinedh.X@meta.data$labels = combined.X.hesc.mm$labels
saveRDS(object = combinedh.X, file = "./Result_for_X/04.Cell_Subtype.Result/SingleR_result/RDS/combinedh.X.aftersingleR.RDS")
DimPlot(combinedh.X, group.by = "labels",reduction = "umap")
ggsave('./Result_for_X/03.Annotation.Result/Plot/X_Annotation.pdf', plot = DimPlot(combinedh.X, group.by = "labels",reduction = "umap") , height = 5, width = 7)
DimPlot(combinedh.X, group.by = "labels",reduction = "umap")  + ggtitle("")
DimPlot(combinedh.X, group.by = "labels",reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(12, "Paired")) + NoLegend() + ggtitle("")
new.cluster.ids <- c("Granulocytes","Granulocytes","Granulocytes","Granulocytes","Monocytes",
                     "B cells","Granulocytes","Monocytes","Granulocytes","Monocytes",
                     "T cells","B cells","B cells","Monocytes","Erythrocytes",
                     "Monocytes","Granulocytes","Granulocytes","Granulocytes","Granulocytes",
                     "NK cells","Granulocytes","Monocytes","Granulocytes","Granulocytes",
                     "Monocytes","B cells","Granulocytes","B cells","Granulocytes",
                     "Granulocytes","Granulocytes")
names(new.cluster.ids) <- levels(combinedh.X)

combinedh.X <- RenameIdents(combinedh.X, new.cluster.ids)
combinedh.X@meta.data$predicted_celltype <- combinedh.X@active.ident
DimPlot(combinedh.X, reduction = "umap") + scale_color_manual(values = brewer.pal(6, "Paired"))
ggsave('./Result_for_X/03.Annotation.Result/Plot/X_Annotation.pdf', plot = DimPlot(combinedh.X, reduction = "umap") + scale_color_manual(values = brewer.pal(6, "Paired")), height = 5, width = 7)
DimPlot(combinedh.X, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Set2")) + NoLegend()
ggsave('./Result_for_X/03.Annotation.Result/Plot/X_Annotation_split_nolegend.pdf', plot = DimPlot(combinedh.X, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Paired")) + NoLegend(), height = 5, width = 5)
DimPlot(combinedh.X, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Set2"))
ggsave('./Result_for_X/03.Annotation.Result/Plot/X_Annotation_split.pdf', plot = DimPlot(combinedh.X, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Paired")) , height = 5, width = 7)


saveRDS(combinedh.X,file = './Result_for_X/03.Annotation.Result/RDS/2M_combinedh_afterAnnotation.RDS')

