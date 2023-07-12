library(SingleR)
library(RColorBrewer)
library(ggplot2)
library(Seurat)
Mouse.se = MouseRNAseqData()
Imm.se = celldex::ImmGenData()

setwd('/data/yihan/Mxra7_2m/')
combinedh.WT = readRDS('/data/yihan/Mxra7_2m/Result_for_WT/02.Merge.Result/RDS/WT_combinedh_afterHMCluster.rds')
meta.WT = combinedh.WT@meta.data
head(meta.WT)
combinedh.WT_for_SingleR = GetAssayData(combinedh.WT, slot="data")

combined.WT.hesc.mm = SingleR(test = combinedh.WT_for_SingleR, ref = Mouse.se, labels = Mouse.se$label.fine)
combined.WT.hesc.mm
saveRDS(object = combined.WT.hesc.mm, file = "./Result_for_WT/04.SingleR/RDS/combined.WT.hesc.mm.RDS")
table(combined.WT.hesc.mm$labels, meta.WT$seurat_clusters)
combinedh.WT@meta.data$labels = combined.WT.hesc.mm$labels
saveRDS(object = combinedh.WT, file = "./Result_for_WT/04.SingleR/RDS/combinedh.WT.aftersingleR.RDS")
DimPlot(combinedh.WT, group.by = "labels",reduction = "umap")
ggsave('./Result_for_WT/03.Annotation.Result/Plot/WT_Annotation.pdf', plot = DimPlot(combinedh.WT, group.by = "labels",reduction = "umap") , height = 5, width = 7)
DimPlot(combinedh.WT, group.by = "labels",reduction = "umap")  + ggtitle("")
DimPlot(combinedh.WT, group.by = "labels",reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(12, "Paired")) + NoLegend() + ggtitle("")
new.cluster.ids <- c("Granulocytes","Granulocytes","B cells","Granulocytes","Granulocytes",
                     "B cells","Granulocytes", "Erythrocytes","Monocytes","B cells",
                     "B cells","T cells","Granulocytes","Monocytes","Monocytes",
                     "Granulocytes","Granulocytes","Granulocytes","Granulocytes","Monocytes",
                     "B cells","B cells","B cells","Monocytes","Granulocytes",
                     "Monocytes","Monocytes","Monocytes","Granulocytes","B cells",
                     "Granulocytes")
names(new.cluster.ids) <- levels(combinedh.WT)

combinedh.WT <- RenameIdents(combinedh.WT, new.cluster.ids)
DimPlot(combinedh.WT, reduction = "umap") + scale_color_manual(values = brewer.pal(5, "Paired"))
ggsave('./Result_for_WT/03.Annotation.Result/Plot/WT_Annotation.pdf', plot = DimPlot(combinedh.WT, reduction = "umap") + scale_color_manual(values = brewer.pal(5, "Paired")), height = 5, width = 7)
DimPlot(combinedh.WT, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(5, "Set2")) + NoLegend()
ggsave('./Result_for_WT/03.Annotation.Result/Plot/WT_Annotation_split_nolegend.pdf', plot = DimPlot(combinedh.WT, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(5, "Paired")) + NoLegend(), height = 5, width = 5)
DimPlot(combinedh.WT, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(5, "Set2"))
ggsave('./Result_for_WT/03.Annotation.Result/Plot/WT_Annotation_split.pdf', plot = DimPlot(combinedh.WT, reduction = "umap", split.by = "split") + scale_color_manual(values = brewer.pal(6, "Paired")) , height = 5, width = 7)
combinedh.WT@meta.data$predicted_celltype <- combinedh.WT@active.ident

saveRDS(combinedh.WT,file = './Result_for_WT/03.Annotation.Result/RDS/2M_combinedh_afterAnnotation.RDS')

