library(Seurat)
library(RColorBrewer)
library(ggplot2)
setwd('/data/yihan/Mxra7_2m/')
combinedh.WT = readRDS('/data/yihan/Mxra7_2m/Result_for_WT/02.Merge.Result/RDS/WT_combinedh_afterHMCluster.rds')
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