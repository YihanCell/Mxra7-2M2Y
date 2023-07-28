library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(SingleR)
Imm.se = ImmGenData()

combined.2m = readRDS("./03.Annotation.Result/RDS/2M_combinedh_afterAnnotation.RDS")
combined.2m.Monocytes =  subset(combined.2m, predicted_celltype == "Monocytes")
combined.2m.Monocytes_PCA = RunPCA(combined.2m.Monocytes, verbose = FALSE, npcs=20)
dev.new()
ElbowPlot(combined.2m.Monocytes_PCA)
dev.off()
combined.2m.Monocytes = RunUMAP(combined.2m.Monocytes_PCA, dims = 1:10,verbose=FALSE)
combined.2m.Monocytes = FindNeighbors(combined.2m.Monocytes, reduction = "pca", dims = 1:10) 
combined.2m.Monocytes = FindClusters(combined.2m.Monocytes)
dev.new()
DimPlot(combined.2m.Monocytes, reduction = "umap",label = T)
dev.off()

combined.2m_for_SingleR_Monocytes = GetAssayData(combined.2m.Monocytes, slot="data")
Monocytes.singleR = SingleR(test = combined.2m_for_SingleR_Monocytes, ref = Imm.se, labels = Imm.se$label.fine)
table(Monocytes.singleR$labels, combined.2m.Monocytes$seurat_clusters)
labels = as.data.frame(Monocytes.singleR$labels, row.names=Monocytes.singleR@rownames)
colnames(labels)[1] = "labels"
scores = as.data.frame(Monocytes.singleR$scores, row.names=Monocytes.singleR@rownames)

new.cluster.ids_Bcell_2m = c("B.FrF","preB.FrD","B.FrE","preB.FrC","proB.FrA",
                             "preB.FrC","proB.FrA", "proB.FrA","proB.FrA","proB.FrA",
                             "B.FrE","proB.FrA","preB.FrC","proB.FrA")
names(new.cluster.ids_Bcell_2m) = levels(combined.2m.Monocytes)
combined.2m.Monocytes = RenameIdents(combined.2m.Monocytes, new.cluster.ids_Bcell_2m)
combined.2m.Monocytes@meta.data$predicted_celltype_bcell = combined.2m.Monocytes@active.ident

sub.2m = table(combined.2m.Monocytes@meta.data$split, combined.2m.Monocytes@meta.data$predicted_celltype_bcell)
sub.2m
prop.2m = prop.table(sub.2m, margin = 1)
prop.2m
cell_count.2m = as.data.frame(table(combined.2m@meta.data$split, combined.2m@meta.data$predicted_celltype))
cell_count.2m = as.data.frame(prop.2m)
head(cell_count.2m)
names(cell_count.2m) = c("Group", "CellType", "Percentage")

#====2y=====

combined.2y = readRDS("/data/yihan/Mxra7_2Y/03.Annotation.Result/RDS/combinedh_afterAnnotation.rds")
combined.2y.Monocytes =  subset(combined.2y, predicted_celltype == "Monocytes")
combined.2y.Monocytes_PCA = RunPCA(combined.2y.Monocytes, verbose = FALSE, npcs=20)
dev.new()
ElbowPlot(combined.2y.Monocytes_PCA)
dev.off()
combined.2y.Monocytes = RunUMAP(combined.2y.Monocytes_PCA, dims = 1:10,verbose=FALSE)
combined.2y.Monocytes = FindNeighbors(combined.2y.Monocytes, reduction = "pca", dims = 1:10) 
combined.2y.Monocytes = FindClusters(combined.2y.Monocytes)
dev.new()
DimPlot(combined.2y.Monocytes, reduction = "umap",label = T)
dev.off()

combined.2y_for_SingleR_Bcell = GetAssayData(combined.2y.Monocytes, slot="data")
Monocytes_2y.singleR = SingleR(test = combined.2y_for_SingleR_Bcell, ref = Imm.se, labels = Imm.se$label.fine)
table(Monocytes_2y.singleR$labels, combined.2y.Monocytes$seurat_clusters)
new.cluster.ids_Bcell_2y = c("preB.FrD","B.FrF","B.FrE","B.FrF","preB.FrC",
                             "B.FrF","B.FrF", "B.FrF","preB.FrC","B.FrF",
                             "proB.FrA","proB.FrA","preB.FrD","proB.FrA","B.FrF")
names(new.cluster.ids_Bcell_2y) = levels(combined.2y.Monocytes)
combined.2y.Monocytes = RenameIdents(combined.2y.Monocytes, new.cluster.ids_Bcell_2y)
combined.2y.Monocytes@meta.data$predicted_celltype_bcell = combined.2y.Monocytes@active.ident

sub.2y = table(combined.2y.Monocytes@meta.data$split, combined.2y.Monocytes@meta.data$predicted_celltype_bcell)
sub.2y
prop.2y = prop.table(sub.2y, margin = 1)
prop.2y
cell_count.2y = as.data.frame(table(combined.2y@meta.data$split, combined.2y@meta.data$predicted_celltype))
cell_count.2y = as.data.frame(prop.2y)
head(cell_count.2y)
names(cell_count.2y) = c("Group", "CellType", "Percentage")