library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(SingleR)
Imm.se = ImmGenData()

combined.2m = readRDS("./03.Annotation.Result/RDS/2M_combinedh_afterAnnotation.RDS")
combined.2m.bcells =  subset(combined.2m, predicted_celltype == "B cells")
combined.2m.bcells_PCA = RunPCA(combined.2m.bcells, verbose = FALSE, npcs=20)
dev.new()
ElbowPlot(combined.2m.bcells_PCA)
dev.off()
combined.2m.bcells = RunUMAP(combined.2m.bcells_PCA, dims = 1:5,verbose=FALSE)
combined.2m.bcells = FindNeighbors(combined.2m.bcells, reduction = "pca", dims = 1:5) 
combined.2m.bcells = FindClusters(combined.2m.bcells)
dev.new()
DimPlot(combined.2m.bcells, reduction = "umap",label = T)
dev.off()
FindAllMarkers(object = combined.2m.bcells)

combined.2m_for_SingleR_Bcell = GetAssayData(combined.2m.bcells, slot="data")
bcells.singleR = SingleR(test = combined.2m_for_SingleR_Bcell, ref = Imm.se, labels = Imm.se$label.fine)
table(bcells.singleR$labels, combined.2m.bcells$seurat_clusters)
labels = as.data.frame(bcells.singleR$labels, row.names=bcells.singleR@rownames)
colnames(labels)[1] = "labels"
scores = as.data.frame(bcells.singleR$scores, row.names=bcells.singleR@rownames)

new.cluster.ids_Bcell_2m = c("B.FrF","preB.FrD","B.FrE","preB.FrC","proB.FrA",
                     "preB.FrC","proB.FrA", "proB.FrA","proB.FrA","proB.FrA",
                     "B.FrE","proB.FrA","preB.FrC","proB.FrA")
names(new.cluster.ids_Bcell_2m) = levels(combined.2m.bcells)
combined.2m.bcells = RenameIdents(combined.2m.bcells, new.cluster.ids_Bcell_2m)
combined.2m.bcells@meta.data$predicted_celltype_bcell = combined.2m.bcells@active.ident

sub.2m = table(combined.2m.bcells@meta.data$split, combined.2m.bcells@meta.data$predicted_celltype_bcell)
sub.2m
prop.2m = prop.table(sub.2m, margin = 1)
prop.2m
cell_count.2m = as.data.frame(table(combined.2m@meta.data$split, combined.2m@meta.data$predicted_celltype))
cell_count.2m = as.data.frame(prop.2m)
head(cell_count.2m)
names(cell_count.2m) = c("Group", "CellType", "Percentage")

#====2y=====

combined.2y = readRDS("/data/yihan/Mxra7_2Y/03.Annotation.Result/RDS/combinedh_afterAnnotation.rds")
combined.2y.bcells =  subset(combined.2y, predicted_celltype == "B cells")
combined.2y.bcells_PCA = RunPCA(combined.2y.bcells, verbose = FALSE, npcs=20)
dev.new()
ElbowPlot(combined.2y.bcells_PCA)
dev.off()
combined.2y.bcells = RunUMAP(combined.2y.bcells_PCA, dims = 1:6,verbose=FALSE)
combined.2y.bcells = FindNeighbors(combined.2y.bcells, reduction = "pca", dims = 1:6) 
combined.2y.bcells = FindClusters(combined.2y.bcells)
dev.new()
DimPlot(combined.2y.bcells, reduction = "umap",label = T)
dev.off()

combined.2y_for_SingleR_Bcell = GetAssayData(combined.2y.bcells, slot="data")
bcells_2y.singleR = SingleR(test = combined.2y_for_SingleR_Bcell, ref = Imm.se, labels = Imm.se$label.fine)
table(bcells_2y.singleR$labels, combined.2y.bcells$seurat_clusters)
new.cluster.ids_Bcell_2y = c("preB.FrD","B.FrF","B.FrE","B.FrF","preB.FrC",
                             "B.FrF","B.FrF", "B.FrF","preB.FrC","B.FrF",
                             "proB.FrA","proB.FrA","preB.FrD","proB.FrA","B.FrF")
names(new.cluster.ids_Bcell_2y) = levels(combined.2y.bcells)
combined.2y.bcells = RenameIdents(combined.2y.bcells, new.cluster.ids_Bcell_2y)
combined.2y.bcells@meta.data$predicted_celltype_bcell = combined.2y.bcells@active.ident

sub.2y = table(combined.2y.bcells@meta.data$split, combined.2y.bcells@meta.data$predicted_celltype_bcell)
sub.2y
prop.2y = prop.table(sub.2y, margin = 1)
prop.2y
cell_count.2y = as.data.frame(table(combined.2y@meta.data$split, combined.2y@meta.data$predicted_celltype))
cell_count.2y = as.data.frame(prop.2y)
head(cell_count.2y)
names(cell_count.2y) = c("Group", "CellType", "Percentage")