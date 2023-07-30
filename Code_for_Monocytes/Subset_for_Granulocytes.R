library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(SingleR)
Imm.se = ImmGenData()

combined.2m = readRDS("./03.Annotation.Result/RDS/2M_combinedh_afterAnnotation.RDS")
combined.2m.Granulocytes =  subset(combined.2m, predicted_celltype == "Granulocytes")
combined.2m.Granulocytes_PCA = RunPCA(combined.2m.Granulocytes, verbose = FALSE, npcs=20)
dev.new()
ElbowPlot(combined.2m.Granulocytes_PCA)
dev.off()
combined.2m.Granulocytes = RunUMAP(combined.2m.Granulocytes_PCA, dims = 1:10,verbose=FALSE)
combined.2m.Granulocytes = FindNeighbors(combined.2m.Granulocytes, reduction = "pca", dims = 1:10) 
combined.2m.Granulocytes = FindClusters(combined.2m.Granulocytes)
dev.new()
DimPlot(combined.2m.Granulocytes, reduction = "umap",label = T)
dev.off()

combined.2m_for_SingleR_Granulocytes = GetAssayData(combined.2m.Granulocytes, slot="data")
Granulocytes.singleR = SingleR(test = combined.2m_for_SingleR_Granulocytes, ref = Imm.se, labels = Imm.se$label.fine)
table(Granulocytes.singleR$labels, combined.2m.Granulocytes$seurat_clusters)
labels = as.data.frame(Granulocytes.singleR$labels, row.names=Granulocytes.singleR@rownames)
colnames(labels)[1] = "labels"
scores = as.data.frame(Granulocytes.singleR$scores, row.names=Granulocytes.singleR@rownames)

new.cluster.ids_Gcell_2m = c("Neutrophils","Neutrophils","Neutrophils","Neutrophils","Neutrophils",
                             "Neutrophils","Neutrophils", "Neutrophils","Neutrophils","Neutrophils",
                             "Neutrophils","Neutrophils","Neutrophils","Neutrophils","Neutrophils",
                             "Basophils","Neutrophils","Neutrophils","Neutrophils")
names(new.cluster.ids_Gcell_2m) = levels(combined.2m.Granulocytes)
combined.2m.Granulocytes = RenameIdents(combined.2m.Granulocytes, new.cluster.ids_Gcell_2m)
combined.2m.Granulocytes@meta.data$predicted_celltype_bcell = combined.2m.Granulocytes@active.ident

sub.2m = table(combined.2m.Granulocytes@meta.data$split, combined.2m.Granulocytes@meta.data$predicted_celltype_bcell)
sub.2m
prop.2m = prop.table(sub.2m, margin = 1)
prop.2m
cell_count.2m = as.data.frame(table(combined.2m@meta.data$split, combined.2m@meta.data$predicted_celltype))
cell_count.2m = as.data.frame(prop.2m)
head(cell_count.2m)
names(cell_count.2m) = c("Group", "CellType", "Percentage")

#====2y=====

combined.2y = readRDS("/data/yihan/Mxra7_2Y/03.Annotation.Result/RDS/combinedh_afterAnnotation.rds")
combined.2y.Granulocytes =  subset(combined.2y, predicted_celltype == "Granulocytes")
combined.2y.Granulocytes_PCA = RunPCA(combined.2y.Granulocytes, verbose = FALSE, npcs=20)
dev.new()
ElbowPlot(combined.2y.Granulocytes_PCA)
dev.off()
combined.2y.Granulocytes = RunUMAP(combined.2y.Granulocytes_PCA, dims = 1:5,verbose=FALSE)
combined.2y.Granulocytes = FindNeighbors(combined.2y.Granulocytes, reduction = "pca", dims = 1:5) 
combined.2y.Granulocytes = FindClusters(combined.2y.Granulocytes)
dev.new()
DimPlot(combined.2y.Granulocytes, reduction = "umap",label = T)
dev.off()

combined.2y_for_SingleR_Gcell = GetAssayData(combined.2y.Granulocytes, slot="data")
Granulocytes_2y.singleR = SingleR(test = combined.2y_for_SingleR_Gcell, ref = Imm.se, labels = Imm.se$label.fine)
table(Granulocytes_2y.singleR$labels, combined.2y.Granulocytes$seurat_clusters)
new.cluster.ids_Gcell_2y = c("Neutrophils","Neutrophils","Neutrophils","Neutrophils","Neutrophils",
                             "Neutrophils","Neutrophils", "Neutrophils","Neutrophils","Neutrophils",
                             "Neutrophils","Neutrophils","Neutrophils","Neutrophils","Neutrophils",
                             "Basophils","Neutrophils","Neutrophils")
names(new.cluster.ids_Gcell_2y) = levels(combined.2y.Granulocytes)
combined.2y.Granulocytes = RenameIdents(combined.2y.Granulocytes, new.cluster.ids_Gcell_2y)
combined.2y.Granulocytes@meta.data$predicted_celltype_bcell = combined.2y.Granulocytes@active.ident

sub.2y = table(combined.2y.Granulocytes@meta.data$split, combined.2y.Granulocytes@meta.data$predicted_celltype_bcell)
sub.2y
prop.2y = prop.table(sub.2y, margin = 1)
prop.2y
cell_count.2y = as.data.frame(table(combined.2y@meta.data$split, combined.2y@meta.data$predicted_celltype))
cell_count.2y = as.data.frame(prop.2y)
head(cell_count.2y)
names(cell_count.2y) = c("Group", "CellType", "Percentage")