library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(SingleR)
Imm.se = ImmGenData()

combined.2m.bcells =  subset(combined.2m, predicted_celltype == "B cells")
combined.2m.bcells_PCA = RunPCA(combined.2m.bcells, verbose = FALSE, npcs=20)
dev.new()
ElbowPlot(combined.2m.bcells_PCA)
dev.off()
combined.2m.bcells = RunUMAP(combined.2m.bcells_PCA, dims = 1:3,verbose=FALSE)
combined.2m.bcells = FindNeighbors(combined.2m.bcells, reduction = "pca", dims = 1:3) 
combined.2m.bcells = FindClusters(combined.2m.bcells)
dev.new()
DimPlot(combined.2m.bcells, reduction = "umap")
dev.off()

FindAllMarkers(object = combined.2m.bcells)

combined.2m_for_SingleR_Bcell = GetAssayData(combined.2m.bcells, slot="data")
bcells.singleR = SingleR(test = combined.2m_for_SingleR_Bcell, ref = Imm.se, labels = Imm.se$label.fine)
table(bcells.singleR$labels, combined.2m.bcells$seurat_clusters)

# 将表格拆分为 WT 类的子表
table_WT <- subset(data.frame(bcells.singleR$labels, combined.2m.bcells$seurat_clusters),
                  combined.2m.bcells$split == "WT2m")

# 将表格拆分为 X 类的子表
table_X <- subset(data.frame(bcells.singleR$labels, combined.2m.bcells$seurat_clusters),
                  combined.2m.bcells$split == "X2m")
