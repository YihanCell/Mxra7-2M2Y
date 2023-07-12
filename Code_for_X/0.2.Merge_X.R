library(harmony)
library(Seurat)
library(ggplot2)


X2m = readRDS(file = '/data/yihan/Mxra7_2m/01.Seurat.Result/RDS/X2m_afterUMAP.rds')
X2y = readRDS(file = '/data/yihan/Mxra7_2m/01.Seurat.Result/RDS/X2y_afterUMAP.rds')
setwd('/data/yihan/Mxra7_2m/')
combined.X = merge(x = X2m, y = c(X2y))
Group = c(rep("X2m",dim(X2m)[2]),rep("X2y",dim(X2y)[2]))
combined.X = AddMetaData(combined.X, Group, col.name = "split")

# merge
combined.X = ScaleData(combined.X, verbose = TRUE)
combined.X = FindVariableFeatures(combined.X, selection.method = "vst", nfeatures = 3000,verbose=FALSE)
combined.X = RunPCA(combined.X, verbose = FALSE, npcs=30)
ElbowPlot(combined.X, ndims= 30)
combined.X = RunUMAP(combined.X,reduction = "pca", dims = 1:30,verbose=FALSE)
DimPlot(combined.X, reduction = "umap", group.by="split")
ggsave('./Result_for_X/02.Merge.Result/Plot/X_Merge_beforeHM.pdf', plot = DimPlot(combined.X, reduction = "umap", group.by="split"), height = 5, width = 7)
saveRDS(combined.X,"./Result_for_X/02.Merge.Result/RDS/X_combined_beforeHM.rds")

# harmony
combinedh.X = RunHarmony(combined.X,dims = 1:20, "split",verbose = FALSE)
combinedh.X = RunUMAP(combinedh.X,reduction = "harmony", dims = 1:20)
DimPlot(combinedh.X, reduction = "harmony",group.by="split")
DimPlot(combinedh.X, reduction = "umap",group.by="split")
ggsave('./Result_for_X/02.Merge.Result/Plot/X_Merge_afterHM.pdf', plot = DimPlot(combined.X, reduction = "umap", group.by="split"), height = 5, width = 7)

# cluster
combinedh.X = FindNeighbors(combinedh.X, reduction = "harmony", dims = 1:20) %>% FindClusters()
DimPlot(combinedh.X, reduction = "umap", label=T, split.by = "split") 
ggsave('./Result_for_X/02.Merge.Result/Plot/X_Merge_afterHMCluster_split.pdf', plot = DimPlot(combinedh.X, reduction = "umap", label=T, split.by = "split") , height = 5, width = 10)
DimPlot(combinedh.X, reduction = "umap", label=T) 
ggsave('./Result_for_X/02.Merge.Result/Plot/X_Merge_afterHMCluster_nosplit.pdf', plot = DimPlot(combinedh.X, reduction = "umap", label=T) , height = 5, width = 8)
saveRDS(combinedh.X,"./Result_for_X/02.Merge.Result/RDS/X_combinedh_afterHMCluster.rds")
