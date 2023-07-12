library(harmony)
library(Seurat)
library(ggplot2)


WT2m = readRDS(file = '/data/yihan/Mxra7_2m/01.Seurat.Result/RDS/WT2m_afterUMAP.rds')
WT2y = readRDS(file = '/data/yihan/Mxra7_2m/01.Seurat.Result/RDS/WT2y_afterUMAP.rds')
setwd('/data/yihan/Mxra7_2m/')
combined.WT = merge(x = WT2m, y = c(WT2y))
Group = c(rep("WT2m",dim(WT2m)[2]),rep("WT2y",dim(WT2y)[2]))
combined.WT = AddMetaData(combined.WT, Group, col.name = "split")

# merge
combined.WT = ScaleData(combined.WT, verbose = TRUE)
combined.WT = FindVariableFeatures(combined.WT, selection.method = "vst", nfeatures = 3000,verbose=FALSE)
combined.WT = RunPCA(combined.WT, verbose = FALSE, npcs=30)
ElbowPlot(combined.WT, ndims= 30)
combined.WT = RunUMAP(combined.WT,reduction = "pca", dims = 1:30,verbose=FALSE)
DimPlot(combined.WT, reduction = "umap", group.by="split")
ggsave('./Result_for_WT/02.Merge.Result/Plot/WT_Merge_beforeHM.pdf', plot = DimPlot(combined.WT, reduction = "umap", group.by="split"), height = 5, width = 7)
saveRDS(combined.WT,"./Result_for_WT/02.Merge.Result/RDS/WT_combined_beforeHM.rds")

# harmony
combinedh.WT = RunHarmony(combined.WT,dims = 1:20, "split",verbose = FALSE)
combinedh.WT = RunUMAP(combinedh.WT,reduction = "harmony", dims = 1:20)
DimPlot(combinedh.WT, reduction = "harmony",group.by="split")
DimPlot(combinedh.WT, reduction = "umap",group.by="split")
ggsave('./Result_for_WT/02.Merge.Result/Plot/WT_Merge_afterHM.pdf', plot = DimPlot(combined.WT, reduction = "umap", group.by="split"), height = 5, width = 7)

# cluster
combinedh.WT = FindNeighbors(combinedh.WT, reduction = "harmony", dims = 1:20) %>% FindClusters()
DimPlot(combinedh.WT, reduction = "umap", label=T, split.by = "split") 
ggsave('./Result_for_WT/02.Merge.Result/Plot/WT_Merge_afterHMCluster_split.pdf', plot = DimPlot(combinedh.WT, reduction = "umap", label=T, split.by = "split") , height = 5, width = 10)
DimPlot(combinedh.WT, reduction = "umap", label=T) 
ggsave('./Result_for_WT/02.Merge.Result/Plot/WT_Merge_afterHMCluster_nosplit.pdf', plot = DimPlot(combinedh.WT, reduction = "umap", label=T) , height = 5, width = 8)
saveRDS(combinedh.WT,"./Result_for_WT/02.Merge.Result/RDS/WT_combinedh_afterHMCluster.rds")
