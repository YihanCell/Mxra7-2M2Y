library(harmony)
library(Seurat)
library(ggplot2)

WT2m = readRDS(file = '/data/yihan/Mxra7_2m/01.Seurat.Result/RDS/WT2m_afterUMAP.rds')
X2m = readRDS(file = '/data/yihan/Mxra7_2m/01.Seurat.Result/RDS/X2m_afterUMAP.rds')
combined.2m = merge(x = WT2m, y = c(X2m))
Group = c(rep("WT2m",dim(WT2m)[2]),rep("X2m",dim(X2m)[2]))
combined.2m = AddMetaData(combined.2m, Group, col.name = "split")

# merge
combined.2m = ScaleData(combined.2m, verbose = TRUE)
combined.2m = FindVariableFeatures(combined.2m, selection.method = "vst", nfeatures = 3000,verbose=FALSE)
combined.2m = RunPCA(combined.2m, verbose = FALSE, npcs=30)
ElbowPlot(combined.2m, ndims= 30)
combined.2m = RunUMAP(combined.2m,reduction = "pca", dims = 1:30,verbose=FALSE)
DimPlot(combined.2m, reduction = "umap", group.by="split")
ggsave('./02.Merge.Result/Plot/New_2M_Merge_beforeHM.pdf', plot = DimPlot(combined.2m, reduction = "umap", group.by="split"), height = 5, width = 7)
saveRDS(combined.2m,"./02.Merge.Result/RDS/New_combined_beforeHM.2M.rds")

# harmony
combinedh.2m = RunHarmony(combined.2m,dims = 1:20, "split",verbose = FALSE)
combinedh.2m = RunUMAP(combinedh.2m,reduction = "harmony", dims = 1:20)
DimPlot(combinedh.2m, reduction = "harmony",group.by="split")
DimPlot(combinedh.2m, reduction = "umap",group.by="split")
ggsave('./02.Merge.Result/Plot/New_2M_Merge_afterHM.pdf', plot = DimPlot(combined.2m, reduction = "umap", group.by="split"), height = 5, width = 7)

# cluster
combinedh.2m = FindNeighbors(combinedh.2m, reduction = "harmony", dims = 1:20) %>% FindClusters()
DimPlot(combinedh.2m, reduction = "umap", label=T, split.by = "split") 
ggsave('./02.Merge.Result/Plot/New_2M_Merge_afterHMCluster_split.pdf', plot = DimPlot(combinedh.2m, reduction = "umap", label=T, split.by = "split") , height = 5, width = 10)
DimPlot(combinedh.2m, reduction = "umap", label=T) 
ggsave('./02.Merge.Result/Plot/New_2M_Merge_afterHMCluster_nosplit.pdf', plot = DimPlot(combinedh.2m, reduction = "umap", label=T) , height = 5, width = 8)
saveRDS(combinedh.2m,"./02.Merge.Result/RDS/New_2M_combinedh_afterHMCluster.rds")
