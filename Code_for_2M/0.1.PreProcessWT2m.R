library(Seurat)
library(ggplot2)

# Read data for WT2m
setwd('/data/yihan/Mxra7_2m/')
WT2m = Read10X('./WT2m')
WT2mobj = CreateSeuratObject(WT2m, project = 'WT2m', min.cells = 50, min.features = 500)
dim(WT2mobj)
hist(colSums(WT2mobj),
     breaks = 150, main = "UMI count per cell",
     xlab = "UMI count per cell")
rm(WT2m)
gc()

# Cell QC
WT2mobj[["percent.mt"]] =PercentageFeatureSet(WT2mobj, pattern = "^Mt.")
head(WT2mobj@meta.data, 5)
VlnPlot(WT2mobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)
ggsave('./01.Seurat.Result/Plot/WT2m_beforeQC.pdf', plot = VlnPlot(WT2mobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2), height = 5, width = 5)
WT2m = subset(WT2mobj, subset = nFeature_RNA < 4000 & percent.mt < 0.5) 
VlnPlot(WT2m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)
ggsave('./01.Seurat.Result/Plot/WT2m_afterQC.pdf', plot = VlnPlot(WT2m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2), height = 5, width = 5)
saveRDS(object = WT2m,file = './01.Seurat.Result/RDS/WT2m_afterQC.rds')
saveRDS(object = WT2mobj,file = './01.Seurat.Result/RDS/WT2m_beforeQC.rds')

# Normalization 
WT2m = NormalizeData(WT2m, normalization.method = "LogNormalize", scale.factor = 10000)
WT2m[["RNA"]]@data[1:4, 1:4]
WT2m[["RNA"]]@counts[1:4,1:4]
hist(colSums(WT2m$RNA@data),
     breaks = 150,
     main = "UMI count per cell after normalization",
     xlab = "UMI count per cell") 

# Feature Selection
WT2m = FindVariableFeatures(WT2m, selection.method = "vst", nfeatures = 2000,verbose=FALSE)
top10 = head(VariableFeatures(WT2m), 10)
VariableFeaturePlot(WT2m)+ theme(legend.position = "none")
LabelPoints(plot = VariableFeaturePlot(WT2m) + theme(legend.position = "none"), points = top10, repel = F)
ggsave('./01.Seurat.Result/Plot/WT2m_Top10Feature.pdf', plot = LabelPoints(plot = VariableFeaturePlot(WT2m), points = top10, repel = F), height = 5, width = 5)

# Scaling the Data
WT2m = ScaleData(WT2m)
hist(colSums(WT2m$RNA@scale.data),
     breaks = 150,
     main = "UMI count per cell after scale",
     xlab = "UMI count per cell")

# PCA
hvg.WT2m = VariableFeatures(object = WT2m)
WT2m = RunPCA(WT2m,features = hvg.WT2m)
ElbowPlot(WT2m, ndims= 40)
ggsave('./01.Seurat.Result/Plot/WT2m_Elbowplot.pdf', plot = ElbowPlot(WT2m, ndims= 40), height = 5, width = 5)
DimHeatmap(WT2m, dims = 1:20, cells = 500, balanced = TRUE)

# UMAP
WT2m = RunUMAP(WT2m, dims = 1:20, n.neighbors = 30L, min.dist = 0.3,verbose=FALSE)
DimPlot(WT2m, reduction = "umap") + NoLegend()
ggsave('./01.Seurat.Result/Plot/WT2m_UMAP.pdf', plot = DimPlot(WT2m, reduction = "umap")+NoLegend(), height = 5, width = 5)
saveRDS(object = WT2m, file = "./01.Seurat.Result/RDS/WT2m_afterUMAP.rds")
