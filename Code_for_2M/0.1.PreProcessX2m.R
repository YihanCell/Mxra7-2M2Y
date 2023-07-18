library(Seurat)
library(ggplot2)

# Read data for X2m
setwd('/data/yihan/Mxra7_2m/')
X2m = Read10X('./X2m')
X2mobj = CreateSeuratObject(X2m, project = 'X2m', min.cells = 50, min.features = 500)
dim(X2mobj)
hist(colSums(X2mobj),
     breaks = 150, main = "UMI count per cell",
     xlab = "UMI count per cell")
rm(X2m)
gc()

# Cell QC
X2mobj[["percent.mt"]] =PercentageFeatureSet(X2mobj, pattern = "^Mt.")
head(X2mobj@meta.data, 5)
VlnPlot(X2mobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)
ggsave('./01.Seurat.Result/Plot/X2m_beforeQC.pdf', plot = VlnPlot(X2mobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2), height = 5, width = 5)
X2m = subset(X2mobj, subset = nFeature_RNA < 4000 & percent.mt < 0.5) 
VlnPlot(X2m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)
ggsave('./01.Seurat.Result/Plot/X2m_afterQC.pdf', plot = VlnPlot(X2m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2), height = 5, width = 5)
saveRDS(object = X2m,file = './01.Seurat.Result/RDS/X2m_afterQC.rds')
saveRDS(object = X2mobj,file = './01.Seurat.Result/RDS/X2m_beforeQC.rds')

# Normalization 
X2m = NormalizeData(X2m, normalization.method = "LogNormalize", scale.factor = 10000)
X2m[["RNA"]]@data[1:4, 1:4]
X2m[["RNA"]]@counts[1:4,1:4]
hist(colSums(X2m$RNA@data),
     breaks = 150,
     main = "UMI count per cell after normalization",
     xlab = "UMI count per cell") 

# Feature Selection
X2m = FindVariableFeatures(X2m, selection.method = "vst", nfeatures = 2000,verbose=FALSE)
top10 = head(VariableFeatures(X2m), 10)
VariableFeaturePlot(X2m)+ theme(legend.position = "none")
LabelPoints(plot = VariableFeaturePlot(X2m) + theme(legend.position = "none"), points = top10, repel = F)
ggsave('./01.Seurat.Result/Plot/X2m_Top10Feature.pdf', plot = LabelPoints(plot = VariableFeaturePlot(X2m), points = top10, repel = F), height = 5, width = 5)

# Scaling the Data
X2m = ScaleData(X2m)
hist(colSums(X2m$RNA@scale.data),
     breaks = 150,
     main = "UMI count per cell after scale",
     xlab = "UMI count per cell")

# PCA
hvg.X2m = VariableFeatures(object = X2m)
X2m = RunPCA(X2m,features = hvg.X2m)
ElbowPlot(X2m, ndims= 40)
ggsave('./01.Seurat.Result/Plot/X2m_Elbowplot.pdf', plot = ElbowPlot(X2m, ndims= 40), height = 5, width = 5)
DimHeatmap(X2m, dims = 1:20, cells = 500, balanced = TRUE)

# UMAP
X2m = RunUMAP(X2m, dims = 1:20, n.neighbors = 30L, min.dist = 0.3,verbose=FALSE)
DimPlot(X2m, reduction = "umap") + NoLegend()
ggsave('./01.Seurat.Result/Plot/X2m_UMAP.pdf', plot = DimPlot(X2m, reduction = "umap")+NoLegend(), height = 5, width = 5)
saveRDS(object = X2m, file = "./01.Seurat.Result/RDS/X2m_afterUMAP.rds")
