library(Seurat)
library(ggplot2)

# Read data for X2y
setwd('/data/yihan/Mxra7_2m/')
X2y = Read10X('/data/yihan/S6=2W2Y骨髓单细胞测序/BHSC190212-1-1 补测后/0_Cellranger/X2Y-homo')
X2yobj = CreateSeuratObject(X2y, project = 'X2y', min.cells = 50, min.features = 500)
dim(X2yobj)
hist(colSums(X2yobj),
     breaks = 150, main = "UMI count per cell",
     xlab = "UMI count per cell")
rm(X2y)
gc()

# Cell QC
X2yobj[["percent.mt"]] =PercentageFeatureSet(X2yobj, pattern = "^Mt.")
head(X2yobj@meta.data, 5)
VlnPlot(X2yobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)
ggsave('./01.Seurat.Result/Plot/X2y_beforeQC.pdf', plot = VlnPlot(X2yobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2), height = 5, width = 5)
X2y = subset(X2yobj, subset = nFeature_RNA < 5000 & percent.mt < 0.5) 
VlnPlot(X2y, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)
ggsave('./01.Seurat.Result/Plot/X2y_afterQC.pdf', plot = VlnPlot(X2y, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2), height = 5, width = 5)
saveRDS(object = X2y,file = './01.Seurat.Result/RDS/X2y_afterQC.rds')
saveRDS(object = X2yobj,file = './01.Seurat.Result/RDS/X2y_beforeQC.rds')

# Normalization 
X2y = NormalizeData(X2y, normalization.method = "LogNormalize", scale.factor = 10000)
X2y[["RNA"]]@data[1:4, 1:4]
X2y[["RNA"]]@counts[1:4,1:4]
hist(colSums(X2y$RNA@data),
     breaks = 150,
     main = "UMI count per cell after normalization",
     xlab = "UMI count per cell") 

# Feature Selection
X2y = FindVariableFeatures(X2y, selection.method = "vst", nfeatures = 2000,verbose=FALSE)
top10 = head(VariableFeatures(X2y), 10)
VariableFeaturePlot(X2y)+ theme(legend.position = "none")
LabelPoints(plot = VariableFeaturePlot(X2y) + theme(legend.position = "none"), points = top10, repel = F)
ggsave('./01.Seurat.Result/Plot/X2y_Top10Feature.pdf', plot = LabelPoints(plot = VariableFeaturePlot(X2y), points = top10, repel = F), height = 5, width = 5)

# Scaling the Data
X2y = ScaleData(X2y)
hist(colSums(X2y$RNA@scale.data),
     breaks = 150,
     main = "UMI count per cell after scale",
     xlab = "UMI count per cell")

# PCA
hvg.X2y = VariableFeatures(object = X2y)
X2y = RunPCA(X2y,features = hvg.X2y)
ElbowPlot(X2y, ndims= 40)
ggsave('./01.Seurat.Result/Plot/X2y_Elbowplot.pdf', plot = ElbowPlot(X2y, ndims= 40), height = 5, width = 5)
DimHeatmap(X2y, dims = 1:20, cells = 500, balanced = TRUE)

# UMAP
X2y = RunUMAP(X2y, dims = 1:20, n.neighbors = 30L, min.dist = 0.3,verbose=FALSE)
DimPlot(X2y, reduction = "umap") + NoLegend()
ggsave('./01.Seurat.Result/Plot/X2y_UMAP.pdf', plot = DimPlot(X2y, reduction = "umap")+NoLegend(), height = 5, width = 5)
saveRDS(object = X2y, file = "./01.Seurat.Result/RDS/X2y_afterUMAP.rds")
