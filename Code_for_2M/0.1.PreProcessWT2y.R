library(Seurat)
library(ggplot2)

# Read data for WT2y
setwd('/data/yihan/')
WT2y= Read10X('/data/yihan/S6=2W2Y骨髓单细胞测序/BHSC190212-1-1 补测后/0_Cellranger/WT2Y')
WT2yobj = CreateSeuratObject(WT2y, project = 'WT2y', min.cells = 50, min.features = 500)
dim(WT2yobj)
hist(colSums(WT2yobj),
     breaks = 150, main = "UMI count per cell",
     xlab = "UMI count per cell")
rm(WT2y)
gc()

# Cell QC
WT2yobj[["percent.mt"]] =PercentageFeatureSet(WT2yobj, pattern = "^Mt.")
head(WT2yobj@meta.data, 5)
VlnPlot(WT2yobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)
setwd('/data/yihan/Mxra7_2m/')
ggsave('./01.Seurat.Result/Plot/WT2y_beforeQC.pdf', plot = VlnPlot(WT2yobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2), height = 5, width = 5)
WT2y= subset(WT2yobj, subset = nFeature_RNA < 4000 & percent.mt < 0.5) 
VlnPlot(WT2y, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)
ggsave('./01.Seurat.Result/Plot/WT2y_afterQC.pdf', plot = VlnPlot(WT2y, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2), height = 5, width = 5)
saveRDS(object = WT2y,file = './01.Seurat.Result/RDS/WT2y_afterQC.rds')
saveRDS(object = WT2yobj,file = './01.Seurat.Result/RDS/WT2y_beforeQC.rds')

# Normalization 
WT2y= NormalizeData(WT2y, normalization.method = "LogNormalize", scale.factor = 10000)
WT2y[["RNA"]]@data[1:4, 1:4]
WT2y[["RNA"]]@counts[1:4,1:4]
hist(colSums(WT2y$RNA@data),
     breaks = 150,
     main = "UMI count per cell after normalization",
     xlab = "UMI count per cell") 

# Feature Selection
WT2y= FindVariableFeatures(WT2y, selection.method = "vst", nfeatures = 2000,verbose=FALSE)
top10 = head(VariableFeatures(WT2y), 10)
VariableFeaturePlot(WT2y)+ theme(legend.position = "none")
LabelPoints(plot = VariableFeaturePlot(WT2y) + theme(legend.position = "none"), points = top10, repel = F)
ggsave('./01.Seurat.Result/Plot/WT2y_Top10Feature.pdf', plot = LabelPoints(plot = VariableFeaturePlot(WT2y), points = top10, repel = F), height = 5, width = 5)

# Scaling the Data
WT2y= ScaleData(WT2y)
hist(colSums(WT2y$RNA@scale.data),
     breaks = 150,
     main = "UMI count per cell after scale",
     xlab = "UMI count per cell")

# PCA
hvg.WT2y= VariableFeatures(object = WT2y)
WT2y= RunPCA(WT2y,features = hvg.WT2y)
ElbowPlot(WT2y, ndims= 40)
ggsave('./01.Seurat.Result/Plot/WT2y_Elbowplot.pdf', plot = ElbowPlot(WT2y, ndims= 40), height = 5, width = 5)
DimHeatmap(WT2y, dims = 1:20, cells = 500, balanced = TRUE)

# UMAP
WT2y= RunUMAP(WT2y, dims = 1:20, n.neighbors = 30L, min.dist = 0.3,verbose=FALSE)
DimPlot(WT2y, reduction = "umap") + NoLegend()
ggsave('./01.Seurat.Result/Plot/WT2y_UMAP.pdf', plot = DimPlot(WT2y, reduction = "umap")+NoLegend(), height = 5, width = 5)
saveRDS(object = WT2y, file = "./01.Seurat.Result/RDS/WT2y_afterUMAP.rds")
