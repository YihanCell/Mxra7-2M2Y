library(Seurat)
library(SingleR)
library(ggplot2)
Imm.se = ImmGenData()
setwd('/data/yihan/Mxra7_2m/')

cell_types = c("B cells",
               "Erythrocytes",
               "Granulocytes",
               "NK cells",
               "Monocytes",
               "T cells")

combinedh.X = readRDS("./Result_for_X/03.Annotation.Result/RDS/X_combinedh_afterAnnotation.RDS")

SingleR_Subset <- function(cell_type){

    print(paste0("Subsetting for ", cell_type, "..."))
    combined.X.subset =  subset(combinedh.X, predicted_celltype == cell_type)
    
    print(paste0("RunPCA for ", cell_type, "..."))
    combined.X.subset = RunPCA(combined.X.subset, verbose = FALSE, npcs=30)
    
    print(paste0("RunUMAP for ", cell_type, "..."))
    combined.X.subset = RunUMAP(combined.X.subset, dims = 1:20,verbose=FALSE)
    
    print(paste0("FindNeighbor for ", cell_type, "..."))
    combined.X.subset = FindNeighbors(combined.X.subset, reduction = "pca", dims = 1:20) 
    
    print(paste0("FindCluster for ", cell_type, "..."))
    combined.X.subset = FindClusters(combined.X.subset)
    
    saveRDS(combined.X.subset,file = (paste0("./Result_for_X/04.Cell_Subtype.Result/RDS/",gsub(" ", "_", cell_type),"_subset.rds")))
    
    print(paste0("DimPlot for", cell_type, "..."))
    DimPlot(combined.X.subset, reduction = "umap", label=T,split.by = "split") 
    ggsave(paste0("./Result_for_X/04.Cell_Subtype.Result/Plot/",gsub(" ", "_", cell_type),"UMAP_split.pdf"), plot = DimPlot(combined.X.subset, reduction = "umap", label=T, split.by = "split") , height = 5, width = 5)
    DimPlot(combined.X.subset, reduction = "umap", label=T) 
    ggsave(paste0("./Result_for_X/04.Cell_Subtype.Result/Plot/",gsub(" ", "_", cell_type),"UMAP_nosplit.pdf"), plot = DimPlot(combined.X.subset, reduction = "umap", label=T) , height = 5, width = 5)
    
    print(paste0("singleR for", cell_type, "..."))
    meta.X = combined.X.subset@meta.data
    combined.X.subset_for_SingleR = GetAssayData(combined.X.subset, slot="data")
    combined.X.subset.hesc.mm = SingleR(test = combined.X.subset_for_SingleR, ref = Imm.se, labels = Imm.se$label.fine)
    
    print(paste0("Write singleR table for", cell_type, "..."))
    table1 = table(combined.X.subset.hesc.mm$labels, meta.X$seurat_clusters)
    write.csv(x = table1, file = paste0("./Result_for_X/04.Cell_Subtype.Result/RDS/",gsub(" ", "_", cell_type),"SingleR.csv"))}

for (cell_type in cell_types){
  SingleR_Subset(cell_type)
}
