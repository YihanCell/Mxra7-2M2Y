library(Seurat)
library(SingleR)
library(ggplot2)
Imm.se = ImmGenData()
setwd('/data/yihan/Mxra7_2m/')

cell_types = c("B cells",
               "Erythrocytes",
               "Granulocytes",
               "Monocytes",
               "T cells")

combinedh.2m = readRDS("./03.Annotation.Result/RDS/2M_combinedh_afterAnnotation.RDS")

SingleR_Subset <- function(cell_type){

    print(paste0("Subsetting for ", cell_type, "..."))
    combined.2m.subset =  subset(combinedh.2m, predicted_celltype == cell_type)
    
    print(paste0("RunPCA for ", cell_type, "..."))
    combined.2m.subset = RunPCA(combined.2m.subset, verbose = FALSE, npcs=30)
    
    print(paste0("RunUMAP for ", cell_type, "..."))
    combined.2m.subset = RunUMAP(combined.2m.subset, dims = 1:20,verbose=FALSE)
    
    print(paste0("FindNeighbor for ", cell_type, "..."))
    combined.2m.subset = FindNeighbors(combined.2m.subset, reduction = "pca", dims = 1:20) 
    
    print(paste0("FindCluster for ", cell_type, "..."))
    combined.2m.subset = FindClusters(combined.2m.subset)
    
    saveRDS(combined.2m.subset,file = (paste0("./04.Cell_Subtype.Result/RDS/",gsub(" ", "_", cell_type),"_subset.rds")))
    
    print(paste0("DimPlot for", cell_type, "..."))
    DimPlot(combined.2m.subset, reduction = "umap", label=T,split.by = "split") 
    ggsave(paste0("./04.Cell_Subtype.Result/Plot/",gsub(" ", "_", cell_type),"UMAP_split.pdf"), plot = DimPlot(combined.2m.subset, reduction = "umap", label=T, split.by = "split") , height = 5, width = 5)
    DimPlot(combined.2m.subset, reduction = "umap", label=T) 
    ggsave(paste0("./04.Cell_Subtype.Result/Plot/",gsub(" ", "_", cell_type),"UMAP_nosplit.pdf"), plot = DimPlot(combined.2m.subset, reduction = "umap", label=T) , height = 5, width = 5)
    
    print(paste0("singleR for", cell_type, "..."))
    meta.2m = combined.2m.subset@meta.data
    combined.2m.subset_for_SingleR = GetAssayData(combined.2m.subset, slot="data")
    combined.2m.subset.hesc.mm = SingleR(test = combined.2m.subset_for_SingleR, ref = Imm.se, labels = Imm.se$label.fine)
    
    print(paste0("Write singleR table for", cell_type, "..."))
    table1 = table(combined.2m.subset.hesc.mm$labels, meta.2m$seurat_clusters)
    write.csv(x = table1, file = paste0("./04.Cell_Subtype.Result/RDS/",gsub(" ", "_", cell_type),"SingleR.csv"))}

for (cell_type in cell_types){
  SingleR_Subset(cell_type)
}
