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

combinedh.WT = readRDS("./Result_for_WT/03.Annotation.Result/RDS/WT_combinedh_afterAnnotation.RDS")

SingleR_Subset <- function(cell_type){

    print(paste0("Subsetting for ", cell_type, "..."))
    combined.WT.subset =  subset(combinedh.WT, predicted_celltype == cell_type)
    
    print(paste0("RunPCA for ", cell_type, "..."))
    combined.WT.subset = RunPCA(combined.WT.subset, verbose = FALSE, npcs=30)
    
    print(paste0("RunUMAP for ", cell_type, "..."))
    combined.WT.subset = RunUMAP(combined.WT.subset, dims = 1:20,verbose=FALSE)
    
    print(paste0("FindNeighbor for ", cell_type, "..."))
    combined.WT.subset = FindNeighbors(combined.WT.subset, reduction = "pca", dims = 1:20) 
    
    print(paste0("FindCluster for ", cell_type, "..."))
    combined.WT.subset = FindClusters(combined.WT.subset)
    
    saveRDS(combined.WT.subset,file = (paste0("./Result_for_WT/04.Cell_Subtype.Result/RDS/",gsub(" ", "_", cell_type),"_subset.rds")))
    
    print(paste0("DimPlot for", cell_type, "..."))
    DimPlot(combined.WT.subset, reduction = "umap", label=T,split.by = "split") 
    ggsave(paste0("./Result_for_WT/04.Cell_Subtype.Result/Plot/",gsub(" ", "_", cell_type),"UMAP_split.pdf"), plot = DimPlot(combined.WT.subset, reduction = "umap", label=T, split.by = "split") , height = 5, width = 5)
    DimPlot(combined.WT.subset, reduction = "umap", label=T) 
    ggsave(paste0("./Result_for_WT/04.Cell_Subtype.Result/Plot/",gsub(" ", "_", cell_type),"UMAP_nosplit.pdf"), plot = DimPlot(combined.WT.subset, reduction = "umap", label=T) , height = 5, width = 5)
    
    print(paste0("singleR for", cell_type, "..."))
    meta.WT = combined.WT.subset@meta.data
    combined.WT.subset_for_SingleR = GetAssayData(combined.WT.subset, slot="data")
    combined.WT.subset.hesc.mm = SingleR(test = combined.WT.subset_for_SingleR, ref = Imm.se, labels = Imm.se$label.fine)
    
    print(paste0("Write singleR table for", cell_type, "..."))
    table1 = table(combined.WT.subset.hesc.mm$labels, meta.WT$seurat_clusters)
    write.csv(x = table1, file = paste0("./Result_for_WT/04.Cell_Subtype.Result/RDS/",gsub(" ", "_", cell_type),"SingleR.csv"))}

for (cell_type in cell_types){
  SingleR_Subset(cell_type)
}
