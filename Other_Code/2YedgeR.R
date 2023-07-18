library(Seurat)
library(edgeR)
library(limma)
setwd('/data/yihan/Mxra7_2m/')
cell_types = c("B cells",
               "Erythrocytes",
               "Granulocytes",
               "Monocytes",
               "NK cells",
               "T cells")

combined.2y = readRDS("/data/yihan/Mxra7_2Y/03.Annotation.Result/RDS/combinedh_afterAnnotation.rds")
cell_type_estimateDisp <- function(cell_type){
  
  print(paste0("Subsetting for ", cell_type, "..."))
  combined.2y.subset =  subset(combined.2y, predicted_celltype == cell_type)
  counts.subset = as.matrix(GetAssayData(combined.2y.subset, assay = "RNA"))
  rownames(counts.subset) = rownames(combined.2y.subset@assays$RNA@data)
  group.subset = as.factor(combined.2y.subset@meta.data$split)
  
  print("Generating DGEList...")
  dge.list = DGEList(counts = counts.subset, group = group.subset)
  
  print("Calculating normalization factors...")
  dge = calcNormFactors(dge.list)
  design = model.matrix(~group.subset)
  
  print("Estimating dispersion...")
  dge.subset = estimateDisp(dge, design, robust = TRUE)
  saveRDS(dge.subset, paste0("./KEGGResult/2Y/", gsub(" ", "_", cell_type), "_Disp.rds"))
  
  print("Building model...")
  fit = glmFit(dge.subset, design, robust = TRUE)
  lrt = topTags(glmLRT(fit), n = nrow(dge.list$counts))
  write.table(lrt, paste0("./KEGGResult/2Y/",gsub(" ", "_", cell_type),"control_treat.glmLRT.txt"), sep = '\t', col.names = NA, quote = FALSE)
  gene_diff = read.delim(paste0("./KEGGResult/2Y/",gsub(" ", "_", cell_type),"control_treat.glmLRT.txt"), row.names = 1, sep = '\t', check.names = FALSE)
  
  print("Setting genes order")
  gene_diff <- gene_diff[order(gene_diff$FDR, gene_diff$logFC, decreasing = c(FALSE, TRUE)), ]
  
  #log2FC≥0.4 & FDR<0.01 标识 up，代表显著上调的基因
  #log2FC≤-0.4 & FDR<0.01 标识 down，代表显著下调的基因
  #其余标识 none，代表非差异的基因
  print("Find up&down marker genes")
  gene_diff[which(gene_diff$logFC >= 0.4 & gene_diff$FDR < 0.05),'sig'] <- 'up'
  gene_diff[which(gene_diff$logFC <= -0.4 & gene_diff$FDR < 0.05),'sig'] <- 'down'
  gene_diff[which(abs(gene_diff$logFC) <= 0.4 | gene_diff$FDR >= 0.05),'sig'] <- 'none'
  write.table(gene_diff, file = paste0("./KEGGResult/2Y/", gsub(" ", "_", cell_type), ".control_treat.glmLRT.plot.txt"), sep = '\t', col.names = NA, quote = FALSE)
  
  #输出选择的差异基因总表
  print("Write marker genes")
  gene_diff_select <- subset(gene_diff, sig %in% c('up', 'down'))
  write.table(gene_diff_select, file = paste0("./KEGGResult/2Y/", gsub(" ", "_", cell_type), ".control_treat.glmLRT.select.txt"), sep = '\t', col.names = NA, quote = FALSE)
  
  #根据 up 和 down 分开输出
  gene_diff_up <- subset(gene_diff, sig == 'up')
  gene_diff_down <- subset(gene_diff, sig == 'down')
  
  write.table(gene_diff_up, file = paste0("./KEGGResult/2Y/", gsub(" ", "_", cell_type), ".control_treat.glmLRT.up.txt"), sep = '\t', col.names = NA, quote = FALSE)
  write.table(gene_diff_down, file = paste0("./KEGGResult/2Y/", gsub(" ", "_", cell_type), ".control_treat.glmLRT.down.txt"), sep = '\t', col.names = NA, quote = FALSE)}


for (cell_type in cell_types){
  cell_type_estimateDisp(cell_type)
}
