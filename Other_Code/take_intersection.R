setwd('/data/yihan/Mxra7_2m/')
cell_types = c("B cells",
               "Erythrocytes",
               "Granulocytes",
               "Monocytes",
               "T cells")

SingleR_Subset <- function(cell_type){
  file1 = read.table(paste0("./Result_for_WT/05.edgeRMarkgenes.Result/", gsub(" ", "_", cell_type), ".control_treat.glmLRT.up.txt"))
  genes1 = rownames(file1)
  file2 = read.table(paste0("./Result_for_X/05.edgeRMarkgenes.Result/", gsub(" ", "_", cell_type), ".control_treat.glmLRT.up.txt"))
  genes2 = rownames(file2)
  intersection = intersect(genes1, genes2)
  writeLines(intersection, paste0("./Result_for_combine/", gsub(" ", "_", cell_type), ".up.txt"))
  
  file3 = read.table(paste0("./Result_for_WT/05.edgeRMarkgenes.Result/", gsub(" ", "_", cell_type), ".control_treat.glmLRT.down.txt"))
  genes3 = rownames(file3)
  file4 = read.table(paste0("./Result_for_X/05.edgeRMarkgenes.Result/", gsub(" ", "_", cell_type), ".control_treat.glmLRT.down.txt"))
  genes4 = rownames(file4)
  intersection = intersect(genes3, genes4)
  writeLines(intersection, paste0("./Result_for_combine/", gsub(" ", "_", cell_type), ".down.txt"))}

for (cell_type in cell_types){
  SingleR_Subset(cell_type)
}
