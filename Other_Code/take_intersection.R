setwd('/data/yihan/Mxra7_2m/')
cell_types = c("B_cells",
               "Erythrocytes",
               "Granulocytes",
               "Monocytes",
               "T_cells")

SingleR_Subset <- function(cell_type){
  markersX = read.delim(paste0("./Result_for_X/05.edgeRMarkgenes.Result/", cell_type, ".control_treat.glmLRT.plot.txt"), row.names = 1, sep = '\t', check.names = FALSE)
  gene_upX = markers[which(markersX$sig=='up'),]
  gene_up_nameX = rownames(gene_upX)
  markersWT = read.delim(paste0("./Result_for_WT/05.edgeRMarkgenes.Result/", cell_type, ".control_treat.glmLRT.plot.txt"), row.names = 1, sep = '\t', check.names = FALSE)
  gene_upWT = markers[which(markersWT$sig=='up'),]
  gene_up_nameWT = rownames(gene_upWT)
  intersectionUP = intersect(gene_up_nameX, gene_up_nameWT)
  writeLines(intersectionUP, paste0("./Result_for_interstection/", gsub(" ", "_", cell_type), ".up.txt"))
  
  gene_downX = markers[which(markersX$sig=='down'),]
  gene_down_nameX = rownames(gene_downX)
  gene_downWT = markers[which(markersWT$sig=='down'),]
  gene_down_nameWT = rownames(gene_downWT)
  intersectionDOWN = intersect(gene_down_nameX, gene_down_nameWT)
  writeLines(intersectionDOWN, paste0("./Result_for_interstection/", gsub(" ", "_", cell_type), ".down.txt"))}
  
for (cell_type in cell_types){
  SingleR_Subset(cell_type)
}
