library(CellChat)
library(patchwork)

combined.2m.ko =  subset(combined.2m, split == "X2m")
combined.2m.wt =  subset(combined.2m, split == "WT2m")
combined.2y.ko =  subset(combined.2y, split == "x2y")
combined.2y.wt =  subset(combined.2y, split == "wt2y")

# Setting up the Ligand-Receptor Interaction Database
CellChatDB = CellChatDB.mouse
showDatabaseCategory(CellChatDB.mouse)
CellChatDB.use = subsetDB(CellChatDB.mouse, search = "Secreted Signaling") 

# Build CellChat Object

## wt2y
labels = combined.2y.wt@meta.data$predicted_celltype
data.input.wt2y = GetAssayData(combined.2y.wt, assay = "RNA", slot = "data") # normalized data matrix
cellchat.wt2y = createCellChat(object = data.input.wt2y)
meta1 = data.frame(group = labels, row.names = names(labels)) 
cellchat.wt2y = addMeta(cellchat.wt2y, meta = meta1, meta.name = "labels")
cellchat.wt2y = setIdent(cellchat.wt2y, ident.use = "labels")
cellchat.wt2y@DB = CellChatDB.use
cellchat.wt2y = subsetData(cellchat.wt2y)
cellchat.wt2y = identifyOverExpressedGenes(cellchat.wt2y)
cellchat.wt2y = identifyOverExpressedInteractions(cellchat.wt2y)
cellchat.wt2y = computeCommunProb(cellchat.wt2y)
cellchat.wt2y = filterCommunication(cellchat.wt2y, min.cells = 10)
cellchat.wt2y = computeCommunProbPathway(cellchat.wt2y)
cellchat.wt2y = aggregateNet(cellchat.wt2y)
groupSize.wt2y = as.numeric(table(cellchat.wt2y@idents))
cellchat.wt2y = netAnalysis_computeCentrality(cellchat.wt2y, slot.name = "netP")


## ko2y
labels = combined.2y.ko@meta.data$predicted_celltype
data.input.ko2y = GetAssayData(combined.2y.ko, assay = "RNA", slot = "data") # normalized data matrix
cellchat.ko2y = createCellChat(object = data.input.ko2y)
meta2 = data.frame(group = labels, row.names = names(labels)) 
cellchat.ko2y = addMeta(cellchat.ko2y, meta = meta2, meta.name = "labels")
cellchat.ko2y = setIdent(cellchat.ko2y, ident.use = "labels")
cellchat.ko2y@DB = CellChatDB.use
cellchat.ko2y = subsetData(cellchat.ko2y)
cellchat.ko2y = identifyOverExpressedGenes(cellchat.ko2y)
cellchat.ko2y = identifyOverExpressedInteractions(cellchat.ko2y)
cellchat.ko2y = computeCommunProb(cellchat.ko2y)
cellchat.ko2y = filterCommunication(cellchat.ko2y, min.cells = 10)
cellchat.ko2y = computeCommunProbPathway(cellchat.ko2y)
cellchat.ko2y = aggregateNet(cellchat.ko2y)
groupSize.ko2y = as.numeric(table(cellchat.ko2y@idents))
cellchat.ko2y = netAnalysis_computeCentrality(cellchat.ko2y, slot.name = "netP")

## wt2m
labels = combined.2m.wt@meta.data$predicted_celltype
data.input.wt2m = GetAssayData(combined.2m.wt, assay = "RNA", slot = "data") # normalized data matrix
cellchat.wt2m = createCellChat(object = data.input.wt2m)
meta3 = data.frame(group = labels, row.names = names(labels)) 
cellchat.wt2m = addMeta(cellchat.wt2m, meta = meta3, meta.name = "labels")
cellchat.wt2m = setIdent(cellchat.wt2m, ident.use = "labels")
cellchat.wt2m@DB = CellChatDB.use
cellchat.wt2m = subsetData(cellchat.wt2m)
cellchat.wt2m = identifyOverExpressedGenes(cellchat.wt2m)
cellchat.wt2m = identifyOverExpressedInteractions(cellchat.wt2m)
cellchat.wt2m = computeCommunProb(cellchat.wt2m)
cellchat.wt2m = filterCommunication(cellchat.wt2m, min.cells = 10)
cellchat.wt2m = computeCommunProbPathway(cellchat.wt2m)
cellchat.wt2m = aggregateNet(cellchat.wt2m)
groupSize.wt2m = as.numeric(table(cellchat.wt2m@idents))
cellchat.wt2m = netAnalysis_computeCentrality(cellchat.wt2m, slot.name = "netP")

## ko2m
labels = combined.2m.ko@meta.data$predicted_celltype
data.input.ko2m = GetAssayData(combined.2m.ko, assay = "RNA", slot = "data") # normalized data matrix
cellchat.ko2m = createCellChat(object = data.input.ko2m)
meta4 = data.frame(group = labels, row.names = names(labels)) 
cellchat.ko2m = addMeta(cellchat.ko2m, meta = meta4, meta.name = "labels")
cellchat.ko2m = setIdent(cellchat.ko2m, ident.use = "labels")
cellchat.ko2m@DB = CellChatDB.use
cellchat.ko2m = subsetData(cellchat.ko2m)
cellchat.ko2m = identifyOverExpressedGenes(cellchat.ko2m)
cellchat.ko2m = identifyOverExpressedInteractions(cellchat.ko2m)
cellchat.ko2m = computeCommunProb(cellchat.ko2m)
cellchat.ko2m = filterCommunication(cellchat.ko2m, min.cells = 10)
cellchat.ko2m = computeCommunProbPathway(cellchat.ko2m)
cellchat.ko2m = aggregateNet(cellchat.ko2m)
groupSize.ko2m = as.numeric(table(cellchat.ko2m@idents))
cellchat.ko2m = netAnalysis_computeCentrality(cellchat.ko2m, slot.name = "netP")


saveRDS(object = cellchat.ko2m,file = "./Result_for_CellChat/ko2m.rds")
saveRDS(object = cellchat.wt2y,file = "./Result_for_CellChat/wt2y.rds")
saveRDS(object = cellchat.wt2m,file = "./Result_for_CellChat/wt2m.rds")
saveRDS(object = cellchat.ko2y,file = "./Result_for_CellChat/ko2y.rds")