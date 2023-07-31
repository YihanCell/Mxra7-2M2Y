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

## ko2y
labels = combined.2y.ko@meta.data$predicted_celltype
data.input.ko2y = GetAssayData(combined.2y.ko, assay = "RNA", slot = "data") # normalized data matrix
cellchat.ko2y = createCellChat(object = data.input.ko2y)
meta2 = data.frame(group = labels, row.names = names(labels)) 
cellchat.ko2y = addMeta(cellchat.ko2y, meta = meta2, meta.name = "labels")
cellchat.ko2y = setIdent(cellchat.ko2y, ident.use = "labels")
cellchat.ko2y@DB = CellChatDB.use

## wt2m
labels = combined.2m.wt@meta.data$predicted_celltype
data.input.wt2m = GetAssayData(combined.2m.wt, assay = "RNA", slot = "data") # normalized data matrix
cellchat.wt2m = createCellChat(object = data.input.wt2m)
meta3 = data.frame(group = labels, row.names = names(labels)) 
cellchat.wt2m = addMeta(cellchat.wt2m, meta = meta1, meta.name = "labels")
cellchat.wt2m = setIdent(cellchat.wt2m, ident.use = "labels")
cellchat.wt2m@DB = CellChatDB.use

## ko2m
labels = combined.2m.ko@meta.data$predicted_celltype
data.input.ko2m = GetAssayData(combined.2m.ko, assay = "RNA", slot = "data") # normalized data matrix
cellchat.ko2m = createCellChat(object = data.input.ko2m)
meta4 = data.frame(group = labels, row.names = names(labels)) 
cellchat.ko2m = addMeta(cellchat.ko2m, meta = meta1, meta.name = "labels")
cellchat.ko2m = setIdent(cellchat.ko2m, ident.use = "labels")
cellchat.wt2m@DB = CellChatDB.use
