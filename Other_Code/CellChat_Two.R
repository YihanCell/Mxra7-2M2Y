library(CellChat)
library(patchwork)

object.list.all = list(wt2m = cellchat.wt2m, ko2m = cellchat.ko2m, wt2y = cellchat.wt2y, ko2y = cellchat.ko2y)
cellchat.all = mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
cellchat.all

object.list.2m = list(wt2m = cellchat.wt2m, ko2m = cellchat.ko2m)
cellchat.2m = mergeCellChat(object.list.2m, add.names = names(object.list.2m))
cellchat.2m

object.list.2y = list(wt2y = cellchat.wt2y, ko2y = cellchat.ko2y)
cellchat.2y = mergeCellChat(object.list.2y, add.names = names(object.list.2y))
cellchat.2y

dev.new()
gg1 = compareInteractions(cellchat.all, show.legend = F, group = c(1,2))
gg2 = compareInteractions(cellchat.all, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

dev.new()
netVisual_heatmap(cellchat.2m)
dev.off()
dev.new()
netVisual_heatmap(cellchat.2m, measure = "weight")
dev.off()

custom_variable_order <- c("Granulocytes", "T cells", "B cells", "Monocytes", "Erythrocytes")
dev.new()
netVisual_heatmap(cellchat.2y)
dev.off()
dev.new()
netVisual_heatmap(cellchat.2y, measure = "weight")
dev.off()

dev.new()
netVisual_heatmap(cellchat)
dev.off()

#===========
dev.new() 
gg1 <- rankNet(cellchat.2y, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat.2y, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()
