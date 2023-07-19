library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

combined.2m = readRDS("./03.Annotation.Result/RDS/2M_combinedh_afterAnnotation.RDS")
combined.2m.wt =  subset(combined.2m, split == "WT2m")
combined.2m.ko =  subset(combined.2m, split == "X2m")
unique(combined.2m$split)

object.list = list(WT2Y = cellchat.wt, X2Y = cellchat.ko)
cellchat = mergeCellChat(object.list, add.names = names(object.list))
##=== 比较交互总数和交互强度
gg1 = compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 = compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg3 = gg1 + gg2
gg3
ggsave(filename = "./08.CellChat/Plot/Total_Interactions_Num_Interaction_Intensity.pdf",plot = gg3,width = 5,height = 5)

##=== Comparing the number and strength of interactions between different cell populations
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
ggsave(filename = "./08.CellChat/Plot/interactions.pdf",plot = gg3,width = 5,height = 5)

gg4 = netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg5 = netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg6 = gg4 + gg5
gg6

weight.max = getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


# ===
cellchat = computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat = netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat = netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
#> 


# ===
gg7 = rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg8 = rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg9 = gg7 + gg8
gg9

# ====
pathways.show = c("MIF") 
weight.max = getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
