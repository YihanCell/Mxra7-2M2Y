setwd('/data/yihan/M2mra7_2m/')

markers2m = read.delim("./05.edgeRMarkgenes.Result/B_cells.control_treat.glmLRT.plot.txt", row.names = 1, sep = '\t', check.names = FALSE)
markers2y = read.delim("/data/yihan/Mxra7_2Y/05.edgeRMarkgenes.Result/B_cells.control_treat.glmLRT.plot.txt", row.names = 1, sep = '\t', check.names = FALSE)
gene_up2m = markers2m[which(markers2m$sig=='up'),]
gene_up2y = markers2y[which(markers2y$sig=='up'),]
gene_down2m = markers2m[which(markers2m$sig=='down'),]
gene_down2y = markers2y[which(markers2y$sig=='down'),]
intersectionUP = intersect(gene_up2m, gene_up2y)
intersectionDOWN = intersect(gene_down2m, gene_down2y)
writeLines(intersectionUP, "./Result_for_interstection/Bcell2yvs2m.up.txt")
writeLines(intersectionUP, "./Result_for_interstection/Bcell2yvs2m.down.txt")

