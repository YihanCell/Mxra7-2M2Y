library(ggplot2)
library(SingleR)
library(RColorBrewer)
library(plot1cell)

###Check and see the meta data info on your Seurat object
colnames(combined.2m@meta.data)  

###Prepare data for ploting
circ_data <- prepare_circlize_data(combined.2m, scale = 0.7 )
set.seed(1234)
cluster_colors<-c("#B2DF8A","#FB9A99","#A6CEE3","#33A02C","#1F78B4")
circ_data$PCA = circ_data$seurat_clusters
group_colors<-rand_color(length(names(table(combined.2m$split))))
rep_colors<-rand_color(length(names(table(combined.2m$seurat_clusters))))

###plot and save figures
pdf(file = "./03.Annotation.Result/Plot/Cell1Plot.pdf",width = 6,height = 6)
ciecle_plot = plot_circlize(circ_data,do.label = T, pt.size = 0.05, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.8) +
  add_track(circ_data, group = "split", colors = group_colors, track_num = 2) +
  add_track(circ_data, group = "PCA",colors = rep_colors, track_num = 3) 
dev.off()
png(file = "./03.Annotation.Result/Plot/Cell1Plot.png",width = 6,height = 6,res = 350)
par(mar = c(0, 0, 0, 0))
ciecle_plot = plot_circlize(circ_data,do.label = T, pt.size = 0.05, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.8) +
  add_track(circ_data, group = "split", colors = group_colors, track_num = 2) +
  add_track(circ_data, group = "PCA",colors = rep_colors, track_num = 3) 
dev.off()


