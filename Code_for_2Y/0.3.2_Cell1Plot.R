library(ggplot2)
library(SingleR)
library(RColorBrewer)
library(plot1cell)

###Check and see the meta data info on your Seurat object
colnames(combinedh.2y@meta.data)  

###Prepare data for ploting
circ_data <- prepare_circlize_data(combinedh.2y, scale = 0.7 )
set.seed(1234)
cluster_colors<-c("#33A02C","#FB9A99","#A6CEE3","#1F78B4","#E31A1C","#B2DF8A")
circ_data$PCA = circ_data$seurat_clusters
group_colors<-rand_color(length(names(table(combinedh.2y$split))))
rep_colors<-rand_color(length(names(table(combinedh.2y$seurat_clusters))))

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


