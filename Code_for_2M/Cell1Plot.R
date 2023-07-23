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
colors_new = c("#B2DF8A", "#FB9A99", "#A6CEE3", "#33A02C", "#1F78B4", "#FDBF6F", "#FF33A4", "#CAB2D6", "#6A3D9A", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#FFFF99", "#B15928", "#8DD3C7", "#BEBADA", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#C0C0C0", "#8DD3C7", "#FFFFB3", "#80B1D3", "#B3DE69", "#FCCDE5", "#D9D9D9")
circ_data$PCA = circ_data$seurat_clusters
group_colors<-rand_color(length(names(table(combined.2m$split))))
rep_colors<-rand_color(length(names(table(combined.2m$seurat_clusters))))

###plot and save figures
pdf(file = "./03.Annotation.Result/Plot/Cell1Plot_New.pdf",width = 6,height = 6)
ciecle_plot = plot_circlize(circ_data,do.label = T, pt.size = 0.05, col.use = colors_new,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.8) +
  add_track(circ_data, group = "split", colors = group_colors, track_num = 2) +
  add_track(circ_data, group = "seurat_clusters",colors = rep_colors, track_num = 3) 
dev.off()


