library(SingleR)
library(RColorBrewer)
library(ggplot2)
library(Seurat)
setwd('/data/yihan/Mxra7/')

color <- c('lightgrey', 'blue','seagreen2')

pdf(file = './03.Annotation.Result/Gene_map1.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'Mmp9',cols = color, pt.size = 1) +  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

pdf(file = './03.Annotation.Result/Gene_map2.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'Wfdc21',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

pdf(file = './03.Annotation.Result/Gene_map3.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'S100a4',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

pdf(file = './03.Annotation.Result/Gene_map4.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'Ms4a6c',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

pdf(file = './03.Annotation.Result/Gene_map5.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'Ccl5',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

pdf(file = './03.Annotation.Result/Gene_map6.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'Ms4a4b',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

pdf(file = './03.Annotation.Result/Gene_map7.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'Cd79a',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

pdf(file = './03.Annotation.Result/Gene_map8.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'Ebf1',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

pdf(file = './03.Annotation.Result/Gene_map9.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'Car2',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

pdf(file = './03.Annotation.Result/Gene_map10.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'Slc4a1',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

pdf(file = './03.Annotation.Result/Gene_map11.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'Klrk1',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

pdf(file = './03.Annotation.Result/Gene_map12.pdf',width = 5,height = 5)
FeaturePlot(combinedh.2y, features = 'Ctla2a',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect( ),         plot.title = element_blank(),legend.position="none",         axis.line=element_blank(),         axis.text.x=element_blank(),         axis.text.y=element_blank(),         axis.ticks=element_blank(),         axis.title.x=element_blank(),         axis.title.y=element_blank())
dev.off()

