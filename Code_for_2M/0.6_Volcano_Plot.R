library(edgeR)
library(statmod)
library(ggplot2)
library(svglite)
library(ggpubr)
setwd('/data/yihan/Mxra7_2m/')

cell_types = c("B cells",
               "Erythrocytes",
               "Granulocytes",
               "Monocytes",
               "T cells")

markergene<- function(cell_type){
  
  gene_diff <- read.delim(paste0("./05.edgeRMarkgenes.Result/", gsub(" ", "_", cell_type), ".control_treat.glmLRT.plot.txt"), row.names = 1, sep = '\t', check.names = FALSE)
  #横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 FDR
  p.subset.cell_type = ggplot(data = gene_diff, aes(x = logFC, y = -log10(FDR), color = sig)) +
    geom_point(size = 1) +  #绘制散点图
    scale_color_manual(values = c("#2f5688","#BBBBBB","#C00000"), limits = c('up', 'none', 'down')) +  #自定义点的颜色
    labs(x = 'log2 Fold Change', y = '-log10 FDR', title = paste0(gsub(" ", "_", cell_type), " ","DEGs"), color = '') +  #坐标轴标题
    theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
          panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = c(-0.4, 0.4), lty = 3, color = 'black') +  #添加阈值线
    geom_hline(yintercept = 1.5, lty = 3, color = 'black') +
    xlim(-1.4, 1.4) + ylim(0, 35) + #定义刻度边
    coord_flip()
  
  ggsave(paste0("./06.VolcanoPlot.Result/", gsub(" ", "_", cell_type), ".gene.pdf"), plot = p.subset.cell_type, width = 5,height = 5)}

for (cell_type in cell_types){
  markergene(cell_type)
}
