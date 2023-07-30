gene_diff_2Y = read.delim("/data/yihan/Mxra7_2Y/05.edgeRMarkgenes.Result/Granulocytes.control_treat.glmLRT.plot.txt", row.names = 1, sep = '\t', check.names = FALSE)

#横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 FDR
dev.new()
ggplot(data = gene_diff_2Y, aes(x = logFC, y = -log10(FDR), color = sig)) +
  geom_point(size = 1) +  #绘制散点图
  scale_color_manual(values = c("#2f5688","#BBBBBB","#C00000"), limits = c('up', 'none', 'down')) +  #自定义点的颜色
  labs(x = 'log2 Fold Change', y = '-log10 FDR',  color = '') +  #坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.4, 0.4), lty = 3, color = 'black') +  #添加阈值线
  geom_hline(yintercept = 1.5, lty = 3, color = 'black') +
  xlim(-1.4, 1.4) + ylim(0, 100)

dev.off()

gene_diff_2m = read.delim("/data/yihan/Mxra7_2m/05.edgeRMarkgenes.Result/Granulocytes.control_treat.glmLRT.plot.txt", row.names = 1, sep = '\t', check.names = FALSE)

#横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 FDR
dev.new()
ggplot(data = gene_diff_2m, aes(x = logFC, y = -log10(FDR), color = sig)) +
  geom_point(size = 1) +  #绘制散点图
  scale_color_manual(values = c("#2f5688","#BBBBBB","#C00000"), limits = c('up', 'none', 'down')) +  #自定义点的颜色
  labs(x = 'log2 Fold Change', y = '-log10 FDR',  color = '') +  #坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.4, 0.4), lty = 3, color = 'black') +  #添加阈值线
  geom_hline(yintercept = 1.5, lty = 3, color = 'black') +
  xlim(-1.4, 1.4) + ylim(0, 100)

dev.off()
