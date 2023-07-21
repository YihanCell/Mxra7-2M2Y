library(ggplot2)
library(ggpubr)
library(tidyr)
library(rstatix)

cellper.new = read.csv(file = "./cellper.csv")
# Convert wide format to long format
cellper_long <- cellper.new %>% 
  gather(key = "Cell_type", value = "Percent", -group)

# Create a violin plot with jittered data points
ggplot(cellper_long, aes(x = Cell_type, y = Percent, fill = group)) +
  geom_violin(scale = "width", alpha = 0.5) +
  geom_jitter(position = position_jitter(width = 0.1), size = 1, shape = 21) +
  stat_compare_means(comparisons = list(c("x2y", "wt2y")), method = "t.test") +
  scale_fill_manual(values = c("grey", "black")) +
  labs(x = "Cell type", y = "Percent", fill = "Group")

ggplot(cellper_long, aes(x = Cell_type, y = Percent, fill = group)) +
  geom_violin(scale = "width", alpha = 0.5) +
  geom_jitter(position = position_jitter(width = 0.1), size = 1, shape = 21) +
  stat_compare_means(comparisons = list(c("x2y", "wt2y")), method = "t.test") +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF")) + # 修改颜色
  scale_color_manual(values = c("#440154FF", "#FDE725FF")) +
  labs(x = "Cell type", y = "Percent", fill = "Group",
       title = "Cell type percentages in different groups", 
       subtitle = "Comparison between x2y and wt2y") + # 添加标题和子标题
  annotate("text", x = 1, y = 0.9, label = "P value = 0.015", size = 4, color = "black") + # 添加注释
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 22, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 18, color = "black"),
        legend.position = "bottom") + # 修改字体和标签
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + # 修改 y 轴范围和显示格式
  scale_x_discrete(labels = c("Granulocytes", "Monocytes", "T cells", "B cells", "Erythrocytes", "NK cells")) # 修改 x 轴标签
