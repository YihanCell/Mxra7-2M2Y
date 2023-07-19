setwd('/data/yihan/Mxra7_2m/')
combinedh.WT = readRDS("./Result_for_WT/03.Annotation.Result/RDS/WT_combinedh_afterAnnotation.RDS")

sub = table(combinedh.WT@meta.data$split, combinedh.WT@meta.data$predicted_celltype)
sub
prop = prop.table(sub, margin = 1)
prop
cell_count = as.data.frame(table(combinedh.WT@meta.data$split, combinedh.WT@meta.data$predicted_celltype))
cell_count = as.data.frame(prop)
head(cell_count)
names(cell_count) = c("Group", "CellType", "Percentage")
write.csv(x = cell_count,file = "./Result_for_WT/09.Proportation/cell_count.csv")

Per1 = ggplot(cell_count, aes(x=Group, y=Percentage, fill=CellType)) +
  geom_col(width = 0.6) +
  scale_fill_brewer(palette = "Paired") +
  xlab(NULL) + ylab(NULL) +
  theme_minimal() +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
Per1
ggsave('./Result_for_WT/09.Proportation/X_CelltypePercentage_nolegend.pdf', plot = Per1, height = 5, width = 5)

Per2 = ggplot(cell_count, aes(x=Group, y=Percentage, fill=CellType)) +
  geom_col() +
  scale_fill_brewer(palette = "Paired") +
  ggtitle("Cell type distribution in two groups") +
  xlab("Group") + ylab("Cell Percentage") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        legend.title = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
Per2
ggsave('./Result_for_WT/09.Proportation/X_CelltypePercentage.pdf', plot = Per2, height = 5, width = 8)


#my_col = brewer.pal(9, "Oranges")[3:6]

Per3 = ggplot(cell_count, aes(x=CellType, y=Group, fill=Percentage)) +
  geom_tile() +
  #scale_fill_gradientn(colours = my_col) +
  ggtitle("Percentage of Cells by Type") +
  xlab("Cell Type") + ylab("Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        legend.position = "right", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
Per3
ggsave('./Result_for_WT/09.Proportation/X_CelltypePercentage_Compare.pdf', plot = Per3, height = 5, width = 6)

