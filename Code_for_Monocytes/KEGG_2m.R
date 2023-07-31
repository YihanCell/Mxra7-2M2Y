library(org.Mm.eg.db)
library(R.utils)
library(MAST)
library(clusterProfiler)
library(ggrepel)
library(ggplot2)

setwd('/data/yihan/Mxra7_2m/')

R.utils::setOption("clusterProfiler.download.method",'auto')  
KEGG_database = 'mmu'

markers = read.delim("/data/yihan/Mxra7_2m/05.edgeRMarkgenes.Result/Granulocytes.control_treat.glmLRT.plot.txt", row.names = 1, sep = '\t', check.names = FALSE)

gene_all = markers[which(markers$sig %in% c('down', 'up')), ]
gene_counts = table(gene_all$sig)
print(gene_counts)
gene_all = rownames(gene_all)
gene_all = as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                      keys = gene_all,
                                                      columns = 'ENTREZID',
                                                      keytype = 'SYMBOL')[,2]))
geneset = gene_all

kegg.histogram.down = as.data.frame(enrichKEGG(gene = geneset,organism= KEGG_database, qvalueCutoff = 1, pvalueCutoff= 1))
write.csv(kegg.histogram.down, "/data/yihan/Mxra7_2m/Result_for_Granulocytes/2m_all_gene_KEGG.csv")

rownames(kegg.histogram.down) = 1:nrow(kegg.histogram.down)
kegg.histogram.down$order=factor(rev(as.integer(rownames(kegg.histogram.down))),labels = rev(kegg.histogram.down$Description))
kegg.histogram.down$Description = gsub(' - Mus musculus \\(house mouse\\)','',kegg.histogram.down$Description)
kegg.histogram.down$order=factor(rev(as.integer(rownames(kegg.histogram.down))),labels = rev(kegg.histogram.down$Description))

selected_pathways = kegg.histogram.down[kegg.histogram.down$pvalue < 0.08, ]
selected_pathways = selected_pathways[selected_pathways$ID != "mmu05171", ]


top_10_pathways = head(kegg.histogram.down, 10)
top_20_pathways = head(selected_pathways, 20)


kegg.plot.down = ggplot(top_20_pathways,aes(y=reorder(order,Count),x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+
  scale_fill_gradient(low = "red",high ="blue" )+
  labs(title = "B cells KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 8),
        axis.title.y = element_text(face = "bold",size = 8),
        legend.title = element_text(face = "bold",size = 8),
        axis.text.y = element_text(size=8)) 

dev.new()
ggplot(top_20_pathways, aes(y = reorder(order, Count), x = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(limits = c(0,0.1),low = "steelblue", high = "orange") +  # Using a different color palette
  labs(title = "Enriched KEGG Pathways in Granulocytes",
       subtitle = "Top pathways based on gene count and adjusted p-value",
       x = "Gene Count",
       y = "Pathways") +  # Add a caption with data source information
  theme_minimal() +  # Using a clean and minimalist theme
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(face = "bold", size = 10),
    axis.title.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8),
    legend.text = element_text(size = 8)
  )
dev.off()
