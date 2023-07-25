library(org.Mm.eg.db)
library(R.utils)
library(MAST)
library(clusterProfiler)
library(ggrepel)
library(ggplot2)

setwd('/data/yihan/Mxra7_2m/')

R.utils::setOption("clusterProfiler.download.method",'auto')  
KEGG_database = 'mmu'

markers = read.delim("/data/yihan/Mxra7_2Y/05.edgeRMarkgenes.Result/B_cells.control_treat.glmLRT.plot.txt", row.names = 1, sep = '\t', check.names = FALSE)

gene_all = markers[which(markers$sig %in% c('down', 'up')), ]
gene_all = rownames(gene_all)

gene_all = as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                      keys = gene_all,
                                                      columns = 'ENTREZID',
                                                      keytype = 'SYMBOL')[,2]))
geneset = gene_all

kegg.histogram.down = as.data.frame(enrichKEGG(gene = geneset,organism= KEGG_database, qvalueCutoff = 1, pvalueCutoff= 1))
write.csv(kegg.histogram.down, "/data/yihan/Mxra7_2m/Result_for_Bcells/2Y_all_gene_KEGG.csv")

rownames(kegg.histogram.down) = 1:nrow(kegg.histogram.down)
kegg.histogram.down$order=factor(rev(as.integer(rownames(kegg.histogram.down))),labels = rev(kegg.histogram.down$Description))
kegg.histogram.down$Description = gsub(' - Mus musculus \\(house mouse\\)','',kegg.histogram.down$Description)
kegg.histogram.down$order=factor(rev(as.integer(rownames(kegg.histogram.down))),labels = rev(kegg.histogram.down$Description))

top_10_pathways = head(kegg.histogram.down, 10)
kegg.plot.down = ggplot(top_10_pathways,aes(y=reorder(order,Count),x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+
  scale_fill_gradient(low = "red",high ="blue" )+
  labs(title = "B cells KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 8),
        axis.title.y = element_text(face = "bold",size = 8),
        legend.title = element_text(face = "bold",size = 8),
        axis.text.y = element_text(size=8)) 

ggsave("/data/yihan/Mxra7_2m/Result_for_Bcells/2Y_all_gene.KEGG_new.pdf", plot = kegg.plot.down)
