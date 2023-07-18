library(org.Mm.eg.db)
library(R.utils)
library(MAST)
library(clusterProfiler)
library(ggrepel)
library(svglite)
library(ggplot2)

setwd('/data/yihan/Mxra7_2m/')

R.utils::setOption("clusterProfiler.download.method",'auto')  
cell_types = c("B_cells","T_cells","Monocytes","Granulocytes","Erythrocytes")
KEGG_database = 'mmu'

goenrichment = function(cell_type){
  print(paste0("Enriching for ", cell_type, "..."))
  markers = read.delim(paste0("./Result_for_X/06/", cell_type, ".control_treat.glmLRT.plot.txt"), row.names = 1, sep = '\t', check.names = FALSE)
  insection = read.delim(paste0("/data/yihan/Mxra7_2m/Result_for_interstection/more/",cell_type,".up.txt"), header = FALSE, sep = '\t', stringsAsFactors = FALSE)[, 1]
  print("Drawing up genes...")
  gene_up = markers[which(markers$sig=='up'),]
  gene_up = rownames(gene_up)
  gene_up = setdiff(gene_up, insection)
  
  gene_up = as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                       keys = gene_up,
                                                       columns = 'ENTREZID',
                                                       keytype = 'SYMBOL')[,2]))
  geneset = gene_up
  
  kegg.histogram.up = as.data.frame(enrichKEGG(gene = geneset,organism= KEGG_database, qvalueCutoff = 1, pvalueCutoff= 1))
  write.csv(kegg.histogram.up, paste0("/data/yihan/Mxra7_2m/Result_for_X/07.Enrichment.Result.KEGG/",cell_type,".up_gene.KEGG.csv"))
  
  rownames(kegg.histogram.up) = 1:nrow(kegg.histogram.up)
  kegg.histogram.up$order=factor(rev(as.integer(rownames(kegg.histogram.up))),labels = rev(kegg.histogram.up$Description))
  kegg.histogram.up$Description = gsub(' - Mus musculus \\(house mouse\\)','',kegg.histogram.up$Description)
  kegg.histogram.up$order=factor(rev(as.integer(rownames(kegg.histogram.up))),labels = rev(kegg.histogram.up$Description))
  
  top_10_pathways = head(kegg.histogram.up, 10)
  kegg.plot.up = ggplot(top_10_pathways,aes(y=reorder(order,Count),x=Count,fill=p.adjust))+
    geom_bar(stat = "identity",width=0.7)+
    scale_fill_gradient(low = "red",high ="blue" )+
    labs(title = paste0("KEGG Pathways Enrichment",cell_type),
         x = "Gene numbers", 
         y = "Pathways")+
    theme(axis.title.x = element_text(face = "bold",size = 8),
          axis.title.y = element_text(face = "bold",size = 8),
          legend.title = element_text(face = "bold",size = 8),
          axis.text.y = element_text(size=8)) 
  kegg.plot.up
  ggsave(paste0("./Result_for_X/07.Enrichment.Result.KEGG/", cell_type, ".up_gene.KEGG.pdf"), plot = kegg.plot.up)
}

for (cell_type in cell_types){
  goenrichment(cell_type)
}
