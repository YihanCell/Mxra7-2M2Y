library(org.Mm.eg.db)
library(R.utils)
library(MAST)
library(clusterProfiler)
library(ggrepel)
library(svglite)
library(ggplot2)

setwd('/data/yihan/Mxra7_2m/')

R.utils::setOption("clusterProfiler.download.method",'auto')
cell_types = c("Monocytes")

R.utils::setOption("clusterProfiler.download.method",'auto')
KEGG_database = 'mmu'


goenrichment <- function(cell_type){
  print(paste0("Enriching for ", cell_type, "..."))
  markers <- read.delim(paste0("./Result_for_X/05.edgeRMarkgenes.Result/", cell_type, ".control_treat.glmLRT.plot.txt"), row.names = 1, sep = '\t', check.names = FALSE)
  print("Drawing up genes...")
  gene_up = markers[which(markers$sig=='down'),]
  gene_up = rownames(gene_up)
  print("Rownamed process...")
  gene_up = as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                       keys = gene_up,
                                                       columns = 'ENTREZID',
                                                       keytype = 'SYMBOL')[,2]))
  geneset = gene_up
  markers.screen_type = "up"
  ego_CC = enrichGO(gene = geneset,
                    OrgDb=org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "CC",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
  ego_BP = enrichGO(gene = geneset,
                    OrgDb=org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
  ego_MF = enrichGO(gene = geneset,
                    OrgDb=org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "MF",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
  BP.Go.up = dotplot(ego_BP)
  print("Prepare for save...")
  ggsave(paste0("./Result_for_X/05.edgeRMarkgenes.Result/Enrich_GO/", cell_type, ".down_gene.Go_BP.pdf"), plot = BP.Go.up)
  
  # GO term 柱形图
  display_number = c(10,10,10)
  ego_result_BP = as.data.frame(ego_BP)[1:display_number[1], ]
  ego_result_CC = as.data.frame(ego_CC)[1:display_number[2], ] 
  ego_result_MF = as.data.frame(ego_MF)[1:display_number[3], ]
  
  ##摘取的部分通路重新组合成数据框
  go_enrich_df = data.frame(
    ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID), Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
    GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
    type=factor(c(rep("biological process", display_number[1]), 
                  rep("cellular component", display_number[2]),
                  rep("molecular function", display_number[3])), 
                levels=c("biological process", "cellular component","molecular function" )))
  
  ##通路的名字太长，选取通路的前10个单词作为通路的名字
  for(i in 1:nrow(go_enrich_df)){
    description_splite=strsplit(go_enrich_df$Description[i],split = " ")
    description_collapse=paste(description_splite[[1]][1:10],collapse = " ") 
    go_enrich_df$Description[i]=description_collapse
    go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
  }
  
  ##绘制GO柱状图
  go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
  COLS = c("#66C3A5", "#8DA1CB", "#FD8D62")
  p.Go.up = ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
    geom_bar(stat="identity", width=0.8) + 
    scale_fill_manual(values = COLS) + 
    coord_flip() + 
    xlab("GO term") + 
    ylab("Gene_Number") + 
    labs(title = paste0("The Most Enriched GO Terms ",markers.screen_type))+
    theme_bw() + theme(axis.text.y=element_text(size=8, angle=45, hjust=1))
  ggsave(paste0("./Result_for_X/05.edgeRMarkgenes.Result/Enrich_GO/", cell_type, ".down_gene.Go.pdf"), plot = p.Go.up)}


for (cell_type in cell_types){
  goenrichment(cell_type)
}
