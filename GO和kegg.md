rm(list = ls())
setwd("/data/nas1/xujiajie/13-87-BJTC-731-12/04.GO_KEGG/")
if(!dir.exists("01_GO")){dir.create("01_GO")}
setwd("01_GO")

#GO_KEGG------
library(clusterProfiler)
library(stringr)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(enrichplot)
library(GOplot)

options(stringsAsFactors = F)
####GO####
genes <- read.csv('/data/nas1/xujiajie/13-87-BJTC-731-12/03.venn/DEG_Candidate.csv')
#genes$x <- substr(genes$x, 3, nchar(genes$x))
GO <- enrichGO(gene = genes$x,
               OrgDb  ="org.Hs.eg.db",
               keyType = "SYMBOL",
               pAdjustMethod = "BH",
               pvalueCutoff =1,
               minGSSize = 5,
               ont="all",
               readable =T)

GO_result <- GO@result 
go <- as.data.frame(GO_result)
rownames(go) <- 1:nrow(go)
go <- go[go$pvalue <= 0.05,]
go <- go[order(go$Count,decreasing = T),]
dim(go) # 223  10
table(go$ONTOLOGY)
# BP  CC  MF 
# 187  11  25
write.table(go,file = "00.GO.txt",sep = "\t",quote = F,row.names = F)

save(GO,file = 'GO.Rdata')

go_res <- go

Go_bp <- go_res[go_res$ONTOLOGY == 'BP', ]
Go_bp <- arrange(Go_bp, desc(Count))
Go_bp <- Go_bp[1:5,]

# 提取CC的行并排序
Go_cc <- go_res[go_res$ONTOLOGY == 'CC', ]
Go_cc <- arrange(Go_cc, desc(Count))
Go_cc <- Go_cc[1:5,]

# 提取MF的行并排序
Go_mf <- go_res[go_res$ONTOLOGY == 'MF', ]
Go_mf <- arrange(Go_mf, desc(Count))
Go_mf <- Go_mf[1:5,]

go_res2 <- rbind(Go_bp, Go_cc, Go_mf)
write.csv(go_res2, file = "./GO_res2.csv")


DEGs <- read.csv("/data/nas1/xujiajie/13-87-BJTC-731-12/02.DEG/02.DEG_sig_res.csv",row.names = 1)
id.fc <- DEGs
id.fc$X <- rownames(id.fc)
genelist <- data.frame(ID = id.fc$X, logFC = id.fc$logFC)#列名

go <-  data.frame(Category = go_res2[,'ONTOLOGY'],
                  ID = go_res2[,'ID'],
                  Term = go_res2[,'Description'], 
                  Genes = gsub("/", ", ", go_res2[,'geneID']), 
                  adj_pval = go_res2[,"pvalue"])

go_circle <- circle_dat(go, genelist)

pdf('GO_ChordalChart.pdf',w=18,h=8,family = 'Times')
GOCircle(data=go_circle,
         nsub=15, ###指定显示通路个数
         rad1 = 2.5, rad2 = 3.5, ##rad1和rad2分别代表内圈和外圈的大小
         zsc.col=c("#96C37D",'white','#5F97D2'),  ##内圈Z-score颜色设置
         lfc.col = c('#96C37D','#5F97D2'),  ##logFC颜色设置
         label.size = 3.2,  ##标签字号设置
         label.fontface='plain',   ##标签字体设置
         table.legend = T   ##是否显示右侧表格，若为F则不显示
)
dev.off()

png('GO_ChordalChart.png',w=1800,h=800,family = 'Times')
GOCircle(data=go_circle,
         nsub=15, ###指定显示通路个数
         rad1 = 2.5, rad2 = 3.5, ##rad1和rad2分别代表内圈和外圈的大小
         zsc.col=c("#96C37D",'white','#5F97D2'),  ##内圈Z-score颜色设置
         lfc.col = c('#96C37D','#5F97D2'),  ##logFC颜色设置
         label.size = 3.3,  ##标签字号设置
         label.fontface='plain',   ##标签字体设置
         table.legend = T   ##是否显示右侧表格，若为F则不显示
)
dev.off()

####KEGG####
setwd("/data/nas1/xujiajie/13-87-BJTC-731-12/04.GO_KEGG/")
if(!dir.exists("02_KEGG")){dir.create("02_KEGG")}
setwd("02_KEGG")
gene_transform <- bitr(genes$x,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 organism = "hsa",
                 keyType = "kegg",
                 pvalueCutoff =0.05,
                 qvalueCutoff =1)

kk <- setReadable(kk, #前面分析的结果
                  OrgDb = "org.Hs.eg.db", #人类数据库
                  keyType = "ENTREZID") #要转换的基因类型

kegg_result <- kk@result 
hh <- as.data.frame(kegg_result)
rownames(hh) <- 1:nrow(hh)
hh <- hh[hh$pvalue <= 0.05,]
hh <- hh[order(hh$Count,decreasing = T),]
dim(hh) #4 9
write.table(hh,file = "KEGG.txt",sep = "\t",quote = F,row.names = F)

pdf('KEGG_ChordalChart.pdf',w=6,h=5,family = 'Times')
ggplot(hh[1:4,], aes(pvalue, Description)) +
  geom_point(aes(y=reorder(Description,pvalue),color=Count))+
  labs(y=(""),x="pvalue") +
  theme_bw()+
  theme(axis.text.x=element_text(vjust=0.5,
                                 size = 10),
        axis.text.y=element_text(size = 10),
        axis.text = element_text(color = 'black', size = 12)
  )+
  scale_color_gradient(low = '#d90424', high = '#374a89')+
  guides(size=guide_legend(order=3))
dev.off()

png('KEGG_ChordalChart.png', w=6, h=5, units = 'in',res = 600)
ggplot(hh[1:4,], aes(pvalue, Description)) +
  geom_point(aes(y=reorder(Description,pvalue),color=Count))+
  labs(y=(""),x="pvalue") +
  theme_bw()+
  theme(axis.text.x=element_text(vjust=0.5,
                                 size = 10),
        axis.text.y=element_text(size = 10),
        axis.text = element_text(color = 'black', size = 12)
  )+
  scale_color_gradient(low = '#d90424', high = '#374a89')+
  guides(size=guide_legend(order=3))
dev.off()


# rownames(hh) <- 1:nrow(hh)
# hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
# kk <- hh[1:10,]
# paste0(kk$Description,collapse = "（）；")
# kegg <-  data.frame(Category = hh[,1],ID = hh[,'ID'],Term = hh[,'Description'], 
#                     Genes = gsub("/", ", ", hh[,'geneID']), adj_pval = hh[,'pvalue'])
# library(tibble)
# library(GOplot)
# 
# id.fc <- DEGs
# id.fc$X <- rownames(id.fc)
# #基因变化倍数
# head(id.fc)
# genelist <- data.frame(ID = id.fc$X, logFC = id.fc$log2FoldChange)          
# 
# #把富集分析和倍数整合在一起
# circ <- circle_dat(kegg, genelist)
# head(circ)
# circ.gsym <- circ
# 
# # KEGG--八卦图--
# pdf('KEGG_ChordalChart.pdf',w=12,h=6,family = 'Times')
# GOCircle(data=circ,
#          nsub=10, ###指定显示通路个数
#          rad1 = 2, rad2 = 3, ##rad1和rad2分别代表内圈和外圈的大小
#          zsc.col=c("#96C37D",'white','#5F97D2'),  ##内圈Z-score颜色设置
#          lfc.col = c('#96C37D','#5F97D2'),  ##logFC颜色设置
#          label.size = 5,  ##标签字号设置
#          label.fontface='plain',   ##标签字体设置
#          table.legend = T   ##是否显示右侧表格，若为F则不显示
# )
# dev.off()
# 
# png('KEGG_ChordalChart.png',w=1200,h=600,family = 'Times')
# GOCircle(data=circ,
#          nsub=10, ###指定显示通路个数
#          rad1 = 2, rad2 = 3, ##rad1和rad2分别代表内圈和外圈的大小
#          zsc.col=c("#96C37D",'white','#5F97D2'),  ##内圈Z-score颜色设置
#          lfc.col = c('#96C37D','#5F97D2'),  ##logFC颜色设置
#          label.size = 5,  ##标签字号设置
#          label.fontface='plain',   ##标签字体设置
#          table.legend = T   ##是否显示右侧表格，若为F则不显示
# )
# dev.off()
