library("gplots")
library("limma")
library("edgeR")
setwd("E:\\daily_work\\926\\gai")
expr=read.table("combined_sine_GSE98422.txt",sep="\t",header=T,check.names=F) 
expr=as.matrix(expr)  
rownames(expr)=expr[,1]  
exp=expr[,2:ncol(expr)] 
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data) 
data=data[rowMeans(data)>1,] 

group=c(rep("health",19),rep("LN",6)) 
design <- model.matrix(~group) 
y <- DGEList(counts=data,group=group) 

y <- calcNormFactors(y)

y <- estimateCommonDisp(y)

y <- estimateTagwiseDisp(y)

et <- exactTest(y,pair = c("health","LN"))

topTags(et)

ordered_tags <- topTags(et, n=100000) 



allDiff=ordered_tags$table

allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]

diff=allDiff

newData=y$pseudo.counts


padj = 0.05 # ?Զ???
foldChange= 1 # ?Զ???

write.table(diff,file="list_normal_squam_counts_DE_sine_GSE98422.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]

write.table(diffSig, file="list_normal_squam_counts_diffSig_sine_GSE98422.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up_list_normal_squam_counts_sine_GSE98422.xls",sep="\t",quote=F)
diffDown = diff[(diff$PValue < padj & (diff$logFC<(-foldChange))),]

write.table(diffDown, file="down_list_normal_squam_counts_sine_GSE98422.xls",sep="\t",quote=F)



下方为诗博注释


#cell
library("gplots")
library("limma")
library("edgeR")
setwd("/bio/cell_and_tissue/cell-cancer-normal/featurecount_cell_DEG_final/LTR/A549-GSE88943/")
expr=read.table("A549-GSE88943_LTR_featureCounts.txt",sep="\t",header=T,check.names=F) 
expr=as.matrix(expr)  #矩阵化
rownames(expr)=expr[,1]  #第一列作为行名#
exp=expr[,2:ncol(expr)] #第二列到最后一列是表达的数据
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)  #这两行将带引号的数据转换成数值
data=avereps(data) #有的基因出现过多行，把出现多行的gene取平均值
data=data[rowMeans(data)>1,] #去除低表达的数据#

group=c(rep("CT",2),rep("IPF",3))   
design <- model.matrix(~group)  #把group设置成一个model matrix#
y <- DGEList(counts=data,group=group) #group哪些是正常，哪些是样本，让edgeR可以识别#

y <- calcNormFactors(y) #对因子矫正#

y <- estimateCommonDisp(y)#25,26估计变异系数，即估计方差；估计内部差异程度，看组间差异是否比内部差异大，如果大，可选为差异基因#

y <- estimateTagwiseDisp(y)

et <- exactTest(y,pair = c("CT","IPF"))

topTags(et)

ordered_tags <- topTags(et, n=100000) #显示排名前十万的#



allDiff=ordered_tags$table

allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]

diff=allDiff

newData=y$pseudo.counts


padj = 0.05 # 自定义
foldChange= 1 # 自定义

write.table(diff,file="/bio/cell_and_tissue/cell-cancer-normal_XLS/LTR/A549-GSE88943/A549-GSE88943_DE_LTR.xls",sep="\t",quote=F)#先把所有差异表达，输入这个叫edgerOut的文件#
diffSig = diff[(diff$PValue < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]#筛选有显著差异的#

write.table(diffSig, file="/bio/cell_and_tissue/cell-cancer-normal_XLS/LTR/A549-GSE88943/A549-GSE88943_diffSig_LTR.xls",sep="\t",quote=F)#输出有显著差异表达的到diffSig这个文件#
diffUp = diff[(diff$PValue < padj & (diff$logFC>foldChange)),]#foldchange>0是上调，foldchange<0是下调#

write.table(diffUp, file="/bio/cell_and_tissue/cell-cancer-normal_XLS/LTR/A549-GSE88943/up_A549-GSE88943_LTR.xls",sep="\t",quote=F)#39-42把上调和下调分别输入up和down两个文件#

diffDown = diff[(diff$PValue < padj & (diff$logFC<(-foldChange))),]

write.table(diffDown, file="/bio/cell_and_tissue/cell-cancer-normal_XLS/LTR/A549-GSE88943/down_A549-GSE88943_LTR.xls",sep="\t",quote=F)