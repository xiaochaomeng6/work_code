library(Seurat) 
dir="E:\\data1\\dalunwen\\single_cell\\filtered_feature_bc_matrix" list.files(dir) 
counts <- Read10X(data.dir = dir) 
class(counts) 
scRNA <- CreateSeuratObject(counts = counts) 
scRNA sce <-CreateSeuratObject(counts,min.cells = 10,min.features = 200)

sce

counts[c("CD3D","TCL1A","MS4A1"),1:30] sce[['percent.mt']]<-PercentageFeatureSet(sce,pattern='^MT-') head(sce@meta.data,5) vin=VlnPlot(sce,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol =3) library(ggplot2) ggsave(filename="D:\\R-4.2.2\\tupian\\w1\\vin.pdf",plot=vin) plot1 <- FeatureScatter(sce,feature1 ="nCount_RNA",feature2 ="percent.mt") plot2 <- FeatureScatter(sce,feature1 ="nCount_RNA",feature2 ="nFeature_RNA") ggsave(filename="D:\\R-4.2.2\\tupian\\w1\\1.pdf",plot=plot1) ggsave(filename="D:\\R-4.2.2\\tupian\\w1\\\2.pdf",plot=plot2)

scesce <-subset(sce,subset=nFeature_RNA>200 & nFeature_RNA< 7500& percent.mt <10) sce <-NormalizeData(sce,normalization.method= 'LogNormalize',scale.factor= 10000)

sce <-FindVariableFeatures(sce,selection.method ="vst",nfeatures =2000) top10 <-head(VariableFeatures(sce),10) plot3<-VariableFeaturePlot(sce) plot4 <-LabelPoints(plot= plot3,points=top10,repel=TRUE) var_lab= plot3 + plot4 ggsave(filename="D:\\R-4.2.2\\tupian\\w1\\3.pdf",plot=plot3) ggsave(filename="D:\\R-4.2.2\\tupian\\w1\\4.pdf",plot=plot4) ggsave(filename="D:\\R-4.2.2\\tupian\\w1\\ var_lab.pdf",plot= var_lab)

all.genes <-rownames(sce) sce <-ScaleData(sce,features =all.genes) sce <-RunPCA(sce, features = VariableFeatures(object = sce)) print(sce[["pca"]],dims= 1:5,nfeatures= 5)

VizDimLoadings(sce,dims= 1:2,reduction= "pca") DimPlot(sce,reduction= "pca") ggsave(filename="D:\\R-4.2.2\\tupian\\w1\\dim.pdf",plot=DimPlot)

ggsave(filename="D:\\R-4.2.2\\tupian\\w1\\dim.pdf",DimPlot(sce,reduction= "pca")) DimHeatmap(sce,dims=1,cells=500,balanced=TRUE)#1个PC500个细胞

ggsave(filename="D:\\R-4.2.2\\tupian\\w1\\dimheat.pdf",DimHeatmap(sce, dims = 1:2, cells = 500, balanced = TRUE))

DimHeatmap(sce,dims=1:15,cells=500,balanced=TRUE)

sce <- JackStraw(sce, num.replicate = 100) sce <-ScoreJackStraw(sce,dims =1:20) JackStrawPlot(sce,dims =1:20) ggsave(filename="/data3/hxl/A445-single-cell/jieguo1/JSplot.pdf",JackStrawPlot(sce,dims =1:20)) ElbowPlot(sce) ggsave(filename="/data3/hxl/A445-single-cell/jieguo1/elbowplot.pdf",ElbowPlot(sce))

sce<-FindNeighbors(sce,dims=1:20)#选取前20个主成分来分类细胞 sce<-FindClusters(sce,resolution =0.6) #查看前5个细胞的分类ID head(Idents(sce), 5) table(sce@active.ident) sce <-RunUMAP(sce,dims= 1:20)

DimPlot(sce,reduction ='umap')

#显示在聚类标签 DimPlot(sce,reduction= 'umap',label= TRUE) ggsave(filename="/data3/hxl/A445-single-cell/jieguo1/umap1.pdf",DimPlot(sce,reduction= 'umap',label= TRUE))

cluster1.markers <-FindMarkers(sce,ident.1 =1,min.pct= 0.25) head(cluster1.markers,n =5) #找出区分c1 uster2与cluster0和cluster3的所有标记(结果) #MT-ATP6 1.633143e-28 1.771912 1 0.934 5.045595e-24 #MT-CO3 4.109575e-28 1.770980 1 0.917 1.269653e-23 #MT-CO1 2.824917e-27 1.692531 1 0.950 8.727580e-23 #MT-CO2 3.259550e-27 1.881263 1 0.934 1.007038e-22 #MT-CYB 5.083406e-27 1.733684 1 0.917 1.570518e-22

cluster5.markers <-FindMarkers(sce,ident.1=5,ident.2 = c(0,3),min.pct =0.25) head(cluster5.markers,n= 5) #找出每个cluster的标记与所有剩余的细胞相比较，只报告阳性细胞 sce.markers <-FindAllMarkers(sce,only.pos =TRUE,min.pct =0.25,logfc.threshold= 0.25) sce.markers %>% group_by(cluster)%>%top_n(n= 2,wt= avg_log2FC)

FeaturePlot(sce,features= c('NCF2','PLXDC2','SLC7A7'),reduction='umap')

FeaturePlot(sce,features= c('MOG','PIP1','AQP4'),reduction='umap')

ggsave(filename="/data3/hxl/A445-single-cell/jieguo1/fplot1.pdf",FeaturePlot(sce,features= c('CK19','CA125','HBME-1'),reduction=umap'))

new.cluster.ids <-c("Naive CD4 T","CD14+Mono","Memory CD4 T","B","CD8 T","FCGR3A+Mono", "NK","DC","Platelet") names(new.cluster.ids)<-levels(sce) sce <-Renameldents(pbmc,new.cluster.ids) DimPlot(sce,reduction ="umap",label= TRUE,pt.size= 0.5)+NoLegend() 保存分析的结果