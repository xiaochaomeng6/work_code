```bash
library(dplyr)
library(Seurat)

rhesus.data <- Read10X(data.dir = "E:\\data1\\dalunwen\\single_cell\\filtered_feature_bc_matrix") ##加载cellranger cout分析后的数据        
rhesus.data <- Read10X(data.dir = "/bio/hxl/06/cell/SRR18042389/outs/filtered_feature_bc_matrix/")
##创建对象，min.feature为基因；counts为rawdata（10x genomic）或者TPM
rhesus <- CreateSeuratObject(counts = rhesus.data, project = "rhesus_singlecell", min.cells = 3, min.features = 200)  
##如果报错则options(Seurat.object.assay.version = 'v3')
head(rhesus)
pig

pig[["percent.mt"]] <- PercentageFeatureSet(pig, pattern = "^MT-") ##添加线粒体基因表达属性

png("vlnplot.png",width=1000,height=600)

##绘制单个细胞的基因表达个数、基因表达量以及线粒体基因百分比之间的关系，用于筛选过滤细胞（死细胞等离群细胞）
VlnPlot(pig, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) ##nFeature_RNA为单个细胞的基因表达个数，nCount_RNA基因表达量，percent.mt为线粒体基因百分比；
dev.off()

##绘制nCount_RNA、percent.mt以及nFeature_RNA的相关性散点图，用于筛选过滤细胞
png("boxplot.png",width=1000,height=600)
plot1 <- FeatureScatter(pig, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pig, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

##根据上面分析结果确定筛选条件
pig <- subset(pig, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & percent.mt < 5) ##过滤：基因数目大于200，不大于2500，MT-percent

##取Log进行数据标准化处理，
pig <- NormalizeData(pig, normalization.method = "LogNormalize", scale.factor = 10000)

##寻找变异系数较大的基因进行分析
pig <- FindVariableFeatures(pig, selection.method = "vst", nfeatures = 1500)

##可视化变异系数较大的细胞
png("feature.png",width=1000,height=600)
par(mfrow=c(1,2))
top10 <- head(VariableFeatures(pig), 10)
plot1 <- VariableFeaturePlot(pig)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

##
all.genes <- rownames(pig)
pig <- ScaleData(pig, features = all.genes)

pig <- RunPCA(pig, features = VariableFeatures(object = pig))

print(pig[["pca"]], dims = 1:5, nfeatures = 5)

png("dim.png")
VizDimLoadings(pig, dims = 1:2, reduction = "pca")
DimPlot(pig, reduction = "pca")
dev.off()


png("dimheatmap.png")
DimHeatmap(pig, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()


pig <- JackStraw(pig, num.replicate = 100)
pig <- ScoreJackStraw(pig, dims = 1:20)
png("pc_jack.png")
JackStrawPlot(pig, dims = 1:20)
dev.off()
png("elbow.png")
ElbowPlot(pig)
dev.off()



pig <- FindNeighbors(pig, dims = 1:15)
pig <- FindClusters(pig, resolution = 0.4)

pig <- RunUMAP(pig, dims = 1:15)

png("umap.png")
DimPlot(pig, reduction = "umap")
dev.off()
saveRDS(pig, file = "./pig_tutorial.rds")

cluster1.markers <- FindMarkers(pig, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster5.markers <- FindMarkers(pig, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

pig.markers <- FindAllMarkers(pig, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pig.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.table(file='tsne.txt' ,Embeddings(pig@reductions$umap),sep="\t",col.names=F,quote=F)
write.table(file='tsne.class' ,pig@active.ident,sep="\t",col.names=F,quote=F)
write.table(pig.markers,sep="\t",file="markers.xlsx",row.names=F)

cluster1.markers <- FindMarkers(pig, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

png("de1.png")
VlnPlot(pig, features = c("DES", "HSPB7"))
dev.off()
png("de2.png")
VlnPlot(pig, features = c("UBE2C", "CKS2"), slot = "counts", log = TRUE)
dev.off()

#FeaturePlot(pig, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
png("gene.png")
FeaturePlot(pig, features = c("ACTA2", "TAGLN", "THBS1", "TPM2", "DES", "CRYAB", "TAGLN", "MYL9","CNN1"))
dev.off()

png("heat.png")
top10 <- pig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pig, features = top10$gene) + NoLegend()
dev.off()
```