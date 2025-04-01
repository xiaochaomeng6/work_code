rm(list = ls())
library("gplots")
library("limma")
library("edgeR")
library("pheatmap")
library("ggplot2")


setwd("D:/desktop/大论文/第四章数据/表达量数据/gene")

# 读取数据并保留所有列
expr <- read.table("D:/desktop/大论文/第四章数据/表达量数据/gene/combined_TPM.txt", sep="\t", header=TRUE, check.names=FALSE)

# 检查数据维度和前几行，确保所有列都被读取
print(dim(expr)) # 打印数据维度
print(head(expr)) # 打印前几行数据，检查列名

# 保留原始行名（假设第一列是基因名）
rownames(expr) <- expr[,1]

# 保留所有数据列
exp <- expr[, -1]

# 检查表达数据维度
print(dim(exp)) # 打印表达数据维度

# 将基因表达数据转换成数值矩阵
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)

# 检查处理后的数据维度
print(dim(data)) # 打印处理后的数据维度

# 有的基因出现过多行，把出现多行的gene取平均值
data <- avereps(data)

# 去除低表达的数据
data <- data[rowMeans(data) > 1, ]

# 检查处理后的数据维度
print(dim(data)) # 打印处理后的数据维度

# 设置分组
group <- c(rep("si-171", 3),rep("si-331", 3), rep("nc", 3))

# 检查group的长度
print(length(group)) # 打印group的长度

# 检查data的列数
print(ncol(data)) # 打印data的列数

if (length(group) != ncol(data)) {
  stop("The length of 'group' must equal the number of columns in 'counts'")
}

design <- model.matrix(~group)  # 把group设置成一个model matrix
data <- data[complete.cases(data), ]
# 创建DGEList对象
y <- DGEList(counts=data, group=group)

# 对因子矫正
y <- calcNormFactors(y)

# 估计变异系数
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

# 进行差异表达分析
#et <- exactTest(y, pair=c("nc","si-171"))
et <- exactTest(y, pair=c("nc","si-331"))
topTags(et)

# 获取排名前十万的基因
ordered_tags <- topTags(et, n=100000)
allDiff <- ordered_tags$table

# 去除NA值
allDiff <- allDiff[is.na(allDiff$FDR) == FALSE, ]

diff <- allDiff

newData <- y$pseudo.counts

# 自定义参数
padj <- 0.05
foldChange <- 1

# 筛选出差异表达的基因
diff$abs_logFC <- abs(diff$logFC)  # 计算绝对值
diff <- diff[diff$FDR < padj & diff$abs_logFC > log2(foldChange), ]  # 筛选满足条件的基因

# 提取筛选出的基因名称
selected_genes <- rownames(diff)

# 筛选出前15个基因
selected_genes <- head(rownames(diff), 15)

# 提取这些基因在原始数据中的表达值
selected_data <- newData[selected_genes, ]

# 转换为矩阵
selected_matrix <- as.matrix(selected_data)

# 对每个值取 log2(TPM + 1)
selected_matrix_log2 <- log2(selected_matrix + 1)
# 将 selected_matrix_log2 输出为制表符分割的 TXT 文件
#write.table(selected_matrix_log2, 
#            file = "D:/desktop/大论文/第三章数据/表达量数据/gene/selected_matrix_log2.txt", 
#            sep = "\t", 
#            row.names = TRUE, 
#            col.names = TRUE)

# 定义文件路径
pdf_file_path <- "D:/desktop/大论文/第四章数据/表达量数据/gene/331前十五差异基因热图.pdf"

# 打开 PDF 图形设备
pdf(pdf_file_path)

# 绘制热图
pheatmap(selected_matrix_log2, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         show_rownames = TRUE,   # 显示基因名称
         show_colnames = TRUE,   # 显示样本名称
         angle_col = 0,          # 旋转列标签
         cluster_rows = FALSE,   # 禁用行聚类，去除左侧树状图
         cluster_cols = FALSE)   # 禁用列聚类，去除顶部树状图

# 关闭图形设备
dev.off()





