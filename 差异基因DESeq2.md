rm(list = ls())
setwd('/data/nas1/gongjingyang_OD/project/01_Program094/01_DEGs/')
if(!dir.exists('./01_DEGs')){
  dir.create('./01_DEGs')
}

library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)


# 读取基因表达数据（raw counts）
expr_data <- read.csv("/data/nas1/gongjingyang_OD/project/01_Program094/00_rawdata/01.dat.GSE199939_count.csv", row.names = 1, check.names = FALSE)

# 读取样本分组信息
sample_info <- read.csv("/data/nas1/gongjingyang_OD/project/01_Program094/00_rawdata/01.group.GSE199939.csv", row.names = 1)

# 确保样本顺序一致
expr_data <- expr_data[, rownames(sample_info)]

# 检查数据格式
dim(expr_data)  # 行: 基因, 列: 样本
head(expr_data[, 1:5])  # 查看前几列
head(sample_info)  # 查看分组信息
table(sample_info$group)  # 确保组别正确

# 创建 DESeq2 数据对象
dds <- DESeqDataSetFromMatrix(
  countData = expr_data, 
  colData = sample_info, 
  design = ~ group
)

# 过滤低表达基因
keep <- rowSums(counts(dds) > 10) >= 10
dds <- dds[keep, ]

# 设置 "Control" 作为对照组
dds$group <- relevel(dds$group, ref = "Control")


# 运行差异分析
dds <- DESeq(dds)

# 获取差异分析结果
res <- results(dds)

# 查看差异分析结果摘要
summary(res)

# 按照 adj.P.Val 排序
res <- res[order(res$padj), ]

# 查看前几行
head(res)

# 过滤显著差异基因
DEGs <- res %>%
  as.data.frame() %>%
  filter(abs(log2FoldChange) > 2 & padj < 0.05)

# 查看前10个差异基因
head(DEGs, 10)

# 保存所有差异基因结果
write.csv(res, "DEGs_all_results.csv")

# 保存过滤后的显著差异基因
write.csv(DEGs, "DEGs_res_results.csv")

# 1. 标记显著差异基因
res$threshold <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, 
                        ifelse(res$log2FoldChange > 0, "Upregulated", "Downregulated"), "Stable")

# 2. 绘制火山图，标注上调、下调基因及稳定基因个数
# 按 |log2FoldChange| 降序排列并取前 10 个基因
top_DEGs <- res %>%
  as.data.frame() %>%
  filter(abs(log2FoldChange) > 1 & padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(10)

# 计算上调、下调和稳定基因的数量
upregulated <- sum(res$log2FoldChange > 1 & res$padj < 0.05)  # 上调基因数量
downregulated <- sum(res$log2FoldChange < -1 & res$padj < 0.05)  # 下调基因数量
stable <- sum(res$padj >= 0.05 | abs(res$log2FoldChange) <= 1)  # 稳定基因数量

# 绘制火山