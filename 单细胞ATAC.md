使用cellranger-atac (v1.2.0)对scATAC-seq测序数据进行预处理。除“--force-cells”外，运行参数均使用默认。其中肝脏、肺脏和结肠的“--force-cell”设置为10000，脾脏设置为8000，其余器官使用默认。后续scATAC-seq数据分析使用ArchR (v1.0.1) [33](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9279386/#CR33)。具体而言，分别使用createGenomeAnnotation和createGeneAnnotation函数构建和注释_M.fascicularis_基因组。然后使用createArrowFiles函数使用默认参数创建arrow文件。使用addDoubletScores函数推断潜在的双联体，使用filterDoublets函数去除潜在双联体，参数为“filterRatio = 1.0”。使用ArchRProject函数使用默认参数创建ArchR项目。为了降维，我们使用了 ArchR 中的 addIterativeLSI 函数，其参数如下：“iterations = 4、clusterParams = list（resolution = c（0.2、0.4、0.6）、sampleCells = 10,000、n.start = 10、maxClusters = 6）、varFeatures = 20,000、dimsToUse = 1:50、scaleDims = FALSE”。接下来，利用 Harmony 方法通过 addHarmony 函数[32](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9279386/#CR32)消除批次效应。AddClusters 函数通过其默认参数对细胞进行聚类。对于单细胞嵌入，我们选择了具有 harmony 的 ReducedDims 对象，并使用带有参数“perplexity = 30”的 addTSNE 函数进行可视化


curl -o cellranger-atac-2.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-atac/cellranger-atac-2.1.0.tar.gz?Expires=1723035804&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hdGFjL2NlbGxyYW5nZXItYXRhYy0yLjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE3MjMwMzU4MDR9fX1dfQ__&Signature=RRac9wQDU6u0yCKVRmYdeLiOINmlEgNpPzEscfRalHbZy6l4G7M0E-BZcxg0QuhrxisRGxHDL07eDVtTlzKg0W93Qkr6KJHXoNYdq2tVuOPqat2GmAF8F2gT~8RgA9EVIyXtUfG1pQrEDfY2M-SHZ829NHKAz2c6evwHfeauNtTJKkHNReDg09G0lrOL5yjMyNlWjzYLsSG75ADLsbotZZpfvhfDKa5yno4bZP9AqP78zcm-eKATKw7kgDl2DEtk25Zpp2iDRGbT04kvt84GZpu5PmdmVwdQsBTiKgVcG2uXEIBhF~fI9jwVTa6NWkzzHYb4hlLutEqLoN8gXZB5-A__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"



{
    organism: "macaque"        #物种明后才能
    genome: ["06macaque"]       #输出文件夹名称
    input_fasta: ["/bio/hxl/06/06-macFas5.fasta"]
    input_gtf: ["/bio/hxl/06/cell/06_all_sorted_0gai1.gtf"]
    
}


cellranger-atac count --id=06SRR18042389 --reference=/bio/hxl/06/06macaque --fastqs=/bio/hxl/06/gai/ --sample=SRR18042389 --localcores=8 --localmem=64




nohup /bio/hxl/cellranger-atac-2.1.0/bin/cellranger-atac count --id=06SRR18042389 --reference=/bio/hxl/06/06macaque --fastqs=/bio/hxl/06/gai/ --sample=SRR18042389 --localcores=60 --localmem=120 >/bio/hxl/06/123.log 2>&1 &