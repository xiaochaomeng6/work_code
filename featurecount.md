重复序列
#!/bin/bash

# 输入和输出目录
input_dir="/bio/GJY/oocyte-KO/Bam/GSE139966"
output_dir="/bio/GJY/oocyte-KO/featurecounts/GSE139966"
gtf_files=(
    "/data/genome/Mouse-GRCm39/Mouse-GRCm39-LTR_num.gtf"
    "/data/genome/Mouse-GRCm39/Mouse-GRCm39-LINE_num.gtf"
    "/data/genome/Mouse-GRCm39/Mouse-GRCm39-SINE_num.gtf"
)

# 设置线程数
threads=10

# 确保输出目录存在
mkdir -p ${output_dir}

# 循环处理所有排序后的 BAM 文件
for bamfile in ${input_dir}/*sorted*.bam; do
    # 提取文件名（去掉路径和扩展名）
    base=$(basename ${bamfile} .bam)
    
    # 对每个 GTF 文件进行计数
    for gtf in "${gtf_files[@]}"; do
        # 提取 GTF 文件的基名（去掉路径和扩展名）
        gtf_base=$(basename ${gtf} .gtf)
        
        # 定义输出文件名
        output=${output_dir}/${base}_${gtf_base}_featureCounts.txt
        
        # 运行 featureCounts
        echo "Counting ${bamfile} with ${gtf}..."
        /data/software/subread/subread-2.0.1-Linux-x86_64/bin/featureCounts \
            -s 0 -f -p -C -T ${threads} \
            -t exon -g gene_id \
            -a ${gtf} \
            -o ${output} \
            ${bamfile}
        
        if [ $? -ne 0 ]; then
            echo "Error processing ${bamfile} with ${gtf}"
            exit 1
        fi
    done
    
    echo "${bamfile} has been processed successfully."
done

echo "All files have been processed."









全基因组
#!/bin/bash

# 设置featureCounts工具的路径
featurecounts_path="/data/software/subread/subread-2.0.1-Linux-x86_64/bin/featureCounts"

# 设置GTF文件路径
gtf_file="/bio/GJY/GRCm39vM27.gtf"

# 设置输入和输出路径
input_path="/bio/GJY/oocyte-KO/Bam/GSE198414"
output_path="/bio/GJY/oocyte-KO/featurecounts/GSE198414/gene"

# 设置线程数
threads=35

# 匹配输入文件
input_files=$(ls ${input_path}/SRR*_sorted.bam)

# 遍历所有匹配的文件
for input_file in ${input_files}; do
    # 提取文件名中的SRR编号
    base=$(basename ${input_file} _sorted.bam)
    
    # 设置输出文件路径
    output_file="${output_path}/${base}_gene_featureCounts.txt"
    
    # 运行featureCounts
    ${featurecounts_path} -C -T -p ${threads} -t exon -g gene_name -a ${gtf_file} -o ${output_file} ${input_file}
    
    echo "Processed ${base}"
done

echo "All files processed."



#!/bin/bash

# 输入和输出目录
input_dir="/bio/hxl/GSE175839_srr/fa/sam"
output_dir="/bio/hxl/GSE175839_srr/fa/sam/featurecount/gene/"
gtf_files=(
    "/bio/hxl/01-human.gtf")

# 设置线程数
threads=35

# 确保输出目录存在
mkdir -p ${output_dir}

# 循环处理所有排序后的 BAM 文件
for bamfile in ${input_dir}/*sort*.bam; do
    # 提取文件名（去掉路径和扩展名）
    base=$(basename ${bamfile} .bam)
    
    # 对每个 GTF 文件进行计数
    for gtf in "${gtf_files[@]}"; do
        # 提取 GTF 文件的基名（去掉路径和扩展名）
        gtf_base=$(basename ${gtf} .gtf)
        
        # 定义输出文件名
        output=${output_dir}/${base}_${gtf_base}_featureCounts.txt
        
        # 运行 featureCounts
        echo "Counting ${bamfile} with ${gtf}..."
        /data/software/subread/subread-2.0.1-Linux-x86_64/bin/featureCounts \
            -p  -C  -T ${threads} \
            -t  exon  -g  gene_name \
            -a ${gtf} \
            -o ${output} \
            ${bamfile}
        
        if [ $? -ne 0 ]; then
            echo "Error processing ${bamfile} with ${gtf}"
            exit 1
        fi
    done
    
    echo "${bamfile} has been processed successfully."
done

echo "All files have been processed."



基因组
 /data/software/subread/subread-2.0.1-Linux-x86_64/bin/featureCounts  -p  -C  -T  35  -t  exon  -g  gene_name -a  gencode_v36_primary_assembly_annotation.gtf  -o  gene_featureCounts.txt  sort.bam

 /data/software/subread/subread-2.0.1-Linux-x86_64/bin/featureCounts  -p  -C  -T  60  -t  exon  -g  gene_name -a  /bio/hxl/01-human.gtf  -o  /bio/hxl/GSE98422_srr/fa/sam/featurecount/gene/SRR5494037.sam.bam.sort_01-human_featureCounts.txt  /bio/hxl/GSE98422_srr/fa/sam/SRR5494037.sort.bam





