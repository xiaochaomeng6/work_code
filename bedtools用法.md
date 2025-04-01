获取基因间基因的原理，用bedtools将LTR的gtf文件和基因组的gtf文件相连接得到A文件，然后再求A文件对于LTR的gtf文件的补集（即用LTR的gtf文件-A文件）


bedtools intersect -a 07-rhesusLTR.txt -b /zhlab/date/vert/gtf/07-rheMac10.gtf -wa -wb >07LTR.txt

b后面跟的是原本的gtf文件，即为基因组文件


a后面跟的是LTR,SINE,LINE等文件


a后面的文件是，从已经收集到的基因的名字的文件中，返回到LTR，line，sine的文件中提取到的位置的信息的文件，比如chr, start,end ,-+ ,基因名


/zhlab/date/vert/gtf/06-macFas5.gtf
/zhlab/date/vert/gtf/07-rheMac10.gtf



bedtools intersect -a /bio/F24A610000313_MACuvbdO_0621/sam/0607/gai/07-rhesusLTR_gtf_gai.txt -b /bio/F24A610000313_MACuvbdO_0621/sam/0607/gai/07_gtf_gai.txt -wa -wb >07LTR_gai.txt






