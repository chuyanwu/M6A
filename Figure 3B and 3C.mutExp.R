

#install.packages("ggpubr")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


#引用包
library(limma)
library(ggpubr)
expFile="m6aGeneExp.txt"     #表达数据文件
mutFile="mutMatrix.txt"      #突变的矩阵文件
mutGene="ZC3H13"             #选择突变分组基因
setwd("D:\\biowolf\\m6aTME\\20.mutExp")      #设置工作目录

#读取表达数据文件
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(exp)=gsub("(.*?)\\_(.*?)", "\\2", colnames(exp))
exp=t(exp)

#读取突变数据文件
mut=read.table(mutFile, header=T, sep="\t", check.names=F, row.names=1)
mut=t(mut[mutGene,,drop=F])
colnames(mut)=c("Type")

#合并数据
sameSample=intersect(row.names(mut), row.names(exp))
mut=mut[sameSample,,drop=F]
exp=exp[sameSample,,drop=F]
data=cbind(as.data.frame(exp), as.data.frame(mut))
data$Type=paste0(mutGene, " " , data$Type)

#设置比较组
data$Type=factor(data$Type, levels=c(paste0(mutGene, " Wild"), paste0(mutGene, " Mutation")) )
group=levels(factor(data$Type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

for(gene in colnames(data)[1:(ncol(data)-1)]){
	data1=data[,c(gene, "Type")]
	colnames(data1)=c("expression", "Type")
	#绘制箱线图
	boxplot=ggboxplot(data1, x="Type", y="expression", fill="Type",
				      xlab="",
				      ylab=paste0(gene, " expression"),
				      legend.title="",
				      palette=c("#0066FF", "#FF0000") )+ 
		stat_compare_means(comparisons = my_comparisons)
	#输出图片
	pdf(file=paste0(mutGene, "(mut)_", gene,".pdf"), width=5, height=4.5)
	print(boxplot)
	dev.off()
}

