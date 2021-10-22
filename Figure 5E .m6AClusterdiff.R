#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)
expFile="m6aGeneExp.txt"          #表达输入文件
geneCluFile="geneCluster.txt"     #基因分型文件
setwd("D:\\biowolf\\m6aTME\\36.m6AClusterDiff")       #设置工作目录

#读取表达输入文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

#读取基因分型文件
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(data), row.names(geneClu))
expClu=cbind(data[sameSample,,drop=F], geneClu[sameSample,,drop=F])

#把数据转换成ggplot2输入文件
data=melt(expClu, id.vars=c("geneCluster"))
colnames(data)=c("geneCluster", "Gene", "Expression")

#设置颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"geneCluster"])))]

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color = "geneCluster", 
	     ylab="Gene expression",
	     xlab="",
	     legend.title="geneCluster",
	     palette = bioCol,
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=geneCluster),
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#输出箱线图
pdf(file="boxplot.pdf", width=7, height=5)
print(p1)
dev.off()

