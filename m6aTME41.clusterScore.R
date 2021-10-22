######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)
m6aCluFile="m6aCluster.txt"        #m6A分型文件
geneCluFile="geneCluster.txt"      #基因分型文件
scoreFile="m6Ascore.txt"           #m6A打分文件
setwd("D:\\biowolf\\m6aTME\\41.clusterScore")     #设置工作目录

#读取输入文件
m6aClu=read.table(m6aCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
twoCluster=cbind(m6aClu, geneClu)
sameSample=intersect(row.names(twoCluster), row.names(score))
data=cbind(score[sameSample,,drop=F], twoCluster[sameSample,,drop=F])

#######m6A分型与打分相关性########
#设置比较组
data$m6Acluster=factor(data$m6Acluster, levels=levels(factor(data$m6Acluster)))
group=levels(factor(data$m6Acluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#定义颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$m6Acluster)))]
	
#绘制boxplot
boxplot=ggboxplot(data, x="m6Acluster", y="m6Ascore", color="m6Acluster",
			      xlab="m6Acluster",
			      ylab="m6Ascore",
			      legend.title="m6Acluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
#输出图片
pdf(file="m6Acluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
#######m6A分型与打分相关性########


#######基因分型与打分相关性########
#设置比较组
data$geneCluster=factor(data$geneCluster, levels=levels(factor(data$geneCluster)))
group=levels(factor(data$geneCluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#定义颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$geneCluster)))]
	
#绘制boxplot
boxplot=ggboxplot(data, x="geneCluster", y="m6Ascore", color="geneCluster",
			      xlab="geneCluster",
			      ylab="m6Ascore",
			      legend.title="geneCluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
#输出图片
pdf(file="geneCluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
#######基因分型与打分相关性########


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
