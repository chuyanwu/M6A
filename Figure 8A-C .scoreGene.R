
#install.packages("ggpubr")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


#引用包
library(limma)
library(ggpubr)
gene="CD274"            #基因的标准名字
showName="PD-L1"        #图形里面显示的基因名称
expFile="merge.txt"                #表达数据文件
scoreFile="m6Ascore.group.txt"     #m6A打分分组文件
setwd("D:\\biowolf\\m6aTME\\48.scoreGene")      #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
colnames(data)=gsub("(.*?)\\_(.*?)", "\\2", colnames(data))

#提取目标基因表达量
data=rbind(data, gene=data[gene,])
exp=t(data[c("gene",gene),])
exp=avereps(exp)

#读取m6A打分分组文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
	
#合并数据
sameSample=intersect(row.names(exp), row.names(score))
exp=exp[sameSample,]
exp[exp>quantile(exp,0.975)]=quantile(exp,0.975)
score=score[sameSample,]
data=cbind(as.data.frame(exp), as.data.frame(score))
	
#设置比较组
data$group=factor(data$group, levels=c("Low", "High"))
group=levels(factor(data$group))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
#绘制boxplot
boxplot=ggboxplot(data, x="group", y="gene", fill="group",
			      xlab="m6Ascore",
			      ylab=paste(showName, "expression"),
			      legend.title="m6Ascore",
			      palette=c("#0066FF","#FF0000"))+ 
	stat_compare_means(comparisons = my_comparisons)
	
#输出图片
pdf(file=paste0(gene, ".pdf"), width=5, height=4.5)
print(boxplot)
dev.off()
