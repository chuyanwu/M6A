######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("ggpubr")


#引用包
library(ggpubr)
library(reshape2)
tmbFile="TMB.txt"     #肿瘤突变负荷文件
scoreFile="m6Ascore.group.txt"     #m6A打分的分组文件
cluFile="geneCluster.txt"          #基因分型文件
setwd("D:\\biowolf\\m6aTME\\42.scoreTMB")       #修改工作目录

#读取输入文件
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)        #读取TMB数据文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)    #读取m6A打分的分组文件
clu=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)        #读取基因分型类文件

#合并数据
tmb=as.matrix(tmb)
tmb[tmb>quantile(tmb,0.975)]=quantile(tmb,0.975)
sameSample=intersect(row.names(tmb), row.names(score))
tmb=tmb[sameSample,,drop=F]
score=score[sameSample,,drop=F]
rownames(clu)=gsub("(.*?)\\_(.*?)", "\\2", rownames(clu))
clu=clu[sameSample,,drop=F]
data=cbind(score, tmb, clu)
data=data[,c("m6Ascore", "group", "geneCluster", "TMB")]

#设置比较组
data$group=factor(data$group, levels=c("Low", "High"))
group=levels(factor(data$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#设置颜色
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]

#绘制箱线图
boxplot=ggboxplot(data, x="group", y="TMB", fill="group",
		          xlab="",
		          ylab="Tumor Burden Mutation",
		          legend.title="m6aScore",
		          palette = bioCol )+ 
	    stat_compare_means(comparisons = my_comparisons)
pdf(file="boxplot.pdf",width=5,height=4.5)
print(boxplot)
dev.off()

#相关性图形
length=length(levels(factor(data$geneCluster)))
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
p1=ggplot(data, aes(m6Ascore, TMB)) + 
		  xlab("m6Ascore")+ylab("Tumor Burden Mutation")+
		  geom_point(aes(colour=geneCluster))+
		  scale_color_manual(values=bioCol[1:length])+ 
		  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =m6Ascore, y =TMB))
#相关性图形
pdf(file="cor.pdf", width=6, height=4.5)
print(p1)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
