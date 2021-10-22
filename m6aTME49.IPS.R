######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("ggpubr")


library(ggpubr)                    #引用包
tciaFile="TCIA.txt"                #免疫治疗打分文件
scoreFile="m6Ascore.group.txt"     #m6A打分分组文件
setwd("D:\\biowolf\\m6aTME\\49.IPS")     #修改工作目录

#读取免疫治疗打分文件
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

#读取m6A打分分组文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(ips), row.names(score))
ips=ips[sameSample, , drop=F]
score=score[sameSample, "group", drop=F]
data=cbind(ips, score)

#设置比较组
data$group=factor(data$group, levels=c("Low", "High"))
group=levels(factor(data$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#对免疫治疗打分进行循环,分别绘制小提琴图
for(i in colnames(data)[1:(ncol(data)-1)]){
	rt=data[,c(i, "group")]
	colnames(rt)=c("IPS", "group")
	gg1=ggviolin(rt, x="group", y="IPS", fill = "group", 
	         xlab="m6Ascore", ylab=i,
	         legend.title="m6Ascore",
	         palette=c("#0066FF", "#FF0000"),
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons)
	         #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	
	pdf(file=paste0(i, ".pdf"), width=6, height=5)
	print(gg1)
	dev.off()
}


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
