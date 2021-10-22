######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("corrplot")


#引用包
library(corrplot)
scoreFile="m6Ascore.txt"         #m6A打分文件
immFile="ssGSEA.result.txt"      #ssGSEA结果文件
setwd("D:\\biowolf\\m6aTME\\40.scoreImm")     #设置工作目录

#读取m6A打分文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

#读取免疫细胞文件
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=t(immune)

#数据合并
sameSample=intersect(row.names(score), row.names(immune))
data=cbind(score[sameSample,,drop=F], immune[sameSample,,drop=F])

#相关性矩阵
M=cor(data)
res1=cor.mtest(data, conf.level = 0.95)

#绘制相关性图形
pdf(file="cor.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.8, pch=T,
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
