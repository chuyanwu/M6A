######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)       #引用包
setwd("D:\\biowolf\\m6aTME\\45.maftools")      #设置工作目录

#读取m6A评分的分组文件
score=read.table("m6Ascore.group.txt", header=T, sep="\t", check.names=F)
outTab=score[,c(1, ncol(score))]
colnames(outTab)=c("Tumor_Sample_Barcode", "m6Ascore")
write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)

#读取基因突变文件
geneNum=20
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]

#颜色
ann_colors=list()
col=c("#0066FF","#FF0000")
names(col)=c("Low", "High")
ann_colors[["m6Ascore"]]=col

#低评分组瀑布图
pdf(file="low.pdf", width=6, height=6)
maf=read.maf(maf="low.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="m6Ascore", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

#高评分组瀑布图
pdf(file="high.pdf", width=6, height=6)
maf=read.maf(maf="high.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="m6Ascore", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
