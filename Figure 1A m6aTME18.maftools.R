######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)       #引用包
setwd("D:\\biowolf\\m6aTME\\18.maftools")      #设置工作目录

#读取突变基因文件
geneRT=read.table("gene.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneRT)

#绘制瀑布图
pdf(file="oncoplot.pdf", width=6.5, height=6)
maf=read.maf(maf="input.maf")
oncoplot(maf=maf, genes=gene, draw_titv=T)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
