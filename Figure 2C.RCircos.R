######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("RCircos")


library("RCircos")       #引用包
setwd("D:\\biowolf\\m6aTME\\14.Rcircos")    #设置工作目录

#初始化圈图
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t")
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 5
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

#设置圈图参数
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=1
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)

#输出文件
pdf(file="RCircos.pdf", width=8, height=8)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

#散点图
RCircos.Scatter.Data=read.table("Rcircos.scatter.txt", header=T, sep="\t", check.names=F)
data.col <- 4
track.num <- 1
side <- "in"
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, side, by.fold=0.1)

#加上基因名称
RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
name.col <- 4
side <- "in"
track.num <- 2
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 3
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
