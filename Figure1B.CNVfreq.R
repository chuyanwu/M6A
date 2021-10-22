######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

inputFile="cnvMatrix.txt"     #输入文件
setwd("D:\\biowolf\\m6aTME\\12.CNVfreq")      #设置工作目录
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
GAIN=rowSums(rt> 0)       #拷贝数增加的样品数目
LOSS=rowSums(rt< 0)       #拷贝数缺失的样品数目
GAIN=GAIN/ncol(rt)*100      #拷贝数增加的百分率
LOSS=LOSS/ncol(rt)*100      #拷贝数缺失的百分率
data=cbind(GAIN, LOSS)
data=data[order(data[,"GAIN"],decreasing = T),]

#绘制图形
data.max = apply(data, 1, max)
pdf(file="CNVfreq.pdf", width=9, height=6)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(data.max, col="grey80", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(data.max)))
points(bar,data[,"GAIN"], pch=20, col=2, cex=3)
points(bar,data[,"LOSS"], pch=20, col=3, cex=3)
legend("top", legend=c('GAIN','LOSS'), col=2:3, pch=20, bty="n", cex=2, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(data), adj=1)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
