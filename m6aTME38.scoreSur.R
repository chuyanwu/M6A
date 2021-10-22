######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)
scoreFile="m6Ascore.txt"     #m6A打分文件
cliFile="time.txt"           #生存数据文件
setwd("D:\\biowolf\\m6aTME\\38.scoreSur")      #设置工作目录

#读取输入文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
sampleType=gsub("(.*?)\\_.*", "\\1", row.names(score))
score=cbind(score, sampleType)
rownames(score)=gsub("(.*?)\\_(.*?)", "\\2", rownames(score))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(score), row.names(cli))
data=cbind(cli[sameSample,], score[sameSample,])

#获取最优cutoff
res.cut=surv_cutpoint(data, time="futime", event="fustat", variables=c("m6Ascore"))
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(data[,"m6Ascore"]<=cutoff, "Low", "High")
data$group=Type
outTab=rbind(id=colnames(data), data)
write.table(outTab, file="m6Ascore.group.txt", sep="\t", quote=F, col.names=F)

#计算高低风险组生存差异
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data)
#print(surv_median(fit))
	
#绘制生存曲线
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
			       data=data,
			       conf.int=F,
			       pval=pValue,
			       pval.size=6,
			       legend.title="m6Ascore",
			       legend.labs=levels(factor(data[,"group"])),
			       legend = c(0.8, 0.8),
			       font.legend=12,
			       xlab="Time(years)",
			       break.time.by = 1,
			       palette = bioCol,
			       surv.median.line = "hv",
			       risk.table=T,
			       cumevents=F,
			       risk.table.height=.25)

#保存图片
pdf(file="survival.pdf", onefile = FALSE, width=7, height=5.5)
print(surPlot)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
