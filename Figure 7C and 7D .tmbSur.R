
#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)
tmbFile="TMB.txt"                  #肿瘤突变负荷文件
scoreFile="m6Ascore.group.txt"     #m6A打分的分组文件
setwd("D:\\biowolf\\m6aTME\\43.tmbSur")       #修改工作目录

#读取输入文件
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)    #读取m6A打分的分组文件
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)        #读取TMB数据文件

#合并数据
sameSample=intersect(row.names(tmb), row.names(score))
tmb=tmb[sameSample,,drop=F]
score=score[sameSample,,drop=F]
data=cbind(score, tmb)

#获取最优cutoff
res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
scoreType=ifelse(data$group=="Low", "L-m6Ascore", "H-m6Ascore")
mergeType=paste0(tmbType, "+", scoreType)

#生存曲线函数
bioSurvival=function(surData=null, outFile=null){
	diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
	length=length(levels(factor(surData[,"group"])))
	pValue=1-pchisq(diff$chisq, df=length-1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
	#print(surv_median(fit))
	
	#绘制生存曲线
	width=6.5
	height=5.5
	if(length(levels(factor(surData[,"group"])))>2){
		width=8
		height=6.5
	}
	bioCol=c("#FF0000","#0066FF","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
	bioCol=bioCol[1:length]
	surPlot=ggsurvplot(fit, 
			           data=surData,
			           conf.int=F,
			           pval=pValue,
			           pval.size=6,
			           legend.title="",
			           legend.labs=levels(factor(surData[,"group"])),
			           font.legend=10,
			           legend = c(0.8, 0.8),
			           xlab="Time(years)",
			           break.time.by = 1,
			           palette = bioCol,
			           surv.median.line = "hv",
			           risk.table=T,
			           cumevents=F,
			           risk.table.height=.25)
	#输出图形
	pdf(file=outFile, onefile = FALSE, width=width, height=height)
	print(surPlot)
	dev.off()
}

#绘制TMB的生存曲线
data$group=tmbType
bioSurvival(surData=data, outFile="TMB.survival.pdf")

#绘制TMB联合m6A打分的生存曲线
data$group=mergeType
bioSurvival(surData=data, outFile="TMB-score.survival.pdf")
