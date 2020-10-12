load("../settings.rda")
# load diffmethylated enhancer region DNA meth data #
load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
`%ni%` <- Negate(`%in%`)
selected=read.delim(paste("../step3/hyper.output/hyper.G+.output/", "hyper.G+.link.zscore.perm.txt", sep=""), header=T)
if (dim(selected)[1]<1){
print("no links found from step3, so cannot run step4")}
selected$mean.expN=paste(NA)
selected$mean.expT=paste(NA)
selected$wilcox.expNcTc=paste(NA)
selected$hypermeth.tumor.length=paste(NA)
selected$meanHypermethT=paste(NA)
selected$hypermeth.lower.tumor.length=paste(NA)
selected$max.metTc=paste(NA)
selected$mean.expN.high.expT=paste(NA)
metDataTcat=ifelse(hypermethDataT<hypercutoff,1,0)
expDataT0=ifelse(expDataT==0,1,0)
expDataN0=ifelse(expDataN==0,1,0)
getmeanexpN=function(geneID){
mean(expDataN[as.character(geneID),], na.rm=T)}
getmeanexpT=function(geneID){
mean(expDataT[as.character(geneID),], na.rm=T)}
getWilcoxexpNcTc=function(probe, geneID){
metT=metDataTcat[as.character(probe),]
metTc=metT[which(metT==0)]
expN=expDataN[as.character(geneID),]
expN0=expDataN0[as.character(geneID),]
expN0=expN0[which(expN0==1)]
expNc=expN[which(names(expN)%ni%names(expN0))]
expT=expDataT[as.character(geneID),]
expT0=expDataT0[as.character(geneID),]
expT0=expT0[which(expT0==1)]
expTc=expT[match(names(metTc),names(expT))]
expTc=expTc[which(names(expTc)%ni%names(expT0))]
if (length(expNc)>0 & length(expTc)>0){
wilcoxPval=paste(wilcox.test(expNc, expTc)$p.value)
}
else{
wilcoxPval=paste(NA)
}
}
gethypermethTlength=function(probe){
sum(hypermethDataT[as.character(probe),]>hypercutoff, na.rm=T)}
getmeanhypermethT=function(geneID, probe){
metT=metDataTcat[as.character(probe),]
metTc=metT[which(metT==0)]
expT=expDataT[as.character(geneID),]
expT0=expDataT0[as.character(geneID),]
expT0=expT0[which(expT0==1)]
expTc=expT[match(names(metTc),names(expT))]
expTc=expTc[which(names(expTc)%ni%names(expT0))]
m=mean(expTc, na.rm=T)
print(m)
}
gethypermethlowerTlength=function(geneID, probe){
metT=metDataTcat[as.character(probe),]
metTc=metT[which(metT==0)]
expT=expDataT[as.character(geneID),]
expT0=expDataT0[as.character(geneID),]
expT0=expT0[which(expT0==1)]
expTc=expT[match(names(metTc),names(expT))]
expTc=expTc[which(names(expTc)%ni%names(expT0))]
mean.expT=mean(expT, na.rm=T)
m=sum(expTc<mean.expT, na.rm=T)
print(m)
}
getmaxmetTc=function(probe){
metT=metDataTcat[as.character(probe),]
metTc=metT[which(metT==0)]
m=paste(max(as.numeric(hypermethDataT[as.character(probe),match(names(metTc), colnames(hypermethDataT))]),na.rm=T))
print(m)
}
top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe)) 
library(parallel)
selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
selected$hypermeth.tumor.length=mcmapply(gethypermethTlength, top$probe, mc.cores=cores)
selected$meanHypermethT=mcmapply(getmeanhypermethT, top$geneID, top$probe, mc.cores=cores)
selected$hypermeth.lower.tumor.length=mcmapply(gethypermethlowerTlength, top$geneID, top$probe, mc.cores=cores)
selected$max.metTc=mcmapply(getmaxmetTc, top$probe, mc.cores=cores)
selected$mean.expN.high.expT=ifelse(as.numeric(selected$mean.expN)>as.numeric(selected$mean.expT),1,0)
selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.high.expT)==1),]
selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypermeth.tumor.length)>minTumor & as.numeric(selectedc$hypermeth.lower.tumor.length)>minTumor & as.numeric(selectedc$max.metTc)>hyper.stringency),]
dir.create("hyper.G+.output")
setwd("./hyper.G+.output")
write.table(selectedc2, file=paste("hyper.G+.link.zscore.perm.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
setwd("../")
