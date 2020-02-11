load("../settings.rda")
#################
library(parallel)
load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
expDataT_subC=expDataT
##### for hypo.G- links ####
links=("hypo.G-")
LS=list.files("../step4/hypo.G-.output/", pattern="all.optimized.links.txt")
x=read.delim(paste("../step4/hypo.G-.output/", LS[1], sep=""), header=T)
metDataT_subC=ifelse(metDataT<hypocutoff,1,0)
datM=t(metDataT_subC[match(x$probe, rownames(metDataT_subC)),])
datE=t(expDataT_subC[match(x$geneID, rownames(expDataT_subC)),])
dir.create("hypo.G-.output.states")
setwd("./hypo.G-.output.states")
if(is.na(match("mean.expT", colnames(x)))==TRUE){
getmeanexpT=function(geneID){
mean(expDataT[as.character(geneID),], na.rm=T)}
top=list(geneID=as.character(x$geneID))
x$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
}
for (i in 1:dim(x)[1]){
datc=data.frame(meth=datM[,i], exp=datE[,i])
datc$exp.changed=ifelse(datc$exp<x$mean.expT[i],1,0)
rownames(datc)=rownames(datM)
datc$FIN=ifelse(datc$meth==1 & datc$exp.changed==1,1,0)
colnames(datc)=c(as.character(x$probe[i]), as.character(x$geneID[i]), "exp.changed","FIN")
write.table(datc, paste(x$probe[i], x$geneID[i], x$geneSymbol[i], prefix, links, "states.txt", sep="."), row.names=T, col.names=T, sep="\t", quote=F)
print(i)
}
TUMOR_Call=matrix("NA", nrow=dim(metDataT)[2], ncol=dim(x)[1])
rownames(TUMOR_Call)=rownames(datc)
colnames(TUMOR_Call)=paste(x$probe, x$geneID, x$geneSymbol,sep=".")
ID=colnames(TUMOR_Call)
findState=function(ID){
LS=list.files(getwd(), pattern=as.character(ID))
LF=read.table(LS, header=T)
LF$FIN
}
TUMOR_Call=mcmapply(ID=ID,findState,mc.cores=cores)
rownames(TUMOR_Call)=rownames(datc)
write.table(TUMOR_Call, paste(prefix, links, "links.states.table.txt", sep="."), row.names=T, col.names=T, sep="\t", quote=F)
# in order to remove NA - in vim change NA to 0 :%s/^INA/^I0/g
