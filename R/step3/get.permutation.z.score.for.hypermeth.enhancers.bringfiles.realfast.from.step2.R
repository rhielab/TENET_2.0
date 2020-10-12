load("../../settings.rda")
library(parallel)
# let's do for hyper.G+ #
PL=read.delim("../../step2/hyper.output/hyper.G+.link.zscore.txt", header=F)
colnames(PL)=c("probe", "geneSymbol", "geneID", "Z.real")
dir.create("hyper.G+.output")
setwd("./hyper.G+.output")
top=list(probe=as.character(PL$probe), geneID=as.character(PL$geneID), Z.real=as.numeric(PL$Z.real))
getPosEpval=function(probe, geneID, Z.real){
LS=list.files("../", pattern=paste("^", as.character(geneID), ".hyper", sep=""))
TESTSR=read.delim(paste("../",LS[1], sep=""), header=F)
colnames(TESTSR)=c("probe", "geneSymbol", "geneID", "Z.real")
P=TESTSR[which(TESTSR$geneID==as.character(geneID)),]
P2=P[order(as.numeric(P$Z.real), decreasing=F),]
P2$rank=c(1:dim(TESTSR)[1])
epval=as.numeric(P2[which(P2$probe==as.character(probe)),5])/as.numeric(length(as.numeric(Z.real)[which(is.finite(as.numeric(P2$Z.real))==TRUE)]))
print(epval)
}
pvals=mcmapply(getPosEpval, top$probe, top$geneID, top$Z.real, mc.cores=cores)
PL$pval=pvals
selected=PL[which(PL$pval<permutation.cutoff),]
if (dim(selected)[1]<1){
print("did you run steps correctly? no significant hyper.G+ links are found from step3")
}
write.table(selected, file=paste("hyper.G+.link.zscore.perm.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
setwd("../")
## let's do for hyper.G- #
PL=read.delim("../../step2/hyper.output/hyper.G-.link.zscore.txt", header=F)
colnames(PL)=c("probe", "geneSymbol", "geneID", "Z.real")
dir.create("hyper.G-.output")
setwd("./hyper.G-.output")
top=list(probe=as.character(PL$probe), geneID=as.character(PL$geneID), Z.real=as.numeric(PL$Z.real))
getNegEpval=function(probe, geneID, Z.real){
LS=list.files("../", pattern=paste("^", as.character(geneID), ".hyper", sep=""))
TESTSR=read.delim(paste("../",LS[1], sep=""), header=F)
colnames(TESTSR)=c("probe", "geneSymbol", "geneID", "Z.real")
P=TESTSR[which(TESTSR$geneID==as.character(geneID)),]
P2=P[order(as.numeric(P$Z.real), decreasing=T),]
P2$rank=c(1:dim(TESTSR)[1])
epval=as.numeric(P2[which(P2$probe==as.character(probe)),5])/as.numeric(length(as.numeric(Z.real)[which(is.finite(as.numeric(P2$Z.real))==TRUE)]))
print(epval)
}
pvals=mcmapply(getNegEpval, top$probe, top$geneID, top$Z.real, mc.cores=cores)
PL$pval=pvals
selected=PL[which(PL$pval<permutation.cutoff),]
if (dim(selected)[1]<1){
print("did you run steps correctly? no significant hyper.G- links are found from step3")
}
write.table(selected, file=paste("hyper.G-.link.zscore.perm.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
setwd("../")
