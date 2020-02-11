load("../../settings.rda")
library(parallel)
## let's do for hypo.G- #
PL=read.delim("../../step2/hypo.output/hypo.G-.link.zscore.txt", header=F)
colnames(PL)=c("probe", "geneSymbol", "geneID", "Z.real")
dir.create("hypo.G-.output")
setwd("./hypo.G-.output")
top=list(probe=as.character(PL$probe), geneID=as.character(PL$geneID), Z.real=as.numeric(PL$Z.real))
getNegEpval=function(probe, geneID, Z.real){
LS=list.files("../", pattern=paste("^", as.character(geneID), ".hypo", sep=""))
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
print("did you run steps correctly? no significant hypo.G- links are found from step3")
}
write.table(selected, file=paste("hypo.G-.link.zscore.perm.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
setwd("../")
