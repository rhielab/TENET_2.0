load("../settings.rda")
library(parallel)
GL=list.files("../external.data/otherinfo", "probe")
probelist=read.delim(paste("../external.data/otherinfo/",GL[1],sep=""),header=F)
links=paste(as.character(strsplit(GL[1], "\\.")[[1]][1]), as.character(strsplit(GL[1], "\\.")[[1]][2]), sep=".")
top=list(hm450probe=as.character(probelist$V1), links=as.character(links), prefix=as.character(prefix))
load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
library(ggplot2)
DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
metDataF_subC=cbind(metDataT, metDataN)
expDataF_subC=cbind(expDataT, expDataN)
makeScatterP=function(hm450probe=c(), links=c(), prefix=c()){
if (links=="hyper.G+"){
TESTSR=read.delim(paste("../step4/hyper.G+.output/", "hyper.G+.link.zscore.perm.all.optimized.links.txt", sep=""), header=T)
TESTSR=TESTSR[which(TESTSR$probe==hm450probe),]
dir.create("hyper.G+.output.scatterplot")
setwd("./hyper.G+.output.scatterplot")
for (i in 1:length(TESTSR[,1])){
EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=16)+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "scatterplots.pdf", sep="."),  useDingbats=FALSE)
print(i)
}
}
if (links=="hyper.G-"){
TESTSR=read.delim(paste("../step4/hyper.G-.output/", "hyper.G-.link.zscore.perm.all.optimized.links.txt", sep=""), header=T)
TESTSR=TESTSR[which(TESTSR$probe==hm450probe),]
dir.create("hyper.G-.output.scatterplot")
setwd("./hyper.G-.output.scatterplot")
for (i in 1:length(TESTSR[,1])){
EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=16)+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "scatterplots.pdf", sep="."),  useDingbats=FALSE)
print(i)
}
}
if (links=="hypo.G+"){
TESTSR=read.delim(paste("../step4/hypo.G+.output/", "hypo.G+.link.zscore.perm.all.optimized.links.txt", sep=""), header=T)
TESTSR=TESTSR[which(TESTSR$probe==hm450probe),]
dir.create("hypo.G+.output.scatterplot")
setwd("./hypo.G+.output.scatterplot")
for (i in 1:length(TESTSR[,1])){
EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=16)+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "scatterplots.pdf", sep="."),  useDingbats=FALSE)
print(i)
}
}
if (links=="hypo.G-"){
TESTSR=read.delim(paste("../step4/hypo.G-.output/", "hypo.G-.link.zscore.perm.all.optimized.links.txt", sep=""), header=T)
TESTSR=TESTSR[which(TESTSR$probe==hm450probe),]
dir.create("hypo.G-.output.scatterplot")
setwd("./hypo.G-.output.scatterplot")
for (i in 1:length(TESTSR[,1])){
EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=16)+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "scatterplots.pdf", sep="."),  useDingbats=FALSE)
print(i)
}
}
}
mcmapply(makeScatterP, top$hm450probe, top$links, top$prefix, mc.cores=cores)
