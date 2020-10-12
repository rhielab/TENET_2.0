load("../settings.rda")
##### let's obtain methylation and expression data if step1 was not performed #####
load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
library(ggplot2)
DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
metDataF_subC=cbind(metDataT, metDataN)
expDataF_subC=cbind(expDataT, expDataN)
LS=list.files("../step4/hyper.G-.output/", pattern="all.optimized.links.txt")
TESTSR=read.delim(paste("../step4/hyper.G-.output/", LS[1], sep=""), header=T)
LC=list.files("../external.data/otherinfo", "CNV")
CNV=read.delim(paste("../external.data/otherinfo/", LC[1], sep=""), sep="\t", header=T, check.names=F)
dir.create("hyper.G-.output.complex.scatterplot")
setwd("./hyper.G-.output.complex.scatterplot")
M1=match(unique(TESTSR$geneSymbol),rownames(CNV))
CNV.c=CNV[na.omit(M1),]
M3=match(colnames(metDataT), colnames(CNV.c))
CNV.c2=as.matrix(CNV.c)[,M3]
for (i in 1:length(TESTSR[,1])){
EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
CNV.sh2=rep(16,dim(metDataF_subC)[2])
if (is.na(match(TESTSR[i,2],rownames(CNV.c2)))==FALSE){
CNV.t=CNV.c2[as.character(TESTSR[i,2]),]
CNV.t2=t(CNV.t)
CNV.t3=c(CNV.t2, rep(0,dim(metDataN)[2]))
CNV.sh=CNV.t3
CNV.sh2=CNV.sh
CNV.sh2[which(CNV.sh==0)]=paste(16)
CNV.sh2[which(CNV.sh==-1)]=paste(17)
CNV.sh2[which(CNV.sh==-2)]=paste(24)
CNV.sh2[which(CNV.sh==1)]=paste(15)
CNV.sh2[which(CNV.sh==2)]=paste(22)
}
p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
print(i)
}
setwd("../")
