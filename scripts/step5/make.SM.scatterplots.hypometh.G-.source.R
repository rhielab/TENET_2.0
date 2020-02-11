load("../settings.rda")
##### let's obtain methylation and expression data if step1 was not performed #####
load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
library(ggplot2)
DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
metDataF_subC=cbind(metDataT, metDataN)
expDataF_subC=cbind(expDataT, expDataN)
LS=list.files("../step4/hypo.G-.output/", pattern="all.optimized.links.txt")
TESTSR=read.delim(paste("../step4/hypo.G-.output/", LS[1], sep=""), header=T)
LM=list.files("../external.data/otherinfo", "SM")
SM=read.delim(paste("../external.data/otherinfo/", LM[1], sep=""), header=T, sep="\t")
dir.create("hypo.G-.output.complex.scatterplot")
setwd("./hypo.G-.output.complex.scatterplot")
SM=t(SM)
M2=match(unique(TESTSR$geneSymbol),rownames(SM))
SM.c=SM[na.omit(M2),]
M4=match(colnames(metDataT), colnames(SM.c))
SM.c2=as.matrix(SM.c)[,M4]
for (i in 1:length(TESTSR[,1])){
EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
CNV.sh2=rep(16,dim(metDataF_subC)[2])
SM.sh=rep(16,dim(metDataF_subC)[2])
if (is.na(match(TESTSR[i,2],rownames(SM.c2)))==FALSE){
SM.t=SM.c2[as.character(TESTSR[i,2]),]
SM.t2=t(SM.t)
SM.t3=c(SM.t2, rep(0,dim(metDataN)[2]))
SM.sh=SM.t3
SM.sh2=SM.sh
SM.sh2[which(SM.sh==0)]=paste(16)
SM.sh2[which(SM.sh==1)]=paste(8)
CNV.sh2[which(SM.sh2==8)]=paste(8)
}
p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
print(i)
}
setwd("../")
