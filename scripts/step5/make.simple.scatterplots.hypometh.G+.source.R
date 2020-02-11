load("../settings.rda")
links=c("hypo.G+")
##### let's obtain methylation and expression data if step1 was not performed #####
load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
library(ggplot2)
DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
metDataF_subC=cbind(metDataT, metDataN)
expDataF_subC=cbind(expDataT, expDataN)
if (links=="hypo.G+"){
LS=list.files("../step4/hypo.G+.output/", pattern="all.optimized.links.txt")
TESTSR=read.delim(paste("../step4/hypo.G+.output/", LS[1], sep=""), header=T)
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
setwd("../")
}
