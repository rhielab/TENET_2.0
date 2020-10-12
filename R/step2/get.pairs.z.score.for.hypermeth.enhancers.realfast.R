library(parallel)
load("../settings.rda")

# let's load diffmethylated enhancer region DNA meth data #
load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
## for hypermethylated probes ##
### let's make a table containing Z score ## 
if (usecaseonly==T){
expData=expDataT
metData=ifelse(hypermethDataT<hypercutoff,1,0)}
if (usecaseonly==F){
expData=cbind(expDataT, expDataN)
metData=ifelse(cbind(hypermethDataT, hypermethDataN)<hypercutoff,1,0)}

## Load the hg38 genes
geneanno <- read.delim("../scripts/data/gene.anno.hg38.txpt.bed", header=T)

## Truncate geneanno to genes contained in expression datasets:
rownames(geneanno) <- geneanno$geneID
geneanno <- geneanno[rownames(expData),]

dir.create("hyper.output")
setwd("./hyper.output")
getZscore=function(hm450probe=hm450probe, metData=metData, expData=expData, geneanno=geneanno){
TESTSR=data.frame(probe=hm450probe, geneSymbol=geneanno$geneName, geneID=geneanno$geneID)
TEST=TESTSR
TEST$Z.real=paste(NA)
P_probe=subset(metData, rownames(metData)==hm450probe)
T_P_probe=t(P_probe)
T_P_probe_1=subset(T_P_probe, T_P_probe[,1]=="1")
T_P_probe_0=subset(T_P_probe, T_P_probe[,1]=="0")
if(length(T_P_probe_1)>minTumor & length(T_P_probe_0)>0){
for (i in 1:length(TEST[,1])){
P_gene=subset(expData, rownames(expData)==TEST[i,3])
P_gene_1=P_gene[which(colnames(P_gene) %in% rownames(T_P_probe_1))]
P_gene_0=P_gene[which(colnames(P_gene) %in% rownames(T_P_probe_0))]
TEST[i,4]=paste((mean(P_gene_0, na.rm=T)-mean(P_gene_1, na.rm=T))/sd(P_gene_1, na.rm=T))
}
assign(hm450probe, TEST)
write.table(TEST, file=paste(hm450probe, "zscore.all.genes.rda.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
}
}
sel=rownames(metData)
mclapply(sel, getZscore, metData=metData, expData=expData, geneanno=geneanno, mc.cores=cores)
