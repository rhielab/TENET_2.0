load("../../settings.rda")
ET=list.files("txt")
annoHM450=read.delim("../../scripts/data/hm450cg_GEO.GR.bed", header=F)
eD=read.delim(ET[1], header=T)
extData=eD[match(as.character(annoHM450$V4), as.character(rownames(eD))),]
rownames(extData)=as.character(annoHM450$V4)
extData=as.matrix(extData)
save(extData,file=paste("othercells.methylation.dataset.rda", sep=""))
