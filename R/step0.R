make.rdafiles.from.textfiles <- function(settings) {
    # load("../../settings.rda")
    annoGenes=read.delim("data/gene.anno.hg38.txt.bed", header=T)
    annoHM450=read.delim("data/hm450cg_GEO.GR.bed", header=F)
    MT=list.files("./methylation/", "case")
    mT=read.delim(paste("./methylation/", MT[1], sep=""), header=T, check.names=F)
    metDataT=mT[match(as.character(annoHM450$V4), as.character(rownames(mT))),]
    rownames(metDataT)=as.character(annoHM450$V4)
    MN=list.files("./methylation/", "ctrl")
    mN=read.delim(paste("./methylation/", MN[1], sep=""), header=T, check.names=F)
    metDataN=mN[match(as.character(annoHM450$V4), as.character(rownames(mN))),]
    rownames(metDataN)=as.character(annoHM450$V4)
    ET=list.files("./expression/", "case")
    eT=read.delim(paste("./expression/", ET[1], sep=""), header=T, check.names=F)
    expDataT=eT[match(as.character(annoGenes$geneID), as.character(rownames(eT))),]
    rownames(expDataT)=as.character(annoGenes$geneID)
    EN=list.files("./expression/", "ctrl")
    eN=read.delim(paste("./expression/", EN[1], sep=""), header=T, check.names=F)
    expDataN=eN[match(as.character(annoGenes$geneID), as.character(rownames(eN))),]
    rownames(expDataN)=as.character(annoGenes$geneID)
    if(length(setdiff(colnames(metDataT), colnames(expDataT)))>0){
        print("sample IDs of case datasets are not matched")}
    if(length(setdiff(colnames(metDataN), colnames(expDataN)))>0){
        print("sample IDs of ctrl datasets are not matched")}
    metDataT=as.matrix(metDataT)
    metDataN=as.matrix(metDataN)
    expDataT=as.matrix(expDataT)
    expDataN=as.matrix(expDataN)
    save(metDataT,metDataN,expDataT,expDataN, file=paste(prefix, ".methylation.expression.dataset.rda", sep=""))
}
