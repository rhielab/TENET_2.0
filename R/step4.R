gethypermethhigherTlength=function(geneID, probe){
    metT=metDataTcat[as.character(probe),]
    metTc=metT[which(metT==0)]
    expT=expDataT[as.character(geneID),]
    expT0=expDataT0[as.character(geneID),]
    expT0=expT0[which(expT0==1)]
    expTc=expT[match(names(metTc),names(expT))]
    expTc=expTc[which(names(expTc)%ni%names(expT0))]
    mean.expT=mean(expT, na.rm=T)
    m=sum(expTc>mean.expT, na.rm=T)
    print(m)
}
gethypermethlowerTlength=function(geneID, probe){
    metT=metDataTcat[as.character(probe),]
    metTc=metT[which(metT==0)]
    expT=expDataT[as.character(geneID),]
    expT0=expDataT0[as.character(geneID),]
    expT0=expT0[which(expT0==1)]
    expTc=expT[match(names(metTc),names(expT))]
    expTc=expTc[which(names(expTc)%ni%names(expT0))]
    mean.expT=mean(expT, na.rm=T)
    m=sum(expTc<mean.expT, na.rm=T)
    print(m)
}
gethypermethTlength=function(probe){
    sum(hypermethDataT[as.character(probe),]>hypercutoff, na.rm=T)}
gethypomethhigherTlength=function(geneID, probe){
    metT=metDataTcat[as.character(probe),]
    metTc=metT[which(metT==0)]
    expT=expDataT[as.character(geneID),]
    expT0=expDataT0[as.character(geneID),]
    expT0=expT0[which(expT0==1)]
    expTc=expT[match(names(metTc),names(expT))]
    expTc=expTc[which(names(expTc)%ni%names(expT0))]
    mean.expT=mean(expT, na.rm=T)
    m=sum(expTc>mean.expT, na.rm=T)
    print(m)
}
gethypomethlowerTlength=function(geneID, probe){
    metT=metDataTcat[as.character(probe),]
    metTc=metT[which(metT==0)]
    expT=expDataT[as.character(geneID),]
    expT0=expDataT0[as.character(geneID),]
    expT0=expT0[which(expT0==1)]
    expTc=expT[match(names(metTc),names(expT))]
    expTc=expTc[which(names(expTc)%ni%names(expT0))]
    mean.expT=mean(expT, na.rm=T)
    m=sum(expTc<mean.expT, na.rm=T)
    print(m)
}
gethypomethTlength=function(probe){
    sum(hypomethDataT[as.character(probe),]<hypocutoff, na.rm=T)}
getmaxmetTc=function(probe){
    metT=metDataTcat[as.character(probe),]
    metTc=metT[which(metT==0)]
    m=paste(max(as.numeric(hypermethDataT[as.character(probe),match(names(metTc), colnames(hypermethDataT))]),na.rm=T))
    print(m)
}
getmeanexpN=function(geneID){
    mean(expDataN[as.character(geneID),], na.rm=T)}
getmeanexpT=function(geneID){
    mean(expDataT[as.character(geneID),], na.rm=T)}
getmeanhypermethT=function(geneID, probe){
    metT=metDataTcat[as.character(probe),]
    metTc=metT[which(metT==0)]
    expT=expDataT[as.character(geneID),]
    expT0=expDataT0[as.character(geneID),]
    expT0=expT0[which(expT0==1)]
    expTc=expT[match(names(metTc),names(expT))]
    expTc=expTc[which(names(expTc)%ni%names(expT0))]
    m=mean(expTc, na.rm=T)
    print(m)
}
getmeanhypomethT=function(geneID, probe){
    metT=metDataTcat[as.character(probe),]
    metTc=metT[which(metT==0)]
    expT=expDataT[as.character(geneID),]
    expT0=expDataT0[as.character(geneID),]
    expT0=expT0[which(expT0==1)]
    expTc=expT[match(names(metTc),names(expT))]
    expTc=expTc[which(names(expTc)%ni%names(expT0))]
    m=mean(expTc, na.rm=T)
    print(m)
}
getminmetTc=function(probe){
    metT=metDataTcat[as.character(probe),]
    metTc=metT[which(metT==0)]
    m=paste(min(as.numeric(hypomethDataT[as.character(probe),match(names(metTc), colnames(hypomethDataT))]),na.rm=T))
    print(m)
}
getWilcoxexpNcTc=function(probe, geneID){
    metT=metDataTcat[as.character(probe),]
    metTc=metT[which(metT==0)]
    expN=expDataN[as.character(geneID),]
    expN0=expDataN0[as.character(geneID),]
    expN0=expN0[which(expN0==1)]
    expNc=expN[which(names(expN)%ni%names(expN0))]
    expT=expDataT[as.character(geneID),]
    expT0=expDataT0[as.character(geneID),]
    expT0=expT0[which(expT0==1)]
    expTc=expT[match(names(metTc),names(expT))]
    expTc=expTc[which(names(expTc)%ni%names(expT0))]
    if (length(expNc)>0 & length(expTc)>0){
        wilcoxPval=paste(wilcox.test(expNc, expTc)$p.value)
    }
    else{
        wilcoxPval=paste(NA)
    }
}


get.optimized.hypermeth.G_pos.pairs.realfast.from.step1 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    geneanno=read.delim("../scripts/data/gene.anno.txt", header=T)
    hm450probelist=rep(rownames(hypermethDataT),each=dim(geneanno)[1])
    selected=data.frame(probe=hm450probelist, geneSymbol=rep(geneanno$geneName,length(rownames(hypermethDataT))), geneID=rep(geneanno$geneID, length(rownames(hypermethDataT))))
    selected$Z.real=paste(NA)
    selected$pval=paste(NA)
    if (dim(selected)[1]<1){
        print("no data found from step1, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypermeth.tumor.length=paste(NA)
    selected$meanHypermethT=paste(NA)
    selected$hypermeth.lower.tumor.length=paste(NA)
    selected$max.metTc=paste(NA)
    selected$mean.expN.high.expT=paste(NA)
    metDataTcat=ifelse(hypermethDataT<hypercutoff,1,0)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypermeth.tumor.length=mcmapply(gethypermethTlength, top$probe, mc.cores=cores)
    selected$meanHypermethT=mcmapply(getmeanhypermethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypermeth.lower.tumor.length=mcmapply(gethypermethlowerTlength, top$geneID, top$probe, mc.cores=cores)
    selected$max.metTc=mcmapply(getmaxmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.high.expT=ifelse(as.numeric(selected$mean.expN)>as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.high.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypermeth.tumor.length)>minTumor & as.numeric(selectedc$hypermeth.lower.tumor.length)>minTumor & as.numeric(selectedc$max.metTc)>hyper.stringency),]
    dir.create("hyper.G+.output")
    setwd("./hyper.G+.output")
    write.table(selectedc2, file=paste("hyper.G+.link.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.optimized.hypermeth.G_pos.pairs.realfast.from.step2 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    selected=read.delim(paste("../step2/hyper.output/", "hyper.G+.link.zscore.txt", sep=""), header=F)
    colnames(selected)=c("probe", "geneSymbol", "geneID", "Z.real")
    selected$pval=paste(NA)
    if (dim(selected)[1]<1){
        print("no links found from step2, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypermeth.tumor.length=paste(NA)
    selected$meanHypermethT=paste(NA)
    selected$hypermeth.lower.tumor.length=paste(NA)
    selected$max.metTc=paste(NA)
    selected$mean.expN.high.expT=paste(NA)
    metDataTcat=ifelse(hypermethDataT<hypercutoff,1,0)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypermeth.tumor.length=mcmapply(gethypermethTlength, top$probe, mc.cores=cores)
    selected$meanHypermethT=mcmapply(getmeanhypermethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypermeth.lower.tumor.length=mcmapply(gethypermethlowerTlength, top$geneID, top$probe, mc.cores=cores)
    selected$max.metTc=mcmapply(getmaxmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.high.expT=ifelse(as.numeric(selected$mean.expN)>as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.high.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypermeth.tumor.length)>minTumor & as.numeric(selectedc$hypermeth.lower.tumor.length)>minTumor & as.numeric(selectedc$max.metTc)>hyper.stringency),]
    dir.create("hyper.G+.output")
    setwd("./hyper.G+.output")
    write.table(selectedc2, file=paste("hyper.G+.link.zscore.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.optimized.hypermeth.G_pos.pairs.realfast.from.step3 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    selected=read.delim(paste("../step3/hyper.output/hyper.G+.output/", "hyper.G+.link.zscore.perm.txt", sep=""), header=T)
    if (dim(selected)[1]<1){
        print("no links found from step3, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypermeth.tumor.length=paste(NA)
    selected$meanHypermethT=paste(NA)
    selected$hypermeth.lower.tumor.length=paste(NA)
    selected$max.metTc=paste(NA)
    selected$mean.expN.high.expT=paste(NA)
    metDataTcat=ifelse(hypermethDataT<hypercutoff,1,0)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypermeth.tumor.length=mcmapply(gethypermethTlength, top$probe, mc.cores=cores)
    selected$meanHypermethT=mcmapply(getmeanhypermethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypermeth.lower.tumor.length=mcmapply(gethypermethlowerTlength, top$geneID, top$probe, mc.cores=cores)
    selected$max.metTc=mcmapply(getmaxmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.high.expT=ifelse(as.numeric(selected$mean.expN)>as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.high.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypermeth.tumor.length)>minTumor & as.numeric(selectedc$hypermeth.lower.tumor.length)>minTumor & as.numeric(selectedc$max.metTc)>hyper.stringency),]
    dir.create("hyper.G+.output")
    setwd("./hyper.G+.output")
    write.table(selectedc2, file=paste("hyper.G+.link.zscore.perm.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.optimized.hypermeth.G_neg.pairs.realfast.from.step1 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    geneanno=read.delim("../scripts/data/gene.anno.txt", header=T)
    hm450probelist=rep(rownames(hypermethDataT),each=dim(geneanno)[1])
    selected=data.frame(probe=hm450probelist, geneSymbol=rep(geneanno$geneName,length(rownames(hypermethDataT))), geneID=rep(geneanno$geneID, length(rownames(hypermethDataT))))
    selected$Z.real=paste(NA)
    selected$pval=paste(NA)
    if (dim(selected)[1]<1){
        print("no data found from step1, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypermeth.tumor.length=paste(NA)
    selected$meanHypomethT=paste(NA)
    selected$hypermeth.higher.tumor.length=paste(NA)
    selected$max.metTc=paste(NA)
    selected$mean.expN.low.expT=paste(NA)
    metDataTcat=ifelse(hypermethDataT<hypercutoff,1,0)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypermeth.tumor.length=mcmapply(gethypermethTlength, top$probe, mc.cores=cores)
    selected$meanHypermethT=mcmapply(getmeanhypermethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypermeth.higher.tumor.length=mcmapply(gethypermethhigherTlength, top$geneID, top$probe, mc.cores=cores)
    selected$max.metTc=mcmapply(getmaxmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.low.expT=ifelse(as.numeric(selected$mean.expN)<as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.low.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypermeth.tumor.length)>minTumor & as.numeric(selectedc$hypermeth.higher.tumor.length)>minTumor & as.numeric(selectedc$max.metTc)>hyper.stringency),]
    dir.create("hyper.G-.output")
    setwd("./hyper.G-.output")
    write.table(selectedc2, file=paste("hyper.G-.link.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.optimized.hypermeth.G_neg.pairs.realfast.from.step2 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    selected=read.delim(paste("../step2/hyper.output/","hyper.G-.link.zscore.txt", sep=""), header=F)
    colnames(selected)=c("probe", "geneSymbol", "geneID", "Z.real")
    selected$pval=paste(NA)
    if (dim(selected)[1]<1){
        print("no links found from step2, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypermeth.tumor.length=paste(NA)
    selected$meanHypomethT=paste(NA)
    selected$hypermeth.higher.tumor.length=paste(NA)
    selected$max.metTc=paste(NA)
    selected$mean.expN.low.expT=paste(NA)
    metDataTcat=ifelse(hypermethDataT<hypercutoff,1,0)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypermeth.tumor.length=mcmapply(gethypermethTlength, top$probe, mc.cores=cores)
    selected$meanHypermethT=mcmapply(getmeanhypermethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypermeth.higher.tumor.length=mcmapply(gethypermethhigherTlength, top$geneID, top$probe, mc.cores=cores)
    selected$max.metTc=mcmapply(getmaxmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.low.expT=ifelse(as.numeric(selected$mean.expN)<as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.low.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypermeth.tumor.length)>minTumor & as.numeric(selectedc$hypermeth.higher.tumor.length)>minTumor & as.numeric(selectedc$max.metTc)>hyper.stringency),]
    dir.create("hyper.G-.output")
    setwd("./hyper.G-.output")
    write.table(selectedc2, file=paste("hyper.G-.link.zscore.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.optimized.hypermeth.G_neg.pairs.realfast.from.step3 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    selected=read.delim(paste("../step3/hyper.output/hyper.G-.output/", "hyper.G-.link.zscore.perm.txt", sep=""), header=T)
    if (dim(selected)[1]<1){
        print("no links found from step3, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypermeth.tumor.length=paste(NA)
    selected$meanHypomethT=paste(NA)
    selected$hypermeth.higher.tumor.length=paste(NA)
    selected$max.metTc=paste(NA)
    selected$mean.expN.low.expT=paste(NA)
    metDataTcat=ifelse(hypermethDataT<hypercutoff,1,0)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypermeth.tumor.length=mcmapply(gethypermethTlength, top$probe, mc.cores=cores)
    selected$meanHypermethT=mcmapply(getmeanhypermethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypermeth.higher.tumor.length=mcmapply(gethypermethhigherTlength, top$geneID, top$probe, mc.cores=cores)
    selected$max.metTc=mcmapply(getmaxmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.low.expT=ifelse(as.numeric(selected$mean.expN)<as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.low.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypermeth.tumor.length)>minTumor & as.numeric(selectedc$hypermeth.higher.tumor.length)>minTumor & as.numeric(selectedc$max.metTc)>hyper.stringency),]
    dir.create("hyper.G-.output")
    setwd("./hyper.G-.output")
    write.table(selectedc2, file=paste("hyper.G-.link.zscore.perm.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.optimized.hypometh.G_pos.pairs.realfast.from.step1 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    geneanno=read.delim("../scripts/data/gene.anno.txt", header=T)
    hm450probelist=rep(rownames(hypomethDataT),each=dim(geneanno)[1])
    selected=data.frame(probe=hm450probelist, geneSymbol=rep(geneanno$geneName,length(rownames(hypomethDataT))), geneID=rep(geneanno$geneID, length(rownames(hypomethDataT))))
    selected$Z.real=paste(NA)
    selected$pval=paste(NA)
    if (dim(selected)[1]<1){
        print("no data found from step1, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypometh.tumor.length=paste(NA)
    selected$meanHypomethT=paste(NA)
    selected$hypometh.higher.tumor.length=paste(NA)
    selected$min.metTc=paste(NA)
    selected$mean.expN.low.expT=paste(NA)
    metDataTcat=ifelse(hypomethDataT<hypocutoff,0,1)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypometh.tumor.length=mcmapply(gethypomethTlength, top$probe, mc.cores=cores)
    selected$meanHypomethT=mcmapply(getmeanhypomethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypometh.higher.tumor.length=mcmapply(gethypomethhigherTlength, top$geneID, top$probe, mc.cores=cores)
    selected$min.metTc=mcmapply(getminmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.low.expT=ifelse(as.numeric(selected$mean.expN)<as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.low.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypometh.tumor.length)>minTumor & as.numeric(selectedc$hypometh.higher.tumor.length)>minTumor & as.numeric(selectedc$min.metTc)<hypo.stringency),]
    dir.create("hypo.G+.output")
    setwd("./hypo.G+.output")
    write.table(selectedc2, file=paste("hypo.G+.link.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.optimized.hypometh.G_pos.pairs.realfast.from.step2 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    selected=read.delim(paste("../step2/hypo.output/", "hypo.G+.link.zscore.txt", sep=""), header=F)
    colnames(selected)=c("probe", "geneSymbol", "geneID", "Z.real")
    selected$pval=paste(NA)
    if (dim(selected)[1]<1){
        print("no links found from step2, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypometh.tumor.length=paste(NA)
    selected$meanHypomethT=paste(NA)
    selected$hypometh.higher.tumor.length=paste(NA)
    selected$min.metTc=paste(NA)
    selected$mean.expN.low.expT=paste(NA)
    metDataTcat=ifelse(hypomethDataT<hypocutoff,0,1)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypometh.tumor.length=mcmapply(gethypomethTlength, top$probe, mc.cores=cores)
    selected$meanHypomethT=mcmapply(getmeanhypomethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypometh.higher.tumor.length=mcmapply(gethypomethhigherTlength, top$geneID, top$probe, mc.cores=cores)
    selected$min.metTc=mcmapply(getminmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.low.expT=ifelse(as.numeric(selected$mean.expN)<as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.low.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypometh.tumor.length)>minTumor & as.numeric(selectedc$hypometh.higher.tumor.length)>minTumor & as.numeric(selectedc$min.metTc)<hypo.stringency),]
    dir.create("hypo.G+.output")
    setwd("./hypo.G+.output")
    write.table(selectedc2, file=paste("hypo.G+.link.zscore.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.optimized.hypometh.G_pos.pairs.realfast.from.step3 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    selected=read.delim(paste("../step3/hypo.output/hypo.G+.output/","hypo.G+.link.zscore.perm.txt", sep=""), header=T)
    if (dim(selected)[1]<1){
        print("no links found from step3, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypometh.tumor.length=paste(NA)
    selected$meanHypomethT=paste(NA)
    selected$hypometh.higher.tumor.length=paste(NA)
    selected$min.metTc=paste(NA)
    selected$mean.expN.low.expT=paste(NA)
    metDataTcat=ifelse(hypomethDataT<hypocutoff,0,1)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypometh.tumor.length=mcmapply(gethypomethTlength, top$probe, mc.cores=cores)
    selected$meanHypomethT=mcmapply(getmeanhypomethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypometh.higher.tumor.length=mcmapply(gethypomethhigherTlength, top$geneID, top$probe, mc.cores=cores)
    selected$min.metTc=mcmapply(getminmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.low.expT=ifelse(as.numeric(selected$mean.expN)<as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.low.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypometh.tumor.length)>minTumor & as.numeric(selectedc$hypometh.higher.tumor.length)>minTumor & as.numeric(selectedc$min.metTc)<hypo.stringency),]
    dir.create("hypo.G+.output")
    setwd("./hypo.G+.output")
    write.table(selectedc2, file=paste("hypo.G+.link.zscore.perm.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.optimized.hypometh.G_neg.pairs.realfast.from.step1 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    geneanno=read.delim("../scripts/data/gene.anno.txt", header=T)
    hm450probelist=rep(rownames(hypomethDataT),each=dim(geneanno)[1])
    selected=data.frame(probe=hm450probelist, geneSymbol=rep(geneanno$geneName,length(rownames(hypomethDataT))), geneID=rep(geneanno$geneID, length(rownames(hypomethDataT))))
    selected$Z.real=paste(NA)
    selected$pval=paste(NA)
    if (dim(selected)[1]<1){
        print("no data found from step1, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypometh.tumor.length=paste(NA)
    selected$meanHypomethT=paste(NA)
    selected$hypometh.lower.tumor.length=paste(NA)
    selected$min.metTc=paste(NA)
    selected$mean.expN.high.expT=paste(NA)
    metDataTcat=ifelse(hypomethDataT<hypocutoff,0,1)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypometh.tumor.length=mcmapply(gethypomethTlength, top$probe, mc.cores=cores)
    selected$meanHypomethT=mcmapply(getmeanhypomethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypometh.lower.tumor.length=mcmapply(gethypomethlowerTlength, top$geneID, top$probe, mc.cores=cores)
    selected$min.metTc=mcmapply(getminmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.high.expT=ifelse(as.numeric(selected$mean.expN)>as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.high.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypometh.tumor.length)>minTumor & as.numeric(selectedc$hypometh.lower.tumor.length)>minTumor & as.numeric(selectedc$min.metTc)<hypo.stringency),]
    dir.create("hypo.G-.output")
    setwd("./hypo.G-.output")
    write.table(selectedc2, file=paste("hypo.G-.link.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.optimized.hypometh.G_neg.pairs.realfast.from.step3 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    selected=read.delim(paste("../step3/hypo.output/hypo.G-.output/", "hypo.G-.link.zscore.perm.txt", sep=""), header=T)
    if (dim(selected)[1]<1){
        print("no links found from step3, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypometh.tumor.length=paste(NA)
    selected$meanHypomethT=paste(NA)
    selected$hypometh.lower.tumor.length=paste(NA)
    selected$min.metTc=paste(NA)
    selected$mean.expN.high.expT=paste(NA)
    metDataTcat=ifelse(hypomethDataT<hypocutoff,0,1)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypometh.tumor.length=mcmapply(gethypomethTlength, top$probe, mc.cores=cores)
    selected$meanHypomethT=mcmapply(getmeanhypomethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypometh.lower.tumor.length=mcmapply(gethypomethlowerTlength, top$geneID, top$probe, mc.cores=cores)
    selected$min.metTc=mcmapply(getminmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.high.expT=ifelse(as.numeric(selected$mean.expN)>as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.high.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypometh.tumor.length)>minTumor & as.numeric(selectedc$hypometh.lower.tumor.length)>minTumor & as.numeric(selectedc$min.metTc)<hypo.stringency),]
    dir.create("hypo.G-.output")
    setwd("./hypo.G-.output")
    write.table(selectedc2, file=paste("hypo.G-.link.zscore.perm.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.optimized.hypometh.G_neg.pairs.realfast.from.step2 <- function() {
    load("../settings.rda")
    # load diffmethylated enhancer region DNA meth data #
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    `%ni%` <- Negate(`%in%`)
    selected=read.delim(paste("../step2/hypo.output/", "hypo.G-.link.zscore.txt", sep=""), header=F)
    colnames(selected)=c("probe", "geneSymbol", "geneID", "Z.real")
    selected$pval=paste(NA)
    if (dim(selected)[1]<1){
        print("no links found from step2, so cannot run step4")}
    selected$mean.expN=paste(NA)
    selected$mean.expT=paste(NA)
    selected$wilcox.expNcTc=paste(NA)
    selected$hypometh.tumor.length=paste(NA)
    selected$meanHypomethT=paste(NA)
    selected$hypometh.lower.tumor.length=paste(NA)
    selected$min.metTc=paste(NA)
    selected$mean.expN.high.expT=paste(NA)
    metDataTcat=ifelse(hypomethDataT<hypocutoff,0,1)
    expDataT0=ifelse(expDataT==0,1,0)
    expDataN0=ifelse(expDataN==0,1,0)
    top=list(geneID=as.character(selected$geneID),probe=as.character(selected$probe))
    # library(parallel)
    selected$mean.expN=mcmapply(getmeanexpN, top$geneID, mc.cores=cores)
    selected$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    selected$wilcox.expNcTc=mcmapply(getWilcoxexpNcTc, top$probe,top$geneID, mc.cores=cores)
    selected$hypometh.tumor.length=mcmapply(gethypomethTlength, top$probe, mc.cores=cores)
    selected$meanHypomethT=mcmapply(getmeanhypomethT, top$geneID, top$probe, mc.cores=cores)
    selected$hypometh.lower.tumor.length=mcmapply(gethypomethlowerTlength, top$geneID, top$probe, mc.cores=cores)
    selected$min.metTc=mcmapply(getminmetTc, top$probe, mc.cores=cores)
    selected$mean.expN.high.expT=ifelse(as.numeric(selected$mean.expN)>as.numeric(selected$mean.expT),1,0)
    selected$wilcox.expNcTc.adj.pval=p.adjust(as.numeric(selected$wilcox.expNcTc), "BH")
    selectedc=selected[which(as.numeric(selected$mean.expN)!=0 & as.numeric(selected$mean.expT)!=0 & as.numeric(selected$mean.expN.high.expT)==1),]
    selectedc2=selectedc[which(selectedc$wilcox.expNcTc.adj.pval<adj.pval.cutoff & as.numeric(selectedc$hypometh.tumor.length)>minTumor & as.numeric(selectedc$hypometh.lower.tumor.length)>minTumor & as.numeric(selectedc$min.metTc)<hypo.stringency),]
    dir.create("hypo.G-.output")
    setwd("./hypo.G-.output")
    write.table(selectedc2, file=paste("hypo.G-.link.zscore.all.optimized.links.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}
