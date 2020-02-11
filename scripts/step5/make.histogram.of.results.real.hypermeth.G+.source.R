load("../settings.rda")
links=c("hyper.G+")
dir.create("hyper.G+.output.histogram")
setwd("./hyper.G+.output.histogram")
all=read.delim(file=paste("../hyper.G+.output/", prefix, ".", links, ".", "links.hg38.all.anno.txt", sep=""), header=T)
AP=data.frame(table(all$probe))
AP=AP[order(AP$Freq, decreasing=T),]
colnames(AP)=c("probe", "Freq")
AP=AP[which(AP$Freq>0),]
write.table(AP, file=paste(prefix, links, "links.all.hm450.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.all.hm450.freq.hist.pdf", sep="."))
hist(AP$Freq, xlab=c("Number of linked genes per enhancer probe"), main=c(paste("Histogram of", prefix, links, "linked genes per enhancer probe"), sep=" "), xlim=c(1,max(AP$Freq)),breaks=c(seq(1,max(AP$Freq),by=1)),col=histcol)
dev.off()
AG=data.frame(table(all$gene.nameID))
AG=AG[order(AG$Freq, decreasing=T),]
AG$geneSymbol=all[match(AG$Var1, all$gene.nameID),]$geneSymbol
AG$geneID=all[match(AG$Var1, all$gene.nameID),]$geneID
colnames(AG)=c("gene.nameID", "Freq", "geneSymbol", "geneID")
AG=AG[which(AG$Freq>0),]
write.table(AG, file=paste(prefix, links, "links.all.gene.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.all.gene.freq.hist.pdf", sep="."))
hist(AG$Freq, xlab=c("Number of linked enhancer probes per gene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per gene"), sep=" "),xlim=c(1,max(AG$Freq)),breaks=c(seq(1,max(AG$Freq),by=1)),col=histcol)
dev.off()
AT=AG[which(AG$geneID %in% all[which(all$TF==1),]$geneID),]
write.table(AT, file=paste(prefix, links, "links.all.tf.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.all.tf.freq.hist.pdf", sep="."))
hist(AT$Freq, xlab=c("Number of linked enhancer probes per tf"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per tf"), sep=" "),xlim=c(1,max(AT$Freq)),breaks=c(seq(1,max(AT$Freq),by=1)),col=histcol)
dev.off()
AT1=AG[which(AG$geneID %in% all[which(all$OC==1),]$geneID),]
write.table(AT1, file=paste(prefix, links, "links.all.known.oncogene.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.all.known.oncogene.freq.hist.pdf", sep="."))
hist(AT1$Freq, xlab=c("Number of linked enhancer probes per known oncogene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per known oncogene"), sep=" "),xlim=c(1,max(AT1$Freq)),breaks=c(seq(1,max(AT1$Freq),by=1)),col=histcol)
dev.off()
AT2=AG[which(AG$geneID %in% all[which(all$TS==1),]$geneID),]
write.table(AT2, file=paste(prefix, links, "links.all.known.tumorsuppressorgene.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.all.known.tumorsuppressorgene.freq.hist.pdf", sep="."))
hist(AT2$Freq, xlab=c("Number of linked enhancer probes per known tumor suppressor gene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per known tumor suppressor gene"), sep=" "),xlim=c(1,max(AT2$Freq)),breaks=c(seq(1,max(AT2$Freq),by=1)),col=histcol)
dev.off()
# let's do for inter chromosome interactions 
trans=all[which(all$intrachrom==FALSE),]
TP=data.frame(table(trans$probe))
TP=TP[order(TP$Freq, decreasing=T),]
colnames(TP)=c("probe", "Freq")
TP=TP[which(TP$Freq>0),]
write.table(TP, file=paste(prefix, links, "links.trans.hm450.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.trans.hm450.freq.hist.pdf", sep="."))
hist(TP$Freq, xlab=c("Number of linked genes per enhancer probe"), main=c(paste("Histogram of", prefix, links, "linked genes per enhancer probe"), sep=" "),xlim=c(1,max(TP$Freq)),breaks=c(seq(1,max(TP$Freq),by=1)),col=histcol)
dev.off()
TG=data.frame(table(trans$gene.nameID))
TG=TG[order(TG$Freq, decreasing=T),]
TG$geneSymbol=trans[match(TG$Var1, trans$gene.nameID),]$geneSymbol
TG$geneID=trans[match(TG$Var1, trans$gene.nameID),]$geneID
colnames(TG)=c("gene.nameID", "Freq", "geneSymbol", "geneID")
TG=TG[which(TG$Freq>0),]
write.table(TG, file=paste(prefix, links, "links.trans.gene.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.trans.gene.freq.hist.pdf", sep="."))
hist(TG$Freq, xlab=c("Number of linked enhancer probes per gene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per gene"), sep=" "),xlim=c(1,max(TG$Freq)),breaks=c(seq(1,max(TG$Freq),by=1)),col=histcol)
dev.off()
TT=TG[which(TG$geneID %in% trans[which(trans$TF==1),]$geneID),]
write.table(TT, file=paste(prefix, links, "links.trans.tf.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.trans.tf.freq.hist.pdf", sep="."))
hist(TT$Freq, xlab=c("Number of linked enhancer probes per tf"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per tf"), sep=" "),xlim=c(1,max(TT$Freq)),breaks=c(seq(1,max(TT$Freq),by=1)),col=histcol)
dev.off()
TT1=TG[which(TG$geneID %in% trans[which(trans$OC==1),]$geneID),]
write.table(TT1, file=paste(prefix, links, "links.trans.known.oncogene.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.trans.known.oncogene.freq.hist.pdf", sep="."))
hist(TT1$Freq, xlab=c("Number of linked enhancer probes per knwon oncogene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per known oncogene"), sep=" "),xlim=c(1,max(TT1$Freq)),breaks=c(seq(1,max(TT1$Freq),by=1)),col=histcol)
dev.off()
TT2=TG[which(TG$geneID %in% trans[which(trans$TS==1),]$geneID),]
write.table(TT2, file=paste(prefix, links, "links.trans.known.tumorsuppressorgene.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.trans.known.tumorsuppressorgene.freq.hist.pdf", sep="."))
hist(TT2$Freq, xlab=c("Number of linked enhancer probes per known tumor suppressor gene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per known tumor suppressor gene"), sep=" "),xlim=c(1,max(TT2$Freq)),breaks=c(seq(1,max(TT2$Freq),by=1)),col=histcol)
dev.off()
# let's do for cis (intrachromosome) interactions
cis=all[which(all$intrachrom==TRUE),]
cis$downstream=ifelse(cis$probe.start<cis$gene.start,1,0)
CP=data.frame(table(cis$probe))
CP=CP[order(CP$Freq, decreasing=T),]
colnames(CP)=c("probe", "Freq")
CP=CP[which(CP$Freq>0),]
write.table(CP, file=paste(prefix, links, "links.cis.hm450.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.cis.hm450.freq.hist.pdf", sep="."))
hist(CP$Freq, xlab=c("Number of linked genes per enhancer probe"), main=c(paste("Histogram of", prefix, links, "linked genes per enhancer probe"), sep=" "),xlim=c(1,max(CP$Freq)),breaks=c(seq(1,max(CP$Freq),by=1)),col=histcol)
dev.off()
CG=data.frame(table(cis$gene.nameID))
CG=CG[order(CG$Freq, decreasing=T),]
CG$geneSymbol=cis[match(CG$Var1, cis$gene.nameID),]$geneSymbol
CG$geneID=cis[match(CG$Var1, cis$gene.nameID),]$geneID
colnames(CG)=c("gene.nameID", "Freq", "geneSymbol", "geneID")
CG=CG[which(CG$Freq>0),]
write.table(CG, file=paste(prefix, links, "links.cis.gene.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.cis.gene.freq.hist.pdf", sep="."))
hist(CG$Freq, xlab=c("Number of linked enhancer probes per gene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per gene"), sep=" "),xlim=c(1,max(CG$Freq)),breaks=c(seq(1,max(CG$Freq),by=1)),col=histcol)
dev.off()
CT=CG[which(CG$geneID %in% cis[which(cis$TF==1),]$geneID),]
write.table(CT, file=paste(prefix, links, "links.cis.tf.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.cis.tf.freq.hist.pdf", sep="."))
hist(CT$Freq, xlab=c("Number of linked enhancer probes per tf"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per tf"), sep=" "),xlim=c(1,max(CT$Freq)),breaks=c(seq(1,max(CT$Freq),by=1)),col=histcol)
dev.off()
CT1=CG[which(CG$geneID %in% cis[which(cis$OC==1),]$geneID),]
write.table(CT1, file=paste(prefix, links, "links.cis.known.oncogene.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.cis.known.oncogene.freq.hist.pdf", sep="."))
hist(CT1$Freq, xlab=c("Number of linked enhancer probes per known oncogene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per known oncogene"), sep=" "),xlim=c(1,max(CT1$Freq)),breaks=c(seq(1,max(CT1$Freq),by=1)),col=histcol)
dev.off()
CT2=CG[which(CG$geneID %in% cis[which(cis$TS==1),]$geneID),]
write.table(CT2, file=paste(prefix, links, "links.cis.known.tumorsuppressorgene.freq.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
pdf(paste(prefix, links, "links.cis.known.tumorsuppressorgene.freq.hist.pdf", sep="."))
hist(CT2$Freq, xlab=c("Number of linked enhancer probes per known tumor suppressor gene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per known tumor suppreosr gene"), sep=" "),xlim=c(1,max(CT2$Freq)),breaks=c(seq(1,max(CT2$Freq),by=1)),col=histcol)
dev.off()
## For cis, let's check distance from enhancer probe to gene promoter region (do for absolute distance)
cis$dist=cis$gene.start-cis$probe.start
cis$absdist=paste(NA)
cis[which(cis$downstream==1),]$absdist=as.numeric(cis[which(cis$downstream==1),]$gene.start-cis[which(cis$downstream==1),]$probe.start)
cis[which(cis$downstream==0),]$absdist=as.numeric(cis[which(cis$downstream==0),]$probe.start-cis[which(cis$downstream==0),]$gene.start)
cis$absdistcat=paste("NA")
if(dim(cis[which(as.numeric(cis$absdist)<10000),])[1]>0){
cis[which(as.numeric(cis$absdist)<10000),]$absdistcat=paste("1.10kb")}
if(dim(cis[which(as.numeric(cis$absdist)>=10000 & as.numeric(cis$absdist)<100000),])[1]>0){
cis[which(as.numeric(cis$absdist)>=10000 & as.numeric(cis$absdist)<100000),]$absdistcat=paste("2.100kb")}
if(dim(cis[which(as.numeric(cis$absdist)>=100000 & as.numeric(cis$absdist)<1000000),])[1]>0){
cis[which(as.numeric(cis$absdist)>=100000 & as.numeric(cis$absdist)<1000000),]$absdistcat=paste("3.1Mb")}
if(dim(cis[which(as.numeric(cis$absdist)>=1000000 & as.numeric(cis$absdist)<10000000),])[1]>0){
cis[which(as.numeric(cis$absdist)>=1000000 & as.numeric(cis$absdist)<10000000),]$absdistcat=paste("4.10Mb")}
if(dim(cis[which(as.numeric(cis$absdist)>=10000000 & as.numeric(cis$absdist)<100000000),])[1]>0){
cis[which(as.numeric(cis$absdist)>=10000000 & as.numeric(cis$absdist)<100000000),]$absdistcat=paste("5.100Mb")}
if(dim(cis[which(as.numeric(cis$absdist)>=100000000),])[1]>0){
cis[which(as.numeric(cis$absdist)>=100000000),]$absdistcat=paste("6.>100Mb")}
cis.table=table(cis$absdistcat)
cis.tablef=c(0,0,0,0,0,0)
names(cis.tablef)=c("1.10kb", "2.100kb", "3.1Mb", "4.10Mb", "5.100Mb", "6.>100Mb")
cis.tablef=cis.table[match(names(cis.tablef), names(cis.table))]
names(cis.tablef)=c("1.10kb", "2.100kb", "3.1Mb", "4.10Mb", "5.100Mb", "6.>100Mb")
write.table(cis.tablef, file=paste(prefix, links, "links.cis.gene.distance.table.txt", sep="."), row.names=T, col.names=F, quote=F, sep="\t")
pdf(paste(prefix, links, "links.cis.gene.distance.barplot.pdf", sep="."))
barplot(cis.tablef, names.arg=c("10kb", "100kb", "1Mb", "10Mb", "100Mb", ">100Mb"), col="red", ylab=c("Number of enhancer probe:gene links"), xlab=c("Distance between enhancer probe and gene"))
dev.off()
cistf=cis[which(cis$TF==1),]
cistf.table=table(cistf$absdistcat)
cistf.tablef=c(0,0,0,0,0,0)
names(cistf.tablef)=c("1.10kb", "2.100kb", "3.1Mb", "4.10Mb", "5.100Mb", "6.>100Mb")
cistf.tablef=cistf.table[match(names(cistf.tablef), names(cistf.table))]
names(cistf.tablef)=c("1.10kb", "2.100kb", "3.1Mb", "4.10Mb", "5.100Mb", "6.>100Mb")
write.table(cistf.tablef, file=paste(prefix, links, "links.cis.tf.distance.table.txt", sep="."), row.names=T, col.names=F, quote=F, sep="\t")
pdf(paste(prefix, links, "links.cis.tf.distance.barplot.pdf", sep="."))
barplot(cistf.tablef, names.arg=c("10kb", "100kb", "1Mb", "10Mb", "100Mb", ">100Mb"), col="red", ylab=c("Number of enhancer probe:gene (TF) links"), xlab=c("Distance between enhancer probe and gene (TF)"))
dev.off()
cisoc=cis[which(cis$OC==1),]
cisoc.table=table(cisoc$absdistcat)
cisoc.tablef=c(0,0,0,0,0,0)
names(cisoc.tablef)=c("1.10kb", "2.100kb", "3.1Mb", "4.10Mb", "5.100Mb", "6.>100Mb")
cisoc.tablef=cisoc.table[match(names(cisoc.tablef), names(cisoc.table))]
names(cisoc.tablef)=c("1.10kb", "2.100kb", "3.1Mb", "4.10Mb", "5.100Mb", "6.>100Mb")
write.table(cisoc.tablef, file=paste(prefix, links, "links.cis.known.oncogene.distance.table.txt", sep="."), row.names=T, col.names=F, quote=F, sep="\t")
pdf(paste(prefix, links, "links.cis.known.oncogene.distance.barplot.pdf", sep="."))
barplot(cisoc.tablef, names.arg=c("10kb", "100kb", "1Mb", "10Mb", "100Mb", ">100Mb"), col="red", ylab=c("Number of enhancer probe:gene (known oncogene) links"), xlab=c("Distance between enhancer probe and gene (known oncogene)"))
dev.off()
cists=cis[which(cis$TS==1),]
cists.table=table(cists$absdistcat)
cists.tablef=c(0,0,0,0,0,0)
names(cists.tablef)=c("1.10kb", "2.100kb", "3.1Mb", "4.10Mb", "5.100Mb", "6.>100Mb")
cists.tablef=cists.table[match(names(cists.tablef), names(cists.table))]
names(cists.tablef)=c("1.10kb", "2.100kb", "3.1Mb", "4.10Mb", "5.100Mb", "6.>100Mb")
write.table(cists.tablef, file=paste(prefix, links, "links.cis.known.tumorsuppressorgene.distance.table.txt", sep="."), row.names=T, col.names=F, quote=F, sep="\t")
pdf(paste(prefix, links, "links.cis.known.tumorsuppressorgene.distance.barplot.pdf", sep="."))
barplot(cists.tablef, names.arg=c("10kb", "100kb", "1Mb", "10Mb", "100Mb", ">100Mb"), col="red", ylab=c("Number of enhancer probe:gene (known tumor suppressor) links"), xlab=c("Distance between enhancer probe and gene (known tumor suppressor)"))
dev.off()
## let's do for cis within 1Mb ##
cis1Mb=cis[which(as.numeric(cis$absdist)<1000000),]
CP=data.frame(table(cis1Mb$probe))
CP=CP[order(CP$Freq, decreasing=TRUE),]
colnames(CP)=c("probe", "Freq")
CP=CP[which(CP$Freq>0),]
write.table(CP, file=paste(prefix, links, "links.cis.within.1Mb.hm450.freq.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
if (dim(CP)[1]>0){
pdf(paste(prefix, links, "links.cis.within.1Mb.hm450.freq.hist.pdf", sep="."))
hist(CP$Freq, xlab=c("Number of linked genes per enhancer probe"), main=c(paste("Histogram of", prefix, links, "linked genes per enhancer probe"), sep=" "),xlim=c(1,max(CP$Freq)),breaks=c(seq(1,max(CP$Freq),by=1)),col=histcol)
dev.off()
}
CG=data.frame(table(cis1Mb$gene.nameID))
CG=CG[order(CG$Freq, decreasing=TRUE),]
CG$geneSymbol=cis1Mb[match(CG$Var1, cis1Mb$gene.nameID),]$geneSymbol
CG$geneID=cis1Mb[match(CG$Var1, cis1Mb$gene.nameID),]$geneID
colnames(CG)=c("gene.nameID", "Freq", "geneSymbol", "geneID")
CG=CG[which(CG$Freq>0),]
write.table(CG, file=paste(prefix, links, "links.cis.within.1Mb.gene.freq.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
if (dim(CG)[1]>0){
pdf(paste(prefix, links, "links.cis.within.1Mb.gene.freq.hist.pdf", sep="."))
hist(CG$Freq, xlab=c("Number of linked enhancer probes per gene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per gene"), sep=" "),xlim=c(1,max(CG$Freq)),breaks=c(seq(1,max(CG$Freq),by=1)),col=histcol)
dev.off()
}
CT=CG[which(CG$geneID %in% cis1Mb[which(cis1Mb$TF==1),]$geneID),]
write.table(CT, file=paste(prefix, links, "links.cis.within.1Mb.tf.freq.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
if (dim(CT)[1]>0){
pdf(paste(prefix, links, "links.cis.within.1Mb.tf.freq.hist.pdf", sep="."))
hist(CT$Freq, xlab=c("Number of linked enhancer probes per tf"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per tf"), sep=" "),xlim=c(1,max(CT$Freq)),breaks=c(seq(1,max(CT$Freq),by=1)),col=histcol)
dev.off()
}
CT1=CG[which(CG$geneID %in% cis1Mb[which(cis1Mb$OC==1),]$geneID),]
write.table(CT1, file=paste(prefix, links, "links.cis.within.1Mb.known.oncogene.freq.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
if (dim(CT1)[1]>0){
pdf(paste(prefix, links, "links.cis.within.1Mb.known.oncogene.freq.hist.pdf", sep="."))
hist(CT1$Freq, xlab=c("Number of linked enhancer probes per known oncogene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per known oncogene"), sep=" "),xlim=c(1,max(CT1$Freq)),breaks=c(seq(1,max(CT1$Freq),by=1)),col=histcol)
dev.off()
}
CT2=CG[which(CG$geneID %in% cis1Mb[which(cis1Mb$TS==1),]$geneID),]
write.table(CT2, file=paste(prefix, links, "links.cis.within.1Mb.known.tumorsuppressorgene.freq.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
if (dim(CT2)[1]>0){
pdf(paste(prefix, links, "links.cis.within.1Mb.known.tumorsuppressorgene.freq.hist.pdf", sep="."))
hist(CT2$Freq, xlab=c("Number of linked enhancer probes per known tumor suppressor gene"), main=c(paste("Histogram of", prefix, links, "linked enhancer probes per known tumor suppreosr gene"), sep=" "),xlim=c(1,max(CT2$Freq)),breaks=c(seq(1,max(CT2$Freq),by=1)),col=histcol)
dev.off()
}
setwd("../")
