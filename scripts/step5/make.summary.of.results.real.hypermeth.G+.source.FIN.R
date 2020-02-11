load("../settings.rda")
links=c("hyper.G+")
histcol=c("red")
################
LS=list.files("../step4/hyper.G+.output/", pattern="all.optimized.links.txt")
all=read.delim(paste("../step4/hyper.G+.output/", LS[1], sep=""), header=T)
gene=read.delim("../scripts/data/gene.anno.hg38.txpt.bed", header=T)
TF=read.delim("../scripts/data/human.TF.geneID.from.ELMERV2.txt", header=T)
OC=read.delim("../scripts/data/human.527.known.cancer.gene.geneID.from.TheCancerGeneCensus.txt", header=T)
TS=read.delim("../scripts/data/human.637.proteincoding.TSGs.geneID.from.TSGene.txt", header=T)
probe=read.delim("../scripts/data/hm450cg_GEO.GR.bed", header=F)
dir.create("hyper.G+.output")
setwd("./hyper.G+.output")
all$chr=probe[match(all$probe, probe$V4),]$V1
all$probe.start=probe[match(all$probe, probe$V4),]$V2
all$probe.end=probe[match(all$probe, probe$V4),]$V3
all$gene.chr=gene[match(all$geneID, gene$geneID),]$chr
all$gene.start=gene[match(all$geneID, gene$geneID),]$start
all$gene.end=all$gene.start+2
all$gene.nameID=gene[match(all$geneID, gene$geneID),]$Feature
all$TF=ifelse(all$geneID %in% TF$GeneID,1,0)
all$OC=ifelse(all$geneID %in% OC$GeneID,1,0)
all$TS=ifelse(all$geneID %in% TS$GeneID,1,0)
all$intrachrom=as.character(all$chr)==as.character(all$gene.chr)
# let's do for all interactions
write.table(unique(all$probe), file=paste(prefix, links, "links.all.hm450.probelist.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(all$gene.nameID), file=paste(prefix, links, "links.all.gene.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(all[which(all$TF==1),]$gene.nameID), file=paste(prefix, links, "links.all.tf.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(all[which(all$OC==1),]$gene.nameID), file=paste(prefix, links, "links.all.known.oncogene.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(all[which(all$TS==1),]$gene.nameID), file=paste(prefix, links, "links.all.known.tumorsuppressorgene.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(all, file=paste(prefix, links, "links.hg38.all.anno.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
# let's do for inter chromosome interactions 
trans=all[which(all$intrachrom==FALSE),]
write.table(unique(trans$probe), file=paste(prefix, links, "links.trans.hm450.probelist.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(trans$gene.nameID), file=paste(prefix, links, "links.trans.gene.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(trans[which(trans$TF==1),]$gene.nameID), file=paste(prefix, links, "links.trans.tf.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(trans[which(trans$OC==1),]$gene.nameID), file=paste(prefix, links, "links.trans.known.oncogene.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(trans[which(trans$TS==1),]$gene.nameID), file=paste(prefix, links, "links.trans.known.tumorsuppressorgene.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(trans, file=paste(prefix, links, "links.hg38.trans.anno.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
# let's do for cis (intrachromosome) interactions
cis=all[which(all$intrachrom==TRUE),]
cis$downstream=ifelse(cis$probe.start<cis$gene.start,1,0)
write.table(unique(cis$probe), file=paste(prefix, links, "links.cis.hm450.probelist.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(cis$gene.nameID), file=paste(prefix, links, "links.cis.gene.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(cis[which(cis$TF==1),]$gene.nameID), file=paste(prefix, links, "links.cis.tf.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(cis[which(cis$OC==1),]$gene.nameID), file=paste(prefix, links, "links.cis.known.oncogene.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
write.table(unique(cis[which(cis$TS==1),]$gene.nameID), file=paste(prefix, links, "links.cis.known.tumorsuppressorgene.list.txt", sep="."), row.names=F, col.names=F, quote=F, sep="\t")
## For cis, let's check distance from enhancer probe to gene promoter region (do for absolute distance)
cis$dist=cis$gene.start-cis$probe.start
cis$absdist=paste(NA)
cis[which(cis$downstream==1),]$absdist=cis[which(cis$downstream==1),]$gene.start-cis[which(cis$downstream==1),]$probe.start
cis[which(cis$downstream==0),]$absdist=cis[which(cis$downstream==0),]$probe.start-cis[which(cis$downstream==0),]$gene.start
write.table(cis, file=paste(prefix, links, "links.hg38.cis.anno.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
cis1Mb=cis[which(as.numeric(cis$absdist)<1000000),]
write.table(cis1Mb, file=paste(prefix, links, "links.hg38.cis.within.1Mb.anno.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
setwd("../")
