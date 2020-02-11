load("../settings.rda")
links=c("hyper.G-")
LS=list.files("../step4/hyper.G-.output/", pattern="all.optimized.links.txt")
all=read.delim(paste("../step4/hyper.G-.output/", LS[1], sep=""), header=T)
gene=read.delim("../scripts/data/gene.anno.hg38.txpt.bed", header=T)
probe=read.delim("../scripts/data/hm450cg_GEO.GR.bed", header=F)
dir.create("hyper.G-.output.tracks")
setwd("./hyper.G-.output.tracks")
# let's annotate hm450 and gene TSS start site #
all$chr=probe[match(all$probe, probe$V4),]$V1
all$probe.start=probe[match(all$probe, probe$V4),]$V2
all$probe.end=probe[match(all$probe, probe$V4),]$V3
all$gene.chr=gene[match(all$geneID, gene$geneID),]$chr
all$gene.start=gene[match(all$geneID, gene$geneID),]$start
all$gene.end=all$gene.start+2
all$gene.nameID=gene[match(all$geneID, gene$geneID),]$Feature
all$name=paste(all$chr, ":", all$probe.start, "..", all$probe.end, "-", all$gene.chr, ":", all$gene.start, "..", all$gene.end, ",2", sep="")
# Let's select the link which are in intrachromosome 
all$intrachrom=as.character(all$chr)==as.character(all$gene.chr)
write.table(all, file=paste(prefix, links, "links.hg38.anno.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
all=all[which(all$intrachrom==TRUE),]
# here we can use pval empirical pval to calculate score
all$score=round((1-all$pval)*1000)
all$strand=c(".")
all$downstream=ifelse(all$probe.start<all$gene.start,1,0)
all$chromStart=paste(NA)
all[which(all$downstream==1),]$chromStart=all[which(all$downstream==1),]$probe.start
all[which(all$downstream==0),]$chromStart=all[which(all$downstream==0),]$gene.start
all$chromEnd=paste(NA)
all[which(all$downstream==1),]$chromEnd=all[which(all$downstream==1),]$gene.end
all[which(all$downstream==0),]$chromEnd=all[which(all$downstream==0),]$probe.end
# same number for all lines
all$reserved=c(16711680)
all$blockCount=c(2)
all$blockSizes=paste("2,2")
all$blockStarts=paste(NA)
all[which(all$downstream==1),]$blockStarts=paste("0", all[which(all$downstream==1),]$gene.start-all[which(all$downstream==1),]$probe.start, sep=",")
all[which(all$downstream==0),]$blockStarts=paste("0", all[which(all$downstream==0),]$probe.start-all[which(all$downstream==0),]$gene.start, sep=",")
write.table(all, file=paste(prefix, links, "intrachrom.links.hg38.anno.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
selnames=c("chr", "chromStart", "chromEnd", "name", "score", "strand", "chromStart", "chromEnd", "reserved", "blockCount", "blockSizes", "blockStarts")
allc=all[,match(selnames, colnames(all))]
colnames(allc)=c("chr", "chromStart", "chromEnd", "name", "score", "strand", "chromStart.1", "chromEnd.1", "reserved", "blockCount", "blockSizes", "blockStarts")
write.table(allc, file=paste(prefix, links, "intrachrom.links.hg38.bed", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
setwd("../")
