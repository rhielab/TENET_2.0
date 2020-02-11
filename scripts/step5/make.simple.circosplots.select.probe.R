## Note: This is currently depreciated with the update to hg38. This script was not used in the paper and will need to be updated with further versions of TENET

load("../settings.rda")
GL=list.files("../external.data/otherinfo", "probe")
probelist=read.delim(paste("../external.data/otherinfo/",GL[1],sep=""),header=F)
links=paste(as.character(strsplit(GL[1], "\\.")[[1]][1]), as.character(strsplit(GL[1], "\\.")[[1]][2]), sep=".")
hm450probe=as.character(probelist$V1)
library(parallel)
library(ggplot2)
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)
if (links=="hyper.G+"){
LS=list.files("../step4/hyper.G+.output/", pattern="all.optimized.links.txt")
all=read.delim(paste("../step4/hyper.G+.output/", LS[1], sep=""), header=T)
gene=read.delim("../scripts/data/gene.anno.hg38.txpt.bed", header=T)
probe=read.delim("../scripts/data/hm450cg_GEO.GR.bed", header=F)
all$chr=probe[match(all$probe, probe$V4),]$V1
all$start=probe[match(all$probe, probe$V4),]$V2
all$end=probe[match(all$probe, probe$V4),]$V3
all$gene.chr=gene[match(all$geneID, gene$geneID),]$chr
all$gene.start=gene[match(all$geneID, gene$geneID),]$start
all$gene.end=all$gene.start+2
all$gene.nameID=gene[match(all$geneID, gene$geneID),]$Feature
dir.create("hyper.G+.output.circosplot")
setwd("hyper.G+.output.circosplot")
}
if (links=="hyper.G-"){
LS=list.files("../step4/hyper.G-.output/", pattern="all.optimized.links.txt")
all=read.delim(paste("../step4/hyper.G-.output/", LS[1], sep=""), header=T)
gene=read.delim("../scripts/data/gene.anno.hg38.txpt.bed", header=T)
probe=read.delim("../scripts/data/hm450cg_GEO.GR.bed", header=F)
all$chr=probe[match(all$probe, probe$V4),]$V1
all$start=probe[match(all$probe, probe$V4),]$V2
all$end=probe[match(all$probe, probe$V4),]$V3
all$gene.chr=gene[match(all$geneID, gene$geneID),]$chr
all$gene.start=gene[match(all$geneID, gene$geneID),]$start
all$gene.end=all$gene.start+2
all$gene.nameID=gene[match(all$geneID, gene$geneID),]$Feature
dir.create("hyper.G-.output.circosplot")
setwd("hyper.G-.output.circosplot")
}
if (links=="hypo.G+"){
LS=list.files("../step4/hypo.G+.output/", pattern="all.optimized.links.txt")
all=read.delim(paste("../step4/hypo.G+.output/", LS[1], sep=""), header=T)
gene=read.delim("../scripts/data/gene.anno.hg38.txpt.bed", header=T)
probe=read.delim("../scripts/data/hm450cg_GEO.GR.bed", header=F)
all$chr=probe[match(all$probe, probe$V4),]$V1
all$start=probe[match(all$probe, probe$V4),]$V2
all$end=probe[match(all$probe, probe$V4),]$V3
all$gene.chr=gene[match(all$geneID, gene$geneID),]$chr
all$gene.start=gene[match(all$geneID, gene$geneID),]$start
all$gene.end=all$gene.start+2
all$gene.nameID=gene[match(all$geneID, gene$geneID),]$Feature
dir.create("hypo.G+.output.circosplot")
setwd("hypo.G+.output.circosplot")
}
if (links=="hypo.G-"){
LS=list.files("../step4/hypo.G-.output/", pattern="all.optimized.links.txt")
all=read.delim(paste("../step4/hypo.G-.output/", LS[1], sep=""), header=T)
gene=read.delim("../scripts/data/gene.anno.hg38.txpt.bed", header=T)
probe=read.delim("../scripts/data/hm450cg_GEO.GR.bed", header=F)
all$chr=probe[match(all$probe, probe$V4),]$V1
all$start=probe[match(all$probe, probe$V4),]$V2
all$end=probe[match(all$probe, probe$V4),]$V3
all$gene.chr=gene[match(all$geneID, gene$geneID),]$chr
all$gene.start=gene[match(all$geneID, gene$geneID),]$start
all$gene.end=all$gene.start+2
all$gene.nameID=gene[match(all$geneID, gene$geneID),]$Feature
dir.create("hypo.G-.output.circosplot")
setwd("hypo.G-.output.circosplot")
}
makeCircosP=function(hm450probe=hm450probe, all=all, circos.col=circos.col){
allc=all[which(all$probe==hm450probe),]
probes=allc[,c(15:17,1)]
genes=allc[,c(18:20,21)]
colnames(probes)=c("chr", "start", "end", "name")
colnames(genes)=c("chr", "start", "end", "name")
probes.genes=rbind(probes, genes)
probes.genes=unique(probes.genes)
Links=probes
Links$V1=genes$chr
Links$V2=genes$start
Links$V3=genes$end
Links$V4=genes$name
colnames(Links)=c("Chromosome", "chromStart", "chromEnd", "GeneName", "Chromosome.1", "chromStart.1", "chromEnd.1", "GeneName.1")
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
tracks.inside <-5
tracks.outside <-0
chr.exclude=NULL
RCircos.Set.Core.Components(cyto.info,chr.exclude,tracks.inside,tracks.outside)
params <- RCircos.Get.Plot.Parameters()
params$base.per.unit <-4500
params$radius.len <- 1.6
params$text.size <- 0.5
params$heatmap.width <-1300
params$track.hieght <-0.1
params$hist.width <-1000
params$plot.radius <-2
RCircos.Reset.Plot.Parameters(params)
out.file <- paste(prefix, hm450probe, links, "links.circosplot.pdf", sep=".")
pdf(file=out.file, height=8, width=8, compress=TRUE)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Connector.Plot(probes.genes, track.num=1, side="in")
RCircos.Gene.Name.Plot(probes.genes, name.col=4,track.num=2, side="in")
Linksc=Links[,c(1:3,5:7)]
Link.data=Linksc
Link.data$PlotColor=circos.col
RCircos.Link.Plot(Link.data, track.num=5, by.chromosome=FALSE)
dev.off()
}
mclapply(hm450probe, makeCircosP, all=all, circos.col=circos.col, mc.cores=cores)
