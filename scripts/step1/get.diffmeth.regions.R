## Load the settings data:
load("../settings.rda")

## Import Genomic ranges:
library('GenomicRanges')

## Import the gtf data for gencode v22 genes:
gencode_v22_gtf <- read.delim(
  '../scripts/data/gencode_v22_gtf_file.tsv.gz',
  stringsAsFactors = FALSE
)

## Get info for both genes and transcripts 
## This will cause some near duplicate ranges but that is ok
## Since intersecting will be the same no matter how many identical
## ranges there are:
gencode_v22_transcript <- gencode_v22_gtf[
  gencode_v22_gtf$type=='transcript',
  ]

gencode_v22_genes <- gencode_v22_gtf[
  gencode_v22_gtf$type=='gene',
  ]

## Create the dataset of both genes and transcripts:
gencode_v22_TSS <- rbind(gencode_v22_genes, gencode_v22_transcript)

## Determine whether the start or the end coordinate is the 
## TSS based on the strand:
gencode_v22_TSS$TSS <- ifelse(
  gencode_v22_TSS$strand=='+',
  gencode_v22_TSS$start,
  ifelse(
    gencode_v22_TSS$strand=='-',
    gencode_v22_TSS$end,
    ''
  )
)

## Create a dataset with just the TSS and the spacing representing a 1500bp 
## gap from the 
gencode_v22_TSS_df <- data.frame(
  'chr'= gencode_v22_TSS$seqnames,
  'start'= ifelse(
    as.character(gencode_v22_TSS$strand)=='+',
    as.numeric(gencode_v22_TSS$TSS),
    ifelse(
      as.character(gencode_v22_TSS$strand)=='-',
      as.numeric(gencode_v22_TSS$TSS)-1,
      as.numeric(gencode_v22_TSS$TSS)
    )
  )-udist,
  'end'= ifelse(
    as.character(gencode_v22_TSS$strand)=='+',
    as.numeric(gencode_v22_TSS$TSS)+1,
    ifelse(
      as.character(gencode_v22_TSS$strand)=='-',
      as.numeric(gencode_v22_TSS$TSS),
      as.numeric(gencode_v22_TSS$TSS)+1
    )
  )+(udist-2),
  'strand'= rep('*',nrow(gencode_v22_TSS)),
  'name'= paste(
    sub(
      '\\..*', 
      '', 
      gencode_v22_TSS$gene_id
    ),
    gencode_v22_TSS$type,
    sep='_'
  ),
  stringsAsFactors = FALSE
)
rownames(gencode_v22_TSS_df) <- make.unique(gencode_v22_TSS_df$name)

## Create a granges object from the TSS dataframe:
gencode_v22_TSS_granges <- makeGRangesFromDataFrame(
  df= gencode_v22_TSS_df,
  keep.extra.columns = FALSE,
  starts.in.df.are.0based = FALSE
)

## Load in the annotations for the 450k array probes on the hg38 genome:
hg38_450k_probe_info <- read.delim(
  file='../scripts/data/hm450.hg38.manifest.tsv.gz',
  header= TRUE,
  stringsAsFactors = FALSE
)

## Create a file with a single annotation for the probe locations:
hg38_probes <- data.frame(
  'chr'= hg38_450k_probe_info$CpG_chrm,
  'pos'= hg38_450k_probe_info$CpG_beg+1,
  'strand'= hg38_450k_probe_info$probe_strand,
  stringsAsFactors = FALSE
)
rownames(hg38_probes) <- hg38_450k_probe_info$probeID

## Determine probes which did not map to hg38 genome according to http://zwdzwd.github.io/InfiniumAnnotation
hg38_NA_probes <- rownames(hg38_probes[c(which(is.na(hg38_probes[1]))),])

## Remove those probes from the probe database:
hg38_probes <- hg38_probes[!(row.names(hg38_probes) %in% hg38_NA_probes), ]

# make GRanges for all probes of hm450
hg38_probes_granges <- GRanges(
  seqnames = hg38_probes$chr, 
  ranges = IRanges(
    hg38_probes$pos, 
    width = 1, 
    names = rownames(hg38_probes)
  ), 
  strand = hg38_probes$strand, 
  name = rownames(hg38_probes)
)

# find probes not overlap with TSS windows
probe <- hg38_probes[setdiff(1:nrow(hg38_probes), unique(queryHits(findOverlaps(hg38_probes_granges,gencode_v22_TSS_granges)))),]
nonTSS=rownames(probe)

if (elmerENH==T){
  
  # get ELMER enhancers from v1.5.1
  # For future versions, update this to include multiple sets!:
  TENET_1.5.1_file <- read.delim(
    '../scripts/data/ELMERv.1.5.1_LUAD_TENET_probes.tsv',
    stringsAsFactors = FALSE
  )
  
  DistalP <- rownames(TENET_1.5.1_file)
  
}
if (elmerENH==F){
  DistalP=c(NA)
}
# If there's external.data for enhancer
if (extENH==T){
  external.enhancer=read.delim("../scripts/step1/hm450overlap/enhancer/external.data.enhancer.hm450.probelist.txt", header=F)
}
if (extENH==F){
  external.enhancer=data.frame(V1=c(NA))
}
allenh=unique(c(as.character(external.enhancer$V1),as.character(DistalP)))
nonTSS2=intersect(nonTSS, allenh)
if (encodeNDR==T){
  # let's get probes overlapped with NDR
  NDR=read.delim("../scripts/data/wgEncodeAwgDnaseMasterSites.hg19.hm450cg_GEO_GR.bed.probelist.txt", header=T)
}
if (encodeNDR==F){
  NDR=data.frame(V1=c(NA))
}
if (extNDR==T){
  external.NDR=read.delim("../scripts/step1/hm450overlap/NDR/external.data.NDR.hm450.probelist.txt", header=F)
}
if (extNDR==F){
  external.NDR=data.frame(V1=c(NA))
}
NDR2=unique(c(as.character(NDR$V1), as.character(external.NDR$V1)))
nonTSS.NDR=intersect(nonTSS2, NDR2)
if (onlyextFeature==T){
  external.Feature=read.delim("../scripts/step1/hm450overlap/feature/external.data.feature.hm450.probelist.txt", header=F)
  NDR2=unique(c(as.character(external.Feature$V1)))
  nonTSS.NDR=intersect(nonTSS, NDR2)
}
##### let's obtain methylation and expression data #####
LS=list.files("../external.data/data", pattern="rda")
load(paste("../external.data/data/", LS[1], sep=""))
annoGenes=read.delim("../scripts/data/gene.anno.hg19.txpt.bed", header=T)
expDataT=expDataT[match(as.character(annoGenes$geneID), as.character(rownames(expDataT))),]
rownames(expDataT)=as.character(annoGenes$geneID)
expDataN=expDataN[match(as.character(annoGenes$geneID), as.character(rownames(expDataN))),]
rownames(expDataN)=as.character(annoGenes$geneID)
if (dim(expDataT)[1]<1 | dim(expDataN)[1]<1){
  print("Error, please load the expression data")
}
annoHM450=read.delim("../scripts/data/hm450cg_GEO.GR.bed", header=F)
RmetDataT=metDataT[match(as.character(annoHM450$V4), as.character(rownames(metDataT))),]
rownames(RmetDataT)=as.character(annoHM450$V4)
RmetDataN=metDataN[match(as.character(annoHM450$V4), as.character(rownames(metDataN))),]
rownames(RmetDataN)=as.character(annoHM450$V4)
# let's have dataset for probes in nonTSS enhancer NDR after removing NA
naCount=apply(is.na(RmetDataT),1,sum)
useInd2=naCount[which(naCount!=dim(RmetDataT)[2])]
Enhprobe=intersect(nonTSS.NDR, names(useInd2))
metDataT=RmetDataT[match(Enhprobe, rownames(RmetDataT)),]
metDataN=RmetDataN[match(Enhprobe, rownames(RmetDataN)),]
### let's load leukocyte, fibroblast, smoothmuscle data ##
load("../scripts/data/othercells/leukocytes.rda")
load("../scripts/data/othercells/fibroblast.rda")
load("../scripts/data/othercells/smoothmuscle.rda")
metDataL=leukData[match(Enhprobe, rownames(leukData)),]
metDataF=fibroData[match(Enhprobe, rownames(fibroData)),]
metDataS=smData[match(Enhprobe, rownames(smData)),]
Probe=data.frame(probe=Enhprobe)
Probe$mean.leuk=apply(metDataL, 1, mean, na.rm=T)
Probe$mean.fibro=apply(metDataF, 1, mean, na.rm=T)
Probe$mean.sm=apply(metDataS, 1, mean, na.rm=T)
Probe$mean.normal=apply(metDataN, 1, mean, na.rm=T)
Probe$mean.tumor=apply(metDataT, 1, mean, na.rm=T)
if (othercells==T){
  LA=list.files("../external.data/othercells/", pattern="rda")
  load(paste("../external.data/othercells/", LA[1], sep=""))
  metDataO=extData[match(Enhprobe, rownames(extData)),]
  Probe$mean.ext=apply(metDataO, 1, mean, na.rm=T)
}
##### find hypomethylated regions compared to normal but not hypomethylated in leukocytes, fibroblast, smoothmuscle #####
if (leuk==T & fibro==T & sm==T & othercells==F){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.leuk>methcutoff & Probe$mean.fibro>methcutoff & Probe$mean.sm>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.leuk<unmethcutoff & Probe$mean.fibro<unmethcutoff & Probe$mean.sm<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==T & fibro==T & sm==T & othercells==T){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.leuk>methcutoff & Probe$mean.fibro>methcutoff & Probe$mean.sm>methcutoff & Probes$mean.ext>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.leuk<unmethcutoff & Probe$mean.fibro<unmethcutoff & Probe$mean.sm<unmethcutoff & Probe$mean.ext<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==F & fibro==T & sm==F & othercells==F){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.fibro>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.fibro<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==F & fibro==T & sm==F & othercells==T){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.fibro>methcutoff & Probe$mean.ext>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.fibro<unmethcutoff & Probe$mean.ext<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==F & fibro==T & sm==T & othercells==F){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.fibro>methcutoff & Probe$mean.sm>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.fibro<unmethcutoff & Probe$mean.sm<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==F & fibro==T & sm==T & othercells==F){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.fibro>methcutoff & Probe$mean.sm>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.fibro<unmethcutoff & Probe$mean.sm<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==F & fibro==T & sm==T & othercells==T){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.fibro>methcutoff & Probe$mean.sm>methcutoff & Probe$mean.ext>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.fibro<unmethcutoff & Probe$mean.sm<unmethcutoff & Probe$ext<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==T & fibro==T & sm==F & othercells==F){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.leuk>methcutoff & Probe$mean.fibro>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.leuk<unmethcutoff & Probe$mean.fibro<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==T & fibro==T & sm==F & othercells==T){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.leuk>methcutoff & Probe$mean.fibro>methcutoff & Probe$mean.ext>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.leuk<unmethcutoff & Probe$mean.fibro<unmethcutoff & Probe$mean.ext<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==T & fibro==F & sm==F & othercells==F){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.leuk>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.leuk<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==T & fibro==F & sm==F & othercells==T){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.leuk>methcutoff & Probe$mean.ext>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.leuk<unmethcutoff & Probe$mean.leuk<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==T & fibro==F & sm==T & othercells==F){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.leuk>methcutoff & Probe$mean.sm>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.leuk<unmethcutoff & Probe$mean.sm<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==T & fibro==F & sm==T & othercells==T){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.leuk>methcutoff & Probe$mean.sm>methcutoff & Probe$mean.ext>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.leuk<unmethcutoff & Probe$mean.sm<unmethcutoff & Probe$mean.ext<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==F & fibro==F & sm==T & othercells==F){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.sm>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.sm<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==F & fibro==F & sm==T & othercells==T){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.sm>methcutoff & Probe$mean.ext>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.sm<unmethcutoff & Probe$mean.ext<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}
if (leuk==F & fibro==F & sm==F & othercells==T){
  # let's get methylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.ext>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.ext<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}

if (leuk==F & fibro==F & sm==F & othercells==F){
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff & Probe$mean.tumor>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>methcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  methDataT=temp.t1
  methDataN=metDataN[match(rownames(methDataT),rownames(metDataN)),]
  # let's get hypomethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal>methcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<hypocutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypomethDataT=temp.t1
  hypomethDataN=metDataN[match(rownames(hypomethDataT), rownames(metDataN)),]
  # let's get unmethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff & Probe$mean.tumor<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t<unmethcutoff,0,1)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat<minTumor),]
  unmethDataT=temp.t1
  unmethDataN=metDataN[match(rownames(unmethDataT),rownames(metDataN)),]
  # let's get hypermethylated probes in tumor
  temp.t=metDataT[match(Probe[which(Probe$mean.normal<unmethcutoff),]$probe, rownames(metDataT)),]
  temp.t.cat=ifelse(temp.t>hypercutoff,1,0)
  length.temp.t.cat=apply(temp.t.cat, 1, sum, na.rm=T)
  temp.t1=temp.t[which(length.temp.t.cat>minTumor),]
  hypermethDataT=temp.t1
  hypermethDataN=metDataN[match(rownames(hypermethDataT),rownames(metDataN)),]
}

dir.create("output")
setwd("./output")

enhmetDataT=metDataT
enhmetDataN=metDataN
metDataT=RmetDataT
metDataN=RmetDataN

## temporarily rebind expression and methylation data:
temp_met_enh <- cbind(enhmetDataT, enhmetDataN)
temp_met <- cbind(metDataT, metDataN)
temp_exp <- cbind(expDataT, expDataN)

temp_met_enh_NA_Index <- apply(temp_met_enh, 1, function(x) all(is.na(x)))
temp_met_NA_Index <- apply(temp_met, 1, function(x) all(is.na(x)))
temp_exp_NA_Index <- apply(temp_exp, 1, function(x) all(is.na(x)))

temp_met_enh <- temp_met_enh[ !temp_met_enh_NA_Index, ]
temp_met <- temp_met[ !temp_met_NA_Index, ]
temp_exp <- temp_exp[ !temp_exp_NA_Index, ]

expDataN <- expDataN[rownames(temp_exp),]
expDataT <- expDataT[rownames(temp_exp),]

metDataN <- metDataN[rownames(temp_met),]
metDataT <- metDataT[rownames(temp_met),]

enhmetDataN <- enhmetDataN[rownames(temp_met_enh),]
enhmetDataT <- enhmetDataT[rownames(temp_met_enh),]

## Remove non-cg probes:

enhmetDataN <- enhmetDataN[which(substring(rownames(enhmetDataN),1,2)=='cg'),]
enhmetDataT <- enhmetDataT[which(substring(rownames(enhmetDataT),1,2)=='cg'),]

hypermethDataN <- hypermethDataN[which(substring(rownames(hypermethDataN),1,2)=='cg'),]
hypermethDataT <- hypermethDataT[which(substring(rownames(hypermethDataT),1,2)=='cg'),]

methDataN <- methDataN[which(substring(rownames(methDataN),1,2)=='cg'),]
methDataT <- methDataT[which(substring(rownames(methDataT),1,2)=='cg'),]

hypomethDataN <- hypomethDataN[which(substring(rownames(hypomethDataN),1,2)=='cg'),]
hypomethDataT <- hypomethDataT[which(substring(rownames(hypomethDataT),1,2)=='cg'),]

unmethDataN <- unmethDataN[which(substring(rownames(unmethDataN),1,2)=='cg'),]
unmethDataT <- unmethDataT[which(substring(rownames(unmethDataT),1,2)=='cg'),]



names <- c(
  'unmeth_probes', 
  'hypermeth_probes', 
  'meth_probes', 
  'hypometh_probes', 
  'enhmet_probes', 
  'total_genes', 
  'total_probes'
)

values <- c(
  nrow(unmethDataN),
  nrow(hypermethDataN),
  nrow(methDataN),
  nrow(hypomethDataN),
  nrow(enhmetDataN),
  nrow(expDataN),
  nrow(metDataN)
)

overall_data <- data.frame('data'=names,'values'=values, stringsAsFactors = FALSE)

write.table(
  x= overall_data,
  file="TENET_step1_overall_metadata.txt",
  quote=FALSE,
  sep='\t',
  row.names=FALSE
)



if (dim(expDataT)[1]>1 & dim(expDataN)[1]>1){
  save(metDataT, metDataN, expDataT, expDataN, enhmetDataT, enhmetDataN, methDataT, methDataN, hypomethDataT, hypomethDataN, unmethDataT, unmethDataN, hypermethDataT, hypermethDataN, file=paste(prefix, "diff.methylated.datasets.rda", sep="."))
}
if (dim(expDataT)[1]<1 | dim(expDataN)[1]<1){
  save(metDataT, metDataN, enhmetDataT, enhmetDataN, methDataT, methDataN, hypomethDataT, hypomethDataN, unmethDataT, unmethDataN, hypermethDataT, hypermethDataN, file=paste(prefix, "diff.methylated.datasets.rda", sep="."))
}

gene_names <- rownames(expDataN)
gene_names_ordered <- gene_names[order(gene_names)]

write(
  x=gene_names_ordered,
  file=('TENET_complete_gene_list.txt'),
  ncolumns = 1
)

