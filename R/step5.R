find.tumor.states.for.links.hypermeth.G_pos.source.realfast <- function(parallel, prefix, expDataT, states_top_n_genes, metDataT, hypercutoff, cores, i) {
    load("../settings.rda")
    #################
    # library(parallel)
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    expDataT_subC=expDataT
    ##### for hyper.G+ links ####
    links=c("hyper.G+")
    LS=list.files("../step4/hyper.G+.output/", pattern="all.optimized.links.txt")
    x=read.delim(paste("../step4/hyper.G+.output/", LS[1], sep=""), header=T)
    ordered_TFs_by_link_count <- read.delim(file= paste("./hyper.G+.output.histogram/", prefix,".hyper.G+.links.all.tf.freq.txt",sep=''),stringsAsFactors = FALSE)
    top_gene_IDs <- ordered_TFs_by_link_count[1:states_top_n_genes,'geneID']
    x=x[which(x$geneID %in% top_gene_IDs),]
    metDataT_subC=ifelse(metDataT>hypercutoff,1,0)
    datM=t(metDataT_subC[match(x$probe, rownames(metDataT_subC)),])
    datE=t(expDataT_subC[match(x$geneID, rownames(expDataT_subC)),])
    dir.create("hyper.G+.output.states")
    setwd("./hyper.G+.output.states")
    if(is.na(match("mean.expT", colnames(x)))==TRUE){
        getmeanexpT=function(geneID){
            mean(expDataT[as.character(geneID),], na.rm=T)}
        top=list(geneID=as.character(top_gene_IDs))
        x$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    }
    for (i in 1:dim(x)[1]){
        datc=data.frame(meth=datM[,i], exp=datE[,i])
        datc$exp.changed=ifelse(datc$exp>x$mean.expT[i],1,0)
        rownames(datc)=rownames(datM)
        datc$FIN=ifelse(datc$meth==1 & datc$exp.changed==1,1,0)
        colnames(datc)=c(as.character(x$probe[i]), as.character(x$geneID[i]), "exp.changed","FIN")
        write.table(datc, paste(x$probe[i], x$geneID[i], x$geneSymbol[i], prefix, links, "states.txt", sep="."), row.names=T, col.names=T, sep="\t", quote=F)
    }
    TUMOR_Call=matrix("NA", nrow=dim(metDataT)[2], ncol=dim(x)[1])
    rownames(TUMOR_Call)=rownames(datc)
    colnames(TUMOR_Call)=paste(x$probe, x$geneID, x$geneSymbol,sep=".")
    ID=colnames(TUMOR_Call)
    findState=function(ID){
        LS=list.files(getwd(), pattern=as.character(ID))
        LF=read.table(LS, header=T)
        LF$FIN
    }
    TUMOR_Call=mcmapply(ID=ID,findState,mc.cores=cores)
    rownames(TUMOR_Call)=rownames(datc)
    write.table(TUMOR_Call, paste(prefix, links, "links.states.table.txt", sep="."), row.names=T, col.names=T, sep="\t", quote=F)
    # remove unnecessary files #
    notneed=dir(getwd(), pattern="states.txt$")
    file.remove(notneed)
}

find.tumor.states.for.links.hypermeth.G_neg.source.realfast <- function(parallel, prefix, expDataT, states_top_n_genes, metDataT, hypercutoff, cores, i) {
    load("../settings.rda")
    #################
    # library(parallel)
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    expDataT_subC=expDataT
    ##### for hyper.G- links ####
    links=c("hyper.G-")
    LS=list.files("../step4/hyper.G-.output/", pattern="all.optimized.links.txt")
    x=read.delim(paste("../step4/hyper.G-.output/", LS[1], sep=""), header=T)
    ordered_TFs_by_link_count <- read.delim(file= paste("./hyper.G-.output.histogram/", prefix,".hyper.G-.links.all.tf.freq.txt",sep=''),stringsAsFactors = FALSE)
    top_gene_IDs <- ordered_TFs_by_link_count[1:states_top_n_genes,'geneID']
    x=x[which(x$geneID %in% top_gene_IDs),]
    metDataT_subC=ifelse(metDataT>hypercutoff,1,0)
    datM=t(metDataT_subC[match(x$probe, rownames(metDataT_subC)),])
    datE=t(expDataT_subC[match(x$geneID, rownames(expDataT_subC)),])
    dir.create("hyper.G-.output.states")
    setwd("./hyper.G-.output.states")
    if(is.na(match("mean.expT", colnames(x)))==TRUE){
        getmeanexpT=function(geneID){
            mean(expDataT[as.character(geneID),], na.rm=T)}
        top=list(geneID=as.character(top_gene_IDs))
        x$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    }
    for (i in 1:dim(x)[1]){
        datc=data.frame(meth=datM[,i], exp=datE[,i])
        datc$exp.changed=ifelse(datc$exp<x$mean.expT[i],1,0)
        rownames(datc)=rownames(datM)
        datc$FIN=ifelse(datc$meth==1 & datc$exp.changed==1,1,0)
        colnames(datc)=c(as.character(x$probe[i]), as.character(x$geneID[i]), "exp.changed","FIN")
        write.table(datc, paste(x$probe[i], x$geneID[i], x$geneSymbol[i], prefix, links, "states.txt", sep="."), row.names=T, col.names=T, sep="\t", quote=F)
    }
    TUMOR_Call=matrix("NA", nrow=dim(metDataT)[2], ncol=dim(x)[1])
    rownames(TUMOR_Call)=rownames(datc)
    colnames(TUMOR_Call)=paste(x$probe, x$geneID, x$geneSymbol,sep=".")
    ID=colnames(TUMOR_Call)
    findState=function(ID){
        LS=list.files(getwd(), pattern=as.character(ID))
        LF=read.table(LS, header=T)
        LF$FIN
    }
    TUMOR_Call=mcmapply(ID=ID,findState,mc.cores=cores)
    rownames(TUMOR_Call)=rownames(datc)
    write.table(TUMOR_Call, paste(prefix, links, "links.states.table.txt", sep="."), row.names=T, col.names=T, sep="\t", quote=F)
    # remove unnecessary files #
    notneed=dir(getwd(), pattern="states.txt$")
    file.remove(notneed)
}

find.tumor.states.for.links.hypometh.G_pos.source.realfast <- function(parallel, prefix, expDataT, states_top_n_genes, metDataT, hypocutoff, cores, i) {
    load("../settings.rda")
    #################
    # library(parallel)
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    expDataT_subC=expDataT
    ##### for hypo.G+ links ####
    links=c("hypo.G+")
    LS=list.files("../step4/hypo.G+.output/", pattern="all.optimized.links.txt")
    x=read.delim(paste("../step4/hypo.G+.output/", LS[1], sep=""), header=T)
    ordered_TFs_by_link_count <- read.delim(file= paste("./hypo.G+.output.histogram/", prefix,".hypo.G+.links.all.tf.freq.txt",sep=''),stringsAsFactors = FALSE)
    top_gene_IDs <- ordered_TFs_by_link_count[1:states_top_n_genes,'geneID']
    x=x[which(x$geneID %in% top_gene_IDs),]
    metDataT_subC=ifelse(metDataT<hypocutoff,1,0)
    datM=t(metDataT_subC[match(x$probe, rownames(metDataT_subC)),])
    datE=t(expDataT_subC[match(x$geneID, rownames(expDataT_subC)),])
    dir.create("hypo.G+.output.states")
    setwd("./hypo.G+.output.states")
    if(is.na(match("mean.expT", colnames(x)))==TRUE){
        getmeanexpT=function(geneID){
            mean(expDataT[as.character(geneID),], na.rm=T)}
        top=list(geneID=as.character(top_gene_IDs))
        x$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    }
    for (i in 1:dim(x)[1]){
        datc=data.frame(meth=datM[,i], exp=datE[,i])
        datc$exp.changed=ifelse(datc$exp>x$mean.expT[i],1,0)
        rownames(datc)=rownames(datM)
        datc$FIN=ifelse(datc$meth==1 & datc$exp.changed==1,1,0)
        colnames(datc)=c(as.character(x$probe[i]), as.character(x$geneID[i]), "exp.changed","FIN")
        write.table(datc, paste(x$probe[i], x$geneID[i], x$geneSymbol[i], prefix, links, "states.txt", sep="."), row.names=T, col.names=T, sep="\t", quote=F)
    }
    TUMOR_Call=matrix("NA", nrow=dim(metDataT)[2], ncol=dim(x)[1])
    rownames(TUMOR_Call)=rownames(datc)
    colnames(TUMOR_Call)=paste(x$probe, x$geneID, x$geneSymbol,sep=".")
    ID=colnames(TUMOR_Call)
    findState=function(ID){
        LS=list.files(getwd(), pattern=as.character(ID))
        LF=read.table(LS, header=T)
        LF$FIN
    }
    TUMOR_Call=mcmapply(ID=ID,findState,mc.cores=cores)
    rownames(TUMOR_Call)=rownames(datc)
    write.table(TUMOR_Call, paste(prefix, links, "links.states.table.txt", sep="."), row.names=T, col.names=T, sep="\t", quote=F)
    # remove unnecessary files #
    notneed=dir(getwd(), pattern="states.txt$")
    file.remove(notneed)
}

find.tumor.states.for.links.hypometh.G_neg.source.realfast <- function(parallel, prefix, expDataT, states_top_n_genes, metDataT, hypocutoff, cores, i) {
    load("../settings.rda")
    #################
    # library(parallel)
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    expDataT_subC=expDataT
    ##### for hypo.G- links ####
    links=("hypo.G-")
    LS=list.files("../step4/hypo.G-.output/", pattern="all.optimized.links.txt")
    x=read.delim(paste("../step4/hypo.G-.output/", LS[1], sep=""), header=T)
    ordered_TFs_by_link_count <- read.delim(file= paste("./hypo.G-.output.histogram/", prefix,".hypo.G-.links.all.tf.freq.txt",sep=''),stringsAsFactors = FALSE)
    top_gene_IDs <- ordered_TFs_by_link_count[1:states_top_n_genes,'geneID']
    x=x[which(x$geneID %in% top_gene_IDs),]
    metDataT_subC=ifelse(metDataT<hypocutoff,1,0)
    datM=t(metDataT_subC[match(x$probe, rownames(metDataT_subC)),])
    datE=t(expDataT_subC[match(x$geneID, rownames(expDataT_subC)),])
    dir.create("hypo.G-.output.states")
    setwd("./hypo.G-.output.states")
    if(is.na(match("mean.expT", colnames(x)))==TRUE){
        getmeanexpT=function(geneID){
            mean(expDataT[as.character(geneID),], na.rm=T)}
        top=list(geneID=as.character(top_gene_IDs))
        x$mean.expT=mcmapply(getmeanexpT, top$geneID, mc.cores=cores)
    }
    for (i in 1:dim(x)[1]){
        datc=data.frame(meth=datM[,i], exp=datE[,i])
        datc$exp.changed=ifelse(datc$exp<x$mean.expT[i],1,0)
        rownames(datc)=rownames(datM)
        datc$FIN=ifelse(datc$meth==1 & datc$exp.changed==1,1,0)
        colnames(datc)=c(as.character(x$probe[i]), as.character(x$geneID[i]), "exp.changed","FIN")
        write.table(datc, paste(x$probe[i], x$geneID[i], x$geneSymbol[i], prefix, links, "states.txt", sep="."), row.names=T, col.names=T, sep="\t", quote=F)
    }
    TUMOR_Call=matrix("NA", nrow=dim(metDataT)[2], ncol=dim(x)[1])
    rownames(TUMOR_Call)=rownames(datc)
    colnames(TUMOR_Call)=paste(x$probe, x$geneID, x$geneSymbol,sep=".")
    ID=colnames(TUMOR_Call)
    findState=function(ID){
        LS=list.files(getwd(), pattern=as.character(ID))
        LF=read.table(LS, header=T)
        LF$FIN
    }
    TUMOR_Call=mcmapply(ID=ID,findState,mc.cores=cores)
    rownames(TUMOR_Call)=rownames(datc)
    write.table(TUMOR_Call, paste(prefix, links, "links.states.table.txt", sep="."), row.names=T, col.names=T, sep="\t", quote=F)
    # remove unnecessary files #
    notneed=dir(getwd(), pattern="states.txt$")
    file.remove(notneed)
}

make.circos.plots.hypermeth.G_pos <- function(circos_top_n_genes, prefix, cores) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Load BioCircos:
    # library('BioCircos')

    ## Load parallel:
    # library('parallel')

    ## Load htmlwidgets:
    # library('htmlwidgets')

    ## Import the gtf data for gencode v22 genes:
    gencode_v22_gtf <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Get info for the hg38 genes only:
    gencode_v22_genes <- gencode_v22_gtf[
        gencode_v22_gtf$type=='gene',
    ]

    ## Determine whether the start or the end coordinate is the
    ## TSS based on the strand:
    gencode_v22_genes$TSS <- ifelse(
        gencode_v22_genes$strand=='+',
        gencode_v22_genes$start,
        ifelse(
            gencode_v22_genes$strand=='-',
            gencode_v22_genes$end,
            ''
        )
    )

    ## remove the stuff after the period in the ENSG ids:
    gencode_v22_genes$ensembl_ID <- sub(
        '\\..*',
        '',
        gencode_v22_genes$gene_id
    )

    ## Load in the annotations for the 450k array probes on the hg38 genome:
    hg38_450k_probe_info <- read.delim(
        file='../scripts/data/hm450.hg38.manifest.tsv.gz',
        header= TRUE,
        stringsAsFactors = FALSE
    )

    ## Set rownames of hg38_450k_probe_info to be the probe IDs:
    rownames(hg38_450k_probe_info) <- hg38_450k_probe_info$probeID

    ## Set the number of genes to make circos plots for
    ## based on user input:
    circos_top_n_genes <- circos_top_n_genes

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G+.output.histogram/",
            prefix,
            ".hyper.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_names <- ordered_TFs_by_link_count[
        1:circos_top_n_genes,
        'geneID'
    ]

    ## Load the file with all linked CpGs:
    complete_linked_CpGs_to_TFs <- read.delim(
        file= '../step4/hyper.G+.output/hyper.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Load the hg38 chromosome size file:
    chrom_sizes <- read.delim(
        file= '../scripts/data/hg38_chrom_sizes.txt',
        header = FALSE,
        sep='\t',
        stringsAsFactors = FALSE
    )

    ## Ensure dataframe is ordered from largest to smallest chromosome:
    chrom_sizes <- chrom_sizes[
        order(chrom_sizes$V2, decreasing = TRUE),
    ]

    ## Get only the first 24 entries (extras are weird constructs):
    chrom_sizes <- chrom_sizes[1:24,]

    ## Set rownames of chrom sizes to be the chromosome names
    ## then delete the column with the chromosome names
    chrom_names <- chrom_sizes$V1

    ## Get only the first 24 entries (extras are weird constructs):
    ## This converts to vector
    chrom_sizes <- chrom_sizes[1:24,'V2']

    ## Remove the 'chr' from the start of the chrom_names
    ## to match the annotation used by BioCircos
    chrom_names <- sub('.*\\chr', '', chrom_names)

    ## Set chrom_names as the names for the chrom_sizes:
    names(chrom_sizes) <- chrom_names

    ## Remove X and Y chromosomes :
    chrom_sizes <- chrom_sizes[as.character(c(1:22))]

    ## Create a custom BioCircos hg38 genome build:
    hg38_genome <- setNames(
        as.list(chrom_sizes),
        names(chrom_sizes)
    )

    ## Write a function to generate and save a circos plot
    ## for a given TF to each of its CpGs linked:
    TF_linked_CpGs_circos_function <- function(
        TF_ENSG
    ){

        ## Get the name of the gene of interest from
        ## the ENSG:
        gene_name <- gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==TF_ENSG,
            'gene_name'
        ]

        ## Get a list of all the CpGs linked to the given TF:
        linked_CpGs <- unique(
            complete_linked_CpGs_to_TFs[
                complete_linked_CpGs_to_TFs$geneID==TF_ENSG,
                'probe'
            ]
        )

        ## List the promoter chromosome and start for TF
        ## a number of times equal to linked CpGs:
        circos_start_chromosome <- rep(
            as.character(
                gencode_v22_genes[
                    gencode_v22_genes$ensembl_ID==TF_ENSG,
                    'seqnames'
                ]
            ),
            length(linked_CpGs)
        )

        circos_start_position <- rep(
            as.numeric(
                gencode_v22_genes[
                    gencode_v22_genes$ensembl_ID==TF_ENSG,
                    'TSS'
                ]
            ),
            length(linked_CpGs)
        )

        ## List the chromosome and start for probes
        ## linked to given TF as end points on circos:
        circos_end_chromosome <- hg38_450k_probe_info[
            linked_CpGs,
            'CpG_chrm'
        ]

        circos_end_position <- (
            hg38_450k_probe_info[
                linked_CpGs,
                'CpG_beg'
            ]+1
        )

        ## Convert chromosomes to individual chromsome names
        ## without 'chr'
        circos_start_chromosome <- substring(
            circos_start_chromosome,
            4,
            nchar(circos_start_chromosome)
        )

        circos_end_chromosome <- substring(
            circos_end_chromosome,
            4,
            nchar(circos_end_chromosome)
        )

        ## Get the path to the step 5 folder to save
        ## the pdf in:
        path_to_folder <- getwd()

        ## Create a title for the circos plot html:
        circos_plot_html_title <- paste(
            path_to_folder,
            '/hyper.G+.output.circos/',
            gene_name,
            '_circos_plot.html',
            sep=''
        )

        ## Create a BioCircos tracklist:
        tracklist_TF <- BioCircosBackgroundTrack(
            "myBackgroundTrack",
            minRadius = 0,
            maxRadius = 0.95,
            borderSize = 0,
            fillColors = "#FFFFFF"
        )

        tracklist_TF <- tracklist_TF +
            BioCircosLinkTrack(
                'myLinkTrack',
                circos_start_chromosome,
                circos_start_position,
                circos_start_position + 10000,
                circos_end_chromosome,
                circos_end_position,
                circos_end_position + 10000,
                maxRadius = 0.95,
                color='red'
            )

        ## Plot the BioCircos track:
        x <- BioCircos(
            tracklist_TF,
            genomeFillColor = "PuOr",
            chrPad = 0.02,
            displayGenomeBorder = FALSE,
            yChr = TRUE,
            genomeTicksDisplay = FALSE,
            genomeLabelTextSize = "8pt",
            genomeLabelDy = 0,
            genome = hg38_genome
        )

        ## Save the circos plot html file:
        saveWidget(
            x,
            circos_plot_html_title
        )
    }

    ## Create survival plots for each of the genes designated:
    mclapply(
        X= top_gene_names,
        FUN= TF_linked_CpGs_circos_function,
        mc.cores= cores
    )
}

make.circos.plots.hypermeth.G_neg <- function(circos_top_n_genes, prefix, cores) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Load BioCircos:
    # library('BioCircos')

    ## Load parallel:
    # library('parallel')

    ## Load htmlwidgets:
    # library('htmlwidgets')

    ## Import the gtf data for gencode v22 genes:
    gencode_v22_gtf <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Get info for the hg38 genes only:
    gencode_v22_genes <- gencode_v22_gtf[
        gencode_v22_gtf$type=='gene',
    ]

    ## Determine whether the start or the end coordinate is the
    ## TSS based on the strand:
    gencode_v22_genes$TSS <- ifelse(
        gencode_v22_genes$strand=='+',
        gencode_v22_genes$start,
        ifelse(
            gencode_v22_genes$strand=='-',
            gencode_v22_genes$end,
            ''
        )
    )

    ## remove the stuff after the period in the ENSG ids:
    gencode_v22_genes$ensembl_ID <- sub(
        '\\..*',
        '',
        gencode_v22_genes$gene_id
    )

    ## Load in the annotations for the 450k array probes on the hg38 genome:
    hg38_450k_probe_info <- read.delim(
        file='../scripts/data/hm450.hg38.manifest.tsv.gz',
        header= TRUE,
        stringsAsFactors = FALSE
    )

    ## Set rownames of hg38_450k_probe_info to be the probe IDs:
    rownames(hg38_450k_probe_info) <- hg38_450k_probe_info$probeID

    ## Set the number of genes to make circos plots for
    ## based on user input:
    circos_top_n_genes <- circos_top_n_genes

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G-.output.histogram/",
            prefix,
            ".hyper.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_names <- ordered_TFs_by_link_count[
        1:circos_top_n_genes,
        'geneID'
    ]

    ## Load the file with all linked CpGs:
    complete_linked_CpGs_to_TFs <- read.delim(
        file= '../step4/hyper.G-.output/hyper.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Load the hg38 chromosome size file:
    chrom_sizes <- read.delim(
        file= '../scripts/data/hg38_chrom_sizes.txt',
        header = FALSE,
        sep='\t',
        stringsAsFactors = FALSE
    )

    ## Ensure dataframe is ordered from largest to smallest chromosome:
    chrom_sizes <- chrom_sizes[
        order(chrom_sizes$V2, decreasing = TRUE),
    ]

    ## Get only the first 24 entries (extras are weird constructs):
    chrom_sizes <- chrom_sizes[1:24,]

    ## Set rownames of chrom sizes to be the chromosome names
    ## then delete the column with the chromosome names
    chrom_names <- chrom_sizes$V1

    ## Get only the first 24 entries (extras are weird constructs):
    ## This converts to vector
    chrom_sizes <- chrom_sizes[1:24,'V2']

    ## Remove the 'chr' from the start of the chrom_names
    ## to match the annotation used by BioCircos
    chrom_names <- sub('.*\\chr', '', chrom_names)

    ## Set chrom_names as the names for the chrom_sizes:
    names(chrom_sizes) <- chrom_names

    ## Remove X and Y chromosomes :
    chrom_sizes <- chrom_sizes[as.character(c(1:22))]

    ## Create a custom BioCircos hg38 genome build:
    hg38_genome <- setNames(
        as.list(chrom_sizes),
        names(chrom_sizes)
    )

    ## Write a function to generate and save a circos plot
    ## for a given TF to each of its CpGs linked:
    TF_linked_CpGs_circos_function <- function(
        TF_ENSG
    ){

        ## Get the name of the gene of interest from
        ## the ENSG:
        gene_name <- gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==TF_ENSG,
            'gene_name'
        ]

        ## Get a list of all the CpGs linked to the given TF:
        linked_CpGs <- unique(
            complete_linked_CpGs_to_TFs[
                complete_linked_CpGs_to_TFs$geneID==TF_ENSG,
                'probe'
            ]
        )

        ## List the promoter chromosome and start for TF
        ## a number of times equal to linked CpGs:
        circos_start_chromosome <- rep(
            as.character(
                gencode_v22_genes[
                    gencode_v22_genes$ensembl_ID==TF_ENSG,
                    'seqnames'
                ]
            ),
            length(linked_CpGs)
        )

        circos_start_position <- rep(
            as.numeric(
                gencode_v22_genes[
                    gencode_v22_genes$ensembl_ID==TF_ENSG,
                    'TSS'
                ]
            ),
            length(linked_CpGs)
        )

        ## List the chromosome and start for probes
        ## linked to given TF as end points on circos:
        circos_end_chromosome <- hg38_450k_probe_info[
            linked_CpGs,
            'CpG_chrm'
        ]

        circos_end_position <- (
            hg38_450k_probe_info[
                linked_CpGs,
                'CpG_beg'
            ]+1
        )

        ## Convert chromosomes to individual chromsome names
        ## without 'chr'
        circos_start_chromosome <- substring(
            circos_start_chromosome,
            4,
            nchar(circos_start_chromosome)
        )

        circos_end_chromosome <- substring(
            circos_end_chromosome,
            4,
            nchar(circos_end_chromosome)
        )

        ## Get the path to the step 5 folder to save
        ## the pdf in:
        path_to_folder <- getwd()

        ## Create a title for the circos plot html:
        circos_plot_html_title <- paste(
            path_to_folder,
            '/hyper.G-.output.circos/',
            gene_name,
            '_circos_plot.html',
            sep=''
        )

        ## Create a BioCircos tracklist:
        tracklist_TF <- BioCircosBackgroundTrack(
            "myBackgroundTrack",
            minRadius = 0,
            maxRadius = 0.95,
            borderSize = 0,
            fillColors = "#FFFFFF"
        )

        tracklist_TF <- tracklist_TF +
            BioCircosLinkTrack(
                'myLinkTrack',
                circos_start_chromosome,
                circos_start_position,
                circos_start_position + 10000,
                circos_end_chromosome,
                circos_end_position,
                circos_end_position + 10000,
                maxRadius = 0.95,
                color='red'
            )

        ## Plot the BioCircos track:
        x <- BioCircos(
            tracklist_TF,
            genomeFillColor = "PuOr",
            chrPad = 0.02,
            displayGenomeBorder = FALSE,
            yChr = TRUE,
            genomeTicksDisplay = FALSE,
            genomeLabelTextSize = "8pt",
            genomeLabelDy = 0,
            genome = hg38_genome
        )

        ## Save the circos plot html file:
        saveWidget(
            x,
            circos_plot_html_title
        )
    }

    ## Create survival plots for each of the genes designated:
    mclapply(
        X= top_gene_names,
        FUN= TF_linked_CpGs_circos_function,
        mc.cores= cores
    )
}

make.circos.plots.hypometh.G_pos <- function(circos_top_n_genes, prefix, cores) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Load BioCircos:
    # library('BioCircos')

    ## Load parallel:
    # library('parallel')

    ## Load htmlwidgets:
    # library('htmlwidgets')

    ## Import the gtf data for gencode v22 genes:
    gencode_v22_gtf <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Get info for the hg38 genes only:
    gencode_v22_genes <- gencode_v22_gtf[
        gencode_v22_gtf$type=='gene',
    ]

    ## Determine whether the start or the end coordinate is the
    ## TSS based on the strand:
    gencode_v22_genes$TSS <- ifelse(
        gencode_v22_genes$strand=='+',
        gencode_v22_genes$start,
        ifelse(
            gencode_v22_genes$strand=='-',
            gencode_v22_genes$end,
            ''
        )
    )

    ## remove the stuff after the period in the ENSG ids:
    gencode_v22_genes$ensembl_ID <- sub(
        '\\..*',
        '',
        gencode_v22_genes$gene_id
    )

    ## Load in the annotations for the 450k array probes on the hg38 genome:
    hg38_450k_probe_info <- read.delim(
        file='../scripts/data/hm450.hg38.manifest.tsv.gz',
        header= TRUE,
        stringsAsFactors = FALSE
    )

    ## Set rownames of hg38_450k_probe_info to be the probe IDs:
    rownames(hg38_450k_probe_info) <- hg38_450k_probe_info$probeID

    ## Set the number of genes to make circos plots for
    ## based on user input:
    circos_top_n_genes <- circos_top_n_genes

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hypo.G+.output.histogram/",
            prefix,
            ".hypo.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_names <- ordered_TFs_by_link_count[
        1:circos_top_n_genes,
        'geneID'
    ]

    ## Load the file with all linked CpGs:
    complete_linked_CpGs_to_TFs <- read.delim(
        file= '../step4/hypo.G+.output/hypo.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Load the hg38 chromosome size file:
    chrom_sizes <- read.delim(
        file= '../scripts/data/hg38_chrom_sizes.txt',
        header = FALSE,
        sep='\t',
        stringsAsFactors = FALSE
    )

    ## Ensure dataframe is ordered from largest to smallest chromosome:
    chrom_sizes <- chrom_sizes[
        order(chrom_sizes$V2, decreasing = TRUE),
    ]

    ## Get only the first 24 entries (extras are weird constructs):
    chrom_sizes <- chrom_sizes[1:24,]

    ## Set rownames of chrom sizes to be the chromosome names
    ## then delete the column with the chromosome names
    chrom_names <- chrom_sizes$V1

    ## Get only the first 24 entries (extras are weird constructs):
    ## This converts to vector
    chrom_sizes <- chrom_sizes[1:24,'V2']

    ## Remove the 'chr' from the start of the chrom_names
    ## to match the annotation used by BioCircos
    chrom_names <- sub('.*\\chr', '', chrom_names)

    ## Set chrom_names as the names for the chrom_sizes:
    names(chrom_sizes) <- chrom_names

    ## Remove X and Y chromosomes :
    chrom_sizes <- chrom_sizes[as.character(c(1:22))]

    ## Create a custom BioCircos hg38 genome build:
    hg38_genome <- setNames(
        as.list(chrom_sizes),
        names(chrom_sizes)
    )

    ## Write a function to generate and save a circos plot
    ## for a given TF to each of its CpGs linked:
    TF_linked_CpGs_circos_function <- function(
        TF_ENSG
    ){

        ## Get the name of the gene of interest from
        ## the ENSG:
        gene_name <- gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==TF_ENSG,
            'gene_name'
        ]

        ## Get a list of all the CpGs linked to the given TF:
        linked_CpGs <- unique(
            complete_linked_CpGs_to_TFs[
                complete_linked_CpGs_to_TFs$geneID==TF_ENSG,
                'probe'
            ]
        )

        ## List the promoter chromosome and start for TF
        ## a number of times equal to linked CpGs:
        circos_start_chromosome <- rep(
            as.character(
                gencode_v22_genes[
                    gencode_v22_genes$ensembl_ID==TF_ENSG,
                    'seqnames'
                ]
            ),
            length(linked_CpGs)
        )

        circos_start_position <- rep(
            as.numeric(
                gencode_v22_genes[
                    gencode_v22_genes$ensembl_ID==TF_ENSG,
                    'TSS'
                ]
            ),
            length(linked_CpGs)
        )

        ## List the chromosome and start for probes
        ## linked to given TF as end points on circos:
        circos_end_chromosome <- hg38_450k_probe_info[
            linked_CpGs,
            'CpG_chrm'
        ]

        circos_end_position <- (
            hg38_450k_probe_info[
                linked_CpGs,
                'CpG_beg'
            ]+1
        )

        ## Convert chromosomes to individual chromsome names
        ## without 'chr'
        circos_start_chromosome <- substring(
            circos_start_chromosome,
            4,
            nchar(circos_start_chromosome)
        )

        circos_end_chromosome <- substring(
            circos_end_chromosome,
            4,
            nchar(circos_end_chromosome)
        )

        ## Get the path to the step 5 folder to save
        ## the pdf in:
        path_to_folder <- getwd()

        ## Create a title for the circos plot html:
        circos_plot_html_title <- paste(
            path_to_folder,
            '/hypo.G+.output.circos/',
            gene_name,
            '_circos_plot.html',
            sep=''
        )

        ## Create a BioCircos tracklist:
        tracklist_TF <- BioCircosBackgroundTrack(
            "myBackgroundTrack",
            minRadius = 0,
            maxRadius = 0.95,
            borderSize = 0,
            fillColors = "#FFFFFF"
        )

        tracklist_TF <- tracklist_TF +
            BioCircosLinkTrack(
                'myLinkTrack',
                circos_start_chromosome,
                circos_start_position,
                circos_start_position + 10000,
                circos_end_chromosome,
                circos_end_position,
                circos_end_position + 10000,
                maxRadius = 0.95,
                color='red'
            )

        ## Plot the BioCircos track:
        x <- BioCircos(
            tracklist_TF,
            genomeFillColor = "PuOr",
            chrPad = 0.02,
            displayGenomeBorder = FALSE,
            yChr = TRUE,
            genomeTicksDisplay = FALSE,
            genomeLabelTextSize = "8pt",
            genomeLabelDy = 0,
            genome = hg38_genome
        )

        ## Save the circos plot html file:
        saveWidget(
            x,
            circos_plot_html_title
        )
    }

    ## Create survival plots for each of the genes designated:
    mclapply(
        X= top_gene_names,
        FUN= TF_linked_CpGs_circos_function,
        mc.cores= cores
    )
}

make.circos.plots.hypometh.G_neg <- function(circos_top_n_genes, prefix, cores) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Load BioCircos:
    # library('BioCircos')

    ## Load parallel:
    # library('parallel')

    ## Load htmlwidgets:
    # library('htmlwidgets')

    ## Import the gtf data for gencode v22 genes:
    gencode_v22_gtf <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Get info for the hg38 genes only:
    gencode_v22_genes <- gencode_v22_gtf[
        gencode_v22_gtf$type=='gene',
    ]

    ## Determine whether the start or the end coordinate is the
    ## TSS based on the strand:
    gencode_v22_genes$TSS <- ifelse(
        gencode_v22_genes$strand=='+',
        gencode_v22_genes$start,
        ifelse(
            gencode_v22_genes$strand=='-',
            gencode_v22_genes$end,
            ''
        )
    )

    ## remove the stuff after the period in the ENSG ids:
    gencode_v22_genes$ensembl_ID <- sub(
        '\\..*',
        '',
        gencode_v22_genes$gene_id
    )

    ## Load in the annotations for the 450k array probes on the hg38 genome:
    hg38_450k_probe_info <- read.delim(
        file='../scripts/data/hm450.hg38.manifest.tsv.gz',
        header= TRUE,
        stringsAsFactors = FALSE
    )

    ## Set rownames of hg38_450k_probe_info to be the probe IDs:
    rownames(hg38_450k_probe_info) <- hg38_450k_probe_info$probeID

    ## Set the number of genes to make circos plots for
    ## based on user input:
    circos_top_n_genes <- circos_top_n_genes

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hypo.G-.output.histogram/",
            prefix,
            ".hypo.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_names <- ordered_TFs_by_link_count[
        1:circos_top_n_genes,
        'geneID'
    ]

    ## Load the file with all linked CpGs:
    complete_linked_CpGs_to_TFs <- read.delim(
        file= '../step4/hypo.G-.output/hypo.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Load the hg38 chromosome size file:
    chrom_sizes <- read.delim(
        file= '../scripts/data/hg38_chrom_sizes.txt',
        header = FALSE,
        sep='\t',
        stringsAsFactors = FALSE
    )

    ## Ensure dataframe is ordered from largest to smallest chromosome:
    chrom_sizes <- chrom_sizes[
        order(chrom_sizes$V2, decreasing = TRUE),
    ]

    ## Get only the first 24 entries (extras are weird constructs):
    chrom_sizes <- chrom_sizes[1:24,]

    ## Set rownames of chrom sizes to be the chromosome names
    ## then delete the column with the chromosome names
    chrom_names <- chrom_sizes$V1

    ## Get only the first 24 entries (extras are weird constructs):
    ## This converts to vector
    chrom_sizes <- chrom_sizes[1:24,'V2']

    ## Remove the 'chr' from the start of the chrom_names
    ## to match the annotation used by BioCircos
    chrom_names <- sub('.*\\chr', '', chrom_names)

    ## Set chrom_names as the names for the chrom_sizes:
    names(chrom_sizes) <- chrom_names

    ## Remove X and Y chromosomes :
    chrom_sizes <- chrom_sizes[as.character(c(1:22))]

    ## Create a custom BioCircos hg38 genome build:
    hg38_genome <- setNames(
        as.list(chrom_sizes),
        names(chrom_sizes)
    )

    ## Write a function to generate and save a circos plot
    ## for a given TF to each of its CpGs linked:
    TF_linked_CpGs_circos_function <- function(
        TF_ENSG
    ){

        ## Get the name of the gene of interest from
        ## the ENSG:
        gene_name <- gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==TF_ENSG,
            'gene_name'
        ]

        ## Get a list of all the CpGs linked to the given TF:
        linked_CpGs <- unique(
            complete_linked_CpGs_to_TFs[
                complete_linked_CpGs_to_TFs$geneID==TF_ENSG,
                'probe'
            ]
        )

        ## List the promoter chromosome and start for TF
        ## a number of times equal to linked CpGs:
        circos_start_chromosome <- rep(
            as.character(
                gencode_v22_genes[
                    gencode_v22_genes$ensembl_ID==TF_ENSG,
                    'seqnames'
                ]
            ),
            length(linked_CpGs)
        )

        circos_start_position <- rep(
            as.numeric(
                gencode_v22_genes[
                    gencode_v22_genes$ensembl_ID==TF_ENSG,
                    'TSS'
                ]
            ),
            length(linked_CpGs)
        )

        ## List the chromosome and start for probes
        ## linked to given TF as end points on circos:
        circos_end_chromosome <- hg38_450k_probe_info[
            linked_CpGs,
            'CpG_chrm'
        ]

        circos_end_position <- (
            hg38_450k_probe_info[
                linked_CpGs,
                'CpG_beg'
            ]+1
        )

        ## Convert chromosomes to individual chromsome names
        ## without 'chr'
        circos_start_chromosome <- substring(
            circos_start_chromosome,
            4,
            nchar(circos_start_chromosome)
        )

        circos_end_chromosome <- substring(
            circos_end_chromosome,
            4,
            nchar(circos_end_chromosome)
        )

        ## Get the path to the step 5 folder to save
        ## the pdf in:
        path_to_folder <- getwd()

        ## Create a title for the circos plot html:
        circos_plot_html_title <- paste(
            path_to_folder,
            '/hypo.G-.output.circos/',
            gene_name,
            '_circos_plot.html',
            sep=''
        )

        ## Create a BioCircos tracklist:
        tracklist_TF <- BioCircosBackgroundTrack(
            "myBackgroundTrack",
            minRadius = 0,
            maxRadius = 0.95,
            borderSize = 0,
            fillColors = "#FFFFFF"
        )

        tracklist_TF <- tracklist_TF +
            BioCircosLinkTrack(
                'myLinkTrack',
                circos_start_chromosome,
                circos_start_position,
                circos_start_position + 10000,
                circos_end_chromosome,
                circos_end_position,
                circos_end_position + 10000,
                maxRadius = 0.95,
                color='red'
            )

        ## Plot the BioCircos track:
        x <- BioCircos(
            tracklist_TF,
            genomeFillColor = "PuOr",
            chrPad = 0.02,
            displayGenomeBorder = FALSE,
            yChr = TRUE,
            genomeTicksDisplay = FALSE,
            genomeLabelTextSize = "8pt",
            genomeLabelDy = 0,
            genome = hg38_genome
        )

        ## Save the circos plot html file:
        saveWidget(
            x,
            circos_plot_html_title
        )
    }

    ## Create survival plots for each of the genes designated:
    mclapply(
        X= top_gene_names,
        FUN= TF_linked_CpGs_circos_function,
        mc.cores= cores
    )
}

make.CNV.scatterplots.hypermeth.G_pos.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hyper.G+.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hyper.G+.output/", LS[1], sep=""), header=T)
    LC=list.files("../external.data/otherinfo", "CNV")
    CNV=read.delim(paste("../external.data/otherinfo/", LC[1], sep=""), sep="\t", header=T, check.names=F)
    dir.create("hyper.G+.output.complex.scatterplot")
    setwd("./hyper.G+.output.complex.scatterplot")
    M1=match(unique(TESTSR$geneSymbol),rownames(CNV))
    CNV.c=CNV[na.omit(M1),]
    M3=match(colnames(metDataT), colnames(CNV.c))
    CNV.c2=as.matrix(CNV.c)[,M3]
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(CNV.c2)))==FALSE){
            CNV.t=CNV.c2[as.character(TESTSR[i,2]),]
            CNV.t2=t(CNV.t)
            CNV.t3=c(CNV.t2, rep(0,dim(metDataN)[2]))
            CNV.sh=CNV.t3
            CNV.sh2=CNV.sh
            CNV.sh2[which(CNV.sh==0)]=paste(16)
            CNV.sh2[which(CNV.sh==-1)]=paste(17)
            CNV.sh2[which(CNV.sh==-2)]=paste(24)
            CNV.sh2[which(CNV.sh==1)]=paste(15)
            CNV.sh2[which(CNV.sh==2)]=paste(22)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.CNV.scatterplots.hypermeth.G_neg.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hyper.G-.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hyper.G-.output/", LS[1], sep=""), header=T)
    LC=list.files("../external.data/otherinfo", "CNV")
    CNV=read.delim(paste("../external.data/otherinfo/", LC[1], sep=""), sep="\t", header=T, check.names=F)
    dir.create("hyper.G-.output.complex.scatterplot")
    setwd("./hyper.G-.output.complex.scatterplot")
    M1=match(unique(TESTSR$geneSymbol),rownames(CNV))
    CNV.c=CNV[na.omit(M1),]
    M3=match(colnames(metDataT), colnames(CNV.c))
    CNV.c2=as.matrix(CNV.c)[,M3]
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(CNV.c2)))==FALSE){
            CNV.t=CNV.c2[as.character(TESTSR[i,2]),]
            CNV.t2=t(CNV.t)
            CNV.t3=c(CNV.t2, rep(0,dim(metDataN)[2]))
            CNV.sh=CNV.t3
            CNV.sh2=CNV.sh
            CNV.sh2[which(CNV.sh==0)]=paste(16)
            CNV.sh2[which(CNV.sh==-1)]=paste(17)
            CNV.sh2[which(CNV.sh==-2)]=paste(24)
            CNV.sh2[which(CNV.sh==1)]=paste(15)
            CNV.sh2[which(CNV.sh==2)]=paste(22)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.CNV.scatterplots.hypometh.G_pos.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hypo.G+.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hypo.G+.output/", LS[1], sep=""), header=T)
    LC=list.files("../external.data/otherinfo", "CNV")
    CNV=read.delim(paste("../external.data/otherinfo/", LC[1], sep=""), sep="\t", header=T, check.names=F)
    dir.create("hypo.G+.output.complex.scatterplot")
    setwd("./hypo.G+.output.complex.scatterplot")
    M1=match(unique(TESTSR$geneSymbol),rownames(CNV))
    CNV.c=CNV[na.omit(M1),]
    M3=match(colnames(metDataT), colnames(CNV.c))
    CNV.c2=as.matrix(CNV.c)[,M3]
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(CNV.c2)))==FALSE){
            CNV.t=CNV.c2[as.character(TESTSR[i,2]),]
            CNV.t2=t(CNV.t)
            CNV.t3=c(CNV.t2, rep(0,dim(metDataN)[2]))
            CNV.sh=CNV.t3
            CNV.sh2=CNV.sh
            CNV.sh2[which(CNV.sh==0)]=paste(16)
            CNV.sh2[which(CNV.sh==-1)]=paste(17)
            CNV.sh2[which(CNV.sh==-2)]=paste(24)
            CNV.sh2[which(CNV.sh==1)]=paste(15)
            CNV.sh2[which(CNV.sh==2)]=paste(22)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.CNV.scatterplots.hypometh.G_neg.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hypo.G-.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hypo.G-.output/", LS[1], sep=""), header=T)
    LC=list.files("../external.data/otherinfo", "CNV")
    CNV=read.delim(paste("../external.data/otherinfo/", LC[1], sep=""), sep="\t", header=T, check.names=F)
    dir.create("hypo.G-.output.complex.scatterplot")
    setwd("./hypo.G-.output.complex.scatterplot")
    M1=match(unique(TESTSR$geneSymbol),rownames(CNV))
    CNV.c=CNV[na.omit(M1),]
    M3=match(colnames(metDataT), colnames(CNV.c))
    CNV.c2=as.matrix(CNV.c)[,M3]
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(CNV.c2)))==FALSE){
            CNV.t=CNV.c2[as.character(TESTSR[i,2]),]
            CNV.t2=t(CNV.t)
            CNV.t3=c(CNV.t2, rep(0,dim(metDataN)[2]))
            CNV.sh=CNV.t3
            CNV.sh2=CNV.sh
            CNV.sh2[which(CNV.sh==0)]=paste(16)
            CNV.sh2[which(CNV.sh==-1)]=paste(17)
            CNV.sh2[which(CNV.sh==-2)]=paste(24)
            CNV.sh2[which(CNV.sh==1)]=paste(15)
            CNV.sh2[which(CNV.sh==2)]=paste(22)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.complex.scatterplots.hypermeth.G_pos.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hyper.G+.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hyper.G+.output/", LS[1], sep=""), header=T)
    LC=list.files("../external.data/otherinfo", "CNV")
    CNV=read.delim(paste("../external.data/otherinfo/", LC[1], sep=""), sep="\t", header=T, check.names=F)
    LM=list.files("../external.data/otherinfo", "SM")
    SM=read.delim(paste("../external.data/otherinfo/", LM[1], sep=""), header=T, sep="\t")
    LP=list.files("../external.data/otherinfo", "Purity")
    Pur=read.delim(paste("../external.data/otherinfo/", LP[1], sep=""), header=T, sep="\t")
    dir.create("hyper.G+.output.complex.scatterplot")
    setwd("./hyper.G+.output.complex.scatterplot")
    SM=t(SM)
    M1=match(unique(TESTSR$geneSymbol),rownames(CNV))
    CNV.c=CNV[na.omit(M1),]
    M2=match(unique(TESTSR$geneSymbol),rownames(SM))
    SM.c=SM[na.omit(M2),]
    M3=match(colnames(metDataT), colnames(CNV.c))
    M4=match(colnames(metDataT), colnames(SM.c))
    CNV.c2=as.matrix(CNV.c)[,M3]
    SM.c2=as.matrix(SM.c)[,M4]
    M5=match(colnames(metDataT), Pur$ID)
    Purity.c=Pur[M5,]
    Purity.c$score=(1+2*(Purity.c$purity))
    Purity.c$score[which(is.na(Purity.c$score)==TRUE)]=paste(1)
    Purity=c(Purity.c$score, rep(3, dim(metDataN)[2]))
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(CNV.c2)))==FALSE){
            CNV.t=CNV.c2[as.character(TESTSR[i,2]),]
            CNV.t2=t(CNV.t)
            CNV.t3=c(CNV.t2, rep(0,dim(metDataN)[2]))
            CNV.sh=CNV.t3
            CNV.sh2=CNV.sh
            CNV.sh2[which(CNV.sh==0)]=paste(16)
            CNV.sh2[which(CNV.sh==-1)]=paste(17)
            CNV.sh2[which(CNV.sh==-2)]=paste(24)
            CNV.sh2[which(CNV.sh==1)]=paste(15)
            CNV.sh2[which(CNV.sh==2)]=paste(22)
        }
        SM.sh=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(SM.c2)))==FALSE){
            SM.t=SM.c2[as.character(TESTSR[i,2]),]
            SM.t2=t(SM.t)
            SM.t3=c(SM.t2, rep(0,dim(metDataN)[2]))
            SM.sh=SM.t3
            SM.sh2=SM.sh
            SM.sh2[which(SM.sh==0)]=paste(16)
            SM.sh2[which(SM.sh==1)]=paste(8)
            CNV.sh2[which(SM.sh2==8)]=paste(8)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=as.numeric(Purity), shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.complex.scatterplots.hypermeth.G_neg.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hyper.G-.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hyper.G-.output/", LS[1], sep=""), header=T)
    LC=list.files("../external.data/otherinfo", "CNV")
    CNV=read.delim(paste("../external.data/otherinfo/", LC[1], sep=""), sep="\t", header=T, check.names=F)
    LM=list.files("../external.data/otherinfo", "SM")
    SM=read.delim(paste("../external.data/otherinfo/", LM[1], sep=""), header=T, sep="\t")
    LP=list.files("../external.data/otherinfo", "Purity")
    Pur=read.delim(paste("../external.data/otherinfo/", LP[1], sep=""), header=T, sep="\t")
    dir.create("hyper.G-.output.complex.scatterplot")
    setwd("./hyper.G-.output.complex.scatterplot")
    SM=t(SM)
    M1=match(unique(TESTSR$geneSymbol),rownames(CNV))
    CNV.c=CNV[na.omit(M1),]
    M2=match(unique(TESTSR$geneSymbol),rownames(SM))
    SM.c=SM[na.omit(M2),]
    M3=match(colnames(metDataT), colnames(CNV.c))
    M4=match(colnames(metDataT), colnames(SM.c))
    CNV.c2=as.matrix(CNV.c)[,M3]
    SM.c2=as.matrix(SM.c)[,M4]
    M5=match(colnames(metDataT), Pur$ID)
    Purity.c=Pur[M5,]
    Purity.c$score=(1+2*(Purity.c$purity))
    Purity.c$score[which(is.na(Purity.c$score)==TRUE)]=paste(1)
    Purity=c(Purity.c$score, rep(3, dim(metDataN)[2]))
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(CNV.c2)))==FALSE){
            CNV.t=CNV.c2[as.character(TESTSR[i,2]),]
            CNV.t2=t(CNV.t)
            CNV.t3=c(CNV.t2, rep(0,dim(metDataN)[2]))
            CNV.sh=CNV.t3
            CNV.sh2=CNV.sh
            CNV.sh2[which(CNV.sh==0)]=paste(16)
            CNV.sh2[which(CNV.sh==-1)]=paste(17)
            CNV.sh2[which(CNV.sh==-2)]=paste(24)
            CNV.sh2[which(CNV.sh==1)]=paste(15)
            CNV.sh2[which(CNV.sh==2)]=paste(22)
        }
        SM.sh=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(SM.c2)))==FALSE){
            SM.t=SM.c2[as.character(TESTSR[i,2]),]
            SM.t2=t(SM.t)
            SM.t3=c(SM.t2, rep(0,dim(metDataN)[2]))
            SM.sh=SM.t3
            SM.sh2=SM.sh
            SM.sh2[which(SM.sh==0)]=paste(16)
            SM.sh2[which(SM.sh==1)]=paste(8)
            CNV.sh2[which(SM.sh2==8)]=paste(8)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=as.numeric(Purity), shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.complex.scatterplots.hypometh.G_pos.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hypo.G+.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hypo.G+.output/", LS[1], sep=""), header=T)
    LC=list.files("../external.data/otherinfo", "CNV")
    CNV=read.delim(paste("../external.data/otherinfo/", LC[1], sep=""), sep="\t", header=T, check.names=F)
    LM=list.files("../external.data/otherinfo", "SM")
    SM=read.delim(paste("../external.data/otherinfo/", LM[1], sep=""), header=T, sep="\t")
    LP=list.files("../external.data/otherinfo", "Purity")
    Pur=read.delim(paste("../external.data/otherinfo/", LP[1], sep=""), header=T, sep="\t")
    dir.create("hypo.G+.output.complex.scatterplot")
    setwd("./hypo.G+.output.complex.scatterplot")
    SM=t(SM)
    M1=match(unique(TESTSR$geneSymbol),rownames(CNV))
    CNV.c=CNV[na.omit(M1),]
    M2=match(unique(TESTSR$geneSymbol),rownames(SM))
    SM.c=SM[na.omit(M2),]
    M3=match(colnames(metDataT), colnames(CNV.c))
    M4=match(colnames(metDataT), colnames(SM.c))
    CNV.c2=as.matrix(CNV.c)[,M3]
    SM.c2=as.matrix(SM.c)[,M4]
    M5=match(colnames(metDataT), Pur$ID)
    Purity.c=Pur[M5,]
    Purity.c$score=(1+2*(Purity.c$purity))
    Purity.c$score[which(is.na(Purity.c$score)==TRUE)]=paste(1)
    Purity=c(Purity.c$score, rep(3, dim(metDataN)[2]))
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(CNV.c2)))==FALSE){
            CNV.t=CNV.c2[as.character(TESTSR[i,2]),]
            CNV.t2=t(CNV.t)
            CNV.t3=c(CNV.t2, rep(0,dim(metDataN)[2]))
            CNV.sh=CNV.t3
            CNV.sh2=CNV.sh
            CNV.sh2[which(CNV.sh==0)]=paste(16)
            CNV.sh2[which(CNV.sh==-1)]=paste(17)
            CNV.sh2[which(CNV.sh==-2)]=paste(24)
            CNV.sh2[which(CNV.sh==1)]=paste(15)
            CNV.sh2[which(CNV.sh==2)]=paste(22)
        }
        SM.sh=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(SM.c2)))==FALSE){
            SM.t=SM.c2[as.character(TESTSR[i,2]),]
            SM.t2=t(SM.t)
            SM.t3=c(SM.t2, rep(0,dim(metDataN)[2]))
            SM.sh=SM.t3
            SM.sh2=SM.sh
            SM.sh2[which(SM.sh==0)]=paste(16)
            SM.sh2[which(SM.sh==1)]=paste(8)
            CNV.sh2[which(SM.sh2==8)]=paste(8)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=as.numeric(Purity), shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.complex.scatterplots.hypometh.G_neg.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hypo.G-.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hypo.G-.output/", LS[1], sep=""), header=T)
    LC=list.files("../external.data/otherinfo", "CNV")
    CNV=read.delim(paste("../external.data/otherinfo/", LC[1], sep=""), sep="\t", header=T, check.names=F)
    LM=list.files("../external.data/otherinfo", "SM")
    SM=read.delim(paste("../external.data/otherinfo/", LM[1], sep=""), header=T, sep="\t")
    LP=list.files("../external.data/otherinfo", "Purity")
    Pur=read.delim(paste("../external.data/otherinfo/", LP[1], sep=""), header=T, sep="\t")
    dir.create("hypo.G-.output.complex.scatterplot")
    setwd("./hypo.G-.output.complex.scatterplot")
    SM=t(SM)
    M1=match(unique(TESTSR$geneSymbol),rownames(CNV))
    CNV.c=CNV[na.omit(M1),]
    M2=match(unique(TESTSR$geneSymbol),rownames(SM))
    SM.c=SM[na.omit(M2),]
    M3=match(colnames(metDataT), colnames(CNV.c))
    M4=match(colnames(metDataT), colnames(SM.c))
    CNV.c2=as.matrix(CNV.c)[,M3]
    SM.c2=as.matrix(SM.c)[,M4]
    M5=match(colnames(metDataT), Pur$ID)
    Purity.c=Pur[M5,]
    Purity.c$score=(1+2*(Purity.c$purity))
    Purity.c$score[which(is.na(Purity.c$score)==TRUE)]=paste(1)
    Purity=c(Purity.c$score, rep(3, dim(metDataN)[2]))
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(CNV.c2)))==FALSE){
            CNV.t=CNV.c2[as.character(TESTSR[i,2]),]
            CNV.t2=t(CNV.t)
            CNV.t3=c(CNV.t2, rep(0,dim(metDataN)[2]))
            CNV.sh=CNV.t3
            CNV.sh2=CNV.sh
            CNV.sh2[which(CNV.sh==0)]=paste(16)
            CNV.sh2[which(CNV.sh==-1)]=paste(17)
            CNV.sh2[which(CNV.sh==-2)]=paste(24)
            CNV.sh2[which(CNV.sh==1)]=paste(15)
            CNV.sh2[which(CNV.sh==2)]=paste(22)
        }
        SM.sh=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(SM.c2)))==FALSE){
            SM.t=SM.c2[as.character(TESTSR[i,2]),]
            SM.t2=t(SM.t)
            SM.t3=c(SM.t2, rep(0,dim(metDataN)[2]))
            SM.sh=SM.t3
            SM.sh2=SM.sh
            SM.sh2[which(SM.sh==0)]=paste(16)
            SM.sh2[which(SM.sh==1)]=paste(8)
            CNV.sh2[which(SM.sh2==8)]=paste(8)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=as.numeric(Purity), shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.histogram.of.results.real.hypermeth.G_pos.source <- function(prefix, histcol) {
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
}

make.histogram.of.results.real.hypermeth.G_neg.source <- function(prefix, histcol) {
    load("../settings.rda")
    links=c("hyper.G-")
    dir.create("hyper.G-.output.histogram")
    setwd("./hyper.G-.output.histogram")
    all=read.delim(file=paste("../hyper.G-.output/", prefix, ".", links, ".", "links.hg38.all.anno.txt", sep=""), header=T)
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
}

make.histogram.of.results.real.hypometh.G_pos.source <- function(prefix, histcol) {
    load("../settings.rda")
    links=c("hypo.G+")
    dir.create("hypo.G+.output.histogram")
    setwd("./hypo.G+.output.histogram")
    all=read.delim(file=paste("../hypo.G+.output/", prefix, ".", links, ".", "links.hg38.all.anno.txt", sep=""), header=T)
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
}

make.histogram.of.results.real.hypometh.G_neg.source <- function(prefix, histcol) {
    load("../settings.rda")
    links=c("hypo.G-")
    dir.create("hypo.G-.output.histogram")
    setwd("./hypo.G-.output.histogram")
    all=read.delim(file=paste("../hypo.G-.output/", prefix, ".", links, ".", "links.hg38.all.anno.txt", sep=""), header=T)
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
}

make.met.heatmap.hypermeth.G_pos <- function(prefix, probe_heatmap_top_n_genes, metDataT, matlab, jet.colors, expDataT) {
    ## Load heatmap.3 function
    source("../scripts/data/heatmap.3.R")

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G+.output.histogram/",
            prefix,
            ".hyper.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hyper.G+.output/hyper.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_IDs <- ordered_TFs_by_link_count[
        1:probe_heatmap_top_n_genes,
        'geneID'
    ]

    top_gene_names <- gencode_v22_genes_df[
        top_gene_IDs,
        'gene_name'
    ]

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes selected:
    probes_linked_to_significant_genes <- unique(
        CpG_linkage_dataset[
            CpG_linkage_dataset$geneID %in% top_gene_IDs,
            'probe'
        ]
    )

    ## Get the rda file with expression and methylation data:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## Load the RDA files
    load(combined_rda_file)

    ## Create a methylation dataset from the tumor samples of interest:
    probes_methylation <- as.matrix(
        metDataT[
            probes_linked_to_significant_genes,
        ]
    )

    ## Create function to convert numbers to jet color vectors:
    jet_color_function <- function(numeric_values){

        jet_color_numeric_color_grad <- matlab::jet.colors(200)

        color_values <- jet_color_numeric_color_grad[numeric_values]

        return(color_values)
    }

    ## Write function to get expression values for the genes of interest:
    tumor_expression_grabber <- function(gene_id){

        unlist(
            c(
                expDataT[
                    c(gene_id),
                ]
            )
        )

    }

    ## Get expression from tumor samples for all the genes of interest:
    top_gene_tumor_expression <- sapply(
        top_gene_IDs,
        tumor_expression_grabber
    )

    ## Write new Rescale function for use later esp. with SCML2
    rescale_func_zero_ignored <- function(x){

        ## get minimum non-zero, non-NA value:
        non_zero_vec <- x[x!=0]

        non_zero_vec <- non_zero_vec[is.na(non_zero_vec)!=TRUE]

        ## get minimum of non-zero values:
        used_min <- min(non_zero_vec)
        used_max <- max(non_zero_vec)

        place <- ((x-used_min)/(used_max-used_min))*200

        return_value <- ifelse(
            place<=0,
            0.001,
            place
        )

        return(
            ceiling(return_value)
        )
    }

    ## Rescale the expression values of the genes of interest:
    rescaled_top_gene_tumor_expression <- apply(
        top_gene_tumor_expression,
        2,
        rescale_func_zero_ignored
    )

    ## Convert the rescaled values to jetcolor:
    rescaled_top_gene_tumor_expression_jet_color <- apply(
        rescaled_top_gene_tumor_expression,
        2,
        jet_color_function
    )

    rownames(rescaled_top_gene_tumor_expression_jet_color) <- rownames(rescaled_top_gene_tumor_expression)

    ## Create row labels for the histogram
    ## (blank for now):
    Tumor_row_color_labels <- rbind(
        rep('white', nrow(probes_methylation))
    )
    rownames(Tumor_row_color_labels) <- c(
        ""
    )

    ## Create column labels for the histogram
    Tumor_col_color_labels <- rescaled_top_gene_tumor_expression_jet_color[
        colnames(probes_methylation),rev(colnames(rescaled_top_gene_tumor_expression_jet_color))
    ]
    colnames(Tumor_col_color_labels) <- rev(top_gene_names)

    ## Define clustering functions:
    distf=function(d){dist(d,method="euclidean")}
    clustf=function(d){hclust(d,method="ward.D2")}

    ## Create a natural clustering of the tumor sample columns and rows:
    Tumor_col_dist <- distf(
        t(probes_methylation)
    )
    Tumor_col_clust <- clustf(Tumor_col_dist)
    Tumor_col_dend <- as.dendrogram(Tumor_col_clust)

    ## Create a row clustering for the LUAD tumor samples
    Tumor_row_dist <- distf(
        probes_methylation
    )
    Tumor_row_clust <- clustf(Tumor_row_dist)
    Tumor_row_dend <- as.dendrogram(Tumor_row_clust)

    ## Create name for pdf file to put heatmap in:
    path_to_folder <- getwd()

    ## Create a title for the survival plot pdf:
    ## This is a comination of the probe name with the linked gene:
    heatmap_plot_pdf_title <- paste(
        path_to_folder,
        '/hyper.G+.output.probe.heatmap/hyper.G+_top_',
        probe_heatmap_top_n_genes,
        '_linked_genes_heatmap.pdf',
        sep=''
    )

    ## Open a pdf for saving the plot:
    pdf(
        heatmap_plot_pdf_title,
        height= 7,
        width= 10
    )

    ## Create the plot for the tumor samples of interest:
    heatmap.3(
        x= probes_methylation,
        Rowv= Tumor_row_dend,
        Colv= Tumor_col_dend,
        RowSideColors= Tumor_row_color_labels,
        ColSideColors= Tumor_col_color_labels,
        dendrogram= "col",
        labCol= NA,
        labRow= NA,
        lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
        lwid= c(0.1,0.02,2),
        lhei= c(0.4,0.6,2,0.001),
        margins= c(0,0),
        col= matlab::jet.colors(200),
        trace= "none",
        key= FALSE,
        main= NULL,
        ylab= "Enhancer probes linked to key transcriptional regulators",
        xlab= "Samples"
    )

    ## Close the plot:
    dev.off()
}

make.met.heatmap.hypermeth.G_neg <- function(prefix, probe_heatmap_top_n_genes, metDataT, matlab, jet.colors, expDataT) {
    ## Load heatmap.3 function
    source("../scripts/data/heatmap.3.R")

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G-.output.histogram/",
            prefix,
            ".hyper.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hyper.G-.output/hyper.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_IDs <- ordered_TFs_by_link_count[
        1:probe_heatmap_top_n_genes,
        'geneID'
    ]

    top_gene_names <- gencode_v22_genes_df[
        top_gene_IDs,
        'gene_name'
    ]

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes selected:
    probes_linked_to_significant_genes <- unique(
        CpG_linkage_dataset[
            CpG_linkage_dataset$geneID %in% top_gene_IDs,
            'probe'
        ]
    )

    ## Get the rda file with expression and methylation data:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## Load the RDA files
    load(combined_rda_file)

    ## Create a methylation dataset from the tumor samples of interest:
    probes_methylation <- as.matrix(
        metDataT[
            probes_linked_to_significant_genes,
        ]
    )

    ## Create function to convert numbers to jet color vectors:
    jet_color_function <- function(numeric_values){

        jet_color_numeric_color_grad <- matlab::jet.colors(200)

        color_values <- jet_color_numeric_color_grad[numeric_values]

        return(color_values)
    }

    ## Write function to get expression values for the genes of interest:
    tumor_expression_grabber <- function(gene_id){

        unlist(
            c(
                expDataT[
                    c(gene_id),
                ]
            )
        )

    }

    ## Get expression from tumor samples for all the genes of interest:
    top_gene_tumor_expression <- sapply(
        top_gene_IDs,
        tumor_expression_grabber
    )

    ## Write new Rescale function for use later esp. with SCML2
    rescale_func_zero_ignored <- function(x){

        ## get minimum non-zero, non-NA value:
        non_zero_vec <- x[x!=0]

        non_zero_vec <- non_zero_vec[is.na(non_zero_vec)!=TRUE]

        ## get minimum of non-zero values:
        used_min <- min(non_zero_vec)
        used_max <- max(non_zero_vec)

        place <- ((x-used_min)/(used_max-used_min))*200

        return_value <- ifelse(
            place<=0,
            0.001,
            place
        )

        return(
            ceiling(return_value)
        )
    }

    ## Rescale the expression values of the genes of interest:
    rescaled_top_gene_tumor_expression <- apply(
        top_gene_tumor_expression,
        2,
        rescale_func_zero_ignored
    )

    ## Convert the rescaled values to jetcolor:
    rescaled_top_gene_tumor_expression_jet_color <- apply(
        rescaled_top_gene_tumor_expression,
        2,
        jet_color_function
    )

    rownames(rescaled_top_gene_tumor_expression_jet_color) <- rownames(rescaled_top_gene_tumor_expression)

    ## Create row labels for the histogram
    ## (blank for now):
    Tumor_row_color_labels <- rbind(
        rep('white', nrow(probes_methylation))
    )
    rownames(Tumor_row_color_labels) <- c(
        ""
    )

    ## Create column labels for the histogram
    Tumor_col_color_labels <- rescaled_top_gene_tumor_expression_jet_color[
        colnames(probes_methylation),rev(colnames(rescaled_top_gene_tumor_expression_jet_color))
    ]
    colnames(Tumor_col_color_labels) <- rev(top_gene_names)

    ## Define clustering functions:
    distf=function(d){dist(d,method="euclidean")}
    clustf=function(d){hclust(d,method="ward.D2")}

    ## Create a natural clustering of the tumor sample columns and rows:
    Tumor_col_dist <- distf(
        t(probes_methylation)
    )
    Tumor_col_clust <- clustf(Tumor_col_dist)
    Tumor_col_dend <- as.dendrogram(Tumor_col_clust)

    ## Create a row clustering for the LUAD tumor samples
    Tumor_row_dist <- distf(
        probes_methylation
    )
    Tumor_row_clust <- clustf(Tumor_row_dist)
    Tumor_row_dend <- as.dendrogram(Tumor_row_clust)

    ## Create name for pdf file to put heatmap in:
    path_to_folder <- getwd()

    ## Create a title for the survival plot pdf:
    ## This is a comination of the probe name with the linked gene:
    heatmap_plot_pdf_title <- paste(
        path_to_folder,
        '/hyper.G-.output.probe.heatmap/hyper.G-_top_',
        probe_heatmap_top_n_genes,
        '_linked_genes_heatmap.pdf',
        sep=''
    )

    ## Open a pdf for saving the plot:
    pdf(
        heatmap_plot_pdf_title,
        height= 7,
        width= 10
    )

    ## Create the plot for the tumor samples of interest:
    heatmap.3(
        x= probes_methylation,
        Rowv= Tumor_row_dend,
        Colv= Tumor_col_dend,
        RowSideColors= Tumor_row_color_labels,
        ColSideColors= Tumor_col_color_labels,
        dendrogram= "col",
        labCol= NA,
        labRow= NA,
        lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
        lwid= c(0.1,0.02,2),
        lhei= c(0.4,0.6,2,0.001),
        margins= c(0,0),
        col= matlab::jet.colors(200),
        trace= "none",
        key= FALSE,
        main= NULL,
        ylab= "Enhancer probes linked to key transcriptional regulators",
        xlab= "Samples"
    )

    ## Close the plot:
    dev.off()
}

make.met.heatmap.hypometh.G_pos <- function(prefix, probe_heatmap_top_n_genes, metDataT, matlab, jet.colors, expDataT) {
    ## Load heatmap.3 function
    source("../scripts/data/heatmap.3.R")

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hypo.G+.output.histogram/",
            prefix,
            ".hypo.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hypo.G+.output/hypo.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_IDs <- ordered_TFs_by_link_count[
        1:probe_heatmap_top_n_genes,
        'geneID'
    ]

    top_gene_names <- gencode_v22_genes_df[
        top_gene_IDs,
        'gene_name'
    ]

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes selected:
    probes_linked_to_significant_genes <- unique(
        CpG_linkage_dataset[
            CpG_linkage_dataset$geneID %in% top_gene_IDs,
            'probe'
        ]
    )

    ## Get the rda file with expression and methylation data:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## Load the RDA files
    load(combined_rda_file)

    ## Create a methylation dataset from the tumor samples of interest:
    probes_methylation <- as.matrix(
        metDataT[
            probes_linked_to_significant_genes,
        ]
    )

    ## Create function to convert numbers to jet color vectors:
    jet_color_function <- function(numeric_values){

        jet_color_numeric_color_grad <- matlab::jet.colors(200)

        color_values <- jet_color_numeric_color_grad[numeric_values]

        return(color_values)
    }

    ## Write function to get expression values for the genes of interest:
    tumor_expression_grabber <- function(gene_id){

        unlist(
            c(
                expDataT[
                    c(gene_id),
                ]
            )
        )

    }

    ## Get expression from tumor samples for all the genes of interest:
    top_gene_tumor_expression <- sapply(
        top_gene_IDs,
        tumor_expression_grabber
    )

    ## Write new Rescale function for use later esp. with SCML2
    rescale_func_zero_ignored <- function(x){

        ## get minimum non-zero, non-NA value:
        non_zero_vec <- x[x!=0]

        non_zero_vec <- non_zero_vec[is.na(non_zero_vec)!=TRUE]

        ## get minimum of non-zero values:
        used_min <- min(non_zero_vec)
        used_max <- max(non_zero_vec)

        place <- ((x-used_min)/(used_max-used_min))*200

        return_value <- ifelse(
            place<=0,
            0.001,
            place
        )

        return(
            ceiling(return_value)
        )
    }

    ## Rescale the expression values of the genes of interest:
    rescaled_top_gene_tumor_expression <- apply(
        top_gene_tumor_expression,
        2,
        rescale_func_zero_ignored
    )

    ## Convert the rescaled values to jetcolor:
    rescaled_top_gene_tumor_expression_jet_color <- apply(
        rescaled_top_gene_tumor_expression,
        2,
        jet_color_function
    )

    rownames(rescaled_top_gene_tumor_expression_jet_color) <- rownames(rescaled_top_gene_tumor_expression)

    ## Create row labels for the histogram
    ## (blank for now):
    Tumor_row_color_labels <- rbind(
        rep('white', nrow(probes_methylation))
    )
    rownames(Tumor_row_color_labels) <- c(
        ""
    )

    ## Create column labels for the histogram
    Tumor_col_color_labels <- rescaled_top_gene_tumor_expression_jet_color[
        colnames(probes_methylation),rev(colnames(rescaled_top_gene_tumor_expression_jet_color))
    ]
    colnames(Tumor_col_color_labels) <- rev(top_gene_names)

    ## Define clustering functions:
    distf=function(d){dist(d,method="euclidean")}
    clustf=function(d){hclust(d,method="ward.D2")}

    ## Create a natural clustering of the tumor sample columns and rows:
    Tumor_col_dist <- distf(
        t(probes_methylation)
    )
    Tumor_col_clust <- clustf(Tumor_col_dist)
    Tumor_col_dend <- as.dendrogram(Tumor_col_clust)

    ## Create a row clustering for the LUAD tumor samples
    Tumor_row_dist <- distf(
        probes_methylation
    )
    Tumor_row_clust <- clustf(Tumor_row_dist)
    Tumor_row_dend <- as.dendrogram(Tumor_row_clust)

    ## Create name for pdf file to put heatmap in:
    path_to_folder <- getwd()

    ## Create a title for the survival plot pdf:
    ## This is a comination of the probe name with the linked gene:
    heatmap_plot_pdf_title <- paste(
        path_to_folder,
        '/hypo.G+.output.probe.heatmap/hypo.G+_top_',
        probe_heatmap_top_n_genes,
        '_linked_genes_heatmap.pdf',
        sep=''
    )

    ## Open a pdf for saving the plot:
    pdf(
        heatmap_plot_pdf_title,
        height= 7,
        width= 10
    )

    ## Create the plot for the tumor samples of interest:
    heatmap.3(
        x= probes_methylation,
        Rowv= Tumor_row_dend,
        Colv= Tumor_col_dend,
        RowSideColors= Tumor_row_color_labels,
        ColSideColors= Tumor_col_color_labels,
        dendrogram= "col",
        labCol= NA,
        labRow= NA,
        lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
        lwid= c(0.1,0.02,2),
        lhei= c(0.4,0.6,2,0.001),
        margins= c(2,2),
        col= matlab::jet.colors(200),
        trace= "none",
        key= FALSE,
        main= NULL,
        ylab= "Enhancer probes linked to key transcriptional regulators",
        xlab= "Samples"
    )

    ## Close the plot:
    dev.off()
}

make.met.heatmap.hypometh.G_neg <- function(prefix, probe_heatmap_top_n_genes, metDataT, matlab, jet.colors, expDataT) {
    ## Load heatmap.3 function
    source("../scripts/data/heatmap.3.R")

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hypo.G-.output.histogram/",
            prefix,
            ".hypo.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hypo.G-.output/hypo.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_IDs <- ordered_TFs_by_link_count[
        1:probe_heatmap_top_n_genes,
        'geneID'
    ]

    top_gene_names <- gencode_v22_genes_df[
        top_gene_IDs,
        'gene_name'
    ]

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes selected:
    probes_linked_to_significant_genes <- unique(
        CpG_linkage_dataset[
            CpG_linkage_dataset$geneID %in% top_gene_IDs,
            'probe'
        ]
    )

    ## Get the rda file with expression and methylation data:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## Load the RDA files
    load(combined_rda_file)

    ## Create a methylation dataset from the tumor samples of interest:
    probes_methylation <- as.matrix(
        metDataT[
            probes_linked_to_significant_genes,
        ]
    )

    ## Create function to convert numbers to jet color vectors:
    jet_color_function <- function(numeric_values){

        jet_color_numeric_color_grad <- matlab::jet.colors(200)

        color_values <- jet_color_numeric_color_grad[numeric_values]

        return(color_values)
    }

    ## Write function to get expression values for the genes of interest:
    tumor_expression_grabber <- function(gene_id){

        unlist(
            c(
                expDataT[
                    c(gene_id),
                ]
            )
        )

    }

    ## Get expression from tumor samples for all the genes of interest:
    top_gene_tumor_expression <- sapply(
        top_gene_IDs,
        tumor_expression_grabber
    )

    ## Write new Rescale function for use later esp. with SCML2
    rescale_func_zero_ignored <- function(x){

        ## get minimum non-zero, non-NA value:
        non_zero_vec <- x[x!=0]

        non_zero_vec <- non_zero_vec[is.na(non_zero_vec)!=TRUE]

        ## get minimum of non-zero values:
        used_min <- min(non_zero_vec)
        used_max <- max(non_zero_vec)

        place <- ((x-used_min)/(used_max-used_min))*200

        return_value <- ifelse(
            place<=0,
            0.001,
            place
        )

        return(
            ceiling(return_value)
        )
    }

    ## Rescale the expression values of the genes of interest:
    rescaled_top_gene_tumor_expression <- apply(
        top_gene_tumor_expression,
        2,
        rescale_func_zero_ignored
    )

    ## Convert the rescaled values to jetcolor:
    rescaled_top_gene_tumor_expression_jet_color <- apply(
        rescaled_top_gene_tumor_expression,
        2,
        jet_color_function
    )

    rownames(rescaled_top_gene_tumor_expression_jet_color) <- rownames(rescaled_top_gene_tumor_expression)

    ## Create row labels for the histogram
    ## (blank for now):
    Tumor_row_color_labels <- rbind(
        rep('white', nrow(probes_methylation))
    )
    rownames(Tumor_row_color_labels) <- c(
        ""
    )

    ## Create column labels for the histogram
    Tumor_col_color_labels <- rescaled_top_gene_tumor_expression_jet_color[
        colnames(probes_methylation),rev(colnames(rescaled_top_gene_tumor_expression_jet_color))
    ]
    colnames(Tumor_col_color_labels) <- rev(top_gene_names)

    ## Define clustering functions:
    distf=function(d){dist(d,method="euclidean")}
    clustf=function(d){hclust(d,method="ward.D2")}

    ## Create a natural clustering of the tumor sample columns and rows:
    Tumor_col_dist <- distf(
        t(probes_methylation)
    )
    Tumor_col_clust <- clustf(Tumor_col_dist)
    Tumor_col_dend <- as.dendrogram(Tumor_col_clust)

    ## Create a row clustering for the LUAD tumor samples
    Tumor_row_dist <- distf(
        probes_methylation
    )
    Tumor_row_clust <- clustf(Tumor_row_dist)
    Tumor_row_dend <- as.dendrogram(Tumor_row_clust)

    ## Create name for pdf file to put heatmap in:
    path_to_folder <- getwd()

    ## Create a title for the survival plot pdf:
    ## This is a comination of the probe name with the linked gene:
    heatmap_plot_pdf_title <- paste(
        path_to_folder,
        '/hypo.G-.output.probe.heatmap/hypo.G-_top_',
        probe_heatmap_top_n_genes,
        '_linked_genes_heatmap.pdf',
        sep=''
    )

    ## Open a pdf for saving the plot:
    pdf(
        heatmap_plot_pdf_title,
        height= 7,
        width= 10
    )

    ## Create the plot for the tumor samples of interest:
    heatmap.3(
        x= probes_methylation,
        Rowv= Tumor_row_dend,
        Colv= Tumor_col_dend,
        RowSideColors= Tumor_row_color_labels,
        ColSideColors= Tumor_col_color_labels,
        dendrogram= "col",
        labCol= NA,
        labRow= NA,
        lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
        lwid= c(0.1,0.02,2),
        lhei= c(0.4,0.6,2,0.001),
        margins= c(0,0),
        col= matlab::jet.colors(200),
        trace= "none",
        key= FALSE,
        main= NULL,
        ylab= "Enhancer probes linked to key transcriptional regulators",
        xlab= "Samples"
    )

    ## Close the plot:
    dev.off()
}

make.purity.scatterplots.hypermeth.G_pos.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hyper.G+.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hyper.G+.output/", LS[1], sep=""), header=T)
    LP=list.files("../external.data/otherinfo", "Purity")
    Pur=read.delim(paste("../external.data/otherinfo/", LP[1], sep=""), header=T, sep="\t")
    dir.create("hyper.G+.output.complex.scatterplot")
    setwd("./hyper.G+.output.complex.scatterplot")
    M5=match(colnames(metDataT), Pur$ID)
    Purity.c=Pur[M5,]
    Purity.c$score=(1+2*(Purity.c$purity))
    Purity.c$score[which(is.na(Purity.c$score)==TRUE)]=paste(1)
    Purity=c(Purity.c$score, rep(3, dim(metDataN)[2]))
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=as.numeric(Purity), shape=16)+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.purity.scatterplots.hypermeth.G_neg.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hyper.G-.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hyper.G-.output/", LS[1], sep=""), header=T)
    LP=list.files("../external.data/otherinfo", "Purity")
    Pur=read.delim(paste("../external.data/otherinfo/", LP[1], sep=""), header=T, sep="\t")
    dir.create("hyper.G-.output.complex.scatterplot")
    setwd("./hyper.G-.output.complex.scatterplot")
    M5=match(colnames(metDataT), Pur$ID)
    Purity.c=Pur[M5,]
    Purity.c$score=(1+2*(Purity.c$purity))
    Purity.c$score[which(is.na(Purity.c$score)==TRUE)]=paste(1)
    Purity=c(Purity.c$score, rep(3, dim(metDataN)[2]))
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=as.numeric(Purity), shape=16)+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.purity.scatterplots.hypometh.G_pos.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hypo.G+.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hypo.G+.output/", LS[1], sep=""), header=T)
    LP=list.files("../external.data/otherinfo", "Purity")
    Pur=read.delim(paste("../external.data/otherinfo/", LP[1], sep=""), header=T, sep="\t")
    dir.create("hypo.G+.output.complex.scatterplot")
    setwd("./hypo.G+.output.complex.scatterplot")
    M5=match(colnames(metDataT), Pur$ID)
    Purity.c=Pur[M5,]
    Purity.c$score=(1+2*(Purity.c$purity))
    Purity.c$score[which(is.na(Purity.c$score)==TRUE)]=paste(1)
    Purity=c(Purity.c$score, rep(3, dim(metDataN)[2]))
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=as.numeric(Purity), shape=16)+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.purity.scatterplots.hypometh.G_neg.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hypo.G-.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hypo.G-.output/", LS[1], sep=""), header=T)
    LP=list.files("../external.data/otherinfo", "Purity")
    Pur=read.delim(paste("../external.data/otherinfo/", LP[1], sep=""), header=T, sep="\t")
    dir.create("hypo.G-.output.complex.scatterplot")
    setwd("./hypo.G-.output.complex.scatterplot")
    M5=match(colnames(metDataT), Pur$ID)
    Purity.c=Pur[M5,]
    Purity.c$score=(1+2*(Purity.c$purity))
    Purity.c$score[which(is.na(Purity.c$score)==TRUE)]=paste(1)
    Purity=c(Purity.c$score, rep(3, dim(metDataN)[2]))
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=as.numeric(Purity), shape=16)+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.simple.scatterplots.hypermeth.G_pos.source <- function(ggplot2, prefix, i, scatterplot_top_n_genes, metDataT, metDataN, expDataT, expDataN) {
    ## Load ggplot2:
    # library(ggplot2)

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts from step5 histogram function:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G+.output.histogram/",
            prefix,
            ".hyper.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Check to see if a genes_of_interest file has been
    ## provided by the user:
    genes_of_user_interest_file <- list.files("../external.data/otherinfo", "genes_of_interest")

    ## If the genes of interest file has been provided, load it
    ## and add those genes to the list of TFs for the scatterplots:
    if(length(genes_of_user_interest_file)==1){

        ## Load the file:
        genes_of_user_interest <- read.delim(
            paste(
                "../external.data/otherinfo/",
                genes_of_user_interest_file,
                sep=""
            ),
            header= FALSE,
            stringsAsFactors = FALSE
        )

        ## Get the IDs of the genes:
        user_gene_IDs <- genes_of_user_interest$V1

        ## Index an empty vector for the gene ENSG IDs:
        user_gene_ENSG <- character()

        ## Get actual gene name depending on input:
        for(i in c(1:length(user_gene_IDs))){

            if(substring(user_gene_IDs[i],1,4)=='ENSG' & nchar(user_gene_IDs[i])==15){

                ## Input is in ENSG, assume user has provided correctly
                ## (will check later if it is a TR with links)
                user_gene_ENSG[i] <- user_gene_IDs[i]

                ## Get gene ENSG assuming a name is plugged in:
            } else{

                ## Assume what was given was the gene name, get the ENSG:
                gene_ENSG_placeholder <- rownames(
                    gencode_v22_genes_df[gencode_v22_genes_df$gene_name==user_gene_IDs[i],]
                )

                ## If the gene name has an associated ENSG ID, add it to the list:
                if(length(gene_ENSG_placeholder)==1){

                    user_gene_ENSG[i] <- gene_ENSG_placeholder
                }
            }
        }

        ## Remove NA values from user input list:
        user_gene_ENSG <- user_gene_ENSG[
            !is.na(user_gene_ENSG)
        ]

        ## Check that the ENSGs from the user
        ## have links annotated to them:
        user_gene_ENSG_w_links <- intersect(
            unique(user_gene_ENSG),
            ordered_TFs_by_link_count$geneID
        )
    }

    ## Get the name of the top n TRs
    ## based on survival_top_n_genes parameter:
    if(scatterplot_top_n_genes>0){

        top_gene_IDs <- ordered_TFs_by_link_count[
            1:scatterplot_top_n_genes,
            'geneID'
        ]

    } else{

        top_gene_IDs <- character()
    }

    ## Add the user input TFs to the top TRs specified:
    top_gene_IDs <- c(
        top_gene_IDs,
        user_gene_ENSG_w_links
    )

    ## Get the names for the top genes of interest:
    top_gene_names <- gencode_v22_genes_df[
        top_gene_IDs,
        'gene_name'
    ]

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hyper.G+.output/hyper.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the rda file with expression and methylation data:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## Load the RDA files
    load(combined_rda_file)

    ## Combine the methylation and expression data from both tumor and normal samples:
    metDataF_subC <- cbind(metDataT, metDataN)
    expDataF_subC <- cbind(expDataT, expDataN)

    ## Create dataframe with color info for the tumor and normal samples
    ## blue for normal, red for tumor data points:
    DichF <- data.frame(
        group=c(
            colnames(metDataT), colnames(metDataN)
        ),
        cluster=c(
            rep(
                "Tumor",
                dim(metDataT)[2]
            ),
            rep(
                "Normal",
                dim(metDataN)[2]
            )
        ),
        stringsAsFactors = FALSE
    )

    ## Create a new directory for hyper.G+ scatterplot output and
    ## set the working directory ot that new folder:
    dir.create("hyper.G+.output.scatterplot")
    setwd("./hyper.G+.output.scatterplot")

    ## For each TR, create a vector of CpGs linked to that TF:

    ## Index an empty list:
    linked_cpgs_list <- list()

    ## For each TR of interest, get the list of CpGs associated with it:
    for(i in c(1:length(top_gene_IDs))){

        ## Get the TRs ENSG:
        TR_ENSG_placeholder <- top_gene_IDs[i]

        ## Get the probes linked to each TR
        probes_linked_to_significant_gene <- unique(
            CpG_linkage_dataset[
                CpG_linkage_dataset$geneID %in% top_gene_IDs[i],
                'probe'
            ]
        )

        ## Add the probes to the list:
        linked_cpgs_list[[i]] <- probes_linked_to_significant_gene
    }

    ## Add the names of the TRs to the list:
    names(linked_cpgs_list) <- top_gene_IDs

    ## Write a function that, when given a list of CpGs linked to a TF,
    ## will plot scatterplots showing the methylation of the probe on the X
    ## and expression of that TF on the Y across all the tumor and normal smaples:
    scatterplot_function <- function(CpGs_list, gene_ENSG){

        ## Unlist the CpGs linked to each probe:
        unlisted_CpGs <- c(
            unlist(
                CpGs_list
            )
        )

        ## Convert the gene ENSG into the gene name:
        gene_name <- gencode_v22_genes_df[
            gene_ENSG,
            'gene_name'
        ]

        ## Add the expression of the gene of interest to DichF:
        TR_expression <- c(
            unlist(
                expDataF_subC[
                    gene_ENSG,
                    DichF$group
                ]
            )
        )

        # ## REMOVE LATER TESTING ONLY:
        # CpGs_linked_to_TR <- unlisted_CpGs[1]

        ## Now write an internal function that will get each linked CpGs methylation,
        ## and use the methylation and TR expression to create a ggplot2 scatterplot
        ## and save it:
        internal_scatterplot_function <- function(CpGs_linked_to_TR){

            ## Save the CpG name:
            CpG_name_placeholder <- CpGs_linked_to_TR

            ## Get DNA methylation values:
            unlisted_CpG_methylation <- c(
                unlist(
                    metDataF_subC[
                        CpGs_linked_to_TR,
                        DichF$group
                    ]
                )
            )

            ## Manually coloring samples:
            t_v_n_group_colors <- c('Normal'='dodgerblue3', 'Tumor'='red3')

            ## Creating scatter with ggplot2:
            scatter_plot <- qplot(
                x=TR_expression,
                y=unlisted_CpG_methylation,
                geom=c("point"),
                colour=DichF$cluster
            )

            ## Create the plot:
            scatter_plot_updated <- scatter_plot +
                ggtitle(
                    paste(
                        gene_name,
                        ' gene expression vs.\n',
                        CpG_name_placeholder,
                        ' DNA methylation',
                        sep=''
                    )
                ) +
                ylab(
                    paste(
                        CpG_name_placeholder,
                        'DNA methylation',
                        sep=' '
                    )
                ) +
                xlab(
                    paste(
                        gene_name,
                        'gene expression',
                        sep=' '
                    )
                ) +
                theme_bw() +
                scale_color_manual(
                    values=t_v_n_group_colors,
                    name="Sample type"
                ) +
                theme(
                    plot.title = element_text(hjust=0.5, size=18),
                    legend.title = element_text(hjust=0.5, size=12),
                    legend.text = element_text(size=10, colour='black'),
                    panel.border = element_rect(colour = 'black', fill=NA, size=1),
                    axis.title.x = element_text(size=16),
                    axis.title.y = element_text(size=16),
                    axis.text.x = element_text(size=14, colour='black'),
                    axis.text.y = element_text(size=14, colour='black'),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                )

            ## Create a title for the scatter plot pdf:
            ## This is a comination of the probe name with the linked gene:
            scatterplot_pdf_title <- paste(
                './',
                gene_name,
                '_',
                CpG_name_placeholder,
                '_scatterplot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(
                scatterplot_pdf_title,
                height= 7,
                width= 10
            )

            plot(scatter_plot_updated)

            ## Close the plot:
            dev.off()
        }

        lapply(
            X= unlisted_CpGs,
            FUN= internal_scatterplot_function
        )

    }

    ## Generate the plots for all TRs of interest:
    suppressWarnings(
        mapply(
            FUN= scatterplot_function,
            CpGs_list= linked_cpgs_list,
            gene_ENSG= names(linked_cpgs_list)
        )
    )
}

make.simple.scatterplots.hypermeth.G_neg.source <- function(ggplot2, prefix, i, scatterplot_top_n_genes, metDataT, metDataN, expDataT, expDataN) {
    ## Load ggplot2:
    # library(ggplot2)

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts from step5 histogram function:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G-.output.histogram/",
            prefix,
            ".hyper.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Check to see if a genes_of_interest file has been
    ## provided by the user:
    genes_of_user_interest_file <- list.files("../external.data/otherinfo", "genes_of_interest")

    ## If the genes of interest file has been provided, load it
    ## and add those genes to the list of TFs for the scatterplots:
    if(length(genes_of_user_interest_file)==1){

        ## Load the file:
        genes_of_user_interest <- read.delim(
            paste(
                "../external.data/otherinfo/",
                genes_of_user_interest_file,
                sep=""
            ),
            header= FALSE,
            stringsAsFactors = FALSE
        )

        ## Get the IDs of the genes:
        user_gene_IDs <- genes_of_user_interest$V1

        ## Index an empty vector for the gene ENSG IDs:
        user_gene_ENSG <- character()

        ## Get actual gene name depending on input:
        for(i in c(1:length(user_gene_IDs))){

            if(substring(user_gene_IDs[i],1,4)=='ENSG' & nchar(user_gene_IDs[i])==15){

                ## Input is in ENSG, assume user has provided correctly
                ## (will check later if it is a TR with links)
                user_gene_ENSG[i] <- user_gene_IDs[i]

                ## Get gene ENSG assuming a name is plugged in:
            } else{

                ## Assume what was given was the gene name, get the ENSG:
                gene_ENSG_placeholder <- rownames(
                    gencode_v22_genes_df[gencode_v22_genes_df$gene_name==user_gene_IDs[i],]
                )

                ## If the gene name has an associated ENSG ID, add it to the list:
                if(length(gene_ENSG_placeholder)==1){

                    user_gene_ENSG[i] <- gene_ENSG_placeholder
                }
            }
        }

        ## Remove NA values from user input list:
        user_gene_ENSG <- user_gene_ENSG[
            !is.na(user_gene_ENSG)
        ]

        ## Check that the ENSGs from the user
        ## have links annotated to them:
        user_gene_ENSG_w_links <- intersect(
            unique(user_gene_ENSG),
            ordered_TFs_by_link_count$geneID
        )
    }

    ## Get the name of the top n TRs
    ## based on survival_top_n_genes parameter:
    if(scatterplot_top_n_genes>0){

        top_gene_IDs <- ordered_TFs_by_link_count[
            1:scatterplot_top_n_genes,
            'geneID'
        ]

    } else{

        top_gene_IDs <- character()
    }

    ## Add the user input TFs to the top TRs specified:
    top_gene_IDs <- c(
        top_gene_IDs,
        user_gene_ENSG_w_links
    )

    ## Get the names for the top genes of interest:
    top_gene_names <- gencode_v22_genes_df[
        top_gene_IDs,
        'gene_name'
    ]

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hyper.G-.output/hyper.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the rda file with expression and methylation data:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## Load the RDA files
    load(combined_rda_file)

    ## Combine the methylation and expression data from both tumor and normal samples:
    metDataF_subC <- cbind(metDataT, metDataN)
    expDataF_subC <- cbind(expDataT, expDataN)

    ## Create dataframe with color info for the tumor and normal samples
    ## blue for normal, red for tumor data points:
    DichF <- data.frame(
        group=c(
            colnames(metDataT), colnames(metDataN)
        ),
        cluster=c(
            rep(
                "Tumor",
                dim(metDataT)[2]
            ),
            rep(
                "Normal",
                dim(metDataN)[2]
            )
        ),
        stringsAsFactors = FALSE
    )

    ## Create a new directory for hyper.G- scatterplot output and
    ## set the working directory ot that new folder:
    dir.create("hyper.G-.output.scatterplot")
    setwd("./hyper.G-.output.scatterplot")

    ## For each TR, create a vector of CpGs linked to that TF:

    ## Index an empty list:
    linked_cpgs_list <- list()

    ## For each TR of interest, get the list of CpGs associated with it:
    for(i in c(1:length(top_gene_IDs))){

        ## Get the TRs ENSG:
        TR_ENSG_placeholder <- top_gene_IDs[i]

        ## Get the probes linked to each TR
        probes_linked_to_significant_gene <- unique(
            CpG_linkage_dataset[
                CpG_linkage_dataset$geneID %in% top_gene_IDs[i],
                'probe'
            ]
        )

        ## Add the probes to the list:
        linked_cpgs_list[[i]] <- probes_linked_to_significant_gene
    }

    ## Add the names of the TRs to the list:
    names(linked_cpgs_list) <- top_gene_IDs

    ## Write a function that, when given a list of CpGs linked to a TF,
    ## will plot scatterplots showing the methylation of the probe on the X
    ## and expression of that TF on the Y across all the tumor and normal smaples:
    scatterplot_function <- function(CpGs_list, gene_ENSG){

        ## Unlist the CpGs linked to each probe:
        unlisted_CpGs <- c(
            unlist(
                CpGs_list
            )
        )

        ## Convert the gene ENSG into the gene name:
        gene_name <- gencode_v22_genes_df[
            gene_ENSG,
            'gene_name'
        ]

        ## Add the expression of the gene of interest to DichF:
        TR_expression <- c(
            unlist(
                expDataF_subC[
                    gene_ENSG,
                    DichF$group
                ]
            )
        )

        # ## REMOVE LATER TESTING ONLY:
        # CpGs_linked_to_TR <- unlisted_CpGs[1]

        ## Now write an internal function that will get each linked CpGs methylation,
        ## and use the methylation and TR expression to create a ggplot2 scatterplot
        ## and save it:
        internal_scatterplot_function <- function(CpGs_linked_to_TR){

            ## Save the CpG name:
            CpG_name_placeholder <- CpGs_linked_to_TR

            ## Get DNA methylation values:
            unlisted_CpG_methylation <- c(
                unlist(
                    metDataF_subC[
                        CpGs_linked_to_TR,
                        DichF$group
                    ]
                )
            )

            ## Manually coloring samples:
            t_v_n_group_colors <- c('Normal'='dodgerblue3', 'Tumor'='red3')

            ## Creating scatter with ggplot2:
            scatter_plot <- qplot(
                x=TR_expression,
                y=unlisted_CpG_methylation,
                geom=c("point"),
                colour=DichF$cluster
            )

            ## Create the plot:
            scatter_plot_updated <- scatter_plot +
                ggtitle(
                    paste(
                        gene_name,
                        ' gene expression vs.\n',
                        CpG_name_placeholder,
                        ' DNA methylation',
                        sep=''
                    )
                ) +
                ylab(
                    paste(
                        CpG_name_placeholder,
                        'DNA methylation',
                        sep=' '
                    )
                ) +
                xlab(
                    paste(
                        gene_name,
                        'gene expression',
                        sep=' '
                    )
                ) +
                theme_bw() +
                scale_color_manual(
                    values=t_v_n_group_colors,
                    name="Sample type"
                ) +
                theme(
                    plot.title = element_text(hjust=0.5, size=18),
                    legend.title = element_text(hjust=0.5, size=12),
                    legend.text = element_text(size=10, colour='black'),
                    panel.border = element_rect(colour = 'black', fill=NA, size=1),
                    axis.title.x = element_text(size=16),
                    axis.title.y = element_text(size=16),
                    axis.text.x = element_text(size=14, colour='black'),
                    axis.text.y = element_text(size=14, colour='black'),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                )

            ## Create a title for the scatter plot pdf:
            ## This is a comination of the probe name with the linked gene:
            scatterplot_pdf_title <- paste(
                './',
                gene_name,
                '_',
                CpG_name_placeholder,
                '_scatterplot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(
                scatterplot_pdf_title,
                height= 7,
                width= 10
            )

            plot(scatter_plot_updated)

            ## Close the plot:
            dev.off()
        }

        lapply(
            X= unlisted_CpGs,
            FUN= internal_scatterplot_function
        )

    }

    ## Generate the plots for all TRs of interest:
    suppressWarnings(
        mapply(
            FUN= scatterplot_function,
            CpGs_list= linked_cpgs_list,
            gene_ENSG= names(linked_cpgs_list)
        )
    )
}

make.simple.scatterplots.hypometh.G_pos.source <- function(ggplot2, prefix, i, scatterplot_top_n_genes, metDataT, metDataN, expDataT, expDataN) {
    ## Load ggplot2:
    # library(ggplot2)

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts from step5 histogram function:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hypo.G+.output.histogram/",
            prefix,
            ".hypo.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Check to see if a genes_of_interest file has been
    ## provided by the user:
    genes_of_user_interest_file <- list.files("../external.data/otherinfo", "genes_of_interest")

    ## If the genes of interest file has been provided, load it
    ## and add those genes to the list of TFs for the scatterplots:
    if(length(genes_of_user_interest_file)==1){

        ## Load the file:
        genes_of_user_interest <- read.delim(
            paste(
                "../external.data/otherinfo/",
                genes_of_user_interest_file,
                sep=""
            ),
            header= FALSE,
            stringsAsFactors = FALSE
        )

        ## Get the IDs of the genes:
        user_gene_IDs <- genes_of_user_interest$V1

        ## Index an empty vector for the gene ENSG IDs:
        user_gene_ENSG <- character()

        ## Get actual gene name depending on input:
        for(i in c(1:length(user_gene_IDs))){

            if(substring(user_gene_IDs[i],1,4)=='ENSG' & nchar(user_gene_IDs[i])==15){

                ## Input is in ENSG, assume user has provided correctly
                ## (will check later if it is a TR with links)
                user_gene_ENSG[i] <- user_gene_IDs[i]

                ## Get gene ENSG assuming a name is plugged in:
            } else{

                ## Assume what was given was the gene name, get the ENSG:
                gene_ENSG_placeholder <- rownames(
                    gencode_v22_genes_df[gencode_v22_genes_df$gene_name==user_gene_IDs[i],]
                )

                ## If the gene name has an associated ENSG ID, add it to the list:
                if(length(gene_ENSG_placeholder)==1){

                    user_gene_ENSG[i] <- gene_ENSG_placeholder
                }
            }
        }

        ## Remove NA values from user input list:
        user_gene_ENSG <- user_gene_ENSG[
            !is.na(user_gene_ENSG)
        ]

        ## Check that the ENSGs from the user
        ## have links annotated to them:
        user_gene_ENSG_w_links <- intersect(
            unique(user_gene_ENSG),
            ordered_TFs_by_link_count$geneID
        )
    }

    ## Get the name of the top n TRs
    ## based on survival_top_n_genes parameter:
    if(scatterplot_top_n_genes>0){

        top_gene_IDs <- ordered_TFs_by_link_count[
            1:scatterplot_top_n_genes,
            'geneID'
        ]

    } else{

        top_gene_IDs <- character()
    }

    ## Add the user input TFs to the top TRs specified:
    top_gene_IDs <- c(
        top_gene_IDs,
        user_gene_ENSG_w_links
    )

    ## Get the names for the top genes of interest:
    top_gene_names <- gencode_v22_genes_df[
        top_gene_IDs,
        'gene_name'
    ]

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hypo.G+.output/hypo.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the rda file with expression and methylation data:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## Load the RDA files
    load(combined_rda_file)

    ## Combine the methylation and expression data from both tumor and normal samples:
    metDataF_subC <- cbind(metDataT, metDataN)
    expDataF_subC <- cbind(expDataT, expDataN)

    ## Create dataframe with color info for the tumor and normal samples
    ## blue for normal, red for tumor data points:
    DichF <- data.frame(
        group=c(
            colnames(metDataT), colnames(metDataN)
        ),
        cluster=c(
            rep(
                "Tumor",
                dim(metDataT)[2]
            ),
            rep(
                "Normal",
                dim(metDataN)[2]
            )
        ),
        stringsAsFactors = FALSE
    )

    ## Create a new directory for hypo.G+ scatterplot output and
    ## set the working directory ot that new folder:
    dir.create("hypo.G+.output.scatterplot")
    setwd("./hypo.G+.output.scatterplot")

    ## For each TR, create a vector of CpGs linked to that TF:

    ## Index an empty list:
    linked_cpgs_list <- list()

    ## For each TR of interest, get the list of CpGs associated with it:
    for(i in c(1:length(top_gene_IDs))){

        ## Get the TRs ENSG:
        TR_ENSG_placeholder <- top_gene_IDs[i]

        ## Get the probes linked to each TR
        probes_linked_to_significant_gene <- unique(
            CpG_linkage_dataset[
                CpG_linkage_dataset$geneID %in% top_gene_IDs[i],
                'probe'
            ]
        )

        ## Add the probes to the list:
        linked_cpgs_list[[i]] <- probes_linked_to_significant_gene
    }

    ## Add the names of the TRs to the list:
    names(linked_cpgs_list) <- top_gene_IDs

    ## Write a function that, when given a list of CpGs linked to a TF,
    ## will plot scatterplots showing the methylation of the probe on the X
    ## and expression of that TF on the Y across all the tumor and normal smaples:
    scatterplot_function <- function(CpGs_list, gene_ENSG){

        ## Unlist the CpGs linked to each probe:
        unlisted_CpGs <- c(
            unlist(
                CpGs_list
            )
        )

        ## Convert the gene ENSG into the gene name:
        gene_name <- gencode_v22_genes_df[
            gene_ENSG,
            'gene_name'
        ]

        ## Add the expression of the gene of interest to DichF:
        TR_expression <- c(
            unlist(
                expDataF_subC[
                    gene_ENSG,
                    DichF$group
                ]
            )
        )

        # ## REMOVE LATER TESTING ONLY:
        # CpGs_linked_to_TR <- unlisted_CpGs[1]

        ## Now write an internal function that will get each linked CpGs methylation,
        ## and use the methylation and TR expression to create a ggplot2 scatterplot
        ## and save it:
        internal_scatterplot_function <- function(CpGs_linked_to_TR){

            ## Save the CpG name:
            CpG_name_placeholder <- CpGs_linked_to_TR

            ## Get DNA methylation values:
            unlisted_CpG_methylation <- c(
                unlist(
                    metDataF_subC[
                        CpGs_linked_to_TR,
                        DichF$group
                    ]
                )
            )

            ## Manually coloring samples:
            t_v_n_group_colors <- c('Normal'='dodgerblue3', 'Tumor'='red3')

            ## Creating scatter with ggplot2:
            scatter_plot <- qplot(
                x=TR_expression,
                y=unlisted_CpG_methylation,
                geom=c("point"),
                colour=DichF$cluster
            )

            ## Create the plot:
            scatter_plot_updated <- scatter_plot +
                ggtitle(
                    paste(
                        gene_name,
                        ' gene expression vs.\n',
                        CpG_name_placeholder,
                        ' DNA methylation',
                        sep=''
                    )
                ) +
                ylab(
                    paste(
                        CpG_name_placeholder,
                        'DNA methylation',
                        sep=' '
                    )
                ) +
                xlab(
                    paste(
                        gene_name,
                        'gene expression',
                        sep=' '
                    )
                ) +
                theme_bw() +
                scale_color_manual(
                    values=t_v_n_group_colors,
                    name="Sample type"
                ) +
                theme(
                    plot.title = element_text(hjust=0.5, size=18),
                    legend.title = element_text(hjust=0.5, size=12),
                    legend.text = element_text(size=10, colour='black'),
                    panel.border = element_rect(colour = 'black', fill=NA, size=1),
                    axis.title.x = element_text(size=16),
                    axis.title.y = element_text(size=16),
                    axis.text.x = element_text(size=14, colour='black'),
                    axis.text.y = element_text(size=14, colour='black'),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                )

            ## Create a title for the scatter plot pdf:
            ## This is a comination of the probe name with the linked gene:
            scatterplot_pdf_title <- paste(
                './',
                gene_name,
                '_',
                CpG_name_placeholder,
                '_scatterplot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(
                scatterplot_pdf_title,
                height= 7,
                width= 10
            )

            plot(scatter_plot_updated)

            ## Close the plot:
            dev.off()
        }

        lapply(
            X= unlisted_CpGs,
            FUN= internal_scatterplot_function
        )

    }

    ## Generate the plots for all TRs of interest:
    suppressWarnings(
        mapply(
            FUN= scatterplot_function,
            CpGs_list= linked_cpgs_list,
            gene_ENSG= names(linked_cpgs_list)
        )
    )
}

make.simple.scatterplots.hypometh.G_neg.source <- function(ggplot2, prefix, i, scatterplot_top_n_genes, metDataT, metDataN, expDataT, expDataN) {
    ## Load ggplot2:
    # library(ggplot2)

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts from step5 histogram function:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hypo.G-.output.histogram/",
            prefix,
            ".hypo.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Check to see if a genes_of_interest file has been
    ## provided by the user:
    genes_of_user_interest_file <- list.files("../external.data/otherinfo", "genes_of_interest")

    ## If the genes of interest file has been provided, load it
    ## and add those genes to the list of TFs for the scatterplots:
    if(length(genes_of_user_interest_file)==1){

        ## Load the file:
        genes_of_user_interest <- read.delim(
            paste(
                "../external.data/otherinfo/",
                genes_of_user_interest_file,
                sep=""
            ),
            header= FALSE,
            stringsAsFactors = FALSE
        )

        ## Get the IDs of the genes:
        user_gene_IDs <- genes_of_user_interest$V1

        ## Index an empty vector for the gene ENSG IDs:
        user_gene_ENSG <- character()

        ## Get actual gene name depending on input:
        for(i in c(1:length(user_gene_IDs))){

            if(substring(user_gene_IDs[i],1,4)=='ENSG' & nchar(user_gene_IDs[i])==15){

                ## Input is in ENSG, assume user has provided correctly
                ## (will check later if it is a TR with links)
                user_gene_ENSG[i] <- user_gene_IDs[i]

                ## Get gene ENSG assuming a name is plugged in:
            } else{

                ## Assume what was given was the gene name, get the ENSG:
                gene_ENSG_placeholder <- rownames(
                    gencode_v22_genes_df[gencode_v22_genes_df$gene_name==user_gene_IDs[i],]
                )

                ## If the gene name has an associated ENSG ID, add it to the list:
                if(length(gene_ENSG_placeholder)==1){

                    user_gene_ENSG[i] <- gene_ENSG_placeholder
                }
            }
        }

        ## Remove NA values from user input list:
        user_gene_ENSG <- user_gene_ENSG[
            !is.na(user_gene_ENSG)
        ]

        ## Check that the ENSGs from the user
        ## have links annotated to them:
        user_gene_ENSG_w_links <- intersect(
            unique(user_gene_ENSG),
            ordered_TFs_by_link_count$geneID
        )
    }

    ## Get the name of the top n TRs
    ## based on survival_top_n_genes parameter:
    if(scatterplot_top_n_genes>0){

        top_gene_IDs <- ordered_TFs_by_link_count[
            1:scatterplot_top_n_genes,
            'geneID'
        ]

    } else{

        top_gene_IDs <- character()
    }

    ## Add the user input TFs to the top TRs specified:
    top_gene_IDs <- c(
        top_gene_IDs,
        user_gene_ENSG_w_links
    )

    ## Get the names for the top genes of interest:
    top_gene_names <- gencode_v22_genes_df[
        top_gene_IDs,
        'gene_name'
    ]

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hypo.G-.output/hypo.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the rda file with expression and methylation data:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## Load the RDA files
    load(combined_rda_file)

    ## Combine the methylation and expression data from both tumor and normal samples:
    metDataF_subC <- cbind(metDataT, metDataN)
    expDataF_subC <- cbind(expDataT, expDataN)

    ## Create dataframe with color info for the tumor and normal samples
    ## blue for normal, red for tumor data points:
    DichF <- data.frame(
        group=c(
            colnames(metDataT), colnames(metDataN)
        ),
        cluster=c(
            rep(
                "Tumor",
                dim(metDataT)[2]
            ),
            rep(
                "Normal",
                dim(metDataN)[2]
            )
        ),
        stringsAsFactors = FALSE
    )

    ## Create a new directory for hypo.G- scatterplot output and
    ## set the working directory ot that new folder:
    dir.create("hypo.G-.output.scatterplot")
    setwd("./hypo.G-.output.scatterplot")

    ## For each TR, create a vector of CpGs linked to that TF:

    ## Index an empty list:
    linked_cpgs_list <- list()

    ## For each TR of interest, get the list of CpGs associated with it:
    for(i in c(1:length(top_gene_IDs))){

        ## Get the TRs ENSG:
        TR_ENSG_placeholder <- top_gene_IDs[i]

        ## Get the probes linked to each TR
        probes_linked_to_significant_gene <- unique(
            CpG_linkage_dataset[
                CpG_linkage_dataset$geneID %in% top_gene_IDs[i],
                'probe'
            ]
        )

        ## Add the probes to the list:
        linked_cpgs_list[[i]] <- probes_linked_to_significant_gene
    }

    ## Add the names of the TRs to the list:
    names(linked_cpgs_list) <- top_gene_IDs

    ## Write a function that, when given a list of CpGs linked to a TF,
    ## will plot scatterplots showing the methylation of the probe on the X
    ## and expression of that TF on the Y across all the tumor and normal smaples:
    scatterplot_function <- function(CpGs_list, gene_ENSG){

        ## Unlist the CpGs linked to each probe:
        unlisted_CpGs <- c(
            unlist(
                CpGs_list
            )
        )

        ## Convert the gene ENSG into the gene name:
        gene_name <- gencode_v22_genes_df[
            gene_ENSG,
            'gene_name'
        ]

        ## Add the expression of the gene of interest to DichF:
        TR_expression <- c(
            unlist(
                expDataF_subC[
                    gene_ENSG,
                    DichF$group
                ]
            )
        )

        # ## REMOVE LATER TESTING ONLY:
        # CpGs_linked_to_TR <- unlisted_CpGs[1]

        ## Now write an internal function that will get each linked CpGs methylation,
        ## and use the methylation and TR expression to create a ggplot2 scatterplot
        ## and save it:
        internal_scatterplot_function <- function(CpGs_linked_to_TR){

            ## Save the CpG name:
            CpG_name_placeholder <- CpGs_linked_to_TR

            ## Get DNA methylation values:
            unlisted_CpG_methylation <- c(
                unlist(
                    metDataF_subC[
                        CpGs_linked_to_TR,
                        DichF$group
                    ]
                )
            )

            ## Manually coloring samples:
            t_v_n_group_colors <- c('Normal'='dodgerblue3', 'Tumor'='red3')

            ## Creating scatter with ggplot2:
            scatter_plot <- qplot(
                x=TR_expression,
                y=unlisted_CpG_methylation,
                geom=c("point"),
                colour=DichF$cluster
            )

            ## Create the plot:
            scatter_plot_updated <- scatter_plot +
                ggtitle(
                    paste(
                        gene_name,
                        ' gene expression vs.\n',
                        CpG_name_placeholder,
                        ' DNA methylation',
                        sep=''
                    )
                ) +
                ylab(
                    paste(
                        CpG_name_placeholder,
                        'DNA methylation',
                        sep=' '
                    )
                ) +
                xlab(
                    paste(
                        gene_name,
                        'gene expression',
                        sep=' '
                    )
                ) +
                theme_bw() +
                scale_color_manual(
                    values=t_v_n_group_colors,
                    name="Sample type"
                ) +
                theme(
                    plot.title = element_text(hjust=0.5, size=18),
                    legend.title = element_text(hjust=0.5, size=12),
                    legend.text = element_text(size=10, colour='black'),
                    panel.border = element_rect(colour = 'black', fill=NA, size=1),
                    axis.title.x = element_text(size=16),
                    axis.title.y = element_text(size=16),
                    axis.text.x = element_text(size=14, colour='black'),
                    axis.text.y = element_text(size=14, colour='black'),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                )

            ## Create a title for the scatter plot pdf:
            ## This is a comination of the probe name with the linked gene:
            scatterplot_pdf_title <- paste(
                './',
                gene_name,
                '_',
                CpG_name_placeholder,
                '_scatterplot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(
                scatterplot_pdf_title,
                height= 7,
                width= 10
            )

            plot(scatter_plot_updated)

            ## Close the plot:
            dev.off()
        }

        lapply(
            X= unlisted_CpGs,
            FUN= internal_scatterplot_function
        )

    }

    ## Generate the plots for all TRs of interest:
    suppressWarnings(
        mapply(
            FUN= scatterplot_function,
            CpGs_list= linked_cpgs_list,
            gene_ENSG= names(linked_cpgs_list)
        )
    )
}

make.simple.scatterplots.select.probe <- function(ggplot2, metDataT, metDataN, expDataT, expDataN, findhypoGpos, findhypoGneg, findhyperGpos, findhyperGneg, i) {
    ## Load ggplot2:
    # library(ggplot2)

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gTR object and associated data frame:
    rm(gencode_v22_genes)

    ## Check to see if a probes_of_interest file has been
    ## provided by the user:
    probes_of_user_interest_file_name <- list.files("../external.data/otherinfo", "probes_of_interest")

    ## If the genes of interest file has been provided, load it
    ## and get a vector of the probes of interest:
    if(length(probes_of_user_interest_file_name)==1){

        ## Load the file:
        probes_of_user_interest_file <- read.delim(
            paste(
                "../external.data/otherinfo/",
                probes_of_user_interest_file_name,
                sep=""
            ),
            header= FALSE,
            stringsAsFactors = FALSE
        )

        ## Get the probes:
        probes_of_user_interest <- probes_of_user_interest_file$V1

    }

    ## Get the rda file with expression and methylation data:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## Load the RDA files
    load(combined_rda_file)

    ## Combine the methylation and expression data from both tumor and normal samples:
    metDataF_subC <- cbind(metDataT, metDataN)
    expDataF_subC <- cbind(expDataT, expDataN)

    ## Create dataframe with color info for the tumor and normal samples
    ## blue for normal, red for tumor data points:
    DichF <- data.frame(
        group=c(
            colnames(metDataT), colnames(metDataN)
        ),
        cluster=c(
            rep(
                "Tumor",
                dim(metDataT)[2]
            ),
            rep(
                "Normal",
                dim(metDataN)[2]
            )
        ),
        stringsAsFactors = FALSE
    )

    ## Read hypometh Gplus TR link counts from step4:
    if(findhypoGpos==TRUE){

        hypo_Gplus_TR_links <- read.delim(
            file= '../step4/hypo.G+.output/hypo.G+.link.zscore.perm.all.optimized.links.txt',
            stringsAsFactors = FALSE
        )

    }

    ## Read hypometh Gminus TR link counts from step4:
    if(findhypoGneg==TRUE){

        hypo_Gminus_TR_links <- read.delim(
            file= '../step4/hypo.G-.output/hypo.G-.link.zscore.perm.all.optimized.links.txt',
            stringsAsFactors = FALSE
        )

    }

    ## Read hypermeth Gplus TR link counts from step4:
    if(findhyperGpos==TRUE){

        hyper_Gplus_TR_links <- read.delim(
            file= '../step4/hyper.G+.output/hyper.G+.link.zscore.perm.all.optimized.links.txt',
            stringsAsFactors = FALSE
        )

    }

    ## Read hypermeth Gminus TR link counts from step4:
    if(findhyperGneg==TRUE){

        hyper_Gminus_TR_links <- read.delim(
            file= '../step4/hyper.G-.output/hyper.G-.link.zscore.perm.all.optimized.links.txt',
            stringsAsFactors = FALSE
        )

    }

    ## Create a new directory for hypo.G+ scatterplot output and
    ## set the working directory ot that new folder:
    dir.create("probes.of.interest.output.scatterplot")
    setwd("./probes.of.interest.output.scatterplot")

    ## Now for each of the CpGs, we will create a dataframe noting the
    ## TRs that are linked to them, as well as the analysis type the link
    ## was found from:

    ## Index an empty list for the data frames:
    probe_df_list <- list()

    ## For each probe, create an individual dataframe with the linked TR
    ## and analysis type calling them:
    for(i in 1:length(probes_of_user_interest)){

        ## Index an empty vectors of linked TRs
        ## and the analyses they came from:
        linked_TR_vector <- character()
        linked_TR_dataset_vector <- character()

        ## For each of the 4 analysis types, if it was performed
        ## Get a vector of TRs linked to the given probe:
        if(exists('hypo_Gplus_TR_links')){

            hypo_Gplus_TRs <- hypo_Gplus_TR_links[
                hypo_Gplus_TR_links$probe==probes_of_user_interest[i],
                'geneID'
            ]

            ## Then add those TRs to the vector:
            linked_TR_vector <- c(
                linked_TR_vector,
                hypo_Gplus_TRs
            )

            ## And add the analysis type for each to its
            ## vector:
            linked_TR_dataset_vector <- c(
                linked_TR_dataset_vector,
                rep(
                    'Hypo.G+',
                    length(hypo_Gplus_TRs)
                )
            )

        }

        if(exists('hypo_Gminus_TR_links')){

            hypo_Gminus_TRs <- hypo_Gminus_TR_links[
                hypo_Gminus_TR_links$probe==probes_of_user_interest[i],
                'geneID'
            ]

            ## Then add those TRs to the vector:
            linked_TR_vector <- c(
                linked_TR_vector,
                hypo_Gminus_TRs
            )

            ## And add the analysis type for each to its
            ## vector:
            linked_TR_dataset_vector <- c(
                linked_TR_dataset_vector,
                rep(
                    'Hypo.G-',
                    length(hypo_Gminus_TRs)
                )
            )

        }

        if(exists('hyper_Gplus_TR_links')){

            hyper_Gplus_TRs <- hyper_Gplus_TR_links[
                hyper_Gplus_TR_links$probe==probes_of_user_interest[i],
                'geneID'
            ]

            ## Then add those TRs to the vector:
            linked_TR_vector <- c(
                linked_TR_vector,
                hyper_Gplus_TRs
            )

            ## And add the analysis type for each to its
            ## vector:
            linked_TR_dataset_vector <- c(
                linked_TR_dataset_vector,
                rep(
                    'Hyper.G+',
                    length(hyper_Gplus_TRs)
                )
            )

        }

        if(exists('hyper_Gminus_TR_links')){

            hyper_Gminus_TRs <- hyper_Gminus_TR_links[
                hyper_Gminus_TR_links$probe==probes_of_user_interest[i],
                'geneID'
            ]

            ## Then add those TRs to the vector:
            linked_TR_vector <- c(
                linked_TR_vector,
                hyper_Gminus_TRs
            )

            ## And add the analysis type for each to its
            ## vector:
            linked_TR_dataset_vector <- c(
                linked_TR_dataset_vector,
                rep(
                    'Hyper.G-',
                    length(hyper_Gminus_TRs)
                )
            )

        }

        ## Now create a dataframe listing each linked TR and the analysis it came from:
        linked_TR_dataframe <- data.frame(
            'linked_TR_ID'=linked_TR_vector,
            'linked_TR_analysis'= linked_TR_dataset_vector,
            stringsAsFactors = FALSE
        )

        ## Add the data frame to the list for each probe:
        probe_df_list[[i]] <- linked_TR_dataframe
    }

    ## Write a function to unpack
    probe_of_interest_scatterplot_function <- function(probe_dfs, CpG_of_interest){

        ## If the probe is found in the dataset, do graphing
        if((CpG_of_interest %in% rownames(metDataF_subC))==TRUE){

            ## Get the probe methylation values for the samples listed out
            ## in DichF
            unlisted_CpG_methylation <- unlist(
                c(
                    metDataF_subC[
                        CpG_of_interest,
                        DichF$group
                    ]
                )
            )

            ## List the genes of interest and analysis types of interest:
            TRs_of_interest <- probe_dfs$linked_TR_ID
            analysis_types <- probe_dfs$linked_TR_analysis

            ## Write a function to actually create the graph:
            internal_scatterplot_function <- function(TRs_linked_to_CpG, analysis_type_of_TR){

                ## Get the expression of the TR for the samples listed out
                ## in DichF
                TR_expression <- unlist(
                    c(
                        expDataF_subC[
                            TRs_linked_to_CpG,
                            DichF$group
                        ]
                    )
                )

                ## Get the name of the gene from its ID:
                gene_name <- gencode_v22_genes_df[
                    TRs_linked_to_CpG,
                    'gene_name'
                ]

                ## Manually coloring samples:
                t_v_n_group_colors <- c('Normal'='dodgerblue3', 'Tumor'='red3')

                ## Creating scatter with ggplot2:
                scatter_plot <- qplot(
                    x=TR_expression,
                    y=unlisted_CpG_methylation,
                    geom=c("point"),
                    colour=DichF$cluster
                )

                ## Create the plot:
                scatter_plot_updated <- scatter_plot +
                    ggtitle(
                        paste(
                            gene_name,
                            ' gene expression vs.\n',
                            CpG_of_interest,
                            ' DNA methylation',
                            sep=''
                        )
                    ) +
                    ylab(
                        paste(
                            CpG_of_interest,
                            'DNA methylation',
                            sep=' '
                        )
                    ) +
                    xlab(
                        paste(
                            gene_name,
                            'gene expression',
                            sep=' '
                        )
                    ) +
                    theme_bw() +
                    scale_color_manual(
                        values=t_v_n_group_colors,
                        name="Sample type"
                    ) +
                    theme(
                        plot.title = element_text(hjust=0.5, size=18),
                        legend.title = element_text(hjust=0.5, size=12),
                        legend.text = element_text(size=10, colour='black'),
                        panel.border = element_rect(colour = 'black', fill=NA, size=1),
                        axis.title.x = element_text(size=16),
                        axis.title.y = element_text(size=16),
                        axis.text.x = element_text(size=14, colour='black'),
                        axis.text.y = element_text(size=14, colour='black'),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()
                    )

                ## Create a title for the scatter plot pdf:
                ## This is a comination of the probe name with the linked gene:
                scatterplot_pdf_title <- paste(
                    './',
                    gene_name,
                    '_',
                    CpG_of_interest,
                    '_',
                    analysis_type_of_TR,
                    '_scatterplot.pdf',
                    sep=''
                )

                ## Open a pdf for saving the plot:
                pdf(
                    scatterplot_pdf_title,
                    height= 7,
                    width= 10
                )

                plot(scatter_plot_updated)

                ## Close the plot:
                dev.off()
            }

            ## Apply the graphing function:
            mapply(
                FUN= internal_scatterplot_function,
                TRs_linked_to_CpG= TRs_of_interest,
                analysis_type_of_TR= analysis_types
            )

        } else{

        }
    }

    mapply(
        FUN= probe_of_interest_scatterplot_function,
        probe_dfs= probe_df_list,
        CpG_of_interest= probes_of_user_interest
    )
}

make.SM.scatterplots.hypermeth.G_pos.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hyper.G+.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hyper.G+.output/", LS[1], sep=""), header=T)
    LM=list.files("../external.data/otherinfo", "SM")
    SM=read.delim(paste("../external.data/otherinfo/", LM[1], sep=""), header=T, sep="\t")
    dir.create("hyper.G+.output.complex.scatterplot")
    setwd("./hyper.G+.output.complex.scatterplot")
    SM=t(SM)
    M2=match(unique(TESTSR$geneSymbol),rownames(SM))
    SM.c=SM[na.omit(M2),]
    M4=match(colnames(metDataT), colnames(SM.c))
    SM.c2=as.matrix(SM.c)[,M4]
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        SM.sh=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(SM.c2)))==FALSE){
            SM.t=SM.c2[as.character(TESTSR[i,2]),]
            SM.t2=t(SM.t)
            SM.t3=c(SM.t2, rep(0,dim(metDataN)[2]))
            SM.sh=SM.t3
            SM.sh2=SM.sh
            SM.sh2[which(SM.sh==0)]=paste(16)
            SM.sh2[which(SM.sh==1)]=paste(8)
            CNV.sh2[which(SM.sh2==8)]=paste(8)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.SM.scatterplots.hypermeth.G_neg.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hyper.G-.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hyper.G-.output/", LS[1], sep=""), header=T)
    LM=list.files("../external.data/otherinfo", "SM")
    SM=read.delim(paste("../external.data/otherinfo/", LM[1], sep=""), header=T, sep="\t")
    dir.create("hyper.G-.output.complex.scatterplot")
    setwd("./hyper.G-.output.complex.scatterplot")
    SM=t(SM)
    M2=match(unique(TESTSR$geneSymbol),rownames(SM))
    SM.c=SM[na.omit(M2),]
    M4=match(colnames(metDataT), colnames(SM.c))
    SM.c2=as.matrix(SM.c)[,M4]
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        SM.sh=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(SM.c2)))==FALSE){
            SM.t=SM.c2[as.character(TESTSR[i,2]),]
            SM.t2=t(SM.t)
            SM.t3=c(SM.t2, rep(0,dim(metDataN)[2]))
            SM.sh=SM.t3
            SM.sh2=SM.sh
            SM.sh2[which(SM.sh==0)]=paste(16)
            SM.sh2[which(SM.sh==1)]=paste(8)
            CNV.sh2[which(SM.sh2==8)]=paste(8)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.SM.scatterplots.hypometh.G_pos.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hypo.G+.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hypo.G+.output/", LS[1], sep=""), header=T)
    LM=list.files("../external.data/otherinfo", "SM")
    SM=read.delim(paste("../external.data/otherinfo/", LM[1], sep=""), header=T, sep="\t")
    dir.create("hypo.G+.output.complex.scatterplot")
    setwd("./hypo.G+.output.complex.scatterplot")
    SM=t(SM)
    M2=match(unique(TESTSR$geneSymbol),rownames(SM))
    SM.c=SM[na.omit(M2),]
    M4=match(colnames(metDataT), colnames(SM.c))
    SM.c2=as.matrix(SM.c)[,M4]
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        SM.sh=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(SM.c2)))==FALSE){
            SM.t=SM.c2[as.character(TESTSR[i,2]),]
            SM.t2=t(SM.t)
            SM.t3=c(SM.t2, rep(0,dim(metDataN)[2]))
            SM.sh=SM.t3
            SM.sh2=SM.sh
            SM.sh2[which(SM.sh==0)]=paste(16)
            SM.sh2[which(SM.sh==1)]=paste(8)
            CNV.sh2[which(SM.sh2==8)]=paste(8)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.SM.scatterplots.hypometh.G_neg.source <- function(prefix, ggplot2, metDataT, metDataN, expDataT, expDataN, i) {
    load("../settings.rda")
    ##### let's obtain methylation and expression data if step1 was not performed #####
    load(paste("../step1/output/", prefix, ".diff.methylated.datasets.rda", sep=""))
    # library(ggplot2)
    DichF=data.frame(group=c(colnames(metDataT), colnames(metDataN)), cluster=c(rep("red", dim(metDataT)[2]), rep("blue", dim(metDataN)[2])))
    metDataF_subC=cbind(metDataT, metDataN)
    expDataF_subC=cbind(expDataT, expDataN)
    LS=list.files("../step4/hypo.G-.output/", pattern="all.optimized.links.txt")
    TESTSR=read.delim(paste("../step4/hypo.G-.output/", LS[1], sep=""), header=T)
    LM=list.files("../external.data/otherinfo", "SM")
    SM=read.delim(paste("../external.data/otherinfo/", LM[1], sep=""), header=T, sep="\t")
    dir.create("hypo.G-.output.complex.scatterplot")
    setwd("./hypo.G-.output.complex.scatterplot")
    SM=t(SM)
    M2=match(unique(TESTSR$geneSymbol),rownames(SM))
    SM.c=SM[na.omit(M2),]
    M4=match(colnames(metDataT), colnames(SM.c))
    SM.c2=as.matrix(SM.c)[,M4]
    for (i in 1:length(TESTSR[,1])){
        EXP_i=subset(expDataF_subC, rownames(expDataF_subC)==TESTSR[i,3])
        DNAM_i=subset(metDataF_subC, rownames(metDataF_subC)==TESTSR[i,1])
        DATA_i=data.frame(x=t(DNAM_i), y=t(EXP_i), colour=DichF$cluster)
        CNV.sh2=rep(16,dim(metDataF_subC)[2])
        SM.sh=rep(16,dim(metDataF_subC)[2])
        if (is.na(match(TESTSR[i,2],rownames(SM.c2)))==FALSE){
            SM.t=SM.c2[as.character(TESTSR[i,2]),]
            SM.t2=t(SM.t)
            SM.t3=c(SM.t2, rep(0,dim(metDataN)[2]))
            SM.sh=SM.t3
            SM.sh2=SM.sh
            SM.sh2[which(SM.sh==0)]=paste(16)
            SM.sh2[which(SM.sh==1)]=paste(8)
            CNV.sh2[which(SM.sh2==8)]=paste(8)
        }
        p=ggplot(DATA_i, aes(x=DATA_i[,1], y=DATA_i[,2]))+scale_shape_identity()+geom_point(colour=as.character(DATA_i$colour), size=3, shape=as.numeric(CNV.sh2))+xlab(paste("DNA methylation level (", names(DATA_i)[1], ")", sep=""))+ylab(paste("Gene expression level (", TESTSR[i,2], ", ID=",TESTSR[i,3],")", sep=""))+xlim(0,1)+theme_bw()+ggtitle(paste("Scatterplot for ", TESTSR[i,1],":", TESTSR[i,2], " (Pe=", signif(TESTSR[i,5],2), ", cor=",signif(cor(DATA_i[,1], DATA_i[,2], use="complete.obs", method="spearman"),2), ")", sep=""))
        ggsave(p, filename=paste(prefix, TESTSR[i,1], TESTSR[i,2], TESTSR[i,3], "anno.scatterplots.pdf", sep="."),  useDingbats=FALSE)
        print(i)
    }
    setwd("../")
}

make.summary.of.results.real.hypermeth.G_pos.source.FIN <- function(prefix) {
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
}

make.summary.of.results.real.hypermeth.G_neg.source.FIN <- function(prefix) {
    load("../settings.rda")
    links=c("hyper.G-")
    histcol=c("red")
    ################
    LS=list.files("../step4/hyper.G-.output/", pattern="all.optimized.links.txt")
    all=read.delim(paste("../step4/hyper.G-.output/", LS[1], sep=""), header=T)
    gene=read.delim("../scripts/data/gene.anno.hg38.txpt.bed", header=T)
    TF=read.delim("../scripts/data/human.TF.geneID.from.ELMERV2.txt", header=T)
    OC=read.delim("../scripts/data/human.527.known.cancer.gene.geneID.from.TheCancerGeneCensus.txt", header=T)
    TS=read.delim("../scripts/data/human.637.proteincoding.TSGs.geneID.from.TSGene.txt", header=T)
    probe=read.delim("../scripts/data/hm450cg_GEO.GR.bed", header=F)
    dir.create("hyper.G-.output")
    setwd("./hyper.G-.output")
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
}

make.summary.of.results.real.hypometh.G_pos.source.FIN <- function(prefix) {
    load("../settings.rda")
    links=c("hypo.G+")
    histcol=c("red")
    ################
    LS=list.files("../step4/hypo.G+.output/", pattern="all.optimized.links.txt")
    all=read.delim(paste("../step4/hypo.G+.output/", LS[1], sep=""), header=T)
    gene=read.delim("../scripts/data/gene.anno.hg38.txpt.bed", header=T)
    TF=read.delim("../scripts/data/human.TF.geneID.from.ELMERV2.txt", header=T)
    OC=read.delim("../scripts/data/human.527.known.cancer.gene.geneID.from.TheCancerGeneCensus.txt", header=T)
    TS=read.delim("../scripts/data/human.637.proteincoding.TSGs.geneID.from.TSGene.txt", header=T)
    probe=read.delim("../scripts/data/hm450cg_GEO.GR.bed", header=F)
    dir.create("hypo.G+.output")
    setwd("./hypo.G+.output")
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
}

make.summary.of.results.real.hypometh.G_neg.source.FIN <- function(prefix) {
    load("../settings.rda")
    links=c("hypo.G-")
    histcol=c("red")
    ################
    LS=list.files("../step4/hypo.G-.output/", pattern="all.optimized.links.txt")
    all=read.delim(paste("../step4/hypo.G-.output/", LS[1], sep=""), header=T)
    gene=read.delim("../scripts/data/gene.anno.hg38.txpt.bed", header=T)
    TF=read.delim("../scripts/data/human.TF.geneID.from.ELMERV2.txt", header=T)
    OC=read.delim("../scripts/data/human.527.known.cancer.gene.geneID.from.TheCancerGeneCensus.txt", header=T)
    TS=read.delim("../scripts/data/human.637.proteincoding.TSGs.geneID.from.TSGene.txt", header=T)
    probe=read.delim("../scripts/data/hm450cg_GEO.GR.bed", header=F)
    dir.create("hypo.G-.output")
    setwd("./hypo.G-.output")
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
}

make.survival.plots.hypermeth.G_pos <- function(prefix, survival_top_n_genes, clinical_data, expDataN, expDataT, metDataN, metDataT, visualize_survival_plots_genes, high_thresh, low_thresh, cores, visualize_survival_plots_probes) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Load ggplot2:
    # library('ggplot2')

    ## Load survival:
    # library('survival')

    ## Load parallel:
    # library('parallel')

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G+.output.histogram/",
            prefix,
            ".hyper.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hyper.G+.output/hyper.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Set the number of genes to analyze based on user input:
    survival_top_n_genes <- survival_top_n_genes

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_names <- ordered_TFs_by_link_count[
        1:survival_top_n_genes,
        'geneID'
    ]

    ## Load the external data if it is available:

    ## Check if .rda file with expression, methylation, and
    ## perhaps clinical data was placed in TENET folder:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## If the combined rda file exists, load it and see if clinical data is included:
    ## If not, check to see if clinical data is included in the "clinical" folder
    load(combined_rda_file)

    # #### REMOVE THIS LATER:
    # load(
    #     "/Volumes/NO_NAME/TENET_paper_update/tenet_match.matched.dedup.luad.methylation.expression.dataset.rda",
    # )

    ## If data is loaded as clinical_data, rename it to clinical:
    if(exists('clinical_data')){

        clinical <- clinical_data

    }

    ## Check to see if the clinical data is properly loaded
    ## if it isn't, check the clinical folder in external.data/data
    if(exists('clinical')){

        paste('clinical data loaded')

    } else{

        ## Get a path to the rda file containing clinical info
        ## in external.data/data/clinical
        clinical_rda_file <- list.files(
            path='../external.data/data/clinical',
            pattern='\\.rda$',
            full.names= TRUE,
            include.dirs = TRUE
        )

        if(length(clinical_rda_file)==1){

            ## Load the .rda file:
            load(clinical_rda_file)

            ## If data is loaded as clinical_data, rename it to clinical:
            if(exists('clinical_data')){

                clinical <- clinical_data

            }

            ## Check to see if the clinical data is now properly loaded
            if(exists('clinical')){

                paste('clinical data loaded')

            } else{

                ## CLinical data is not found in rda
                ## return an error message and abort script:
                paste('rda file did not contain clinical data or it was not saved as clinical or clinical_data')
                quit(save='no')
            }

        } else{

            ## Clinical rda not found.
            ## Return error message and abort script
            ## Later functionality may be added to try loading clinical data from txt file
            ## See step 0 make.rda.from.textfiles.R for reference:
            paste('rda file with clinical data not found')
            quit(save='no')
        }
    }

    ## Process the expression/methylation/clinical data so the function doesn't have to:

    ## combine datasets:
    expData <- cbind(expDataN, expDataT)
    metData <- cbind(metDataN, metDataT)

    ## Get the list of normal and tumor samples:
    exp_normal_samples <- colnames(expDataN)
    exp_tumor_samples <- colnames(expDataT)

    met_normal_samples <- colnames(metDataN)
    met_tumor_samples <- colnames(metDataT)

    ## Get the relevant columns from the clinical data:
    relevant_clinical <- clinical[
        ,c(
            "bcr_patient_barcode",
            "days_to_death",
            "days_to_last_followup",
            "vital_status"
        )
    ]
    rownames(relevant_clinical) <- relevant_clinical$bcr_patient_barcode

    # If there are subjects that are alive in the
    # dataset, set their days to death info to '-Inf'
    relevant_clinical[
        grep(
            "alive",
            relevant_clinical$vital_status,
            ignore.case = TRUE
        ),
        'days_to_death'
    ] <- '-Inf'

    # If there are subjects that are dead in the
    # dataset, set their days to last followup info to '-Inf'
    relevant_clinical[
        grep(
            "dead",
            relevant_clinical$vital_status,
            ignore.case = TRUE
        ),
        'days_to_last_followup'
    ] <- '-Inf'

    ## Remove subjects that remain unaccounted for:
    ## i.e. NA values remain in previous two columns:
    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_death"]
        ),
    ]

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_last_followup"]
        ),
    ]

    # Changing the days to last followup info in the cancer clinical data
    # to be numeric values through character values and remove any NAs that are induced
    # as a final check:
    relevant_clinical$days_to_death <- as.numeric(
        as.character(relevant_clinical$days_to_death)
    )

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_last_followup"]
        ),
    ]

    # Changing the days to death info in the cancer clinical data
    # to be numeric values through character values and remove any NAs that are induced
    # as a final check:
    relevant_clinical$days_to_last_followup <- as.numeric(
        as.character(relevant_clinical$days_to_last_followup)
    )

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_death"]
        ),
    ]

    # Add the relevant number as a final "time" variable
    # Equal to days to death for dead individuals, and
    # days to last followup for alive individuals
    relevant_clinical$time <- ifelse(
        relevant_clinical$vital_status=='Alive',
        relevant_clinical$days_to_last_followup,
        ifelse(
            relevant_clinical$vital_status=='Dead',
            relevant_clinical$days_to_death,
            NA
        )
    )

    ## Create a function to get the survival p-value or graph
    ## for each gene of interest:
    expression_survival_function_graph <- function(
        gene_of_interest,
        high_cutoff,
        low_cutoff,
        graph
    ){

        ## Get actual gene name depending on input:
        if(substring(gene_of_interest,1,4)=='ENSG' & nchar(gene_of_interest)==15){

            ## Input is in ENSG, get the gene name:
            gene_ENSG <- gene_of_interest

            gene_name <- gencode_v22_genes_df[
                gene_ENSG, 'gene_name'
            ]

            ## Get gene ENSG assuming a name is plugged in:
        } else{

            ## Assume what was given was the gene name, get the ENSG:
            gene_name <- gene_of_interest

            gene_ENSG <- rownames(
                gencode_v22_genes_df[gencode_v22_genes_df$gene_name==gene_name,]
            )
        }

        ## Get expression values for gene of interest:
        expression_values <- unlist(
            expData[
                gene_ENSG,
            ]
        )
        names(expression_values) <- colnames(expData)

        # Split gene expression values into normal and tumor samples
        tumor_expression_values <- expression_values[exp_tumor_samples]
        normal_expression_values <- expression_values[exp_normal_samples]

        # Changing the sample names to match the subject names:
        names(tumor_expression_values)  <- substr(
            names(tumor_expression_values),
            1,
            12
        )

        names(normal_expression_values)  <- substr(
            names(normal_expression_values),
            1,
            12
        )

        ## Calculate some basic data:
        normal_sample_n <- length(normal_expression_values)
        tumor_sample_n <- length(tumor_expression_values)

        ## Count the number of NA samples:
        NA_normal <- as.numeric(
            unname(
                table(
                    is.na(normal_expression_values)
                )['TRUE']
            )
        )

        NA_tumor <- as.numeric(
            unname(
                table(
                    is.na(tumor_expression_values)
                )['TRUE']
            )
        )

        if(is.na(NA_normal)){

            NA_normal <- 0

        }

        if(is.na(NA_tumor)){

            NA_tumor <- 0

        }

        ## Calculate mean expression for all samples present:
        mean_tumor_expression <- mean(
            tumor_expression_values,
            na.rm= TRUE
        )

        mean_normal_expression <- mean(
            normal_expression_values,
            na.rm= TRUE
        )

        ## Calculate the number of samples with present clinical data:
        present_clinical_normal_sample_n <- as.numeric(
            unname(
                table(
                    names(normal_expression_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        present_clinical_tumor_sample_n  <- as.numeric(
            unname(
                table(
                    names(tumor_expression_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        ## Identify the samples which are present in the clinical data:
        present_normal_samples <- names(normal_expression_values)[
            names(normal_expression_values) %in% rownames(relevant_clinical)
        ]

        present_tumor_samples <- names(tumor_expression_values)[
            names(tumor_expression_values) %in% rownames(relevant_clinical)
        ]


        # Subsetting clinical patient data for individuals with
        # gene expression data in the tumor set:
        function_relevant_clinical <- relevant_clinical[
            present_tumor_samples,
        ]

        ## Calculate quantiles:
        high_cutoff_quantile= quantile(
            tumor_expression_values,
            high_cutoff,
            na.rm= TRUE
        )[1]

        low_cutoff_quantile= quantile(
            tumor_expression_values,
            low_cutoff,
            na.rm= TRUE
        )[1]

        ## Determine if a sample is in the high, low, or intermediate quartiles:
        sample_expression_grouping <- ifelse(
            tumor_expression_values >  high_cutoff_quantile,
            "High",
            ifelse(
                tumor_expression_values <= low_cutoff_quantile,
                "Low",
                ifelse(
                    is.na(tumor_expression_values),
                    'Missing',
                    'Intermediate'
                )
            )
        )

        ## Add the expression and relevant clinical information to the relevant clinical dataframe:
        function_relevant_clinical$tumor_expression <- tumor_expression_values[
            rownames(function_relevant_clinical)
        ]

        ## Add the expression grouping to the clinical data:
        function_relevant_clinical$expression_grouping <- sample_expression_grouping[
            rownames(function_relevant_clinical)
        ]

        ## Remove samples missing expression data:
        function_relevant_clinical_complete <- na.omit(
            function_relevant_clinical,
            cols=c(
                "tumor_expression",
                "expression_grouping"
            )
        )

        function_relevant_clinical_complete <- function_relevant_clinical

        ## Count the number of samples of each group remaining:
        present_tumor_sample_high_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='High',
            ]
        )

        present_tumor_sample_intermediate_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Intermediate',
            ]
        )

        present_tumor_sample_low_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Low',
            ]
        )

        present_tumor_sample_missing_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Missing',
            ]
        )

        ## Get the mean expression in each group:
        present_tumor_sample_high_mean <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='High',
                'tumor_expression'
            ]
        )

        present_tumor_sample_intermediate_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Intermediate',
                'tumor_expression'
            ]
        )

        present_tumor_sample_low_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Low',
                'tumor_expression'
            ]
        )

        ## Get clinical results for only samples in the high and low group with info:
        function_relevant_clinical_complete_high_low <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
        ]

        ## Calculate proportion of deceased individuals in each group
        proportion_dead_high <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$expression_grouping=='High' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$expression_grouping=='High'),
                    ]
                )
        )

        proportion_dead_low <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$expression_grouping=='Low' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$expression_grouping=='Low'),
                    ]
                )
        )

        ## Determine which group had greater proportion of fatalities:
        highest_dead_proportion_group <- ifelse(
            proportion_dead_high > proportion_dead_low,
            'high_expression_low_survival',
            ifelse(
                proportion_dead_high < proportion_dead_low,
                'low_expression_low_survival',
                'unclear'
            )
        )

        ## Create a high_low survival object:
        high_low_survival_object <- Surv(
            function_relevant_clinical_complete_high_low$time,
            ifelse(
                function_relevant_clinical_complete_high_low$vital_status=='Alive',
                FALSE,
                TRUE
            )
        )

        # Set the rownames of the survival object to be the
        # patient names for the high/low group smaples:
        rownames(high_low_survival_object) <- rownames(
            function_relevant_clinical_complete[
                !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
            ]
        )

        # Combine the top and down vectors into a single vector
        # in the order top-down
        expression_group <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
            'expression_grouping'
        ]

        # Creating legend names for high/low expression groups:
        legend_name_high <- paste(gene_name,"high")

        legend_name_low <- paste(gene_name,"low")

        # Perform the survival analysis.
        # This uses the combined vector as the x variable
        # and the survival data as the y variable to create
        # a table with information about the test
        # including chi-squared p-value
        survival_table <- survdiff(
            high_low_survival_object  ~ expression_group
        )

        # Get the chi-squared test statistic from the analysis above:
        survival_table_chi_squared <- survival_table$chisq

        # Calculating a p value based on the test statistic to get
        # a precise p-value
        survival_pvalue <- as.numeric(
            1 - pchisq(
                abs(survival_table_chi_squared), df = 1
            )
        )

        ## Round the p-value on the graph to 3 digits:
        survival_pvalue_formatted <- formatC(
            survival_pvalue,
            format='e',
            digits=3
        )

        ## Create the graph if specified:
        if(graph==TRUE){

            # Create the plot title
            # with gene name and p-value included:
            survival_title <- paste(
                gene_name,
                "\nKaplan-Meier Survival analysis\nP = ",
                survival_pvalue_formatted,
                sep=''
            )

            ## Get the path to the step 5 folder to save
            ## the pdf in:
            path_to_folder <- getwd()

            ## Create a title for the survival plot pdf:
            survival_plot_pdf_title <- paste(
                path_to_folder,
                '/hyper.G+.output.survival/',
                gene_name,
                '_survival_plot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(survival_plot_pdf_title)

            ## Now actually create the survival plot
            ## Using a similar structure to what was used to generate the p-value:
            plot(
                survfit(
                    high_low_survival_object ~ expression_group
                ),

                # Color the lines (high in red first!)
                col = c('red', 'black'),

                ## Add thickness to the lines:
                lwd=3,

                # Use the title that was created earlier as the title
                # of the plot:
                main= survival_title,

                # Set titles of the x and y axis:
                # Note: TCGA measures survival in days as noted:
                xlab="Time (days)",
                ylab="Probability",

                # Change axis size:
                cex.axis=1,
                cex.main=1.5,
                cex.lab=1.25
            )

            ## Add a legend to the plot:
            legend(

                # Set X position of legend in graph:
                x= (
                    max(function_relevant_clinical_complete_high_low$days_to_last_followup) - 2500
                ),

                # Set Y position of legend in graph
                y= 1,

                # Use the legend titles that were created earlier
                legend= c(legend_name_high, legend_name_low),

                # As above, use black for low and red for high:
                col= c('red', 'black'),

                # Coloring the text in the legend as well
                text.col= c('red', 'black'),

                # Change the shape of the labels in the legend:
                pch= 15
            )

            ## Close the plot:
            dev.off()

        } else if(graph==FALSE){

            ## create vector of relevant info:
            relevant_vector <- c(
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(NA_normal),
                as.numeric(NA_tumor),
                as.numeric(mean_normal_expression),
                as.numeric(mean_tumor_expression),
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(present_tumor_sample_missing_n),
                as.numeric(present_tumor_sample_low_n),
                as.numeric(present_tumor_sample_intermediate_n),
                as.numeric(present_tumor_sample_high_n),
                as.numeric(present_tumor_sample_low_mean),
                as.numeric(present_tumor_sample_intermediate_mean),
                as.numeric(present_tumor_sample_high_mean),
                as.numeric(proportion_dead_low),
                as.numeric(proportion_dead_high),
                highest_dead_proportion_group,
                as.numeric(survival_pvalue)
            )
            names(relevant_vector) <- c(
                'normal_sample_count',
                'tumor_sample_count',
                'normal_sample_count_missing',
                'tumor_sample_count_missing',
                'mean_normal_expression',
                'mean_tumor_expression',
                'normal_sample_with_clinical_count',
                'tumor_sample_with_clinical_count',
                'tumor_sample_with_clinical_NA_count',
                'tumor_sample_with_clinical_low_count',
                'tumor_sample_with_clinical_intermediate_count',
                'tumor_sample_with_clinical_high_count',
                'mean_tumor_with_clinical_low_expression',
                'mean_tumor_with_clinical_intermediate_expression',
                'mean_tumor_with_clinical_high_expression',
                'proportion_dead_in_low_expression',
                'proportion_dead_in_high_expression',
                'survival_direction_of_effect',
                'survival_p_value'
            )

            ## Return the vector
            return(relevant_vector)
        }

    }

    ## Create a function to get the survival p-value or graph for
    ##beach probe of interest:
    methylation_survival_function_graph <- function(
        probe_of_interest,
        high_cutoff,
        low_cutoff,
        graph
    ){

        ## Get methylation values for probe of interest:
        methylation_values <- unlist(
            metData[
                probe_of_interest,
            ]
        )
        names(methylation_values) <- colnames(metData)

        # Split methylation values into normal and tumor samples
        # If both are present in the dataset
        tumor_methylation_values <- methylation_values[met_tumor_samples]
        normal_methylation_values <- methylation_values[met_normal_samples]

        # Changing the sample names to match the subject names:
        names(tumor_methylation_values)  <- substr(
            names(tumor_methylation_values),
            1,
            12
        )

        names(normal_methylation_values)  <- substr(
            names(normal_methylation_values),
            1,
            12
        )

        ## Calculate some basic data:
        normal_sample_n <- length(normal_methylation_values)
        tumor_sample_n <- length(tumor_methylation_values)

        ## Count the number of NA samples:
        NA_normal <- as.numeric(
            unname(
                table(
                    is.na(normal_methylation_values)
                )['TRUE']
            )
        )

        NA_tumor <- as.numeric(
            unname(
                table(
                    is.na(tumor_methylation_values)
                )['TRUE']
            )
        )

        if(is.na(NA_normal)){

            NA_normal <- 0

        }

        if(is.na(NA_tumor)){

            NA_tumor <- 0

        }

        ## Calculate mean methylation for all samples present:
        mean_tumor_methylation <- mean(
            tumor_methylation_values,
            na.rm= TRUE
        )

        mean_normal_methylation <- mean(
            normal_methylation_values,
            na.rm= TRUE
        )

        ## Calculate the number of samples with present clinical data:
        present_clinical_normal_sample_n <- as.numeric(
            unname(
                table(
                    names(normal_methylation_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        present_clinical_tumor_sample_n  <- as.numeric(
            unname(
                table(
                    names(tumor_methylation_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        ## Identify the samples which are present in the clinical data:
        present_normal_samples <- names(normal_methylation_values)[
            names(normal_methylation_values) %in% rownames(relevant_clinical)
        ]

        present_tumor_samples <- names(tumor_methylation_values)[
            names(tumor_methylation_values) %in% rownames(relevant_clinical)
        ]


        # Subsetting clinical patient data for individuals with
        # DNA methylation data in the tumor set:
        function_relevant_clinical <- relevant_clinical[
            present_tumor_samples,
        ]

        ## Calculate quantiles:
        high_cutoff_quantile= quantile(
            tumor_methylation_values,
            high_cutoff,
            na.rm= TRUE
        )[1]

        low_cutoff_quantile= quantile(
            tumor_methylation_values,
            low_cutoff,
            na.rm= TRUE
        )[1]

        ## Determine if a sample is in the high, low, or intermediate quartiles:
        sample_methylation_grouping <- ifelse(
            tumor_methylation_values >  high_cutoff_quantile,
            "High",
            ifelse(
                tumor_methylation_values <= low_cutoff_quantile,
                "Low",
                ifelse(
                    is.na(tumor_methylation_values),
                    'Missing',
                    'Intermediate'
                )
            )
        )

        ## Add the methylation and relevant clinical information to the relevant clinical dataframe:
        function_relevant_clinical$tumor_methylation <- tumor_methylation_values[
            rownames(function_relevant_clinical)
        ]

        ## Add the methylation grouping to the clinical data:
        function_relevant_clinical$methylation_grouping <- sample_methylation_grouping[
            rownames(function_relevant_clinical)
        ]

        ## Remove samples missing methylation data:
        function_relevant_clinical_complete <- na.omit(
            function_relevant_clinical,
            cols=c(
                "tumor_methylation",
                "methylation_grouping"
            )
        )
        function_relevant_clinical_complete <- function_relevant_clinical

        ## Count the number of samples of each group remaining:
        present_tumor_sample_high_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='High',
            ]
        )

        present_tumor_sample_intermediate_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Intermediate',
            ]
        )

        present_tumor_sample_low_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Low',
            ]
        )

        present_tumor_sample_missing_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Missing',
            ]
        )

        ## Get the mean methylation in each group:
        present_tumor_sample_high_mean <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='High',
                'tumor_methylation'
            ]
        )

        present_tumor_sample_intermediate_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Intermediate',
                'tumor_methylation'
            ]
        )

        present_tumor_sample_low_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Low',
                'tumor_methylation'
            ]
        )

        ## Get clinical results for only samples in the high and low group with info:
        function_relevant_clinical_complete_high_low <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
        ]

        ## Calculate proportion of deceased individuals in each group
        proportion_dead_high <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$methylation_grouping=='High' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$methylation_grouping=='High'),
                    ]
                )
        )

        proportion_dead_low <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$methylation_grouping=='Low' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$methylation_grouping=='Low'),
                    ]
                )
        )

        ## Determine which group had greater proportion of fatalities:
        highest_dead_proportion_group <- ifelse(
            proportion_dead_high > proportion_dead_low,
            'high_methylation_low_survival',
            ifelse(
                proportion_dead_high < proportion_dead_low,
                'low_methylation_low_survival',
                'unclear'
            )
        )

        ## Create a high_low survival object:
        high_low_survival_object <- Surv(
            function_relevant_clinical_complete_high_low$time,
            ifelse(
                function_relevant_clinical_complete_high_low$vital_status=='Alive',
                FALSE,
                TRUE
            )
        )

        # Set the rownames of the survival object to be the
        # patient names for the high/low group smaples:
        rownames(high_low_survival_object) <- rownames(
            function_relevant_clinical_complete[
                !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
            ]
        )

        # Combine the top and down vectors into a single vector
        # in the order top-down
        methylation_group <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
            'methylation_grouping'
        ]

        # Creating legend names for high/low methylation groups:
        legend_name_high <- paste(probe_of_interest,"high")

        legend_name_low <- paste(probe_of_interest,"low")

        # Perform the survival analysis.
        # This uses the combined vector as the x variable
        # and the survival data as the y variable to create
        # a table with information about the test
        # including chi-squared p-value
        survival_table <- survdiff(
            high_low_survival_object  ~ methylation_group
        )

        # Get the chi-squared test statistic from the analysis above:
        survival_table_chi_squared <- survival_table$chisq

        # Calculating a p value based on the test statistic to get
        # a precise p-value
        survival_pvalue <- as.numeric(
            1 - pchisq(
                abs(survival_table_chi_squared), df = 1
            )
        )

        ## Round the p-value on the graph to 3 digits:
        survival_pvalue_formatted <- formatC(
            survival_pvalue,
            format='e',
            digits=3
        )

        ## Create the graph if specified:
        if(graph==TRUE){

            # Create the plot title
            # with gene name and p-value included:
            survival_title <- paste(
                probe_of_interest,
                "\nKaplan-Meier Survival analysis\nP = ",
                survival_pvalue_formatted,
                sep=''
            )

            ## Get the path to the step 5 folder to save
            ## the pdf in:
            path_to_folder <- getwd()

            ## Create a title for the survival plot pdf:
            ## This is a comination of the probe name with the linked gene:
            survival_plot_pdf_title <- paste(
                path_to_folder,
                '/hyper.G+.output.survival/',
                probe_of_interest,
                '_survival_plot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(survival_plot_pdf_title)

            ## Now actually create the survival plot
            ## Using a similar structure to what was used to generate the p-value:
            plot(
                survfit(
                    high_low_survival_object ~ methylation_group
                ),

                # Color the lines (high in red first!)
                col = c('red', 'black'),

                ## Add thickness to the lines:
                lwd=3,

                # Use the title that was created earlier as the title
                # of the plot:
                main= survival_title,

                # Set titles of the x and y axis:
                # Note: TCGA measures survival in days as noted:
                xlab="Time (days)",
                ylab="Probability",

                # Change axis size:
                cex.axis=1,
                cex.main=1.5,
                cex.lab=1.25
            )

            ## Add a legend to the plot:
            legend(

                # Set X position of legend in graph:
                x= (
                    max(function_relevant_clinical_complete_high_low$days_to_last_followup) - 2500
                ),

                # Set Y position of legend in graph
                y= 1,

                # Use the legend titles that were created earlier
                legend= c(legend_name_high, legend_name_low),

                # As above, use black for low and red for high:
                col= c('red', 'black'),

                # Coloring the text in the legend as well
                text.col= c('red', 'black'),

                # Change the shape of the labels in the legend:
                pch= 15
            )

            ## Close the plot:
            dev.off()

        } else if(graph==FALSE){

            ## create vector of relevant info:
            relevant_vector <- c(
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(NA_normal),
                as.numeric(NA_tumor),
                as.numeric(mean_normal_methylation),
                as.numeric(mean_tumor_methylation),
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(present_tumor_sample_missing_n),
                as.numeric(present_tumor_sample_low_n),
                as.numeric(present_tumor_sample_intermediate_n),
                as.numeric(present_tumor_sample_high_n),
                as.numeric(present_tumor_sample_low_mean),
                as.numeric(present_tumor_sample_intermediate_mean),
                as.numeric(present_tumor_sample_high_mean),
                as.numeric(proportion_dead_low),
                as.numeric(proportion_dead_high),
                highest_dead_proportion_group,
                as.numeric(survival_pvalue)
            )
            names(relevant_vector) <- c(
                'normal_sample_count',
                'tumor_sample_count',
                'normal_sample_count_missing',
                'tumor_sample_count_missing',
                'mean_normal_methylation',
                'mean_tumor_methylation',
                'normal_sample_with_clinical_count',
                'tumor_sample_with_clinical_count',
                'tumor_sample_with_clinical_NA_count',
                'tumor_sample_with_clinical_low_count',
                'tumor_sample_with_clinical_intermediate_count',
                'tumor_sample_with_clinical_high_count',
                'mean_tumor_with_clinical_low_methylation',
                'mean_tumor_with_clinical_intermediate_methylation',
                'mean_tumor_with_clinical_high_methylation',
                'proportion_dead_in_low_methylation',
                'proportion_dead_in_high_methylation',
                'survival_direction_of_effect',
                'survival_p_value'
            )

            ## Return the vector
            return(relevant_vector)
        }

    }

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_genes==TRUE){

        ## Create survival plots for each of the genes designated:
        mclapply(
            X= top_gene_names,
            FUN= expression_survival_function_graph,
            high_cutoff= high_thresh,
            low_cutoff= low_thresh,
            graph= TRUE,
            mc.cores= cores
        )

    }

    ## Get survival information for each gene:
    gene_survival_results_list <- mclapply(
        X= top_gene_names,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= FALSE,
        mc.cores= cores
    )

    ## Transform the list into a matrix:
    gene_survival_results_df <- data.frame(
        matrix(
            unlist(gene_survival_results_list),
            nrow=length(gene_survival_results_list),
            byrow=T
        ),
        stringsAsFactors = FALSE
    )

    colnames(gene_survival_results_df) <- c(
        'normal_sample_count',
        'tumor_sample_count',
        'normal_sample_count_missing',
        'tumor_sample_count_missing',
        'mean_normal_expression',
        'mean_tumor_expression',
        'normal_sample_with_clinical_count',
        'tumor_sample_with_clinical_count',
        'tumor_sample_with_clinical_NA_count',
        'tumor_sample_with_clinical_low_count',
        'tumor_sample_with_clinical_intermediate_count',
        'tumor_sample_with_clinical_high_count',
        'mean_tumor_with_clinical_low_expression',
        'mean_tumor_with_clinical_intermediate_expression',
        'mean_tumor_with_clinical_high_expression',
        'proportion_dead_in_low_expression',
        'proportion_dead_in_high_expression',
        'survival_direction_of_effect',
        'survival_p_value'
    )

    rownames(gene_survival_results_df) <- c(
        gene_name <- gencode_v22_genes_df[
            top_gene_names, 'gene_name'
        ]
    )

    ## Save the gene information as a .tsv:
    write.table(
        gene_survival_results_df,
        file='./hyper.G+.output.survival/hyper.G+.top_genes_survival_info.tsv',
        quote= FALSE,
        sep='\t'
    )

    ## Now for each of the probes, we will need to calculate survival info and graphs:

    ## Get the list of probes linked to the significant genes:
    CpGs_linked <- CpG_linkage_dataset[
        CpG_linkage_dataset$geneID %in% top_gene_names,
        'probe'
    ]

    ## Get the unique probes from the list:
    unique_CpGs_linked <- unique(CpGs_linked)

    ## Get the survival info for each probe by mclapplying the
    ## methylation survival function to all of the probes
    ## linked to each gene:
    probe_survival_results_list <- mclapply(
        X= unique_CpGs_linked,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= FALSE,
        mc.cores= cores
    )

    ## Transform the list into a matrix:
    probe_survival_results_df <- data.frame(
        matrix(
            unlist(probe_survival_results_list),
            nrow=length(probe_survival_results_list),
            byrow=T
        ),
        stringsAsFactors = FALSE
    )

    colnames(probe_survival_results_df) <- c(
        'normal_sample_count',
        'tumor_sample_count',
        'normal_sample_count_missing',
        'tumor_sample_count_missing',
        'mean_normal_methylation',
        'mean_tumor_methylation',
        'normal_sample_with_clinical_count',
        'tumor_sample_with_clinical_count',
        'tumor_sample_with_clinical_NA_count',
        'tumor_sample_with_clinical_low_count',
        'tumor_sample_with_clinical_intermediate_count',
        'tumor_sample_with_clinical_high_count',
        'mean_tumor_with_clinical_low_methylation',
        'mean_tumor_with_clinical_intermediate_methylation',
        'mean_tumor_with_clinical_high_methylation',
        'proportion_dead_in_low_methylation',
        'proportion_dead_in_high_methylation',
        'survival_direction_of_effect',
        'survival_p_value'
    )

    rownames(probe_survival_results_df) <- unique_CpGs_linked

    ## For each of the CpGs, list which of the top X genes it was listed to:

    ## Wite a function to get the top X gene names a probe is assigned to:
    top_gene_assignment_function <- function(probe_of_interest){

        all_listed_genes <- CpG_linkage_dataset[
            CpG_linkage_dataset$probe==probe_of_interest,
            'geneSymbol'
        ]

        ## Conver thte top X genes to their gene symbols
        top_gene_symbols <- gencode_v22_genes_df[
            top_gene_names, 'gene_name'
        ]

        ## Get the names of the top genes that were linked to this probe:
        top_gene_symbols_linked_to_probe <- top_gene_symbols[
            top_gene_symbols %in% all_listed_genes
        ]

        ## Return the listed genes:
        return(
            paste(
                top_gene_symbols_linked_to_probe,
                collapse=','
            )
        )
    }

    ## Add the listed top genes to the survival results dataframe
    ## for the probes:
    probe_survival_results_df$top_genes_linked <- unname(
        sapply(
            rownames(probe_survival_results_df),
            top_gene_assignment_function
        )
    )

    ## Create a vector of the CpG p-values and add the CpG entries
    ## as names
    CpGs_linked_p_values_vector <- probe_survival_results_df$survival_p_value
    names(CpGs_linked_p_values_vector) <- rownames(probe_survival_results_df)

    ## Get the CpGs that are nominally significant:
    CpGs_linked_p_values_nominally_significant <- CpGs_linked_p_values_vector[
        CpGs_linked_p_values_vector<0.05
    ]

    ## Order the nominally significant CpGs:
    CpGs_linked_p_values_nominally_significant_ordered <- sort(
        CpGs_linked_p_values_nominally_significant,
        decreasing = FALSE
    )

    ## Get the names of the nominally significant CpGs:
    CpGs_linked_nominally_significant_names <- names(
        CpGs_linked_p_values_nominally_significant_ordered
    )

    ## If visualize_survival_plots_probes is TRUE, generate the plots:
    if(visualize_survival_plots_probes==TRUE){

        #### Save plots for the nominally significant CpGs:
        mclapply(
            X= CpGs_linked_nominally_significant_names,
            FUN= methylation_survival_function_graph,
            high_cutoff= high_thresh,
            low_cutoff= low_thresh,
            graph= TRUE,
            mc.cores= cores
        )

    }

    ## Save the probe information as a .tsv:
    write.table(
        probe_survival_results_df,
        file='./hyper.G+.output.survival/hyper.G+.probes_linked_to_top_genes_survival_info.tsv',
        quote= FALSE,
        sep='\t'
    )
}

make.survival.plots.hypermeth.G_neg <- function(prefix, survival_top_n_genes, clinical_data, expDataN, expDataT, metDataN, metDataT, visualize_survival_plots_genes, high_thresh, low_thresh, cores, visualize_survival_plots_probes) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Load ggplot2:
    # library('ggplot2')

    ## Load survival:
    # library('survival')

    ## Load parallel:
    # library('parallel')

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G-.output.histogram/",
            prefix,
            ".hyper.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hyper.G-.output/hyper.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Set the number of genes to analyze based on user input:
    survival_top_n_genes <- survival_top_n_genes

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_names <- ordered_TFs_by_link_count[
        1:survival_top_n_genes,
        'geneID'
    ]

    ## Load the external data if it is available:

    ## Check if .rda file with expression, methylation, and
    ## perhaps clinical data was placed in TENET folder:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## If the combined rda file exists, load it and see if clinical data is included:
    ## If not, check to see if clinical data is included in the "clinical" folder
    load(combined_rda_file)

    # #### REMOVE THIS LATER:
    # load(
    #     "/Volumes/NO_NAME/TENET_paper_update/tenet_match.matched.dedup.luad.methylation.expression.dataset.rda",
    # )

    ## If data is loaded as clinical_data, rename it to clinical:
    if(exists('clinical_data')){

        clinical <- clinical_data

    }

    ## Check to see if the clinical data is properly loaded
    ## if it isn't, check the clinical folder in external.data/data
    if(exists('clinical')){

        paste('clinical data loaded')

    } else{

        ## Get a path to the rda file containing clinical info
        ## in external.data/data/clinical
        clinical_rda_file <- list.files(
            path='../external.data/data/clinical',
            pattern='\\.rda$',
            full.names= TRUE,
            include.dirs = TRUE
        )

        if(length(clinical_rda_file)==1){

            ## Load the .rda file:
            load(clinical_rda_file)

            ## If data is loaded as clinical_data, rename it to clinical:
            if(exists('clinical_data')){

                clinical <- clinical_data

            }

            ## Check to see if the clinical data is now properly loaded
            if(exists('clinical')){

                paste('clinical data loaded')

            } else{

                ## CLinical data is not found in rda
                ## return an error message and abort script:
                paste('rda file did not contain clinical data or it was not saved as clinical or clinical_data')
                quit(save='no')
            }

        } else{

            ## Clinical rda not found.
            ## Return error message and abort script
            ## Later functionality may be added to try loading clinical data from txt file
            ## See step 0 make.rda.from.textfiles.R for reference:
            paste('rda file with clinical data not found')
            quit(save='no')
        }
    }

    ## Process the expression/methylation/clinical data so the function doesn't have to:

    ## combine datasets:
    expData <- cbind(expDataN, expDataT)
    metData <- cbind(metDataN, metDataT)

    ## Get the list of normal and tumor samples:
    exp_normal_samples <- colnames(expDataN)
    exp_tumor_samples <- colnames(expDataT)

    met_normal_samples <- colnames(metDataN)
    met_tumor_samples <- colnames(metDataT)

    ## Get the relevant columns from the clinical data:
    relevant_clinical <- clinical[
        ,c(
            "bcr_patient_barcode",
            "days_to_death",
            "days_to_last_followup",
            "vital_status"
        )
    ]
    rownames(relevant_clinical) <- relevant_clinical$bcr_patient_barcode

    # If there are subjects that are alive in the
    # dataset, set their days to death info to '-Inf'
    relevant_clinical[
        grep(
            "alive",
            relevant_clinical$vital_status,
            ignore.case = TRUE
        ),
        'days_to_death'
    ] <- '-Inf'

    # If there are subjects that are dead in the
    # dataset, set their days to last followup info to '-Inf'
    relevant_clinical[
        grep(
            "dead",
            relevant_clinical$vital_status,
            ignore.case = TRUE
        ),
        'days_to_last_followup'
    ] <- '-Inf'

    ## Remove subjects that remain unaccounted for:
    ## i.e. NA values remain in previous two columns:
    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_death"]
        ),
    ]

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_last_followup"]
        ),
    ]

    # Changing the days to last followup info in the cancer clinical data
    # to be numeric values through character values and remove any NAs that are induced
    # as a final check:
    relevant_clinical$days_to_death <- as.numeric(
        as.character(relevant_clinical$days_to_death)
    )

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_last_followup"]
        ),
    ]

    # Changing the days to death info in the cancer clinical data
    # to be numeric values through character values and remove any NAs that are induced
    # as a final check:
    relevant_clinical$days_to_last_followup <- as.numeric(
        as.character(relevant_clinical$days_to_last_followup)
    )

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_death"]
        ),
    ]

    # Add the relevant number as a final "time" variable
    # Equal to days to death for dead individuals, and
    # days to last followup for alive individuals
    relevant_clinical$time <- ifelse(
        relevant_clinical$vital_status=='Alive',
        relevant_clinical$days_to_last_followup,
        ifelse(
            relevant_clinical$vital_status=='Dead',
            relevant_clinical$days_to_death,
            NA
        )
    )

    ## Create a function to get the survival p-value or graph
    ## for each gene of interest:
    expression_survival_function_graph <- function(
        gene_of_interest,
        high_cutoff,
        low_cutoff,
        graph
    ){

        ## Get actual gene name depending on input:
        if(substring(gene_of_interest,1,4)=='ENSG' & nchar(gene_of_interest)==15){

            ## Input is in ENSG, get the gene name:
            gene_ENSG <- gene_of_interest

            gene_name <- gencode_v22_genes_df[
                gene_ENSG, 'gene_name'
            ]

            ## Get gene ENSG assuming a name is plugged in:
        } else{

            ## Assume what was given was the gene name, get the ENSG:
            gene_name <- gene_of_interest

            gene_ENSG <- rownames(
                gencode_v22_genes_df[gencode_v22_genes_df$gene_name==gene_name,]
            )
        }

        ## Get expression values for gene of interest:
        expression_values <- unlist(
            expData[
                gene_ENSG,
            ]
        )
        names(expression_values) <- colnames(expData)

        # Split gene expression values into normal and tumor samples
        tumor_expression_values <- expression_values[exp_tumor_samples]
        normal_expression_values <- expression_values[exp_normal_samples]

        # Changing the sample names to match the subject names:
        names(tumor_expression_values)  <- substr(
            names(tumor_expression_values),
            1,
            12
        )

        names(normal_expression_values)  <- substr(
            names(normal_expression_values),
            1,
            12
        )

        ## Calculate some basic data:
        normal_sample_n <- length(normal_expression_values)
        tumor_sample_n <- length(tumor_expression_values)

        ## Count the number of NA samples:
        NA_normal <- as.numeric(
            unname(
                table(
                    is.na(normal_expression_values)
                )['TRUE']
            )
        )

        NA_tumor <- as.numeric(
            unname(
                table(
                    is.na(tumor_expression_values)
                )['TRUE']
            )
        )

        if(is.na(NA_normal)){

            NA_normal <- 0

        }

        if(is.na(NA_tumor)){

            NA_tumor <- 0

        }

        ## Calculate mean expression for all samples present:
        mean_tumor_expression <- mean(
            tumor_expression_values,
            na.rm= TRUE
        )

        mean_normal_expression <- mean(
            normal_expression_values,
            na.rm= TRUE
        )

        ## Calculate the number of samples with present clinical data:
        present_clinical_normal_sample_n <- as.numeric(
            unname(
                table(
                    names(normal_expression_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        present_clinical_tumor_sample_n  <- as.numeric(
            unname(
                table(
                    names(tumor_expression_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        ## Identify the samples which are present in the clinical data:
        present_normal_samples <- names(normal_expression_values)[
            names(normal_expression_values) %in% rownames(relevant_clinical)
        ]

        present_tumor_samples <- names(tumor_expression_values)[
            names(tumor_expression_values) %in% rownames(relevant_clinical)
        ]


        # Subsetting clinical patient data for individuals with
        # gene expression data in the tumor set:
        function_relevant_clinical <- relevant_clinical[
            present_tumor_samples,
        ]

        ## Calculate quantiles:
        high_cutoff_quantile= quantile(
            tumor_expression_values,
            high_cutoff,
            na.rm= TRUE
        )[1]

        low_cutoff_quantile= quantile(
            tumor_expression_values,
            low_cutoff,
            na.rm= TRUE
        )[1]

        ## Determine if a sample is in the high, low, or intermediate quartiles:
        sample_expression_grouping <- ifelse(
            tumor_expression_values >  high_cutoff_quantile,
            "High",
            ifelse(
                tumor_expression_values <= low_cutoff_quantile,
                "Low",
                ifelse(
                    is.na(tumor_expression_values),
                    'Missing',
                    'Intermediate'
                )
            )
        )

        ## Add the expression and relevant clinical information to the relevant clinical dataframe:
        function_relevant_clinical$tumor_expression <- tumor_expression_values[
            rownames(function_relevant_clinical)
        ]

        ## Add the expression grouping to the clinical data:
        function_relevant_clinical$expression_grouping <- sample_expression_grouping[
            rownames(function_relevant_clinical)
        ]

        ## Remove samples missing expression data:
        function_relevant_clinical_complete <- na.omit(
            function_relevant_clinical,
            cols=c(
                "tumor_expression",
                "expression_grouping"
            )
        )

        function_relevant_clinical_complete <- function_relevant_clinical

        ## Count the number of samples of each group remaining:
        present_tumor_sample_high_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='High',
            ]
        )

        present_tumor_sample_intermediate_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Intermediate',
            ]
        )

        present_tumor_sample_low_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Low',
            ]
        )

        present_tumor_sample_missing_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Missing',
            ]
        )

        ## Get the mean expression in each group:
        present_tumor_sample_high_mean <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='High',
                'tumor_expression'
            ]
        )

        present_tumor_sample_intermediate_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Intermediate',
                'tumor_expression'
            ]
        )

        present_tumor_sample_low_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Low',
                'tumor_expression'
            ]
        )

        ## Get clinical results for only samples in the high and low group with info:
        function_relevant_clinical_complete_high_low <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
        ]

        ## Calculate proportion of deceased individuals in each group
        proportion_dead_high <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$expression_grouping=='High' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$expression_grouping=='High'),
                    ]
                )
        )

        proportion_dead_low <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$expression_grouping=='Low' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$expression_grouping=='Low'),
                    ]
                )
        )

        ## Determine which group had greater proportion of fatalities:
        highest_dead_proportion_group <- ifelse(
            proportion_dead_high > proportion_dead_low,
            'high_expression_low_survival',
            ifelse(
                proportion_dead_high < proportion_dead_low,
                'low_expression_low_survival',
                'unclear'
            )
        )

        ## Create a high_low survival object:
        high_low_survival_object <- Surv(
            function_relevant_clinical_complete_high_low$time,
            ifelse(
                function_relevant_clinical_complete_high_low$vital_status=='Alive',
                FALSE,
                TRUE
            )
        )

        # Set the rownames of the survival object to be the
        # patient names for the high/low group smaples:
        rownames(high_low_survival_object) <- rownames(
            function_relevant_clinical_complete[
                !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
            ]
        )

        # Combine the top and down vectors into a single vector
        # in the order top-down
        expression_group <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
            'expression_grouping'
        ]

        # Creating legend names for high/low expression groups:
        legend_name_high <- paste(gene_name,"high")

        legend_name_low <- paste(gene_name,"low")

        # Perform the survival analysis.
        # This uses the combined vector as the x variable
        # and the survival data as the y variable to create
        # a table with information about the test
        # including chi-squared p-value
        survival_table <- survdiff(
            high_low_survival_object  ~ expression_group
        )

        # Get the chi-squared test statistic from the analysis above:
        survival_table_chi_squared <- survival_table$chisq

        # Calculating a p value based on the test statistic to get
        # a precise p-value
        survival_pvalue <- as.numeric(
            1 - pchisq(
                abs(survival_table_chi_squared), df = 1
            )
        )

        ## Round the p-value on the graph to 3 digits:
        survival_pvalue_formatted <- formatC(
            survival_pvalue,
            format='e',
            digits=3
        )

        ## Create the graph if specified:
        if(graph==TRUE){

            # Create the plot title
            # with gene name and p-value included:
            survival_title <- paste(
                gene_name,
                "\nKaplan-Meier Survival analysis\nP = ",
                survival_pvalue_formatted,
                sep=''
            )

            ## Get the path to the step 5 folder to save
            ## the pdf in:
            path_to_folder <- getwd()

            ## Create a title for the survival plot pdf:
            survival_plot_pdf_title <- paste(
                path_to_folder,
                '/hyper.G-.output.survival/',
                gene_name,
                '_survival_plot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(survival_plot_pdf_title)

            ## Now actually create the survival plot
            ## Using a similar structure to what was used to generate the p-value:
            plot(
                survfit(
                    high_low_survival_object ~ expression_group
                ),

                # Color the lines (high in red first!)
                col = c('red', 'black'),

                ## Add thickness to the lines:
                lwd=3,

                # Use the title that was created earlier as the title
                # of the plot:
                main= survival_title,

                # Set titles of the x and y axis:
                # Note: TCGA measures survival in days as noted:
                xlab="Time (days)",
                ylab="Probability",

                # Change axis size:
                cex.axis=1,
                cex.main=1.5,
                cex.lab=1.25
            )

            ## Add a legend to the plot:
            legend(

                # Set X position of legend in graph:
                x= (
                    max(function_relevant_clinical_complete_high_low$days_to_last_followup) - 2500
                ),

                # Set Y position of legend in graph
                y= 1,

                # Use the legend titles that were created earlier
                legend= c(legend_name_high, legend_name_low),

                # As above, use black for low and red for high:
                col= c('red', 'black'),

                # Coloring the text in the legend as well
                text.col= c('red', 'black'),

                # Change the shape of the labels in the legend:
                pch= 15
            )

            ## Close the plot:
            dev.off()

        } else if(graph==FALSE){

            ## create vector of relevant info:
            relevant_vector <- c(
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(NA_normal),
                as.numeric(NA_tumor),
                as.numeric(mean_normal_expression),
                as.numeric(mean_tumor_expression),
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(present_tumor_sample_missing_n),
                as.numeric(present_tumor_sample_low_n),
                as.numeric(present_tumor_sample_intermediate_n),
                as.numeric(present_tumor_sample_high_n),
                as.numeric(present_tumor_sample_low_mean),
                as.numeric(present_tumor_sample_intermediate_mean),
                as.numeric(present_tumor_sample_high_mean),
                as.numeric(proportion_dead_low),
                as.numeric(proportion_dead_high),
                highest_dead_proportion_group,
                as.numeric(survival_pvalue)
            )
            names(relevant_vector) <- c(
                'normal_sample_count',
                'tumor_sample_count',
                'normal_sample_count_missing',
                'tumor_sample_count_missing',
                'mean_normal_expression',
                'mean_tumor_expression',
                'normal_sample_with_clinical_count',
                'tumor_sample_with_clinical_count',
                'tumor_sample_with_clinical_NA_count',
                'tumor_sample_with_clinical_low_count',
                'tumor_sample_with_clinical_intermediate_count',
                'tumor_sample_with_clinical_high_count',
                'mean_tumor_with_clinical_low_expression',
                'mean_tumor_with_clinical_intermediate_expression',
                'mean_tumor_with_clinical_high_expression',
                'proportion_dead_in_low_expression',
                'proportion_dead_in_high_expression',
                'survival_direction_of_effect',
                'survival_p_value'
            )

            ## Return the vector
            return(relevant_vector)
        }

    }

    ## Create a function to get the survival p-value or graph for
    ##beach probe of interest:
    methylation_survival_function_graph <- function(
        probe_of_interest,
        high_cutoff,
        low_cutoff,
        graph
    ){

        ## Get methylation values for probe of interest:
        methylation_values <- unlist(
            metData[
                probe_of_interest,
            ]
        )
        names(methylation_values) <- colnames(metData)

        # Split methylation values into normal and tumor samples
        # If both are present in the dataset
        tumor_methylation_values <- methylation_values[met_tumor_samples]
        normal_methylation_values <- methylation_values[met_normal_samples]

        # Changing the sample names to match the subject names:
        names(tumor_methylation_values)  <- substr(
            names(tumor_methylation_values),
            1,
            12
        )

        names(normal_methylation_values)  <- substr(
            names(normal_methylation_values),
            1,
            12
        )

        ## Calculate some basic data:
        normal_sample_n <- length(normal_methylation_values)
        tumor_sample_n <- length(tumor_methylation_values)

        ## Count the number of NA samples:
        NA_normal <- as.numeric(
            unname(
                table(
                    is.na(normal_methylation_values)
                )['TRUE']
            )
        )

        NA_tumor <- as.numeric(
            unname(
                table(
                    is.na(tumor_methylation_values)
                )['TRUE']
            )
        )

        if(is.na(NA_normal)){

            NA_normal <- 0

        }

        if(is.na(NA_tumor)){

            NA_tumor <- 0

        }

        ## Calculate mean methylation for all samples present:
        mean_tumor_methylation <- mean(
            tumor_methylation_values,
            na.rm= TRUE
        )

        mean_normal_methylation <- mean(
            normal_methylation_values,
            na.rm= TRUE
        )

        ## Calculate the number of samples with present clinical data:
        present_clinical_normal_sample_n <- as.numeric(
            unname(
                table(
                    names(normal_methylation_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        present_clinical_tumor_sample_n  <- as.numeric(
            unname(
                table(
                    names(tumor_methylation_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        ## Identify the samples which are present in the clinical data:
        present_normal_samples <- names(normal_methylation_values)[
            names(normal_methylation_values) %in% rownames(relevant_clinical)
        ]

        present_tumor_samples <- names(tumor_methylation_values)[
            names(tumor_methylation_values) %in% rownames(relevant_clinical)
        ]


        # Subsetting clinical patient data for individuals with
        # DNA methylation data in the tumor set:
        function_relevant_clinical <- relevant_clinical[
            present_tumor_samples,
        ]

        ## Calculate quantiles:
        high_cutoff_quantile= quantile(
            tumor_methylation_values,
            high_cutoff,
            na.rm= TRUE
        )[1]

        low_cutoff_quantile= quantile(
            tumor_methylation_values,
            low_cutoff,
            na.rm= TRUE
        )[1]

        ## Determine if a sample is in the high, low, or intermediate quartiles:
        sample_methylation_grouping <- ifelse(
            tumor_methylation_values >  high_cutoff_quantile,
            "High",
            ifelse(
                tumor_methylation_values <= low_cutoff_quantile,
                "Low",
                ifelse(
                    is.na(tumor_methylation_values),
                    'Missing',
                    'Intermediate'
                )
            )
        )

        ## Add the methylation and relevant clinical information to the relevant clinical dataframe:
        function_relevant_clinical$tumor_methylation <- tumor_methylation_values[
            rownames(function_relevant_clinical)
        ]

        ## Add the methylation grouping to the clinical data:
        function_relevant_clinical$methylation_grouping <- sample_methylation_grouping[
            rownames(function_relevant_clinical)
        ]

        ## Remove samples missing methylation data:
        function_relevant_clinical_complete <- na.omit(
            function_relevant_clinical,
            cols=c(
                "tumor_methylation",
                "methylation_grouping"
            )
        )
        function_relevant_clinical_complete <- function_relevant_clinical

        ## Count the number of samples of each group remaining:
        present_tumor_sample_high_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='High',
            ]
        )

        present_tumor_sample_intermediate_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Intermediate',
            ]
        )

        present_tumor_sample_low_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Low',
            ]
        )

        present_tumor_sample_missing_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Missing',
            ]
        )

        ## Get the mean methylation in each group:
        present_tumor_sample_high_mean <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='High',
                'tumor_methylation'
            ]
        )

        present_tumor_sample_intermediate_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Intermediate',
                'tumor_methylation'
            ]
        )

        present_tumor_sample_low_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Low',
                'tumor_methylation'
            ]
        )

        ## Get clinical results for only samples in the high and low group with info:
        function_relevant_clinical_complete_high_low <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
        ]

        ## Calculate proportion of deceased individuals in each group
        proportion_dead_high <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$methylation_grouping=='High' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$methylation_grouping=='High'),
                    ]
                )
        )

        proportion_dead_low <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$methylation_grouping=='Low' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$methylation_grouping=='Low'),
                    ]
                )
        )

        ## Determine which group had greater proportion of fatalities:
        highest_dead_proportion_group <- ifelse(
            proportion_dead_high > proportion_dead_low,
            'high_methylation_low_survival',
            ifelse(
                proportion_dead_high < proportion_dead_low,
                'low_methylation_low_survival',
                'unclear'
            )
        )

        ## Create a high_low survival object:
        high_low_survival_object <- Surv(
            function_relevant_clinical_complete_high_low$time,
            ifelse(
                function_relevant_clinical_complete_high_low$vital_status=='Alive',
                FALSE,
                TRUE
            )
        )

        # Set the rownames of the survival object to be the
        # patient names for the high/low group smaples:
        rownames(high_low_survival_object) <- rownames(
            function_relevant_clinical_complete[
                !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
            ]
        )

        # Combine the top and down vectors into a single vector
        # in the order top-down
        methylation_group <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
            'methylation_grouping'
        ]

        # Creating legend names for high/low methylation groups:
        legend_name_high <- paste(probe_of_interest,"high")

        legend_name_low <- paste(probe_of_interest,"low")

        # Perform the survival analysis.
        # This uses the combined vector as the x variable
        # and the survival data as the y variable to create
        # a table with information about the test
        # including chi-squared p-value
        survival_table <- survdiff(
            high_low_survival_object  ~ methylation_group
        )

        # Get the chi-squared test statistic from the analysis above:
        survival_table_chi_squared <- survival_table$chisq

        # Calculating a p value based on the test statistic to get
        # a precise p-value
        survival_pvalue <- as.numeric(
            1 - pchisq(
                abs(survival_table_chi_squared), df = 1
            )
        )

        ## Round the p-value on the graph to 3 digits:
        survival_pvalue_formatted <- formatC(
            survival_pvalue,
            format='e',
            digits=3
        )

        ## Create the graph if specified:
        if(graph==TRUE){

            # Create the plot title
            # with gene name and p-value included:
            survival_title <- paste(
                probe_of_interest,
                "\nKaplan-Meier Survival analysis\nP = ",
                survival_pvalue_formatted,
                sep=''
            )

            ## Get the path to the step 5 folder to save
            ## the pdf in:
            path_to_folder <- getwd()

            ## Create a title for the survival plot pdf:
            ## This is a comination of the probe name with the linked gene:
            survival_plot_pdf_title <- paste(
                path_to_folder,
                '/hyper.G-.output.survival/',
                probe_of_interest,
                '_survival_plot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(survival_plot_pdf_title)

            ## Now actually create the survival plot
            ## Using a similar structure to what was used to generate the p-value:
            plot(
                survfit(
                    high_low_survival_object ~ methylation_group
                ),

                # Color the lines (high in red first!)
                col = c('red', 'black'),

                ## Add thickness to the lines:
                lwd=3,

                # Use the title that was created earlier as the title
                # of the plot:
                main= survival_title,

                # Set titles of the x and y axis:
                # Note: TCGA measures survival in days as noted:
                xlab="Time (days)",
                ylab="Probability",

                # Change axis size:
                cex.axis=1,
                cex.main=1.5,
                cex.lab=1.25
            )

            ## Add a legend to the plot:
            legend(

                # Set X position of legend in graph:
                x= (
                    max(function_relevant_clinical_complete_high_low$days_to_last_followup) - 2500
                ),

                # Set Y position of legend in graph
                y= 1,

                # Use the legend titles that were created earlier
                legend= c(legend_name_high, legend_name_low),

                # As above, use black for low and red for high:
                col= c('red', 'black'),

                # Coloring the text in the legend as well
                text.col= c('red', 'black'),

                # Change the shape of the labels in the legend:
                pch= 15
            )

            ## Close the plot:
            dev.off()

        } else if(graph==FALSE){

            ## create vector of relevant info:
            relevant_vector <- c(
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(NA_normal),
                as.numeric(NA_tumor),
                as.numeric(mean_normal_methylation),
                as.numeric(mean_tumor_methylation),
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(present_tumor_sample_missing_n),
                as.numeric(present_tumor_sample_low_n),
                as.numeric(present_tumor_sample_intermediate_n),
                as.numeric(present_tumor_sample_high_n),
                as.numeric(present_tumor_sample_low_mean),
                as.numeric(present_tumor_sample_intermediate_mean),
                as.numeric(present_tumor_sample_high_mean),
                as.numeric(proportion_dead_low),
                as.numeric(proportion_dead_high),
                highest_dead_proportion_group,
                as.numeric(survival_pvalue)
            )
            names(relevant_vector) <- c(
                'normal_sample_count',
                'tumor_sample_count',
                'normal_sample_count_missing',
                'tumor_sample_count_missing',
                'mean_normal_methylation',
                'mean_tumor_methylation',
                'normal_sample_with_clinical_count',
                'tumor_sample_with_clinical_count',
                'tumor_sample_with_clinical_NA_count',
                'tumor_sample_with_clinical_low_count',
                'tumor_sample_with_clinical_intermediate_count',
                'tumor_sample_with_clinical_high_count',
                'mean_tumor_with_clinical_low_methylation',
                'mean_tumor_with_clinical_intermediate_methylation',
                'mean_tumor_with_clinical_high_methylation',
                'proportion_dead_in_low_methylation',
                'proportion_dead_in_high_methylation',
                'survival_direction_of_effect',
                'survival_p_value'
            )

            ## Return the vector
            return(relevant_vector)
        }

    }

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_genes==TRUE){

        ## Create survival plots for each of the genes designated:
        mclapply(
            X= top_gene_names,
            FUN= expression_survival_function_graph,
            high_cutoff= high_thresh,
            low_cutoff= low_thresh,
            graph= TRUE,
            mc.cores= cores
        )

    }

    ## Get survival information for each gene:
    gene_survival_results_list <- mclapply(
        X= top_gene_names,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= FALSE,
        mc.cores= cores
    )

    ## Transform the list into a matrix:
    gene_survival_results_df <- data.frame(
        matrix(
            unlist(gene_survival_results_list),
            nrow=length(gene_survival_results_list),
            byrow=T
        ),
        stringsAsFactors = FALSE
    )

    colnames(gene_survival_results_df) <- c(
        'normal_sample_count',
        'tumor_sample_count',
        'normal_sample_count_missing',
        'tumor_sample_count_missing',
        'mean_normal_expression',
        'mean_tumor_expression',
        'normal_sample_with_clinical_count',
        'tumor_sample_with_clinical_count',
        'tumor_sample_with_clinical_NA_count',
        'tumor_sample_with_clinical_low_count',
        'tumor_sample_with_clinical_intermediate_count',
        'tumor_sample_with_clinical_high_count',
        'mean_tumor_with_clinical_low_expression',
        'mean_tumor_with_clinical_intermediate_expression',
        'mean_tumor_with_clinical_high_expression',
        'proportion_dead_in_low_expression',
        'proportion_dead_in_high_expression',
        'survival_direction_of_effect',
        'survival_p_value'
    )

    rownames(gene_survival_results_df) <- c(
        gene_name <- gencode_v22_genes_df[
            top_gene_names, 'gene_name'
        ]
    )

    ## Save the gene information as a .tsv:
    write.table(
        gene_survival_results_df,
        file='./hyper.G-.output.survival/hyper.G-.top_genes_survival_info.tsv',
        quote= FALSE,
        sep='\t'
    )

    ## Now for each of the probes, we will need to calculate survival info and graphs:

    ## Get the list of probes linked to the significant genes:
    CpGs_linked <- CpG_linkage_dataset[
        CpG_linkage_dataset$geneID %in% top_gene_names,
        'probe'
    ]

    ## Get the unique probes from the list:
    unique_CpGs_linked <- unique(CpGs_linked)

    ## Get the survival info for each probe by mclapplying the
    ## methylation survival function to all of the probes
    ## linked to each gene:
    probe_survival_results_list <- mclapply(
        X= unique_CpGs_linked,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= FALSE,
        mc.cores= cores
    )

    ## Transform the list into a matrix:
    probe_survival_results_df <- data.frame(
        matrix(
            unlist(probe_survival_results_list),
            nrow=length(probe_survival_results_list),
            byrow=T
        ),
        stringsAsFactors = FALSE
    )

    colnames(probe_survival_results_df) <- c(
        'normal_sample_count',
        'tumor_sample_count',
        'normal_sample_count_missing',
        'tumor_sample_count_missing',
        'mean_normal_methylation',
        'mean_tumor_methylation',
        'normal_sample_with_clinical_count',
        'tumor_sample_with_clinical_count',
        'tumor_sample_with_clinical_NA_count',
        'tumor_sample_with_clinical_low_count',
        'tumor_sample_with_clinical_intermediate_count',
        'tumor_sample_with_clinical_high_count',
        'mean_tumor_with_clinical_low_methylation',
        'mean_tumor_with_clinical_intermediate_methylation',
        'mean_tumor_with_clinical_high_methylation',
        'proportion_dead_in_low_methylation',
        'proportion_dead_in_high_methylation',
        'survival_direction_of_effect',
        'survival_p_value'
    )

    rownames(probe_survival_results_df) <- unique_CpGs_linked

    ## For each of the CpGs, list which of the top X genes it was listed to:

    ## Wite a function to get the top X gene names a probe is assigned to:
    top_gene_assignment_function <- function(probe_of_interest){

        all_listed_genes <- CpG_linkage_dataset[
            CpG_linkage_dataset$probe==probe_of_interest,
            'geneSymbol'
        ]

        ## Conver thte top X genes to their gene symbols
        top_gene_symbols <- gencode_v22_genes_df[
            top_gene_names, 'gene_name'
        ]

        ## Get the names of the top genes that were linked to this probe:
        top_gene_symbols_linked_to_probe <- top_gene_symbols[
            top_gene_symbols %in% all_listed_genes
        ]

        ## Return the listed genes:
        return(
            paste(
                top_gene_symbols_linked_to_probe,
                collapse=','
            )
        )
    }

    ## Add the listed top genes to the survival results dataframe
    ## for the probes:
    probe_survival_results_df$top_genes_linked <- unname(
        sapply(
            rownames(probe_survival_results_df),
            top_gene_assignment_function
        )
    )

    ## Create a vector of the CpG p-values and add the CpG entries
    ## as names
    CpGs_linked_p_values_vector <- probe_survival_results_df$survival_p_value
    names(CpGs_linked_p_values_vector) <- rownames(probe_survival_results_df)

    ## Get the CpGs that are nominally significant:
    CpGs_linked_p_values_nominally_significant <- CpGs_linked_p_values_vector[
        CpGs_linked_p_values_vector<0.05
    ]

    ## Order the nominally significant CpGs:
    CpGs_linked_p_values_nominally_significant_ordered <- sort(
        CpGs_linked_p_values_nominally_significant,
        decreasing = FALSE
    )

    ## Get the names of the nominally significant CpGs:
    CpGs_linked_nominally_significant_names <- names(
        CpGs_linked_p_values_nominally_significant_ordered
    )

    ## If visualize_survival_plots_probes is TRUE, generate the plots:
    if(visualize_survival_plots_probes==TRUE){

        #### Save plots for the nominally significant CpGs:
        mclapply(
            X= CpGs_linked_nominally_significant_names,
            FUN= methylation_survival_function_graph,
            high_cutoff= high_thresh,
            low_cutoff= low_thresh,
            graph= TRUE,
            mc.cores= cores
        )

    }

    ## Save the probe information as a .tsv:
    write.table(
        probe_survival_results_df,
        file='./hyper.G-.output.survival/hyper.G-.probes_linked_to_top_genes_survival_info.tsv',
        quote= FALSE,
        sep='\t'
    )
}

make.survival.plots.hypometh.G_pos <- function(prefix, survival_top_n_genes, clinical_data, expDataN, expDataT, metDataN, metDataT, visualize_survival_plots_genes, high_thresh, low_thresh, cores, visualize_survival_plots_probes) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Load ggplot2:
    # library('ggplot2')

    ## Load survival:
    # library('survival')

    ## Load parallel:
    # library('parallel')

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hypo.G+.output.histogram/",
            prefix,
            ".hypo.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hypo.G+.output/hypo.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Set the number of genes to analyze based on user input:
    survival_top_n_genes <- survival_top_n_genes

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_names <- ordered_TFs_by_link_count[
        1:survival_top_n_genes,
        'geneID'
    ]

    ## Load the external data if it is available:

    ## Check if .rda file with expression, methylation, and
    ## perhaps clinical data was placed in TENET folder:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## If the combined rda file exists, load it and see if clinical data is included:
    ## If not, check to see if clinical data is included in the "clinical" folder
    load(combined_rda_file)

    # #### REMOVE THIS LATER:
    # load(
    #     "/Volumes/NO_NAME/TENET_paper_update/tenet_match.matched.dedup.luad.methylation.expression.dataset.rda",
    # )

    ## If data is loaded as clinical_data, rename it to clinical:
    if(exists('clinical_data')){

        clinical <- clinical_data

    }

    ## Check to see if the clinical data is properly loaded
    ## if it isn't, check the clinical folder in external.data/data
    if(exists('clinical')){

        paste('clinical data loaded')

    } else{

        ## Get a path to the rda file containing clinical info
        ## in external.data/data/clinical
        clinical_rda_file <- list.files(
            path='../external.data/data/clinical',
            pattern='\\.rda$',
            full.names= TRUE,
            include.dirs = TRUE
        )

        if(length(clinical_rda_file)==1){

            ## Load the .rda file:
            load(clinical_rda_file)

            ## If data is loaded as clinical_data, rename it to clinical:
            if(exists('clinical_data')){

                clinical <- clinical_data

            }

            ## Check to see if the clinical data is now properly loaded
            if(exists('clinical')){

                paste('clinical data loaded')

            } else{

                ## CLinical data is not found in rda
                ## return an error message and abort script:
                paste('rda file did not contain clinical data or it was not saved as clinical or clinical_data')
                quit(save='no')
            }

        } else{

            ## Clinical rda not found.
            ## Return error message and abort script
            ## Later functionality may be added to try loading clinical data from txt file
            ## See step 0 make.rda.from.textfiles.R for reference:
            paste('rda file with clinical data not found')
            quit(save='no')
        }
    }

    ## Process the expression/methylation/clinical data so the function doesn't have to:

    ## combine datasets:
    expData <- cbind(expDataN, expDataT)
    metData <- cbind(metDataN, metDataT)

    ## Get the list of normal and tumor samples:
    exp_normal_samples <- colnames(expDataN)
    exp_tumor_samples <- colnames(expDataT)

    met_normal_samples <- colnames(metDataN)
    met_tumor_samples <- colnames(metDataT)

    ## Get the relevant columns from the clinical data:
    relevant_clinical <- clinical[
        ,c(
            "bcr_patient_barcode",
            "days_to_death",
            "days_to_last_followup",
            "vital_status"
        )
    ]
    rownames(relevant_clinical) <- relevant_clinical$bcr_patient_barcode

    # If there are subjects that are alive in the
    # dataset, set their days to death info to '-Inf'
    relevant_clinical[
        grep(
            "alive",
            relevant_clinical$vital_status,
            ignore.case = TRUE
        ),
        'days_to_death'
    ] <- '-Inf'

    # If there are subjects that are dead in the
    # dataset, set their days to last followup info to '-Inf'
    relevant_clinical[
        grep(
            "dead",
            relevant_clinical$vital_status,
            ignore.case = TRUE
        ),
        'days_to_last_followup'
    ] <- '-Inf'

    ## Remove subjects that remain unaccounted for:
    ## i.e. NA values remain in previous two columns:
    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_death"]
        ),
    ]

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_last_followup"]
        ),
    ]

    # Changing the days to last followup info in the cancer clinical data
    # to be numeric values through character values and remove any NAs that are induced
    # as a final check:
    relevant_clinical$days_to_death <- as.numeric(
        as.character(relevant_clinical$days_to_death)
    )

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_last_followup"]
        ),
    ]

    # Changing the days to death info in the cancer clinical data
    # to be numeric values through character values and remove any NAs that are induced
    # as a final check:
    relevant_clinical$days_to_last_followup <- as.numeric(
        as.character(relevant_clinical$days_to_last_followup)
    )

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_death"]
        ),
    ]

    # Add the relevant number as a final "time" variable
    # Equal to days to death for dead individuals, and
    # days to last followup for alive individuals
    relevant_clinical$time <- ifelse(
        relevant_clinical$vital_status=='Alive',
        relevant_clinical$days_to_last_followup,
        ifelse(
            relevant_clinical$vital_status=='Dead',
            relevant_clinical$days_to_death,
            NA
        )
    )

    ## Create a function to get the survival p-value or graph
    ## for each gene of interest:
    expression_survival_function_graph <- function(
        gene_of_interest,
        high_cutoff,
        low_cutoff,
        graph
    ){

        ## Get actual gene name depending on input:
        if(substring(gene_of_interest,1,4)=='ENSG' & nchar(gene_of_interest)==15){

            ## Input is in ENSG, get the gene name:
            gene_ENSG <- gene_of_interest

            gene_name <- gencode_v22_genes_df[
                gene_ENSG, 'gene_name'
            ]

            ## Get gene ENSG assuming a name is plugged in:
        } else{

            ## Assume what was given was the gene name, get the ENSG:
            gene_name <- gene_of_interest

            gene_ENSG <- rownames(
                gencode_v22_genes_df[gencode_v22_genes_df$gene_name==gene_name,]
            )
        }

        ## Get expression values for gene of interest:
        expression_values <- unlist(
            expData[
                gene_ENSG,
            ]
        )
        names(expression_values) <- colnames(expData)

        # Split gene expression values into normal and tumor samples
        tumor_expression_values <- expression_values[exp_tumor_samples]
        normal_expression_values <- expression_values[exp_normal_samples]

        # Changing the sample names to match the subject names:
        names(tumor_expression_values)  <- substr(
            names(tumor_expression_values),
            1,
            12
        )

        names(normal_expression_values)  <- substr(
            names(normal_expression_values),
            1,
            12
        )

        ## Calculate some basic data:
        normal_sample_n <- length(normal_expression_values)
        tumor_sample_n <- length(tumor_expression_values)

        ## Count the number of NA samples:
        NA_normal <- as.numeric(
            unname(
                table(
                    is.na(normal_expression_values)
                )['TRUE']
            )
        )

        NA_tumor <- as.numeric(
            unname(
                table(
                    is.na(tumor_expression_values)
                )['TRUE']
            )
        )

        if(is.na(NA_normal)){

            NA_normal <- 0

        }

        if(is.na(NA_tumor)){

            NA_tumor <- 0

        }

        ## Calculate mean expression for all samples present:
        mean_tumor_expression <- mean(
            tumor_expression_values,
            na.rm= TRUE
        )

        mean_normal_expression <- mean(
            normal_expression_values,
            na.rm= TRUE
        )

        ## Calculate the number of samples with present clinical data:
        present_clinical_normal_sample_n <- as.numeric(
            unname(
                table(
                    names(normal_expression_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        present_clinical_tumor_sample_n  <- as.numeric(
            unname(
                table(
                    names(tumor_expression_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        ## Identify the samples which are present in the clinical data:
        present_normal_samples <- names(normal_expression_values)[
            names(normal_expression_values) %in% rownames(relevant_clinical)
        ]

        present_tumor_samples <- names(tumor_expression_values)[
            names(tumor_expression_values) %in% rownames(relevant_clinical)
        ]


        # Subsetting clinical patient data for individuals with
        # gene expression data in the tumor set:
        function_relevant_clinical <- relevant_clinical[
            present_tumor_samples,
        ]

        ## Calculate quantiles:
        high_cutoff_quantile= quantile(
            tumor_expression_values,
            high_cutoff,
            na.rm= TRUE
        )[1]

        low_cutoff_quantile= quantile(
            tumor_expression_values,
            low_cutoff,
            na.rm= TRUE
        )[1]

        ## Determine if a sample is in the high, low, or intermediate quartiles:
        sample_expression_grouping <- ifelse(
            tumor_expression_values >  high_cutoff_quantile,
            "High",
            ifelse(
                tumor_expression_values <= low_cutoff_quantile,
                "Low",
                ifelse(
                    is.na(tumor_expression_values),
                    'Missing',
                    'Intermediate'
                )
            )
        )

        ## Add the expression and relevant clinical information to the relevant clinical dataframe:
        function_relevant_clinical$tumor_expression <- tumor_expression_values[
            rownames(function_relevant_clinical)
        ]

        ## Add the expression grouping to the clinical data:
        function_relevant_clinical$expression_grouping <- sample_expression_grouping[
            rownames(function_relevant_clinical)
        ]

        ## Remove samples missing expression data:
        function_relevant_clinical_complete <- na.omit(
            function_relevant_clinical,
            cols=c(
                "tumor_expression",
                "expression_grouping"
            )
        )

        function_relevant_clinical_complete <- function_relevant_clinical

        ## Count the number of samples of each group remaining:
        present_tumor_sample_high_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='High',
            ]
        )

        present_tumor_sample_intermediate_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Intermediate',
            ]
        )

        present_tumor_sample_low_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Low',
            ]
        )

        present_tumor_sample_missing_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Missing',
            ]
        )

        ## Get the mean expression in each group:
        present_tumor_sample_high_mean <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='High',
                'tumor_expression'
            ]
        )

        present_tumor_sample_intermediate_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Intermediate',
                'tumor_expression'
            ]
        )

        present_tumor_sample_low_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Low',
                'tumor_expression'
            ]
        )

        ## Get clinical results for only samples in the high and low group with info:
        function_relevant_clinical_complete_high_low <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
        ]

        ## Calculate proportion of deceased individuals in each group
        proportion_dead_high <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$expression_grouping=='High' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$expression_grouping=='High'),
                    ]
                )
        )

        proportion_dead_low <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$expression_grouping=='Low' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$expression_grouping=='Low'),
                    ]
                )
        )

        ## Determine which group had greater proportion of fatalities:
        highest_dead_proportion_group <- ifelse(
            proportion_dead_high > proportion_dead_low,
            'high_expression_low_survival',
            ifelse(
                proportion_dead_high < proportion_dead_low,
                'low_expression_low_survival',
                'unclear'
            )
        )

        ## Create a high_low survival object:
        high_low_survival_object <- Surv(
            function_relevant_clinical_complete_high_low$time,
            ifelse(
                function_relevant_clinical_complete_high_low$vital_status=='Alive',
                FALSE,
                TRUE
            )
        )

        # Set the rownames of the survival object to be the
        # patient names for the high/low group smaples:
        rownames(high_low_survival_object) <- rownames(
            function_relevant_clinical_complete[
                !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
            ]
        )

        # Combine the top and down vectors into a single vector
        # in the order top-down
        expression_group <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
            'expression_grouping'
        ]

        # Creating legend names for high/low expression groups:
        legend_name_high <- paste(gene_name,"high")

        legend_name_low <- paste(gene_name,"low")

        # Perform the survival analysis.
        # This uses the combined vector as the x variable
        # and the survival data as the y variable to create
        # a table with information about the test
        # including chi-squared p-value
        survival_table <- survdiff(
            high_low_survival_object  ~ expression_group
        )

        # Get the chi-squared test statistic from the analysis above:
        survival_table_chi_squared <- survival_table$chisq

        # Calculating a p value based on the test statistic to get
        # a precise p-value
        survival_pvalue <- as.numeric(
            1 - pchisq(
                abs(survival_table_chi_squared), df = 1
            )
        )

        ## Round the p-value on the graph to 3 digits:
        survival_pvalue_formatted <- formatC(
            survival_pvalue,
            format='e',
            digits=3
        )

        ## Create the graph if specified:
        if(graph==TRUE){

            # Create the plot title
            # with gene name and p-value included:
            survival_title <- paste(
                gene_name,
                "\nKaplan-Meier Survival analysis\nP = ",
                survival_pvalue_formatted,
                sep=''
            )

            ## Get the path to the step 5 folder to save
            ## the pdf in:
            path_to_folder <- getwd()

            ## Create a title for the survival plot pdf:
            survival_plot_pdf_title <- paste(
                path_to_folder,
                '/hypo.G+.output.survival/',
                gene_name,
                '_survival_plot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(survival_plot_pdf_title)

            ## Now actually create the survival plot
            ## Using a similar structure to what was used to generate the p-value:
            plot(
                survfit(
                    high_low_survival_object ~ expression_group
                ),

                # Color the lines (high in red first!)
                col = c('red', 'black'),

                ## Add thickness to the lines:
                lwd=3,

                # Use the title that was created earlier as the title
                # of the plot:
                main= survival_title,

                # Set titles of the x and y axis:
                # Note: TCGA measures survival in days as noted:
                xlab="Time (days)",
                ylab="Probability",

                # Change axis size:
                cex.axis=1,
                cex.main=1.5,
                cex.lab=1.25
            )

            ## Add a legend to the plot:
            legend(

                # Set X position of legend in graph:
                x= (
                    max(function_relevant_clinical_complete_high_low$days_to_last_followup) - 2500
                ),

                # Set Y position of legend in graph
                y= 1,

                # Use the legend titles that were created earlier
                legend= c(legend_name_high, legend_name_low),

                # As above, use black for low and red for high:
                col= c('red', 'black'),

                # Coloring the text in the legend as well
                text.col= c('red', 'black'),

                # Change the shape of the labels in the legend:
                pch= 15
            )

            ## Close the plot:
            dev.off()

        } else if(graph==FALSE){

            ## create vector of relevant info:
            relevant_vector <- c(
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(NA_normal),
                as.numeric(NA_tumor),
                as.numeric(mean_normal_expression),
                as.numeric(mean_tumor_expression),
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(present_tumor_sample_missing_n),
                as.numeric(present_tumor_sample_low_n),
                as.numeric(present_tumor_sample_intermediate_n),
                as.numeric(present_tumor_sample_high_n),
                as.numeric(present_tumor_sample_low_mean),
                as.numeric(present_tumor_sample_intermediate_mean),
                as.numeric(present_tumor_sample_high_mean),
                as.numeric(proportion_dead_low),
                as.numeric(proportion_dead_high),
                highest_dead_proportion_group,
                as.numeric(survival_pvalue)
            )
            names(relevant_vector) <- c(
                'normal_sample_count',
                'tumor_sample_count',
                'normal_sample_count_missing',
                'tumor_sample_count_missing',
                'mean_normal_expression',
                'mean_tumor_expression',
                'normal_sample_with_clinical_count',
                'tumor_sample_with_clinical_count',
                'tumor_sample_with_clinical_NA_count',
                'tumor_sample_with_clinical_low_count',
                'tumor_sample_with_clinical_intermediate_count',
                'tumor_sample_with_clinical_high_count',
                'mean_tumor_with_clinical_low_expression',
                'mean_tumor_with_clinical_intermediate_expression',
                'mean_tumor_with_clinical_high_expression',
                'proportion_dead_in_low_expression',
                'proportion_dead_in_high_expression',
                'survival_direction_of_effect',
                'survival_p_value'
            )

            ## Return the vector
            return(relevant_vector)
        }

    }

    ## Create a function to get the survival p-value or graph for
    ##beach probe of interest:
    methylation_survival_function_graph <- function(
        probe_of_interest,
        high_cutoff,
        low_cutoff,
        graph
    ){

        ## Get methylation values for probe of interest:
        methylation_values <- unlist(
            metData[
                probe_of_interest,
            ]
        )
        names(methylation_values) <- colnames(metData)

        # Split methylation values into normal and tumor samples
        # If both are present in the dataset
        tumor_methylation_values <- methylation_values[met_tumor_samples]
        normal_methylation_values <- methylation_values[met_normal_samples]

        # Changing the sample names to match the subject names:
        names(tumor_methylation_values)  <- substr(
            names(tumor_methylation_values),
            1,
            12
        )

        names(normal_methylation_values)  <- substr(
            names(normal_methylation_values),
            1,
            12
        )

        ## Calculate some basic data:
        normal_sample_n <- length(normal_methylation_values)
        tumor_sample_n <- length(tumor_methylation_values)

        ## Count the number of NA samples:
        NA_normal <- as.numeric(
            unname(
                table(
                    is.na(normal_methylation_values)
                )['TRUE']
            )
        )

        NA_tumor <- as.numeric(
            unname(
                table(
                    is.na(tumor_methylation_values)
                )['TRUE']
            )
        )

        if(is.na(NA_normal)){

            NA_normal <- 0

        }

        if(is.na(NA_tumor)){

            NA_tumor <- 0

        }

        ## Calculate mean methylation for all samples present:
        mean_tumor_methylation <- mean(
            tumor_methylation_values,
            na.rm= TRUE
        )

        mean_normal_methylation <- mean(
            normal_methylation_values,
            na.rm= TRUE
        )

        ## Calculate the number of samples with present clinical data:
        present_clinical_normal_sample_n <- as.numeric(
            unname(
                table(
                    names(normal_methylation_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        present_clinical_tumor_sample_n  <- as.numeric(
            unname(
                table(
                    names(tumor_methylation_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        ## Identify the samples which are present in the clinical data:
        present_normal_samples <- names(normal_methylation_values)[
            names(normal_methylation_values) %in% rownames(relevant_clinical)
        ]

        present_tumor_samples <- names(tumor_methylation_values)[
            names(tumor_methylation_values) %in% rownames(relevant_clinical)
        ]


        # Subsetting clinical patient data for individuals with
        # DNA methylation data in the tumor set:
        function_relevant_clinical <- relevant_clinical[
            present_tumor_samples,
        ]

        ## Calculate quantiles:
        high_cutoff_quantile= quantile(
            tumor_methylation_values,
            high_cutoff,
            na.rm= TRUE
        )[1]

        low_cutoff_quantile= quantile(
            tumor_methylation_values,
            low_cutoff,
            na.rm= TRUE
        )[1]

        ## Determine if a sample is in the high, low, or intermediate quartiles:
        sample_methylation_grouping <- ifelse(
            tumor_methylation_values >  high_cutoff_quantile,
            "High",
            ifelse(
                tumor_methylation_values <= low_cutoff_quantile,
                "Low",
                ifelse(
                    is.na(tumor_methylation_values),
                    'Missing',
                    'Intermediate'
                )
            )
        )

        ## Add the methylation and relevant clinical information to the relevant clinical dataframe:
        function_relevant_clinical$tumor_methylation <- tumor_methylation_values[
            rownames(function_relevant_clinical)
        ]

        ## Add the methylation grouping to the clinical data:
        function_relevant_clinical$methylation_grouping <- sample_methylation_grouping[
            rownames(function_relevant_clinical)
        ]

        ## Remove samples missing methylation data:
        function_relevant_clinical_complete <- na.omit(
            function_relevant_clinical,
            cols=c(
                "tumor_methylation",
                "methylation_grouping"
            )
        )
        function_relevant_clinical_complete <- function_relevant_clinical

        ## Count the number of samples of each group remaining:
        present_tumor_sample_high_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='High',
            ]
        )

        present_tumor_sample_intermediate_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Intermediate',
            ]
        )

        present_tumor_sample_low_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Low',
            ]
        )

        present_tumor_sample_missing_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Missing',
            ]
        )

        ## Get the mean methylation in each group:
        present_tumor_sample_high_mean <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='High',
                'tumor_methylation'
            ]
        )

        present_tumor_sample_intermediate_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Intermediate',
                'tumor_methylation'
            ]
        )

        present_tumor_sample_low_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Low',
                'tumor_methylation'
            ]
        )

        ## Get clinical results for only samples in the high and low group with info:
        function_relevant_clinical_complete_high_low <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
        ]

        ## Calculate proportion of deceased individuals in each group
        proportion_dead_high <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$methylation_grouping=='High' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$methylation_grouping=='High'),
                    ]
                )
        )

        proportion_dead_low <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$methylation_grouping=='Low' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$methylation_grouping=='Low'),
                    ]
                )
        )

        ## Determine which group had greater proportion of fatalities:
        highest_dead_proportion_group <- ifelse(
            proportion_dead_high > proportion_dead_low,
            'high_methylation_low_survival',
            ifelse(
                proportion_dead_high < proportion_dead_low,
                'low_methylation_low_survival',
                'unclear'
            )
        )

        ## Create a high_low survival object:
        high_low_survival_object <- Surv(
            function_relevant_clinical_complete_high_low$time,
            ifelse(
                function_relevant_clinical_complete_high_low$vital_status=='Alive',
                FALSE,
                TRUE
            )
        )

        # Set the rownames of the survival object to be the
        # patient names for the high/low group smaples:
        rownames(high_low_survival_object) <- rownames(
            function_relevant_clinical_complete[
                !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
            ]
        )

        # Combine the top and down vectors into a single vector
        # in the order top-down
        methylation_group <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
            'methylation_grouping'
        ]

        # Creating legend names for high/low methylation groups:
        legend_name_high <- paste(probe_of_interest,"high")

        legend_name_low <- paste(probe_of_interest,"low")

        # Perform the survival analysis.
        # This uses the combined vector as the x variable
        # and the survival data as the y variable to create
        # a table with information about the test
        # including chi-squared p-value
        survival_table <- survdiff(
            high_low_survival_object  ~ methylation_group
        )

        # Get the chi-squared test statistic from the analysis above:
        survival_table_chi_squared <- survival_table$chisq

        # Calculating a p value based on the test statistic to get
        # a precise p-value
        survival_pvalue <- as.numeric(
            1 - pchisq(
                abs(survival_table_chi_squared), df = 1
            )
        )

        ## Round the p-value on the graph to 3 digits:
        survival_pvalue_formatted <- formatC(
            survival_pvalue,
            format='e',
            digits=3
        )

        ## Create the graph if specified:
        if(graph==TRUE){

            # Create the plot title
            # with gene name and p-value included:
            survival_title <- paste(
                probe_of_interest,
                "\nKaplan-Meier Survival analysis\nP = ",
                survival_pvalue_formatted,
                sep=''
            )

            ## Get the path to the step 5 folder to save
            ## the pdf in:
            path_to_folder <- getwd()

            ## Create a title for the survival plot pdf:
            ## This is a comination of the probe name with the linked gene:
            survival_plot_pdf_title <- paste(
                path_to_folder,
                '/hypo.G+.output.survival/',
                probe_of_interest,
                '_survival_plot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(survival_plot_pdf_title)

            ## Now actually create the survival plot
            ## Using a similar structure to what was used to generate the p-value:
            plot(
                survfit(
                    high_low_survival_object ~ methylation_group
                ),

                # Color the lines (high in red first!)
                col = c('red', 'black'),

                ## Add thickness to the lines:
                lwd=3,

                # Use the title that was created earlier as the title
                # of the plot:
                main= survival_title,

                # Set titles of the x and y axis:
                # Note: TCGA measures survival in days as noted:
                xlab="Time (days)",
                ylab="Probability",

                # Change axis size:
                cex.axis=1,
                cex.main=1.5,
                cex.lab=1.25
            )

            ## Add a legend to the plot:
            legend(

                # Set X position of legend in graph:
                x= (
                    max(function_relevant_clinical_complete_high_low$days_to_last_followup) - 2500
                ),

                # Set Y position of legend in graph
                y= 1,

                # Use the legend titles that were created earlier
                legend= c(legend_name_high, legend_name_low),

                # As above, use black for low and red for high:
                col= c('red', 'black'),

                # Coloring the text in the legend as well
                text.col= c('red', 'black'),

                # Change the shape of the labels in the legend:
                pch= 15
            )

            ## Close the plot:
            dev.off()

        } else if(graph==FALSE){

            ## create vector of relevant info:
            relevant_vector <- c(
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(NA_normal),
                as.numeric(NA_tumor),
                as.numeric(mean_normal_methylation),
                as.numeric(mean_tumor_methylation),
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(present_tumor_sample_missing_n),
                as.numeric(present_tumor_sample_low_n),
                as.numeric(present_tumor_sample_intermediate_n),
                as.numeric(present_tumor_sample_high_n),
                as.numeric(present_tumor_sample_low_mean),
                as.numeric(present_tumor_sample_intermediate_mean),
                as.numeric(present_tumor_sample_high_mean),
                as.numeric(proportion_dead_low),
                as.numeric(proportion_dead_high),
                highest_dead_proportion_group,
                as.numeric(survival_pvalue)
            )
            names(relevant_vector) <- c(
                'normal_sample_count',
                'tumor_sample_count',
                'normal_sample_count_missing',
                'tumor_sample_count_missing',
                'mean_normal_methylation',
                'mean_tumor_methylation',
                'normal_sample_with_clinical_count',
                'tumor_sample_with_clinical_count',
                'tumor_sample_with_clinical_NA_count',
                'tumor_sample_with_clinical_low_count',
                'tumor_sample_with_clinical_intermediate_count',
                'tumor_sample_with_clinical_high_count',
                'mean_tumor_with_clinical_low_methylation',
                'mean_tumor_with_clinical_intermediate_methylation',
                'mean_tumor_with_clinical_high_methylation',
                'proportion_dead_in_low_methylation',
                'proportion_dead_in_high_methylation',
                'survival_direction_of_effect',
                'survival_p_value'
            )

            ## Return the vector
            return(relevant_vector)
        }

    }

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_genes==TRUE){

        ## Create survival plots for each of the genes designated:
        mclapply(
            X= top_gene_names,
            FUN= expression_survival_function_graph,
            high_cutoff= high_thresh,
            low_cutoff= low_thresh,
            graph= TRUE,
            mc.cores= cores
        )

    }

    ## Get survival information for each gene:
    gene_survival_results_list <- mclapply(
        X= top_gene_names,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= FALSE,
        mc.cores= cores
    )

    ## Transform the list into a matrix:
    gene_survival_results_df <- data.frame(
        matrix(
            unlist(gene_survival_results_list),
            nrow=length(gene_survival_results_list),
            byrow=T
        ),
        stringsAsFactors = FALSE
    )

    colnames(gene_survival_results_df) <- c(
        'normal_sample_count',
        'tumor_sample_count',
        'normal_sample_count_missing',
        'tumor_sample_count_missing',
        'mean_normal_expression',
        'mean_tumor_expression',
        'normal_sample_with_clinical_count',
        'tumor_sample_with_clinical_count',
        'tumor_sample_with_clinical_NA_count',
        'tumor_sample_with_clinical_low_count',
        'tumor_sample_with_clinical_intermediate_count',
        'tumor_sample_with_clinical_high_count',
        'mean_tumor_with_clinical_low_expression',
        'mean_tumor_with_clinical_intermediate_expression',
        'mean_tumor_with_clinical_high_expression',
        'proportion_dead_in_low_expression',
        'proportion_dead_in_high_expression',
        'survival_direction_of_effect',
        'survival_p_value'
    )

    rownames(gene_survival_results_df) <- c(
        gene_name <- gencode_v22_genes_df[
            top_gene_names, 'gene_name'
        ]
    )

    ## Save the gene information as a .tsv:
    write.table(
        gene_survival_results_df,
        file='./hypo.G+.output.survival/hypo.G+.top_genes_survival_info.tsv',
        quote= FALSE,
        sep='\t'
    )

    ## Now for each of the probes, we will need to calculate survival info and graphs:

    ## Get the list of probes linked to the significant genes:
    CpGs_linked <- CpG_linkage_dataset[
        CpG_linkage_dataset$geneID %in% top_gene_names,
        'probe'
    ]

    ## Get the unique probes from the list:
    unique_CpGs_linked <- unique(CpGs_linked)

    ## Get the survival info for each probe by mclapplying the
    ## methylation survival function to all of the probes
    ## linked to each gene:
    probe_survival_results_list <- mclapply(
        X= unique_CpGs_linked,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= FALSE,
        mc.cores= cores
    )

    ## Transform the list into a matrix:
    probe_survival_results_df <- data.frame(
        matrix(
            unlist(probe_survival_results_list),
            nrow=length(probe_survival_results_list),
            byrow=T
        ),
        stringsAsFactors = FALSE
    )

    colnames(probe_survival_results_df) <- c(
        'normal_sample_count',
        'tumor_sample_count',
        'normal_sample_count_missing',
        'tumor_sample_count_missing',
        'mean_normal_methylation',
        'mean_tumor_methylation',
        'normal_sample_with_clinical_count',
        'tumor_sample_with_clinical_count',
        'tumor_sample_with_clinical_NA_count',
        'tumor_sample_with_clinical_low_count',
        'tumor_sample_with_clinical_intermediate_count',
        'tumor_sample_with_clinical_high_count',
        'mean_tumor_with_clinical_low_methylation',
        'mean_tumor_with_clinical_intermediate_methylation',
        'mean_tumor_with_clinical_high_methylation',
        'proportion_dead_in_low_methylation',
        'proportion_dead_in_high_methylation',
        'survival_direction_of_effect',
        'survival_p_value'
    )

    rownames(probe_survival_results_df) <- unique_CpGs_linked

    ## For each of the CpGs, list which of the top X genes it was listed to:

    ## Wite a function to get the top X gene names a probe is assigned to:
    top_gene_assignment_function <- function(probe_of_interest){

        all_listed_genes <- CpG_linkage_dataset[
            CpG_linkage_dataset$probe==probe_of_interest,
            'geneSymbol'
        ]

        ## Conver thte top X genes to their gene symbols
        top_gene_symbols <- gencode_v22_genes_df[
            top_gene_names, 'gene_name'
        ]

        ## Get the names of the top genes that were linked to this probe:
        top_gene_symbols_linked_to_probe <- top_gene_symbols[
            top_gene_symbols %in% all_listed_genes
        ]

        ## Return the listed genes:
        return(
            paste(
                top_gene_symbols_linked_to_probe,
                collapse=','
            )
        )
    }

    ## Add the listed top genes to the survival results dataframe
    ## for the probes:
    probe_survival_results_df$top_genes_linked <- unname(
        sapply(
            rownames(probe_survival_results_df),
            top_gene_assignment_function
        )
    )

    ## Create a vector of the CpG p-values and add the CpG entries
    ## as names
    CpGs_linked_p_values_vector <- probe_survival_results_df$survival_p_value
    names(CpGs_linked_p_values_vector) <- rownames(probe_survival_results_df)

    ## Get the CpGs that are nominally significant:
    CpGs_linked_p_values_nominally_significant <- CpGs_linked_p_values_vector[
        CpGs_linked_p_values_vector<0.05
    ]

    ## Order the nominally significant CpGs:
    CpGs_linked_p_values_nominally_significant_ordered <- sort(
        CpGs_linked_p_values_nominally_significant,
        decreasing = FALSE
    )

    ## Get the names of the nominally significant CpGs:
    CpGs_linked_nominally_significant_names <- names(
        CpGs_linked_p_values_nominally_significant_ordered
    )

    ## If visualize_survival_plots_probes is TRUE, generate the plots:
    if(visualize_survival_plots_probes==TRUE){

        #### Save plots for the nominally significant CpGs:
        mclapply(
            X= CpGs_linked_nominally_significant_names,
            FUN= methylation_survival_function_graph,
            high_cutoff= high_thresh,
            low_cutoff= low_thresh,
            graph= TRUE,
            mc.cores= cores
        )

    }

    ## Save the probe information as a .tsv:
    write.table(
        probe_survival_results_df,
        file='./hypo.G+.output.survival/hypo.G+.probes_linked_to_top_genes_survival_info.tsv',
        quote= FALSE,
        sep='\t'
    )
}

make.survival.plots.hypometh.G_neg <- function(prefix, survival_top_n_genes, clinical_data, expDataN, expDataT, metDataN, metDataT, visualize_survival_plots_genes, high_thresh, low_thresh, cores, visualize_survival_plots_probes) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Load ggplot2:
    # library('ggplot2')

    ## Load survival:
    # library('survival')

    ## Load parallel:
    # library('parallel')

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G-.output.histogram/",
            prefix,
            ".hyper.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hyper.G-.output/hyper.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Set the number of genes to analyze based on user input:
    survival_top_n_genes <- survival_top_n_genes

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_names <- ordered_TFs_by_link_count[
        1:survival_top_n_genes,
        'geneID'
    ]

    ## Load the external data if it is available:

    ## Check if .rda file with expression, methylation, and
    ## perhaps clinical data was placed in TENET folder:
    combined_rda_file <- list.files(
        path='../external.data/data',
        pattern='\\.rda$',
        full.names= TRUE,
        include.dirs = TRUE
    )

    ## If the combined rda file exists, load it and see if clinical data is included:
    ## If not, check to see if clinical data is included in the "clinical" folder
    load(combined_rda_file)

    # #### REMOVE THIS LATER:
    # load(
    #     "/Volumes/NO_NAME/TENET_paper_update/tenet_match.matched.dedup.luad.methylation.expression.dataset.rda",
    # )

    ## If data is loaded as clinical_data, rename it to clinical:
    if(exists('clinical_data')){

        clinical <- clinical_data

    }

    ## Check to see if the clinical data is properly loaded
    ## if it isn't, check the clinical folder in external.data/data
    if(exists('clinical')){

        paste('clinical data loaded')

    } else{

        ## Get a path to the rda file containing clinical info
        ## in external.data/data/clinical
        clinical_rda_file <- list.files(
            path='../external.data/data/clinical',
            pattern='\\.rda$',
            full.names= TRUE,
            include.dirs = TRUE
        )

        if(length(clinical_rda_file)==1){

            ## Load the .rda file:
            load(clinical_rda_file)

            ## If data is loaded as clinical_data, rename it to clinical:
            if(exists('clinical_data')){

                clinical <- clinical_data

            }

            ## Check to see if the clinical data is now properly loaded
            if(exists('clinical')){

                paste('clinical data loaded')

            } else{

                ## CLinical data is not found in rda
                ## return an error message and abort script:
                paste('rda file did not contain clinical data or it was not saved as clinical or clinical_data')
                quit(save='no')
            }

        } else{

            ## Clinical rda not found.
            ## Return error message and abort script
            ## Later functionality may be added to try loading clinical data from txt file
            ## See step 0 make.rda.from.textfiles.R for reference:
            paste('rda file with clinical data not found')
            quit(save='no')
        }
    }

    ## Process the expression/methylation/clinical data so the function doesn't have to:

    ## combine datasets:
    expData <- cbind(expDataN, expDataT)
    metData <- cbind(metDataN, metDataT)

    ## Get the list of normal and tumor samples:
    exp_normal_samples <- colnames(expDataN)
    exp_tumor_samples <- colnames(expDataT)

    met_normal_samples <- colnames(metDataN)
    met_tumor_samples <- colnames(metDataT)

    ## Get the relevant columns from the clinical data:
    relevant_clinical <- clinical[
        ,c(
            "bcr_patient_barcode",
            "days_to_death",
            "days_to_last_followup",
            "vital_status"
        )
    ]
    rownames(relevant_clinical) <- relevant_clinical$bcr_patient_barcode

    # If there are subjects that are alive in the
    # dataset, set their days to death info to '-Inf'
    relevant_clinical[
        grep(
            "alive",
            relevant_clinical$vital_status,
            ignore.case = TRUE
        ),
        'days_to_death'
    ] <- '-Inf'

    # If there are subjects that are dead in the
    # dataset, set their days to last followup info to '-Inf'
    relevant_clinical[
        grep(
            "dead",
            relevant_clinical$vital_status,
            ignore.case = TRUE
        ),
        'days_to_last_followup'
    ] <- '-Inf'

    ## Remove subjects that remain unaccounted for:
    ## i.e. NA values remain in previous two columns:
    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_death"]
        ),
    ]

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_last_followup"]
        ),
    ]

    # Changing the days to last followup info in the cancer clinical data
    # to be numeric values through character values and remove any NAs that are induced
    # as a final check:
    relevant_clinical$days_to_death <- as.numeric(
        as.character(relevant_clinical$days_to_death)
    )

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_last_followup"]
        ),
    ]

    # Changing the days to death info in the cancer clinical data
    # to be numeric values through character values and remove any NAs that are induced
    # as a final check:
    relevant_clinical$days_to_last_followup <- as.numeric(
        as.character(relevant_clinical$days_to_last_followup)
    )

    relevant_clinical <- relevant_clinical[
        !is.na(
            relevant_clinical[,"days_to_death"]
        ),
    ]

    # Add the relevant number as a final "time" variable
    # Equal to days to death for dead individuals, and
    # days to last followup for alive individuals
    relevant_clinical$time <- ifelse(
        relevant_clinical$vital_status=='Alive',
        relevant_clinical$days_to_last_followup,
        ifelse(
            relevant_clinical$vital_status=='Dead',
            relevant_clinical$days_to_death,
            NA
        )
    )

    ## Create a function to get the survival p-value or graph
    ## for each gene of interest:
    expression_survival_function_graph <- function(
        gene_of_interest,
        high_cutoff,
        low_cutoff,
        graph
    ){

        ## Get actual gene name depending on input:
        if(substring(gene_of_interest,1,4)=='ENSG' & nchar(gene_of_interest)==15){

            ## Input is in ENSG, get the gene name:
            gene_ENSG <- gene_of_interest

            gene_name <- gencode_v22_genes_df[
                gene_ENSG, 'gene_name'
            ]

            ## Get gene ENSG assuming a name is plugged in:
        } else{

            ## Assume what was given was the gene name, get the ENSG:
            gene_name <- gene_of_interest

            gene_ENSG <- rownames(
                gencode_v22_genes_df[gencode_v22_genes_df$gene_name==gene_name,]
            )
        }

        ## Get expression values for gene of interest:
        expression_values <- unlist(
            expData[
                gene_ENSG,
            ]
        )
        names(expression_values) <- colnames(expData)

        # Split gene expression values into normal and tumor samples
        tumor_expression_values <- expression_values[exp_tumor_samples]
        normal_expression_values <- expression_values[exp_normal_samples]

        # Changing the sample names to match the subject names:
        names(tumor_expression_values)  <- substr(
            names(tumor_expression_values),
            1,
            12
        )

        names(normal_expression_values)  <- substr(
            names(normal_expression_values),
            1,
            12
        )

        ## Calculate some basic data:
        normal_sample_n <- length(normal_expression_values)
        tumor_sample_n <- length(tumor_expression_values)

        ## Count the number of NA samples:
        NA_normal <- as.numeric(
            unname(
                table(
                    is.na(normal_expression_values)
                )['TRUE']
            )
        )

        NA_tumor <- as.numeric(
            unname(
                table(
                    is.na(tumor_expression_values)
                )['TRUE']
            )
        )

        if(is.na(NA_normal)){

            NA_normal <- 0

        }

        if(is.na(NA_tumor)){

            NA_tumor <- 0

        }

        ## Calculate mean expression for all samples present:
        mean_tumor_expression <- mean(
            tumor_expression_values,
            na.rm= TRUE
        )

        mean_normal_expression <- mean(
            normal_expression_values,
            na.rm= TRUE
        )

        ## Calculate the number of samples with present clinical data:
        present_clinical_normal_sample_n <- as.numeric(
            unname(
                table(
                    names(normal_expression_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        present_clinical_tumor_sample_n  <- as.numeric(
            unname(
                table(
                    names(tumor_expression_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        ## Identify the samples which are present in the clinical data:
        present_normal_samples <- names(normal_expression_values)[
            names(normal_expression_values) %in% rownames(relevant_clinical)
        ]

        present_tumor_samples <- names(tumor_expression_values)[
            names(tumor_expression_values) %in% rownames(relevant_clinical)
        ]


        # Subsetting clinical patient data for individuals with
        # gene expression data in the tumor set:
        function_relevant_clinical <- relevant_clinical[
            present_tumor_samples,
        ]

        ## Calculate quantiles:
        high_cutoff_quantile= quantile(
            tumor_expression_values,
            high_cutoff,
            na.rm= TRUE
        )[1]

        low_cutoff_quantile= quantile(
            tumor_expression_values,
            low_cutoff,
            na.rm= TRUE
        )[1]

        ## Determine if a sample is in the high, low, or intermediate quartiles:
        sample_expression_grouping <- ifelse(
            tumor_expression_values >  high_cutoff_quantile,
            "High",
            ifelse(
                tumor_expression_values <= low_cutoff_quantile,
                "Low",
                ifelse(
                    is.na(tumor_expression_values),
                    'Missing',
                    'Intermediate'
                )
            )
        )

        ## Add the expression and relevant clinical information to the relevant clinical dataframe:
        function_relevant_clinical$tumor_expression <- tumor_expression_values[
            rownames(function_relevant_clinical)
        ]

        ## Add the expression grouping to the clinical data:
        function_relevant_clinical$expression_grouping <- sample_expression_grouping[
            rownames(function_relevant_clinical)
        ]

        ## Remove samples missing expression data:
        function_relevant_clinical_complete <- na.omit(
            function_relevant_clinical,
            cols=c(
                "tumor_expression",
                "expression_grouping"
            )
        )

        function_relevant_clinical_complete <- function_relevant_clinical

        ## Count the number of samples of each group remaining:
        present_tumor_sample_high_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='High',
            ]
        )

        present_tumor_sample_intermediate_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Intermediate',
            ]
        )

        present_tumor_sample_low_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Low',
            ]
        )

        present_tumor_sample_missing_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Missing',
            ]
        )

        ## Get the mean expression in each group:
        present_tumor_sample_high_mean <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='High',
                'tumor_expression'
            ]
        )

        present_tumor_sample_intermediate_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Intermediate',
                'tumor_expression'
            ]
        )

        present_tumor_sample_low_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$expression_grouping=='Low',
                'tumor_expression'
            ]
        )

        ## Get clinical results for only samples in the high and low group with info:
        function_relevant_clinical_complete_high_low <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
        ]

        ## Calculate proportion of deceased individuals in each group
        proportion_dead_high <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$expression_grouping=='High' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$expression_grouping=='High'),
                    ]
                )
        )

        proportion_dead_low <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$expression_grouping=='Low' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$expression_grouping=='Low'),
                    ]
                )
        )

        ## Determine which group had greater proportion of fatalities:
        highest_dead_proportion_group <- ifelse(
            proportion_dead_high > proportion_dead_low,
            'high_expression_low_survival',
            ifelse(
                proportion_dead_high < proportion_dead_low,
                'low_expression_low_survival',
                'unclear'
            )
        )

        ## Create a high_low survival object:
        high_low_survival_object <- Surv(
            function_relevant_clinical_complete_high_low$time,
            ifelse(
                function_relevant_clinical_complete_high_low$vital_status=='Alive',
                FALSE,
                TRUE
            )
        )

        # Set the rownames of the survival object to be the
        # patient names for the high/low group smaples:
        rownames(high_low_survival_object) <- rownames(
            function_relevant_clinical_complete[
                !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
            ]
        )

        # Combine the top and down vectors into a single vector
        # in the order top-down
        expression_group <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
            'expression_grouping'
        ]

        # Creating legend names for high/low expression groups:
        legend_name_high <- paste(gene_name,"high")

        legend_name_low <- paste(gene_name,"low")

        # Perform the survival analysis.
        # This uses the combined vector as the x variable
        # and the survival data as the y variable to create
        # a table with information about the test
        # including chi-squared p-value
        survival_table <- survdiff(
            high_low_survival_object  ~ expression_group
        )

        # Get the chi-squared test statistic from the analysis above:
        survival_table_chi_squared <- survival_table$chisq

        # Calculating a p value based on the test statistic to get
        # a precise p-value
        survival_pvalue <- as.numeric(
            1 - pchisq(
                abs(survival_table_chi_squared), df = 1
            )
        )

        ## Round the p-value on the graph to 3 digits:
        survival_pvalue_formatted <- formatC(
            survival_pvalue,
            format='e',
            digits=3
        )

        ## Create the graph if specified:
        if(graph==TRUE){

            # Create the plot title
            # with gene name and p-value included:
            survival_title <- paste(
                gene_name,
                "\nKaplan-Meier Survival analysis\nP = ",
                survival_pvalue_formatted,
                sep=''
            )

            ## Get the path to the step 5 folder to save
            ## the pdf in:
            path_to_folder <- getwd()

            ## Create a title for the survival plot pdf:
            survival_plot_pdf_title <- paste(
                path_to_folder,
                '/hyper.G-.output.survival/',
                gene_name,
                '_survival_plot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(survival_plot_pdf_title)

            ## Now actually create the survival plot
            ## Using a similar structure to what was used to generate the p-value:
            plot(
                survfit(
                    high_low_survival_object ~ expression_group
                ),

                # Color the lines (high in red first!)
                col = c('red', 'black'),

                ## Add thickness to the lines:
                lwd=3,

                # Use the title that was created earlier as the title
                # of the plot:
                main= survival_title,

                # Set titles of the x and y axis:
                # Note: TCGA measures survival in days as noted:
                xlab="Time (days)",
                ylab="Probability",

                # Change axis size:
                cex.axis=1,
                cex.main=1.5,
                cex.lab=1.25
            )

            ## Add a legend to the plot:
            legend(

                # Set X position of legend in graph:
                x= (
                    max(function_relevant_clinical_complete_high_low$days_to_last_followup) - 2500
                ),

                # Set Y position of legend in graph
                y= 1,

                # Use the legend titles that were created earlier
                legend= c(legend_name_high, legend_name_low),

                # As above, use black for low and red for high:
                col= c('red', 'black'),

                # Coloring the text in the legend as well
                text.col= c('red', 'black'),

                # Change the shape of the labels in the legend:
                pch= 15
            )

            ## Close the plot:
            dev.off()

        } else if(graph==FALSE){

            ## create vector of relevant info:
            relevant_vector <- c(
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(NA_normal),
                as.numeric(NA_tumor),
                as.numeric(mean_normal_expression),
                as.numeric(mean_tumor_expression),
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(present_tumor_sample_missing_n),
                as.numeric(present_tumor_sample_low_n),
                as.numeric(present_tumor_sample_intermediate_n),
                as.numeric(present_tumor_sample_high_n),
                as.numeric(present_tumor_sample_low_mean),
                as.numeric(present_tumor_sample_intermediate_mean),
                as.numeric(present_tumor_sample_high_mean),
                as.numeric(proportion_dead_low),
                as.numeric(proportion_dead_high),
                highest_dead_proportion_group,
                as.numeric(survival_pvalue)
            )
            names(relevant_vector) <- c(
                'normal_sample_count',
                'tumor_sample_count',
                'normal_sample_count_missing',
                'tumor_sample_count_missing',
                'mean_normal_expression',
                'mean_tumor_expression',
                'normal_sample_with_clinical_count',
                'tumor_sample_with_clinical_count',
                'tumor_sample_with_clinical_NA_count',
                'tumor_sample_with_clinical_low_count',
                'tumor_sample_with_clinical_intermediate_count',
                'tumor_sample_with_clinical_high_count',
                'mean_tumor_with_clinical_low_expression',
                'mean_tumor_with_clinical_intermediate_expression',
                'mean_tumor_with_clinical_high_expression',
                'proportion_dead_in_low_expression',
                'proportion_dead_in_high_expression',
                'survival_direction_of_effect',
                'survival_p_value'
            )

            ## Return the vector
            return(relevant_vector)
        }

    }

    ## Create a function to get the survival p-value or graph for
    ##beach probe of interest:
    methylation_survival_function_graph <- function(
        probe_of_interest,
        high_cutoff,
        low_cutoff,
        graph
    ){

        ## Get methylation values for probe of interest:
        methylation_values <- unlist(
            metData[
                probe_of_interest,
            ]
        )
        names(methylation_values) <- colnames(metData)

        # Split methylation values into normal and tumor samples
        # If both are present in the dataset
        tumor_methylation_values <- methylation_values[met_tumor_samples]
        normal_methylation_values <- methylation_values[met_normal_samples]

        # Changing the sample names to match the subject names:
        names(tumor_methylation_values)  <- substr(
            names(tumor_methylation_values),
            1,
            12
        )

        names(normal_methylation_values)  <- substr(
            names(normal_methylation_values),
            1,
            12
        )

        ## Calculate some basic data:
        normal_sample_n <- length(normal_methylation_values)
        tumor_sample_n <- length(tumor_methylation_values)

        ## Count the number of NA samples:
        NA_normal <- as.numeric(
            unname(
                table(
                    is.na(normal_methylation_values)
                )['TRUE']
            )
        )

        NA_tumor <- as.numeric(
            unname(
                table(
                    is.na(tumor_methylation_values)
                )['TRUE']
            )
        )

        if(is.na(NA_normal)){

            NA_normal <- 0

        }

        if(is.na(NA_tumor)){

            NA_tumor <- 0

        }

        ## Calculate mean methylation for all samples present:
        mean_tumor_methylation <- mean(
            tumor_methylation_values,
            na.rm= TRUE
        )

        mean_normal_methylation <- mean(
            normal_methylation_values,
            na.rm= TRUE
        )

        ## Calculate the number of samples with present clinical data:
        present_clinical_normal_sample_n <- as.numeric(
            unname(
                table(
                    names(normal_methylation_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        present_clinical_tumor_sample_n  <- as.numeric(
            unname(
                table(
                    names(tumor_methylation_values) %in% rownames(relevant_clinical)
                )['TRUE']
            )
        )

        ## Identify the samples which are present in the clinical data:
        present_normal_samples <- names(normal_methylation_values)[
            names(normal_methylation_values) %in% rownames(relevant_clinical)
        ]

        present_tumor_samples <- names(tumor_methylation_values)[
            names(tumor_methylation_values) %in% rownames(relevant_clinical)
        ]


        # Subsetting clinical patient data for individuals with
        # DNA methylation data in the tumor set:
        function_relevant_clinical <- relevant_clinical[
            present_tumor_samples,
        ]

        ## Calculate quantiles:
        high_cutoff_quantile= quantile(
            tumor_methylation_values,
            high_cutoff,
            na.rm= TRUE
        )[1]

        low_cutoff_quantile= quantile(
            tumor_methylation_values,
            low_cutoff,
            na.rm= TRUE
        )[1]

        ## Determine if a sample is in the high, low, or intermediate quartiles:
        sample_methylation_grouping <- ifelse(
            tumor_methylation_values >  high_cutoff_quantile,
            "High",
            ifelse(
                tumor_methylation_values <= low_cutoff_quantile,
                "Low",
                ifelse(
                    is.na(tumor_methylation_values),
                    'Missing',
                    'Intermediate'
                )
            )
        )

        ## Add the methylation and relevant clinical information to the relevant clinical dataframe:
        function_relevant_clinical$tumor_methylation <- tumor_methylation_values[
            rownames(function_relevant_clinical)
        ]

        ## Add the methylation grouping to the clinical data:
        function_relevant_clinical$methylation_grouping <- sample_methylation_grouping[
            rownames(function_relevant_clinical)
        ]

        ## Remove samples missing methylation data:
        function_relevant_clinical_complete <- na.omit(
            function_relevant_clinical,
            cols=c(
                "tumor_methylation",
                "methylation_grouping"
            )
        )
        function_relevant_clinical_complete <- function_relevant_clinical

        ## Count the number of samples of each group remaining:
        present_tumor_sample_high_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='High',
            ]
        )

        present_tumor_sample_intermediate_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Intermediate',
            ]
        )

        present_tumor_sample_low_n  <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Low',
            ]
        )

        present_tumor_sample_missing_n <- nrow(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Missing',
            ]
        )

        ## Get the mean methylation in each group:
        present_tumor_sample_high_mean <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='High',
                'tumor_methylation'
            ]
        )

        present_tumor_sample_intermediate_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Intermediate',
                'tumor_methylation'
            ]
        )

        present_tumor_sample_low_mean  <- mean(
            function_relevant_clinical_complete[
                function_relevant_clinical_complete$methylation_grouping=='Low',
                'tumor_methylation'
            ]
        )

        ## Get clinical results for only samples in the high and low group with info:
        function_relevant_clinical_complete_high_low <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
        ]

        ## Calculate proportion of deceased individuals in each group
        proportion_dead_high <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$methylation_grouping=='High' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$methylation_grouping=='High'),
                    ]
                )
        )

        proportion_dead_low <- as.numeric(
            nrow(
                function_relevant_clinical_complete_high_low[
                    (function_relevant_clinical_complete_high_low$methylation_grouping=='Low' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
                ]
            ) /
                nrow(
                    function_relevant_clinical_complete_high_low[
                        (function_relevant_clinical_complete_high_low$methylation_grouping=='Low'),
                    ]
                )
        )

        ## Determine which group had greater proportion of fatalities:
        highest_dead_proportion_group <- ifelse(
            proportion_dead_high > proportion_dead_low,
            'high_methylation_low_survival',
            ifelse(
                proportion_dead_high < proportion_dead_low,
                'low_methylation_low_survival',
                'unclear'
            )
        )

        ## Create a high_low survival object:
        high_low_survival_object <- Surv(
            function_relevant_clinical_complete_high_low$time,
            ifelse(
                function_relevant_clinical_complete_high_low$vital_status=='Alive',
                FALSE,
                TRUE
            )
        )

        # Set the rownames of the survival object to be the
        # patient names for the high/low group smaples:
        rownames(high_low_survival_object) <- rownames(
            function_relevant_clinical_complete[
                !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
            ]
        )

        # Combine the top and down vectors into a single vector
        # in the order top-down
        methylation_group <- function_relevant_clinical_complete[
            !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
            'methylation_grouping'
        ]

        # Creating legend names for high/low methylation groups:
        legend_name_high <- paste(probe_of_interest,"high")

        legend_name_low <- paste(probe_of_interest,"low")

        # Perform the survival analysis.
        # This uses the combined vector as the x variable
        # and the survival data as the y variable to create
        # a table with information about the test
        # including chi-squared p-value
        survival_table <- survdiff(
            high_low_survival_object  ~ methylation_group
        )

        # Get the chi-squared test statistic from the analysis above:
        survival_table_chi_squared <- survival_table$chisq

        # Calculating a p value based on the test statistic to get
        # a precise p-value
        survival_pvalue <- as.numeric(
            1 - pchisq(
                abs(survival_table_chi_squared), df = 1
            )
        )

        ## Round the p-value on the graph to 3 digits:
        survival_pvalue_formatted <- formatC(
            survival_pvalue,
            format='e',
            digits=3
        )

        ## Create the graph if specified:
        if(graph==TRUE){

            # Create the plot title
            # with gene name and p-value included:
            survival_title <- paste(
                probe_of_interest,
                "\nKaplan-Meier Survival analysis\nP = ",
                survival_pvalue_formatted,
                sep=''
            )

            ## Get the path to the step 5 folder to save
            ## the pdf in:
            path_to_folder <- getwd()

            ## Create a title for the survival plot pdf:
            ## This is a comination of the probe name with the linked gene:
            survival_plot_pdf_title <- paste(
                path_to_folder,
                '/hyper.G-.output.survival/',
                probe_of_interest,
                '_survival_plot.pdf',
                sep=''
            )

            ## Open a pdf for saving the plot:
            pdf(survival_plot_pdf_title)

            ## Now actually create the survival plot
            ## Using a similar structure to what was used to generate the p-value:
            plot(
                survfit(
                    high_low_survival_object ~ methylation_group
                ),

                # Color the lines (high in red first!)
                col = c('red', 'black'),

                ## Add thickness to the lines:
                lwd=3,

                # Use the title that was created earlier as the title
                # of the plot:
                main= survival_title,

                # Set titles of the x and y axis:
                # Note: TCGA measures survival in days as noted:
                xlab="Time (days)",
                ylab="Probability",

                # Change axis size:
                cex.axis=1,
                cex.main=1.5,
                cex.lab=1.25
            )

            ## Add a legend to the plot:
            legend(

                # Set X position of legend in graph:
                x= (
                    max(function_relevant_clinical_complete_high_low$days_to_last_followup) - 2500
                ),

                # Set Y position of legend in graph
                y= 1,

                # Use the legend titles that were created earlier
                legend= c(legend_name_high, legend_name_low),

                # As above, use black for low and red for high:
                col= c('red', 'black'),

                # Coloring the text in the legend as well
                text.col= c('red', 'black'),

                # Change the shape of the labels in the legend:
                pch= 15
            )

            ## Close the plot:
            dev.off()

        } else if(graph==FALSE){

            ## create vector of relevant info:
            relevant_vector <- c(
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(NA_normal),
                as.numeric(NA_tumor),
                as.numeric(mean_normal_methylation),
                as.numeric(mean_tumor_methylation),
                as.numeric(normal_sample_n),
                as.numeric(tumor_sample_n),
                as.numeric(present_tumor_sample_missing_n),
                as.numeric(present_tumor_sample_low_n),
                as.numeric(present_tumor_sample_intermediate_n),
                as.numeric(present_tumor_sample_high_n),
                as.numeric(present_tumor_sample_low_mean),
                as.numeric(present_tumor_sample_intermediate_mean),
                as.numeric(present_tumor_sample_high_mean),
                as.numeric(proportion_dead_low),
                as.numeric(proportion_dead_high),
                highest_dead_proportion_group,
                as.numeric(survival_pvalue)
            )
            names(relevant_vector) <- c(
                'normal_sample_count',
                'tumor_sample_count',
                'normal_sample_count_missing',
                'tumor_sample_count_missing',
                'mean_normal_methylation',
                'mean_tumor_methylation',
                'normal_sample_with_clinical_count',
                'tumor_sample_with_clinical_count',
                'tumor_sample_with_clinical_NA_count',
                'tumor_sample_with_clinical_low_count',
                'tumor_sample_with_clinical_intermediate_count',
                'tumor_sample_with_clinical_high_count',
                'mean_tumor_with_clinical_low_methylation',
                'mean_tumor_with_clinical_intermediate_methylation',
                'mean_tumor_with_clinical_high_methylation',
                'proportion_dead_in_low_methylation',
                'proportion_dead_in_high_methylation',
                'survival_direction_of_effect',
                'survival_p_value'
            )

            ## Return the vector
            return(relevant_vector)
        }

    }

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_genes==TRUE){

        ## Create survival plots for each of the genes designated:
        mclapply(
            X= top_gene_names,
            FUN= expression_survival_function_graph,
            high_cutoff= high_thresh,
            low_cutoff= low_thresh,
            graph= TRUE,
            mc.cores= cores
        )

    }

    ## Get survival information for each gene:
    gene_survival_results_list <- mclapply(
        X= top_gene_names,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= FALSE,
        mc.cores= cores
    )

    ## Transform the list into a matrix:
    gene_survival_results_df <- data.frame(
        matrix(
            unlist(gene_survival_results_list),
            nrow=length(gene_survival_results_list),
            byrow=T
        ),
        stringsAsFactors = FALSE
    )

    colnames(gene_survival_results_df) <- c(
        'normal_sample_count',
        'tumor_sample_count',
        'normal_sample_count_missing',
        'tumor_sample_count_missing',
        'mean_normal_expression',
        'mean_tumor_expression',
        'normal_sample_with_clinical_count',
        'tumor_sample_with_clinical_count',
        'tumor_sample_with_clinical_NA_count',
        'tumor_sample_with_clinical_low_count',
        'tumor_sample_with_clinical_intermediate_count',
        'tumor_sample_with_clinical_high_count',
        'mean_tumor_with_clinical_low_expression',
        'mean_tumor_with_clinical_intermediate_expression',
        'mean_tumor_with_clinical_high_expression',
        'proportion_dead_in_low_expression',
        'proportion_dead_in_high_expression',
        'survival_direction_of_effect',
        'survival_p_value'
    )

    rownames(gene_survival_results_df) <- c(
        gene_name <- gencode_v22_genes_df[
            top_gene_names, 'gene_name'
        ]
    )

    ## Save the gene information as a .tsv:
    write.table(
        gene_survival_results_df,
        file='./hyper.G-.output.survival/hyper.G-.top_genes_survival_info.tsv',
        quote= FALSE,
        sep='\t'
    )

    ## Now for each of the probes, we will need to calculate survival info and graphs:

    ## Get the list of probes linked to the significant genes:
    CpGs_linked <- CpG_linkage_dataset[
        CpG_linkage_dataset$geneID %in% top_gene_names,
        'probe'
    ]

    ## Get the unique probes from the list:
    unique_CpGs_linked <- unique(CpGs_linked)

    ## Get the survival info for each probe by mclapplying the
    ## methylation survival function to all of the probes
    ## linked to each gene:
    probe_survival_results_list <- mclapply(
        X= unique_CpGs_linked,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= FALSE,
        mc.cores= cores
    )

    ## Transform the list into a matrix:
    probe_survival_results_df <- data.frame(
        matrix(
            unlist(probe_survival_results_list),
            nrow=length(probe_survival_results_list),
            byrow=T
        ),
        stringsAsFactors = FALSE
    )

    colnames(probe_survival_results_df) <- c(
        'normal_sample_count',
        'tumor_sample_count',
        'normal_sample_count_missing',
        'tumor_sample_count_missing',
        'mean_normal_methylation',
        'mean_tumor_methylation',
        'normal_sample_with_clinical_count',
        'tumor_sample_with_clinical_count',
        'tumor_sample_with_clinical_NA_count',
        'tumor_sample_with_clinical_low_count',
        'tumor_sample_with_clinical_intermediate_count',
        'tumor_sample_with_clinical_high_count',
        'mean_tumor_with_clinical_low_methylation',
        'mean_tumor_with_clinical_intermediate_methylation',
        'mean_tumor_with_clinical_high_methylation',
        'proportion_dead_in_low_methylation',
        'proportion_dead_in_high_methylation',
        'survival_direction_of_effect',
        'survival_p_value'
    )

    rownames(probe_survival_results_df) <- unique_CpGs_linked

    ## For each of the CpGs, list which of the top X genes it was listed to:

    ## Wite a function to get the top X gene names a probe is assigned to:
    top_gene_assignment_function <- function(probe_of_interest){

        all_listed_genes <- CpG_linkage_dataset[
            CpG_linkage_dataset$probe==probe_of_interest,
            'geneSymbol'
        ]

        ## Conver thte top X genes to their gene symbols
        top_gene_symbols <- gencode_v22_genes_df[
            top_gene_names, 'gene_name'
        ]

        ## Get the names of the top genes that were linked to this probe:
        top_gene_symbols_linked_to_probe <- top_gene_symbols[
            top_gene_symbols %in% all_listed_genes
        ]

        ## Return the listed genes:
        return(
            paste(
                top_gene_symbols_linked_to_probe,
                collapse=','
            )
        )
    }

    ## Add the listed top genes to the survival results dataframe
    ## for the probes:
    probe_survival_results_df$top_genes_linked <- unname(
        sapply(
            rownames(probe_survival_results_df),
            top_gene_assignment_function
        )
    )

    ## Create a vector of the CpG p-values and add the CpG entries
    ## as names
    CpGs_linked_p_values_vector <- probe_survival_results_df$survival_p_value
    names(CpGs_linked_p_values_vector) <- rownames(probe_survival_results_df)

    ## Get the CpGs that are nominally significant:
    CpGs_linked_p_values_nominally_significant <- CpGs_linked_p_values_vector[
        CpGs_linked_p_values_vector<0.05
    ]

    ## Order the nominally significant CpGs:
    CpGs_linked_p_values_nominally_significant_ordered <- sort(
        CpGs_linked_p_values_nominally_significant,
        decreasing = FALSE
    )

    ## Get the names of the nominally significant CpGs:
    CpGs_linked_nominally_significant_names <- names(
        CpGs_linked_p_values_nominally_significant_ordered
    )

    ## If visualize_survival_plots_probes is TRUE, generate the plots:
    if(visualize_survival_plots_probes==TRUE){

        #### Save plots for the nominally significant CpGs:
        mclapply(
            X= CpGs_linked_nominally_significant_names,
            FUN= methylation_survival_function_graph,
            high_cutoff= high_thresh,
            low_cutoff= low_thresh,
            graph= TRUE,
            mc.cores= cores
        )

    }

    ## Save the probe information as a .tsv:
    write.table(
        probe_survival_results_df,
        file='./hyper.G-.output.survival/hyper.G-.probes_linked_to_top_genes_survival_info.tsv',
        quote= FALSE,
        sep='\t'
    )
}

make.TAD.tables.hypermeth.G_pos <- function(prefix, TAD_top_n_genes, i, gene_ENSG, grep, TAD_file_index) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Load the gtf file for CpG locations:
    hg38_450k_CpG_df <- read.delim(
        '../scripts/data/hm450.hg38.manifest.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Probes are 0-indexed, let's fix this:
    hg38_450k_CpG_df$CpG_beg <- hg38_450k_CpG_df$CpG_beg+1

    ## Set the rownames to be the probe IDs:
    rownames(hg38_450k_CpG_df) <- hg38_450k_CpG_df$probeID

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G+.output.histogram/",
            prefix,
            ".hyper.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hyper.G+.output/hyper.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_IDs <- ordered_TFs_by_link_count[
        1:TAD_top_n_genes,
        'geneID'
    ]

    ## Load TAD files to overlap with:
    TAD_file_paths <- list.files(
        path='../external.data/TAD',
        full.names=TRUE
    )

    ## Create an empty list to load each of the files into:
    TAD_file_list <- vector(
        mode = "list",
        length = length(TAD_file_paths)
    )

    ## Remove the .rtf file from being loaded as a TAD file:
    TAD_file_paths <- TAD_file_paths[
        !(sub('.*\\.', '', TAD_file_paths)=='rtf')
    ]

    ## Prepare and load the TAD files:
    ## This assumes the files are bed3 formatted with the first 3 columns:
    ## Being the chromosome, start, and end coordinates:
    ## And are 0-indexed:
    for(i in c(1:length(TAD_file_paths))){

        ## Load the file:
        placeholder_file <- read.delim(
            TAD_file_paths[i],
            header= FALSE,
            stringsAsFactors = FALSE
        )

        ## Remove extra columns and rename the first three columns:
        ## Which should contain chr, start, end locations:
        if(ncol(placeholder_file)>=3){

            placeholder_file <- placeholder_file[1:3]

            colnames(placeholder_file) <- c(
                'seqnames','start','end'
            )

        } else{

            print("File is improperly formatted")
            break()
        }

        ## Add 1 to the starts to make them 1-indexed:
        placeholder_file$start <- placeholder_file$start+1

        ## Add the final processed TAD file to the list:
        TAD_file_list[[i]] <- placeholder_file
    }

    ## Assign names of the files to the TAD file list:
    TAD_file_names <- sub(
        '\\..*',
        '',
        basename(TAD_file_paths)
    )

    names(TAD_file_list) <- TAD_file_names

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes selected:
    probes_linked_to_significant_genes <- unique(
        CpG_linkage_dataset[
            CpG_linkage_dataset$geneID %in% top_gene_IDs,
            'probe'
        ]
    )

    ## Create a data frame with the probes as the first (and only column for now):
    probe_dataset <- data.frame(
        'probe_ID'=sort(probes_linked_to_significant_genes),
        stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info for hg38:
    probe_dataset$seqnames <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_chrm'
    ]

    probe_dataset$start <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_beg'
    ]

    probe_dataset$end <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_end'
    ]

    ## Create names of columns for each of the probes:
    for(gene_ENSG in top_gene_IDs){

        ## Get the gene name corresponding to the ENSG:
        corresponding_gene_name <- ordered_TFs_by_link_count[
            ordered_TFs_by_link_count$geneID==gene_ENSG,
            'geneSymbol'
        ]

        column_name <- paste(
            corresponding_gene_name,
            '_linked',
            sep=''
        )

        ## Get the probes linked to that gene specifically:
        probes_linked_to_gene <- unique(
            CpG_linkage_dataset[
                CpG_linkage_dataset$geneID %in% gene_ENSG,
                'probe'
            ]
        )

        ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
        ## at least one of the top n TFs is linked to this specific one:
        probe_dataset[column_name] <- probe_dataset$probe_ID %in% probes_linked_to_gene

    }

    ## Write functions to identify genes and probes within each TAD:
    TAD_probe_overlapper <- function(chromosome, start, stop, reference_peaks, buffer){
        paste(
            rownames(
                reference_peaks[
                    (
                        (reference_peaks[['seqnames']] == chromosome) &
                            ((reference_peaks[['start']]-buffer) <= stop) &
                            ((reference_peaks[['end']]+buffer) >= start)
                    ),
                ]
            ), collapse=","
        )
    }

    TAD_gene_overlapper <- function(chromosome, TSS, reference_peaks, buffer){
        paste(
            rownames(
                reference_peaks[
                    (
                        (reference_peaks[['seqnames']] == chromosome) &
                            (((reference_peaks[['start']]-buffer) <= TSS) &
                                 ((reference_peaks[['end']]+buffer) >= TSS))
                    ),
                ]
            ), collapse=","
        )
    }

    ## Now let's find genes in the same TAD(s):
    TAD_gene_finder <- function(TAD_numbers_comma_delimited, target_list){

        TAD_numbers <- unlist(strsplit(TAD_numbers_comma_delimited,','))

        gene_numbers <- sapply(
            TAD_numbers,
            grep,
            target_list
        )

        unique_gene_numbers <- unique(unlist(gene_numbers))

        unique_gene_numbers_combined <- paste(unique_gene_numbers,collapse=',')

        return(unique_gene_numbers_combined)
    }

    ## Write a function to count the genes in a TAD with a given probe:
    gene_count_in_TAD <- function(gene_numbers){

        gene_numbers_separated <- as.numeric(unlist(strsplit(gene_numbers,',')))

        return(length(gene_numbers_separated))
    }

    ## Now let's write a function to convert the gene numbers to the gene name
    ## based on setup of siX_downregulated_gencode_v22_genes_df files:
    gene_name_from_number_lister <- function(gene_numbers, downregulated_gene_df, return_type){

        gene_numbers_separated <- as.numeric(unlist(strsplit(gene_numbers,',')))

        if(return_type=='gene_name'){

            gene_names <- downregulated_gene_df[
                gene_numbers_separated,
                'gene_name'
            ]

        } else if(return_type=='gene_id'){

            gene_names <- downregulated_gene_df[
                gene_numbers_separated,
                'gene_id'
            ]
        }

        unique_gene_names_combined <- paste(unique(gene_names),collapse=',')

        return(unique_gene_names_combined)
    }

    ## Now for each of the probes, let's see what genes are in the same TAD as it:
    ## To do this, we will need to identify the TAD number(s) each probe is in and each gene is in:
    for(TAD_file_index in c(1:length(TAD_file_list))){

        ## Create name for the TAD file overlap columns:
        TAD_file_overlap_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_overlap_number',
            sep=''
        )

        ## Create name for the TAD file gene numbers found in the TAD:
        gene_row_overlap_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_gene_numbers',
            sep=''
        )

        ## Create a name for the total number of genes in a TAD of the probe:
        gene_count_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_gene_count_in_TAD',
            sep=''
        )

        ## Create names for the gene ENSGs and gene names found in TADs:
        gene_ENSG_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_TAD_gene_ENSGs',
            sep=''
        )

        gene_name_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_TAD_gene_names',
            sep=''
        )

        ## Get overlaps for each of the TADs with the probes of interest:
        probe_dataset[TAD_file_overlap_column_name] <- unname(
            mapply(
                TAD_probe_overlapper,
                chromosome= probe_dataset[,2],
                start= probe_dataset[,3],
                stop= probe_dataset[,4],
                MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
            )
        )

        ## Get overlaps for each of the TADs with the genes of interest:
        gencode_v22_genes_df[TAD_file_overlap_column_name] <- unname(
            mapply(
                TAD_gene_overlapper,
                chromosome= gencode_v22_genes_df[,1],
                TSS= gencode_v22_genes_df[,'TSS'],
                MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
            )
        )

        ## Get the numbers of the genes that have overlap TADs with
        ## The probes:
        probe_dataset[gene_row_overlap_column_name] <- unname(sapply(
            probe_dataset[,TAD_file_overlap_column_name],
            TAD_gene_finder,
            target_list=gencode_v22_genes_df[,TAD_file_overlap_column_name]
        ))

        ## Count the number of total genes within a TAD of the probe:
        probe_dataset[gene_count_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_count_in_TAD
        ))

        ## Convert the numbers into gene ENSGs:
        probe_dataset[gene_ENSG_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_name_from_number_lister,
            downregulated_gene_df=gencode_v22_genes_df,
            return_type='gene_id'
        ))

        ## Convert the numbers into gene names:
        probe_dataset[gene_name_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_name_from_number_lister,
            downregulated_gene_df=gencode_v22_genes_df,
            return_type='gene_name'
        ))

        ## If no TAD is found, make a note of that in the gene names/ENSG columns:
        for(i in c(1:nrow(probe_dataset))){

            if(probe_dataset[i,TAD_file_overlap_column_name]==''){

                probe_dataset[i,gene_ENSG_column_name] <- "No_TAD_indentified"
                probe_dataset[i,gene_name_column_name] <- "No_TAD_indentified"
            }

        }

        ## Let's remove unneeded columns now:
        probe_dataset[TAD_file_overlap_column_name] <- NULL
        probe_dataset[gene_row_overlap_column_name] <- NULL

    }

    ## Write the table to file:
    write.table(
        probe_dataset,
        file= paste(
            './hyper.G+.output.TAD/',
            'hyper.G+.TAD_analysis.tsv',
            sep=''
        ),
        quote= FALSE,
        sep='\t'
    )
}

make.TAD.tables.hypermeth.G_neg <- function(prefix, TAD_top_n_genes, i, gene_ENSG, grep, TAD_file_index) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Load the gtf file for CpG locations:
    hg38_450k_CpG_df <- read.delim(
        '../scripts/data/hm450.hg38.manifest.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Probes are 0-indexed, let's fix this:
    hg38_450k_CpG_df$CpG_beg <- hg38_450k_CpG_df$CpG_beg+1

    ## Set the rownames to be the probe IDs:
    rownames(hg38_450k_CpG_df) <- hg38_450k_CpG_df$probeID

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G-.output.histogram/",
            prefix,
            ".hyper.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hyper.G-.output/hyper.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_IDs <- ordered_TFs_by_link_count[
        1:TAD_top_n_genes,
        'geneID'
    ]

    ## Load TAD files to overlap with:
    TAD_file_paths <- list.files(
        path='../external.data/TAD',
        full.names=TRUE
    )

    ## Create an empty list to load each of the files into:
    TAD_file_list <- vector(
        mode = "list",
        length = length(TAD_file_paths)
    )

    ## Remove the .rtf file from being loaded as a TAD file:
    TAD_file_paths <- TAD_file_paths[
        !(sub('.*\\.', '', TAD_file_paths)=='rtf')
    ]

    ## Prepare and load the TAD files:
    ## This assumes the files are bed3 formatted with the first 3 columns:
    ## Being the chromosome, start, and end coordinates:
    ## And are 0-indexed:
    for(i in c(1:length(TAD_file_paths))){

        ## Load the file:
        placeholder_file <- read.delim(
            TAD_file_paths[i],
            header= FALSE,
            stringsAsFactors = FALSE
        )

        ## Remove extra columns and rename the first three columns:
        ## Which should contain chr, start, end locations:
        if(ncol(placeholder_file)>=3){

            placeholder_file <- placeholder_file[1:3]

            colnames(placeholder_file) <- c(
                'seqnames','start','end'
            )

        } else{

            print("File is improperly formatted")
            break()
        }

        ## Add 1 to the starts to make them 1-indexed:
        placeholder_file$start <- placeholder_file$start+1

        ## Add the final processed TAD file to the list:
        TAD_file_list[[i]] <- placeholder_file
    }

    ## Assign names of the files to the TAD file list:
    TAD_file_names <- sub(
        '\\..*',
        '',
        basename(TAD_file_paths)
    )

    names(TAD_file_list) <- TAD_file_names

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes selected:
    probes_linked_to_significant_genes <- unique(
        CpG_linkage_dataset[
            CpG_linkage_dataset$geneID %in% top_gene_IDs,
            'probe'
        ]
    )

    ## Create a data frame with the probes as the first (and only column for now):
    probe_dataset <- data.frame(
        'probe_ID'=sort(probes_linked_to_significant_genes),
        stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info for hg38:
    probe_dataset$seqnames <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_chrm'
    ]

    probe_dataset$start <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_beg'
    ]

    probe_dataset$end <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_end'
    ]

    ## Create names of columns for each of the probes:
    for(gene_ENSG in top_gene_IDs){

        ## Get the gene name corresponding to the ENSG:
        corresponding_gene_name <- ordered_TFs_by_link_count[
            ordered_TFs_by_link_count$geneID==gene_ENSG,
            'geneSymbol'
        ]

        column_name <- paste(
            corresponding_gene_name,
            '_linked',
            sep=''
        )

        ## Get the probes linked to that gene specifically:
        probes_linked_to_gene <- unique(
            CpG_linkage_dataset[
                CpG_linkage_dataset$geneID %in% gene_ENSG,
                'probe'
            ]
        )

        ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
        ## at least one of the top n TFs is linked to this specific one:
        probe_dataset[column_name] <- probe_dataset$probe_ID %in% probes_linked_to_gene

    }

    ## Write functions to identify genes and probes within each TAD:
    TAD_probe_overlapper <- function(chromosome, start, stop, reference_peaks, buffer){
        paste(
            rownames(
                reference_peaks[
                    (
                        (reference_peaks[['seqnames']] == chromosome) &
                            ((reference_peaks[['start']]-buffer) <= stop) &
                            ((reference_peaks[['end']]+buffer) >= start)
                    ),
                ]
            ), collapse=","
        )
    }

    TAD_gene_overlapper <- function(chromosome, TSS, reference_peaks, buffer){
        paste(
            rownames(
                reference_peaks[
                    (
                        (reference_peaks[['seqnames']] == chromosome) &
                            (((reference_peaks[['start']]-buffer) <= TSS) &
                                 ((reference_peaks[['end']]+buffer) >= TSS))
                    ),
                ]
            ), collapse=","
        )
    }

    ## Now let's find genes in the same TAD(s):
    TAD_gene_finder <- function(TAD_numbers_comma_delimited, target_list){

        TAD_numbers <- unlist(strsplit(TAD_numbers_comma_delimited,','))

        gene_numbers <- sapply(
            TAD_numbers,
            grep,
            target_list
        )

        unique_gene_numbers <- unique(unlist(gene_numbers))

        unique_gene_numbers_combined <- paste(unique_gene_numbers,collapse=',')

        return(unique_gene_numbers_combined)
    }

    ## Write a function to count the genes in a TAD with a given probe:
    gene_count_in_TAD <- function(gene_numbers){

        gene_numbers_separated <- as.numeric(unlist(strsplit(gene_numbers,',')))

        return(length(gene_numbers_separated))
    }

    ## Now let's write a function to convert the gene numbers to the gene name
    ## based on setup of siX_downregulated_gencode_v22_genes_df files:
    gene_name_from_number_lister <- function(gene_numbers, downregulated_gene_df, return_type){

        gene_numbers_separated <- as.numeric(unlist(strsplit(gene_numbers,',')))

        if(return_type=='gene_name'){

            gene_names <- downregulated_gene_df[
                gene_numbers_separated,
                'gene_name'
            ]

        } else if(return_type=='gene_id'){

            gene_names <- downregulated_gene_df[
                gene_numbers_separated,
                'gene_id'
            ]
        }

        unique_gene_names_combined <- paste(unique(gene_names),collapse=',')

        return(unique_gene_names_combined)
    }

    ## Now for each of the probes, let's see what genes are in the same TAD as it:
    ## To do this, we will need to identify the TAD number(s) each probe is in and each gene is in:
    for(TAD_file_index in c(1:length(TAD_file_list))){

        ## Create name for the TAD file overlap columns:
        TAD_file_overlap_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_overlap_number',
            sep=''
        )

        ## Create name for the TAD file gene numbers found in the TAD:
        gene_row_overlap_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_gene_numbers',
            sep=''
        )

        ## Create a name for the total number of genes in a TAD of the probe:
        gene_count_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_gene_count_in_TAD',
            sep=''
        )

        ## Create names for the gene ENSGs and gene names found in TADs:
        gene_ENSG_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_TAD_gene_ENSGs',
            sep=''
        )

        gene_name_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_TAD_gene_names',
            sep=''
        )

        ## Get overlaps for each of the TADs with the probes of interest:
        probe_dataset[TAD_file_overlap_column_name] <- unname(
            mapply(
                TAD_probe_overlapper,
                chromosome= probe_dataset[,2],
                start= probe_dataset[,3],
                stop= probe_dataset[,4],
                MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
            )
        )

        ## Get overlaps for each of the TADs with the genes of interest:
        gencode_v22_genes_df[TAD_file_overlap_column_name] <- unname(
            mapply(
                TAD_gene_overlapper,
                chromosome= gencode_v22_genes_df[,1],
                TSS= gencode_v22_genes_df[,'TSS'],
                MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
            )
        )

        ## Get the numbers of the genes that have overlap TADs with
        ## The probes:
        probe_dataset[gene_row_overlap_column_name] <- unname(sapply(
            probe_dataset[,TAD_file_overlap_column_name],
            TAD_gene_finder,
            target_list=gencode_v22_genes_df[,TAD_file_overlap_column_name]
        ))

        ## Count the number of total genes within a TAD of the probe:
        probe_dataset[gene_count_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_count_in_TAD
        ))

        ## Convert the numbers into gene ENSGs:
        probe_dataset[gene_ENSG_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_name_from_number_lister,
            downregulated_gene_df=gencode_v22_genes_df,
            return_type='gene_id'
        ))

        ## Convert the numbers into gene names:
        probe_dataset[gene_name_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_name_from_number_lister,
            downregulated_gene_df=gencode_v22_genes_df,
            return_type='gene_name'
        ))

        ## If no TAD is found, make a note of that in the gene names/ENSG columns:
        for(i in c(1:nrow(probe_dataset))){

            if(probe_dataset[i,TAD_file_overlap_column_name]==''){

                probe_dataset[i,gene_ENSG_column_name] <- "No_TAD_indentified"
                probe_dataset[i,gene_name_column_name] <- "No_TAD_indentified"
            }

        }

        ## Let's remove unneeded columns now:
        probe_dataset[TAD_file_overlap_column_name] <- NULL
        probe_dataset[gene_row_overlap_column_name] <- NULL

    }

    ## Write the table to file:
    write.table(
        probe_dataset,
        file= paste(
            './hyper.G-.output.TAD/',
            'hyper.G-.TAD_analysis.tsv',
            sep=''
        ),
        quote= FALSE,
        sep='\t'
    )
}

make.TAD.tables.hypometh.G_pos <- function(prefix, TAD_top_n_genes, i, gene_ENSG, grep, TAD_file_index) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Load the gtf file for CpG locations:
    hg38_450k_CpG_df <- read.delim(
        '../scripts/data/hm450.hg38.manifest.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Probes are 0-indexed, let's fix this:
    hg38_450k_CpG_df$CpG_beg <- hg38_450k_CpG_df$CpG_beg+1

    ## Set the rownames to be the probe IDs:
    rownames(hg38_450k_CpG_df) <- hg38_450k_CpG_df$probeID

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hypo.G+.output.histogram/",
            prefix,
            ".hypo.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hypo.G+.output/hypo.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_IDs <- ordered_TFs_by_link_count[
        1:TAD_top_n_genes,
        'geneID'
    ]

    ## Load TAD files to overlap with:
    TAD_file_paths <- list.files(
        path='../external.data/TAD',
        full.names=TRUE
    )

    ## Remove the .rtf file from being loaded as a TAD file:
    TAD_file_paths <- TAD_file_paths[
        !(sub('.*\\.', '', TAD_file_paths)=='rtf')
    ]

    ## Create an empty list to load each of the files into:
    TAD_file_list <- vector(
        mode = "list",
        length = length(TAD_file_paths)
    )

    ## Prepare and load the TAD files:
    ## This assumes the files are bed3 formatted with the first 3 columns:
    ## Being the chromosome, start, and end coordinates:
    ## And are 0-indexed:
    for(i in c(1:length(TAD_file_paths))){

        ## Load the file:
        placeholder_file <- read.delim(
            TAD_file_paths[i],
            header= FALSE,
            stringsAsFactors = FALSE
        )

        ## Remove extra columns and rename the first three columns:
        ## Which should contain chr, start, end locations:
        if(ncol(placeholder_file)>=3){

            placeholder_file <- placeholder_file[1:3]

            colnames(placeholder_file) <- c(
                'seqnames','start','end'
            )

        } else{

            print("File is improperly formatted")
            break()
        }

        ## Add 1 to the starts to make them 1-indexed:
        placeholder_file$start <- placeholder_file$start+1

        ## Add the final processed TAD file to the list:
        TAD_file_list[[i]] <- placeholder_file
    }

    ## Assign names of the files to the TAD file list:
    TAD_file_names <- sub(
        '\\..*',
        '',
        basename(TAD_file_paths)
    )

    names(TAD_file_list) <- TAD_file_names

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes selected:
    probes_linked_to_significant_genes <- unique(
        CpG_linkage_dataset[
            CpG_linkage_dataset$geneID %in% top_gene_IDs,
            'probe'
        ]
    )

    ## Create a data frame with the probes as the first (and only column for now):
    probe_dataset <- data.frame(
        'probe_ID'=sort(probes_linked_to_significant_genes),
        stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info for hg38:
    probe_dataset$seqnames <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_chrm'
    ]

    probe_dataset$start <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_beg'
    ]

    probe_dataset$end <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_end'
    ]

    ## Create names of columns for each of the probes:
    for(gene_ENSG in top_gene_IDs){

        ## Get the gene name corresponding to the ENSG:
        corresponding_gene_name <- ordered_TFs_by_link_count[
            ordered_TFs_by_link_count$geneID==gene_ENSG,
            'geneSymbol'
        ]

        column_name <- paste(
            corresponding_gene_name,
            '_linked',
            sep=''
        )

        ## Get the probes linked to that gene specifically:
        probes_linked_to_gene <- unique(
            CpG_linkage_dataset[
                CpG_linkage_dataset$geneID %in% gene_ENSG,
                'probe'
            ]
        )

        ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
        ## at least one of the top n TFs is linked to this specific one:
        probe_dataset[column_name] <- probe_dataset$probe_ID %in% probes_linked_to_gene

    }

    ## Write functions to identify genes and probes within each TAD:
    TAD_probe_overlapper <- function(chromosome, start, stop, reference_peaks, buffer){
        paste(
            rownames(
                reference_peaks[
                    (
                        (reference_peaks[['seqnames']] == chromosome) &
                            ((reference_peaks[['start']]-buffer) <= stop) &
                            ((reference_peaks[['end']]+buffer) >= start)
                    ),
                ]
            ), collapse=","
        )
    }

    TAD_gene_overlapper <- function(chromosome, TSS, reference_peaks, buffer){
        paste(
            rownames(
                reference_peaks[
                    (
                        (reference_peaks[['seqnames']] == chromosome) &
                            (((reference_peaks[['start']]-buffer) <= TSS) &
                                 ((reference_peaks[['end']]+buffer) >= TSS))
                    ),
                ]
            ), collapse=","
        )
    }

    ## Now let's find genes in the same TAD(s):
    TAD_gene_finder <- function(TAD_numbers_comma_delimited, target_list){

        TAD_numbers <- unlist(strsplit(TAD_numbers_comma_delimited,','))

        gene_numbers <- sapply(
            TAD_numbers,
            grep,
            target_list
        )

        unique_gene_numbers <- unique(unlist(gene_numbers))

        unique_gene_numbers_combined <- paste(unique_gene_numbers,collapse=',')

        return(unique_gene_numbers_combined)
    }

    ## Write a function to count the genes in a TAD with a given probe:
    gene_count_in_TAD <- function(gene_numbers){

        gene_numbers_separated <- as.numeric(unlist(strsplit(gene_numbers,',')))

        return(length(gene_numbers_separated))
    }

    ## Now let's write a function to convert the gene numbers to the gene name
    ## based on setup of siX_downregulated_gencode_v22_genes_df files:
    gene_name_from_number_lister <- function(gene_numbers, downregulated_gene_df, return_type){

        gene_numbers_separated <- as.numeric(unlist(strsplit(gene_numbers,',')))

        if(return_type=='gene_name'){

            gene_names <- downregulated_gene_df[
                gene_numbers_separated,
                'gene_name'
            ]

        } else if(return_type=='gene_id'){

            gene_names <- downregulated_gene_df[
                gene_numbers_separated,
                'gene_id'
            ]
        }

        unique_gene_names_combined <- paste(unique(gene_names),collapse=',')

        return(unique_gene_names_combined)
    }

    ## Now for each of the probes, let's see what genes are in the same TAD as it:
    ## To do this, we will need to identify the TAD number(s) each probe is in and each gene is in:
    for(TAD_file_index in c(1:length(TAD_file_list))){

        ## Create name for the TAD file overlap columns:
        TAD_file_overlap_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_overlap_number',
            sep=''
        )

        ## Create name for the TAD file gene numbers found in the TAD:
        gene_row_overlap_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_gene_numbers',
            sep=''
        )

        ## Create a name for the total number of genes in a TAD of the probe:
        gene_count_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_gene_count_in_TAD',
            sep=''
        )

        ## Create names for the gene ENSGs and gene names found in TADs:
        gene_ENSG_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_TAD_gene_ENSGs',
            sep=''
        )

        gene_name_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_TAD_gene_names',
            sep=''
        )

        ## Get overlaps for each of the TADs with the probes of interest:
        probe_dataset[TAD_file_overlap_column_name] <- unname(
            mapply(
                TAD_probe_overlapper,
                chromosome= probe_dataset[,2],
                start= probe_dataset[,3],
                stop= probe_dataset[,4],
                MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
            )
        )

        ## Get overlaps for each of the TADs with the genes of interest:
        gencode_v22_genes_df[TAD_file_overlap_column_name] <- unname(
            mapply(
                TAD_gene_overlapper,
                chromosome= gencode_v22_genes_df[,1],
                TSS= gencode_v22_genes_df[,'TSS'],
                MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
            )
        )

        ## Get the numbers of the genes that have overlap TADs with
        ## The probes:
        probe_dataset[gene_row_overlap_column_name] <- unname(sapply(
            probe_dataset[,TAD_file_overlap_column_name],
            TAD_gene_finder,
            target_list=gencode_v22_genes_df[,TAD_file_overlap_column_name]
        ))

        ## Count the number of total genes within a TAD of the probe:
        probe_dataset[gene_count_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_count_in_TAD
        ))

        ## Convert the numbers into gene ENSGs:
        probe_dataset[gene_ENSG_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_name_from_number_lister,
            downregulated_gene_df=gencode_v22_genes_df,
            return_type='gene_id'
        ))

        ## Convert the numbers into gene names:
        probe_dataset[gene_name_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_name_from_number_lister,
            downregulated_gene_df=gencode_v22_genes_df,
            return_type='gene_name'
        ))

        ## If no TAD is found, make a note of that in the gene names/ENSG columns:
        for(i in c(1:nrow(probe_dataset))){

            if(probe_dataset[i,TAD_file_overlap_column_name]==''){

                probe_dataset[i,gene_ENSG_column_name] <- "No_TAD_indentified"
                probe_dataset[i,gene_name_column_name] <- "No_TAD_indentified"
            }

        }

        ## Let's remove unneeded columns now:
        probe_dataset[TAD_file_overlap_column_name] <- NULL
        probe_dataset[gene_row_overlap_column_name] <- NULL

    }

    ## Write the table to file:
    write.table(
        probe_dataset,
        file= paste(
            './hypo.G+.output.TAD/',
            'hypo.G+.TAD_analysis.tsv',
            sep=''
        ),
        quote= FALSE,
        sep='\t'
    )
}

make.TAD.tables.hypometh.G_neg <- function(prefix, TAD_top_n_genes, i, gene_ENSG, grep, TAD_file_index) {
    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Load the gtf file for CpG locations:
    hg38_450k_CpG_df <- read.delim(
        '../scripts/data/hm450.hg38.manifest.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Probes are 0-indexed, let's fix this:
    hg38_450k_CpG_df$CpG_beg <- hg38_450k_CpG_df$CpG_beg+1

    ## Set the rownames to be the probe IDs:
    rownames(hg38_450k_CpG_df) <- hg38_450k_CpG_df$probeID

    ## Read TF link counts:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hypo.G-.output.histogram/",
            prefix,
            ".hypo.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Load the CpGs linked to each gene:
    CpG_linkage_dataset <- read.delim(
        file= '../step4/hypo.G-.output/hypo.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TFs
    ## based on survival_top_n_genes parameter:
    top_gene_IDs <- ordered_TFs_by_link_count[
        1:TAD_top_n_genes,
        'geneID'
    ]

    ## Load TAD files to overlap with:
    TAD_file_paths <- list.files(
        path='../external.data/TAD',
        full.names=TRUE
    )

    ## Create an empty list to load each of the files into:
    TAD_file_list <- vector(
        mode = "list",
        length = length(TAD_file_paths)
    )

    ## Remove the .rtf file from being loaded as a TAD file:
    TAD_file_paths <- TAD_file_paths[
        !(sub('.*\\.', '', TAD_file_paths)=='rtf')
    ]

    ## Prepare and load the TAD files:
    ## This assumes the files are bed3 formatted with the first 3 columns:
    ## Being the chromosome, start, and end coordinates:
    ## And are 0-indexed:
    for(i in c(1:length(TAD_file_paths))){

        ## Load the file:
        placeholder_file <- read.delim(
            TAD_file_paths[i],
            header= FALSE,
            stringsAsFactors = FALSE
        )

        ## Remove extra columns and rename the first three columns:
        ## Which should contain chr, start, end locations:
        if(ncol(placeholder_file)>=3){

            placeholder_file <- placeholder_file[1:3]

            colnames(placeholder_file) <- c(
                'seqnames','start','end'
            )

        } else{

            print("File is improperly formatted")
            break()
        }

        ## Add 1 to the starts to make them 1-indexed:
        placeholder_file$start <- placeholder_file$start+1

        ## Add the final processed TAD file to the list:
        TAD_file_list[[i]] <- placeholder_file
    }

    ## Assign names of the files to the TAD file list:
    TAD_file_names <- sub(
        '\\..*',
        '',
        basename(TAD_file_paths)
    )

    names(TAD_file_list) <- TAD_file_names

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes selected:
    probes_linked_to_significant_genes <- unique(
        CpG_linkage_dataset[
            CpG_linkage_dataset$geneID %in% top_gene_IDs,
            'probe'
        ]
    )

    ## Create a data frame with the probes as the first (and only column for now):
    probe_dataset <- data.frame(
        'probe_ID'=sort(probes_linked_to_significant_genes),
        stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info for hg38:
    probe_dataset$seqnames <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_chrm'
    ]

    probe_dataset$start <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_beg'
    ]

    probe_dataset$end <- hg38_450k_CpG_df[
        probe_dataset$probe_ID,
        'CpG_end'
    ]

    ## Create names of columns for each of the probes:
    for(gene_ENSG in top_gene_IDs){

        ## Get the gene name corresponding to the ENSG:
        corresponding_gene_name <- ordered_TFs_by_link_count[
            ordered_TFs_by_link_count$geneID==gene_ENSG,
            'geneSymbol'
        ]

        column_name <- paste(
            corresponding_gene_name,
            '_linked',
            sep=''
        )

        ## Get the probes linked to that gene specifically:
        probes_linked_to_gene <- unique(
            CpG_linkage_dataset[
                CpG_linkage_dataset$geneID %in% gene_ENSG,
                'probe'
            ]
        )

        ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
        ## at least one of the top n TFs is linked to this specific one:
        probe_dataset[column_name] <- probe_dataset$probe_ID %in% probes_linked_to_gene

    }

    ## Write functions to identify genes and probes within each TAD:
    TAD_probe_overlapper <- function(chromosome, start, stop, reference_peaks, buffer){
        paste(
            rownames(
                reference_peaks[
                    (
                        (reference_peaks[['seqnames']] == chromosome) &
                            ((reference_peaks[['start']]-buffer) <= stop) &
                            ((reference_peaks[['end']]+buffer) >= start)
                    ),
                ]
            ), collapse=","
        )
    }

    TAD_gene_overlapper <- function(chromosome, TSS, reference_peaks, buffer){
        paste(
            rownames(
                reference_peaks[
                    (
                        (reference_peaks[['seqnames']] == chromosome) &
                            (((reference_peaks[['start']]-buffer) <= TSS) &
                                 ((reference_peaks[['end']]+buffer) >= TSS))
                    ),
                ]
            ), collapse=","
        )
    }

    ## Now let's find genes in the same TAD(s):
    TAD_gene_finder <- function(TAD_numbers_comma_delimited, target_list){

        TAD_numbers <- unlist(strsplit(TAD_numbers_comma_delimited,','))

        gene_numbers <- sapply(
            TAD_numbers,
            grep,
            target_list
        )

        unique_gene_numbers <- unique(unlist(gene_numbers))

        unique_gene_numbers_combined <- paste(unique_gene_numbers,collapse=',')

        return(unique_gene_numbers_combined)
    }

    ## Write a function to count the genes in a TAD with a given probe:
    gene_count_in_TAD <- function(gene_numbers){

        gene_numbers_separated <- as.numeric(unlist(strsplit(gene_numbers,',')))

        return(length(gene_numbers_separated))
    }

    ## Now let's write a function to convert the gene numbers to the gene name
    ## based on setup of siX_downregulated_gencode_v22_genes_df files:
    gene_name_from_number_lister <- function(gene_numbers, downregulated_gene_df, return_type){

        gene_numbers_separated <- as.numeric(unlist(strsplit(gene_numbers,',')))

        if(return_type=='gene_name'){

            gene_names <- downregulated_gene_df[
                gene_numbers_separated,
                'gene_name'
            ]

        } else if(return_type=='gene_id'){

            gene_names <- downregulated_gene_df[
                gene_numbers_separated,
                'gene_id'
            ]
        }

        unique_gene_names_combined <- paste(unique(gene_names),collapse=',')

        return(unique_gene_names_combined)
    }

    ## Now for each of the probes, let's see what genes are in the same TAD as it:
    ## To do this, we will need to identify the TAD number(s) each probe is in and each gene is in:
    for(TAD_file_index in c(1:length(TAD_file_list))){

        ## Create name for the TAD file overlap columns:
        TAD_file_overlap_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_overlap_number',
            sep=''
        )

        ## Create name for the TAD file gene numbers found in the TAD:
        gene_row_overlap_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_gene_numbers',
            sep=''
        )

        ## Create a name for the total number of genes in a TAD of the probe:
        gene_count_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_gene_count_in_TAD',
            sep=''
        )

        ## Create names for the gene ENSGs and gene names found in TADs:
        gene_ENSG_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_TAD_gene_ENSGs',
            sep=''
        )

        gene_name_column_name <- paste(
            names(TAD_file_list[TAD_file_index]),
            '_TAD_gene_names',
            sep=''
        )

        ## Get overlaps for each of the TADs with the probes of interest:
        probe_dataset[TAD_file_overlap_column_name] <- unname(
            mapply(
                TAD_probe_overlapper,
                chromosome= probe_dataset[,2],
                start= probe_dataset[,3],
                stop= probe_dataset[,4],
                MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
            )
        )

        ## Get overlaps for each of the TADs with the genes of interest:
        gencode_v22_genes_df[TAD_file_overlap_column_name] <- unname(
            mapply(
                TAD_gene_overlapper,
                chromosome= gencode_v22_genes_df[,1],
                TSS= gencode_v22_genes_df[,'TSS'],
                MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
            )
        )

        ## Get the numbers of the genes that have overlap TADs with
        ## The probes:
        probe_dataset[gene_row_overlap_column_name] <- unname(sapply(
            probe_dataset[,TAD_file_overlap_column_name],
            TAD_gene_finder,
            target_list=gencode_v22_genes_df[,TAD_file_overlap_column_name]
        ))

        ## Count the number of total genes within a TAD of the probe:
        probe_dataset[gene_count_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_count_in_TAD
        ))

        ## Convert the numbers into gene ENSGs:
        probe_dataset[gene_ENSG_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_name_from_number_lister,
            downregulated_gene_df=gencode_v22_genes_df,
            return_type='gene_id'
        ))

        ## Convert the numbers into gene names:
        probe_dataset[gene_name_column_name] <- unname(sapply(
            probe_dataset[,gene_row_overlap_column_name],
            gene_name_from_number_lister,
            downregulated_gene_df=gencode_v22_genes_df,
            return_type='gene_name'
        ))

        ## If no TAD is found, make a note of that in the gene names/ENSG columns:
        for(i in c(1:nrow(probe_dataset))){

            if(probe_dataset[i,TAD_file_overlap_column_name]==''){

                probe_dataset[i,gene_ENSG_column_name] <- "No_TAD_indentified"
                probe_dataset[i,gene_name_column_name] <- "No_TAD_indentified"
            }

        }

        ## Let's remove unneeded columns now:
        probe_dataset[TAD_file_overlap_column_name] <- NULL
        probe_dataset[gene_row_overlap_column_name] <- NULL

    }

    ## Write the table to file:
    write.table(
        probe_dataset,
        file= paste(
            './hypo.G-.output.TAD/',
            'hypo.G-.TAD_analysis.tsv',
            sep=''
        ),
        quote= FALSE,
        sep='\t'
    )
}

make.ucsc.bed.file.for.hypermeth.G_pos.source <- function(prefix, track_top_n_genes, TR_ENSG) {
    ## Create a directory to place information into:
    dir.create(
        './hyper.G+.output.tracks/'
    )

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Load in the annotations for the 450k array probes on the hg38 genome:
    hg38_450k_probe_info <- read.delim(
        file='../scripts/data/hm450.hg38.manifest.tsv.gz',
        header= TRUE,
        stringsAsFactors = FALSE
    )

    ## Set rownames of hg38_450k_probe_info to be the probe IDs:
    rownames(hg38_450k_probe_info) <- hg38_450k_probe_info$probeID

    ## Read TF link counts from step5 histogram function:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G+.output.histogram/",
            prefix,
            ".hyper.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Read the hypermeth links from step4:
    complete_linked_CpGs_to_TRs <- read.delim(
        file= '../step4/hyper.G+.output/hyper.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TRs
    ## based on survival_top_n_genes parameter:
    top_gene_ids <- ordered_TFs_by_link_count[
        1:track_top_n_genes,
        'geneID'
    ]

    ## Now lets create a dataframe listing the top genes in one column (repeated)
    ## And all the cpgs linked to them in a second column:

    ## Index an empty dataframe:
    intersect_df <- data.frame(
        'TR_ID'=character(),
        'TR_name'=character(),
        'CpG_ID'=character(),
        stringsAsFactors = FALSE
    )

    ## Add the CpGs linked to each TR to the dataframe:
    for(TR_ENSG in top_gene_ids){

        ## Get the name of the gene of interest from
        ## the ENSG:
        gene_name <- gencode_v22_genes_df[
            TR_ENSG,
            'gene_name'
        ]

        ## Get a list of all the CpGs linked to the given TR:
        linked_CpGs <- unique(
            complete_linked_CpGs_to_TRs[
                complete_linked_CpGs_to_TRs$geneID==TR_ENSG,
                'probe'
            ]
        )

        ## Create temporary dataframe with the 3 vectors of info:
        temp_df <- data.frame(
            'TR_ID'= rep(
                TR_ENSG,
                length(linked_CpGs)
            ),
            'TR_name'= rep(
                gene_name,
                length(linked_CpGs)
            ),
            'CpG_ID'= linked_CpGs,
            stringsAsFactors = FALSE
        )

        ## Bind the temp_df to the intersect_df:
        intersect_df <- rbind(
            intersect_df,
            temp_df
        )

    }

    ## For each of the genes in the new intersect df
    ## get the chromosome, "start", and "end" of the gene:
    ## Start and end will both be the TSS
    intersect_df$TR_chr <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'seqnames'
    ]

    intersect_df$TR_start <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'TSS'
    ]

    intersect_df$TR_end <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'TSS'
    ]

    ## For each of the probes in the new intersect df
    ## get the chromosome, "start", and "end" of the probe:
    intersect_df$CpG_chr <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_chrm'
    ]

    intersect_df$CpG_start <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_beg'
    ]

    intersect_df$CpG_end <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_end'
    ]

    ## Now let's use the jetcolor function to set up a gradient of colors equal to
    ## the number of TRs analyzed and assign each TR a color:
    jet_color_numeric_color_grad <- rainbow(track_top_n_genes)

    names(jet_color_numeric_color_grad) <- top_gene_ids

    ### Now that we have the information, let's assemble the real intersect file:
    output_intersect_file <- data.frame(
        'chrom'= intersect_df$TR_chr,
        'chromStart'= ifelse(
            intersect_df$TR_chr==intersect_df$CpG_chr,
            ifelse(
                (intersect_df$CpG_start < intersect_df$TR_start),
                intersect_df$CpG_start,
                intersect_df$TR_start
            ),
            intersect_df$TR_start-1
        ),
        'chromEnd'= ifelse(
            intersect_df$TR_chr==intersect_df$CpG_chr,
            ifelse(
                (intersect_df$CpG_start < intersect_df$TR_start),
                intersect_df$TR_start-1,
                intersect_df$CpG_start
            ),
            intersect_df$TR_end
        ),
        'name'= paste(
            intersect_df$TR_ID,
            intersect_df$CpG_ID,
            'link',
            sep='_'
        ),
        'score'= rep(
            0,
            nrow(intersect_df)
        ),
        'value'= rep(
            0,
            nrow(intersect_df)
        ),
        'exp'= rep(
            '.',
            nrow(intersect_df)
        ),
        'color'= substring(
            jet_color_numeric_color_grad[
                intersect_df$TR_ID
            ],
            1,
            7
        ),
        'sourceChrom'= intersect_df$TR_chr,
        'sourceStart'= intersect_df$TR_start-1,
        'sourceEnd'= intersect_df$TR_end,
        'sourceName'= intersect_df$TR_name,
        'sourceStrand'= rep(
            '.',
            nrow(intersect_df)
        ),
        'targetChrom'= intersect_df$CpG_chr,
        'targetStart'= intersect_df$CpG_start,
        'targetEnd'= intersect_df$CpG_end,
        'targetName'= intersect_df$CpG_ID,
        'targetStrand'= rep(
            '.',
            nrow(intersect_df)
        ),
        stringsAsFactors = FALSE
    )

    ## Create text for the header line:
    header_line_text <- "track type=interact name=\"TENET2.0_hyper.G+_interactions\" description=\"TENET2.0 top TR to enhancer DNA methylation probe links\""

    ## Create a file name for the output bed file:
    bed_file_name <- paste(
        './hyper.G+.output.tracks/',
        'TR.enhancer.probe.links.hg38.bed',
        sep=''
    )

    ## Add the header line to the new bed file:
    cat(header_line_text, "\n", file = bed_file_name)

    ## Write the info to the file:
    write.table(output_intersect_file, file = bed_file_name, append = TRUE, row.names = FALSE, col.names = FALSE, quote= FALSE)
}

make.ucsc.bed.file.for.hypermeth.G_neg.source <- function(prefix, track_top_n_genes, TR_ENSG) {
    ## Create a directory to place information into:
    dir.create(
        './hyper.G-.output.tracks/'
    )

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Load in the annotations for the 450k array probes on the hg38 genome:
    hg38_450k_probe_info <- read.delim(
        file='../scripts/data/hm450.hg38.manifest.tsv.gz',
        header= TRUE,
        stringsAsFactors = FALSE
    )

    ## Set rownames of hg38_450k_probe_info to be the probe IDs:
    rownames(hg38_450k_probe_info) <- hg38_450k_probe_info$probeID

    ## Read TF link counts from step5 histogram function:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hyper.G-.output.histogram/",
            prefix,
            ".hyper.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Read the hypermeth links from step4:
    complete_linked_CpGs_to_TRs <- read.delim(
        file= '../step4/hyper.G-.output/hyper.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TRs
    ## based on survival_top_n_genes parameter:
    top_gene_ids <- ordered_TFs_by_link_count[
        1:track_top_n_genes,
        'geneID'
    ]

    ## Now lets create a dataframe listing the top genes in one column (repeated)
    ## And all the cpgs linked to them in a second column:

    ## Index an empty dataframe:
    intersect_df <- data.frame(
        'TR_ID'=character(),
        'TR_name'=character(),
        'CpG_ID'=character(),
        stringsAsFactors = FALSE
    )

    ## Add the CpGs linked to each TR to the dataframe:
    for(TR_ENSG in top_gene_ids){

        ## Get the name of the gene of interest from
        ## the ENSG:
        gene_name <- gencode_v22_genes_df[
            TR_ENSG,
            'gene_name'
        ]

        ## Get a list of all the CpGs linked to the given TR:
        linked_CpGs <- unique(
            complete_linked_CpGs_to_TRs[
                complete_linked_CpGs_to_TRs$geneID==TR_ENSG,
                'probe'
            ]
        )

        ## Create temporary dataframe with the 3 vectors of info:
        temp_df <- data.frame(
            'TR_ID'= rep(
                TR_ENSG,
                length(linked_CpGs)
            ),
            'TR_name'= rep(
                gene_name,
                length(linked_CpGs)
            ),
            'CpG_ID'= linked_CpGs,
            stringsAsFactors = FALSE
        )

        ## Bind the temp_df to the intersect_df:
        intersect_df <- rbind(
            intersect_df,
            temp_df
        )

    }

    ## For each of the genes in the new intersect df
    ## get the chromosome, "start", and "end" of the gene:
    ## Start and end will both be the TSS
    intersect_df$TR_chr <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'seqnames'
    ]

    intersect_df$TR_start <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'TSS'
    ]

    intersect_df$TR_end <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'TSS'
    ]

    ## For each of the probes in the new intersect df
    ## get the chromosome, "start", and "end" of the probe:
    intersect_df$CpG_chr <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_chrm'
    ]

    intersect_df$CpG_start <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_beg'
    ]

    intersect_df$CpG_end <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_end'
    ]

    ## Now let's use the jetcolor function to set up a gradient of colors equal to
    ## the number of TRs analyzed and assign each TR a color:
    jet_color_numeric_color_grad <- rainbow(track_top_n_genes)

    names(jet_color_numeric_color_grad) <- top_gene_ids

    ### Now that we have the information, let's assemble the real intersect file:
    output_intersect_file <- data.frame(
        'chrom'= intersect_df$TR_chr,
        'chromStart'= ifelse(
            intersect_df$TR_chr==intersect_df$CpG_chr,
            ifelse(
                (intersect_df$CpG_start < intersect_df$TR_start),
                intersect_df$CpG_start,
                intersect_df$TR_start
            ),
            intersect_df$TR_start-1
        ),
        'chromEnd'= ifelse(
            intersect_df$TR_chr==intersect_df$CpG_chr,
            ifelse(
                (intersect_df$CpG_start < intersect_df$TR_start),
                intersect_df$TR_start-1,
                intersect_df$CpG_start
            ),
            intersect_df$TR_end
        ),
        'name'= paste(
            intersect_df$TR_ID,
            intersect_df$CpG_ID,
            'link',
            sep='_'
        ),
        'score'= rep(
            0,
            nrow(intersect_df)
        ),
        'value'= rep(
            0,
            nrow(intersect_df)
        ),
        'exp'= rep(
            '.',
            nrow(intersect_df)
        ),
        'color'= substring(
            jet_color_numeric_color_grad[
                intersect_df$TR_ID
            ],
            1,
            7
        ),
        'sourceChrom'= intersect_df$TR_chr,
        'sourceStart'= intersect_df$TR_start-1,
        'sourceEnd'= intersect_df$TR_end,
        'sourceName'= intersect_df$TR_name,
        'sourceStrand'= rep(
            '.',
            nrow(intersect_df)
        ),
        'targetChrom'= intersect_df$CpG_chr,
        'targetStart'= intersect_df$CpG_start,
        'targetEnd'= intersect_df$CpG_end,
        'targetName'= intersect_df$CpG_ID,
        'targetStrand'= rep(
            '.',
            nrow(intersect_df)
        ),
        stringsAsFactors = FALSE
    )

    ## Create text for the header line:
    header_line_text <- "track type=interact name=\"TENET2.0_hyper.G-_interactions\" description=\"TENET2.0 top TR to enhancer DNA methylation probe links\""

    ## Create a file name for the output bed file:
    bed_file_name <- paste(
        './hyper.G-.output.tracks/',
        'TR.enhancer.probe.links.hg38.bed',
        sep=''
    )

    ## Add the header line to the new bed file:
    cat(header_line_text, "\n", file = bed_file_name)

    ## Write the info to the file:
    write.table(output_intersect_file, file = bed_file_name, append = TRUE, row.names = FALSE, col.names = FALSE, quote= FALSE)
}

make.ucsc.bed.file.for.hypometh.G_pos.source <- function(prefix, track_top_n_genes, TR_ENSG) {
    ## Create a directory to place information into:
    dir.create(
        './hypo.G+.output.tracks/'
    )

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Load in the annotations for the 450k array probes on the hg38 genome:
    hg38_450k_probe_info <- read.delim(
        file='../scripts/data/hm450.hg38.manifest.tsv.gz',
        header= TRUE,
        stringsAsFactors = FALSE
    )

    ## Set rownames of hg38_450k_probe_info to be the probe IDs:
    rownames(hg38_450k_probe_info) <- hg38_450k_probe_info$probeID

    ## Read TF link counts from step5 histogram function:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hypo.G+.output.histogram/",
            prefix,
            ".hypo.G+.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Read the hypometh links from step4:
    complete_linked_CpGs_to_TRs <- read.delim(
        file= '../step4/hypo.G+.output/hypo.G+.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TRs
    ## based on survival_top_n_genes parameter:
    top_gene_ids <- ordered_TFs_by_link_count[
        1:track_top_n_genes,
        'geneID'
    ]

    ## Now lets create a dataframe listing the top genes in one column (repeated)
    ## And all the cpgs linked to them in a second column:

    ## Index an empty dataframe:
    intersect_df <- data.frame(
        'TR_ID'=character(),
        'TR_name'=character(),
        'CpG_ID'=character(),
        stringsAsFactors = FALSE
    )

    ## Add the CpGs linked to each TR to the dataframe:
    for(TR_ENSG in top_gene_ids){

        ## Get the name of the gene of interest from
        ## the ENSG:
        gene_name <- gencode_v22_genes_df[
            TR_ENSG,
            'gene_name'
        ]

        ## Get a list of all the CpGs linked to the given TR:
        linked_CpGs <- unique(
            complete_linked_CpGs_to_TRs[
                complete_linked_CpGs_to_TRs$geneID==TR_ENSG,
                'probe'
            ]
        )

        ## Create temporary dataframe with the 3 vectors of info:
        temp_df <- data.frame(
            'TR_ID'= rep(
                TR_ENSG,
                length(linked_CpGs)
            ),
            'TR_name'= rep(
                gene_name,
                length(linked_CpGs)
            ),
            'CpG_ID'= linked_CpGs,
            stringsAsFactors = FALSE
        )

        ## Bind the temp_df to the intersect_df:
        intersect_df <- rbind(
            intersect_df,
            temp_df
        )

    }

    ## For each of the genes in the new intersect df
    ## get the chromosome, "start", and "end" of the gene:
    ## Start and end will both be the TSS
    intersect_df$TR_chr <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'seqnames'
    ]

    intersect_df$TR_start <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'TSS'
    ]

    intersect_df$TR_end <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'TSS'
    ]

    ## For each of the probes in the new intersect df
    ## get the chromosome, "start", and "end" of the probe:
    intersect_df$CpG_chr <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_chrm'
    ]

    intersect_df$CpG_start <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_beg'
    ]

    intersect_df$CpG_end <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_end'
    ]

    ## Now let's use the jetcolor function to set up a gradient of colors equal to
    ## the number of TRs analyzed and assign each TR a color:
    jet_color_numeric_color_grad <- rainbow(track_top_n_genes)

    names(jet_color_numeric_color_grad) <- top_gene_ids

    ### Now that we have the information, let's assemble the real intersect file:
    output_intersect_file <- data.frame(
        'chrom'= intersect_df$TR_chr,
        'chromStart'= ifelse(
            intersect_df$TR_chr==intersect_df$CpG_chr,
            ifelse(
                (intersect_df$CpG_start < intersect_df$TR_start),
                intersect_df$CpG_start,
                intersect_df$TR_start
            ),
            intersect_df$TR_start-1
        ),
        'chromEnd'= ifelse(
            intersect_df$TR_chr==intersect_df$CpG_chr,
            ifelse(
                (intersect_df$CpG_start < intersect_df$TR_start),
                intersect_df$TR_start-1,
                intersect_df$CpG_start
            ),
            intersect_df$TR_end
        ),
        'name'= paste(
            intersect_df$TR_ID,
            intersect_df$CpG_ID,
            'link',
            sep='_'
        ),
        'score'= rep(
            0,
            nrow(intersect_df)
        ),
        'value'= rep(
            0,
            nrow(intersect_df)
        ),
        'exp'= rep(
            '.',
            nrow(intersect_df)
        ),
        'color'= substring(
            jet_color_numeric_color_grad[
                intersect_df$TR_ID
            ],
            1,
            7
        ),
        'sourceChrom'= intersect_df$TR_chr,
        'sourceStart'= intersect_df$TR_start-1,
        'sourceEnd'= intersect_df$TR_end,
        'sourceName'= intersect_df$TR_name,
        'sourceStrand'= rep(
            '.',
            nrow(intersect_df)
        ),
        'targetChrom'= intersect_df$CpG_chr,
        'targetStart'= intersect_df$CpG_start,
        'targetEnd'= intersect_df$CpG_end,
        'targetName'= intersect_df$CpG_ID,
        'targetStrand'= rep(
            '.',
            nrow(intersect_df)
        ),
        stringsAsFactors = FALSE
    )

    ## Create text for the header line:
    header_line_text <- "track type=interact name=\"TENET2.0_Hypo.G+_interactions\" description=\"TENET2.0 top TR to enhancer DNA methylation probe links\""

    ## Create a file name for the output bed file:
    bed_file_name <- paste(
        './hypo.G+.output.tracks/',
        'TR.enhancer.probe.links.hg38.bed',
        sep=''
    )

    ## Add the header line to the new bed file:
    cat(header_line_text, "\n", file = bed_file_name)

    ## Write the info to the file:
    write.table(output_intersect_file, file = bed_file_name, append = TRUE, row.names = FALSE, col.names = FALSE, quote= FALSE)
}

make.ucsc.bed.file.for.hypometh.G_neg.source <- function(prefix, track_top_n_genes, TR_ENSG) {
    ## Create a directory to place information into:
    dir.create(
        './hypo.G-.output.tracks/'
    )

    ## Load settings information from .rda file:
    load("../settings.rda")

    ## Read in gene info from gtf file:
    gencode_v22_genes <- read.delim(
        '../scripts/data/gencode_v22_gtf_file.tsv.gz',
        stringsAsFactors = FALSE
    )

    ## Subset gene info to only the whole genes themselves:
    gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]

    ## Subset the gene IDs to remove the stuff after the periods (and the periods)
    gencode_v22_genes_df$gene_id <- sub(
        '\\..*',
        '',
        gencode_v22_genes_df$gene_id
    )

    ## Set the rownames equal to the genes ENSG:
    rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id

    ## Identify the TSS of genes based on strand:
    gencode_v22_genes_df$TSS <- ifelse(
        gencode_v22_genes_df$strand=='-',
        gencode_v22_genes_df$end,
        gencode_v22_genes_df$start
    )

    ## Remove the original gtf object and associated data frame:
    rm(gencode_v22_genes)

    ## Load in the annotations for the 450k array probes on the hg38 genome:
    hg38_450k_probe_info <- read.delim(
        file='../scripts/data/hm450.hg38.manifest.tsv.gz',
        header= TRUE,
        stringsAsFactors = FALSE
    )

    ## Set rownames of hg38_450k_probe_info to be the probe IDs:
    rownames(hg38_450k_probe_info) <- hg38_450k_probe_info$probeID

    ## Read TF link counts from step5 histogram function:
    ordered_TFs_by_link_count <- read.delim(
        file= paste(
            "./hypo.G-.output.histogram/",
            prefix,
            ".hypo.G-.links.all.tf.freq.txt",
            sep=''
        ),
        stringsAsFactors = FALSE
    )

    ## Read the hypometh links from step4:
    complete_linked_CpGs_to_TRs <- read.delim(
        file= '../step4/hypo.G-.output/hypo.G-.link.zscore.perm.all.optimized.links.txt',
        stringsAsFactors = FALSE
    )

    ## Get the name of the top n TRs
    ## based on survival_top_n_genes parameter:
    top_gene_ids <- ordered_TFs_by_link_count[
        1:track_top_n_genes,
        'geneID'
    ]

    ## Now lets create a dataframe listing the top genes in one column (repeated)
    ## And all the cpgs linked to them in a second column:

    ## Index an empty dataframe:
    intersect_df <- data.frame(
        'TR_ID'=character(),
        'TR_name'=character(),
        'CpG_ID'=character(),
        stringsAsFactors = FALSE
    )

    ## Add the CpGs linked to each TR to the dataframe:
    for(TR_ENSG in top_gene_ids){

        ## Get the name of the gene of interest from
        ## the ENSG:
        gene_name <- gencode_v22_genes_df[
            TR_ENSG,
            'gene_name'
        ]

        ## Get a list of all the CpGs linked to the given TR:
        linked_CpGs <- unique(
            complete_linked_CpGs_to_TRs[
                complete_linked_CpGs_to_TRs$geneID==TR_ENSG,
                'probe'
            ]
        )

        ## Create temporary dataframe with the 3 vectors of info:
        temp_df <- data.frame(
            'TR_ID'= rep(
                TR_ENSG,
                length(linked_CpGs)
            ),
            'TR_name'= rep(
                gene_name,
                length(linked_CpGs)
            ),
            'CpG_ID'= linked_CpGs,
            stringsAsFactors = FALSE
        )

        ## Bind the temp_df to the intersect_df:
        intersect_df <- rbind(
            intersect_df,
            temp_df
        )

    }

    ## For each of the genes in the new intersect df
    ## get the chromosome, "start", and "end" of the gene:
    ## Start and end will both be the TSS
    intersect_df$TR_chr <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'seqnames'
    ]

    intersect_df$TR_start <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'TSS'
    ]

    intersect_df$TR_end <- gencode_v22_genes_df[
        intersect_df$TR_ID,
        'TSS'
    ]

    ## For each of the probes in the new intersect df
    ## get the chromosome, "start", and "end" of the probe:
    intersect_df$CpG_chr <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_chrm'
    ]

    intersect_df$CpG_start <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_beg'
    ]

    intersect_df$CpG_end <- hg38_450k_probe_info[
        intersect_df$CpG_ID,
        'CpG_end'
    ]

    ## Now let's use the jetcolor function to set up a gradient of colors equal to
    ## the number of TRs analyzed and assign each TR a color:
    jet_color_numeric_color_grad <- rainbow(track_top_n_genes)

    names(jet_color_numeric_color_grad) <- top_gene_ids

    ### Now that we have the information, let's assemble the real intersect file:
    output_intersect_file <- data.frame(
        'chrom'= intersect_df$TR_chr,
        'chromStart'= ifelse(
            intersect_df$TR_chr==intersect_df$CpG_chr,
            ifelse(
                (intersect_df$CpG_start < intersect_df$TR_start),
                intersect_df$CpG_start,
                intersect_df$TR_start
            ),
            intersect_df$TR_start-1
        ),
        'chromEnd'= ifelse(
            intersect_df$TR_chr==intersect_df$CpG_chr,
            ifelse(
                (intersect_df$CpG_start < intersect_df$TR_start),
                intersect_df$TR_start-1,
                intersect_df$CpG_start
            ),
            intersect_df$TR_end
        ),
        'name'= paste(
            intersect_df$TR_ID,
            intersect_df$CpG_ID,
            'link',
            sep='_'
        ),
        'score'= rep(
            0,
            nrow(intersect_df)
        ),
        'value'= rep(
            0,
            nrow(intersect_df)
        ),
        'exp'= rep(
            '.',
            nrow(intersect_df)
        ),
        'color'= substring(
            jet_color_numeric_color_grad[
                intersect_df$TR_ID
            ],
            1,
            7
        ),
        'sourceChrom'= intersect_df$TR_chr,
        'sourceStart'= intersect_df$TR_start-1,
        'sourceEnd'= intersect_df$TR_end,
        'sourceName'= intersect_df$TR_name,
        'sourceStrand'= rep(
            '.',
            nrow(intersect_df)
        ),
        'targetChrom'= intersect_df$CpG_chr,
        'targetStart'= intersect_df$CpG_start,
        'targetEnd'= intersect_df$CpG_end,
        'targetName'= intersect_df$CpG_ID,
        'targetStrand'= rep(
            '.',
            nrow(intersect_df)
        ),
        stringsAsFactors = FALSE
    )

    ## Create text for the header line:
    header_line_text <- "track type=interact name=\"TENET2.0_hypo.G-_interactions\" description=\"TENET2.0 top TR to enhancer DNA methylation probe links\""

    ## Create a file name for the output bed file:
    bed_file_name <- paste(
        './hypo.G-.output.tracks/',
        'TR.enhancer.probe.links.hg38.bed',
        sep=''
    )

    ## Add the header line to the new bed file:
    cat(header_line_text, "\n", file = bed_file_name)

    ## Write the info to the file:
    write.table(output_intersect_file, file = bed_file_name, append = TRUE, row.names = FALSE, col.names = FALSE, quote= FALSE)
}
