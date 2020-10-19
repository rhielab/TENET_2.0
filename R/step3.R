getPosEpval=function(probe, geneID, Z.real){
    LS=list.files("../", pattern=paste("^", as.character(geneID), ".hypo", sep=""))
    TESTSR=read.delim(paste("../",LS[1], sep=""), header=F)
    colnames(TESTSR)=c("probe", "geneSymbol", "geneID", "Z.real")
    P=TESTSR[which(TESTSR$geneID==as.character(geneID)),]
    P2=P[order(as.numeric(P$Z.real), decreasing=F),]
    P2$rank=c(1:dim(TESTSR)[1])
    epval=as.numeric(P2[which(P2$probe==as.character(probe)),5])/as.numeric(length(as.numeric(Z.real)[which(is.finite(as.numeric(P2$Z.real))==TRUE)]))
    print(epval)
}

getNegEpval=function(probe, geneID, Z.real){
    LS=list.files("../", pattern=paste("^", as.character(geneID), ".hyper", sep=""))
    TESTSR=read.delim(paste("../",LS[1], sep=""), header=F)
    colnames(TESTSR)=c("probe", "geneSymbol", "geneID", "Z.real")
    P=TESTSR[which(TESTSR$geneID==as.character(geneID)),]
    P2=P[order(as.numeric(P$Z.real), decreasing=T),]
    P2$rank=c(1:dim(TESTSR)[1])
    epval=as.numeric(P2[which(P2$probe==as.character(probe)),5])/as.numeric(length(as.numeric(Z.real)[which(is.finite(as.numeric(P2$Z.real))==TRUE)]))
    print(epval)
}

get.permutation.z.score.for.hyperG_neg.enhancers.bringfiles.realfast.from.step2 <- function(parallel, cores, permutation.cutoff) {
    load("../../settings.rda")
    # library(parallel)
    ## let's do for hyper.G- #
    PL=read.delim("../../step2/hyper.output/hyper.G-.link.zscore.txt", header=F)
    colnames(PL)=c("probe", "geneSymbol", "geneID", "Z.real")
    dir.create("hyper.G-.output")
    setwd("./hyper.G-.output")
    top=list(probe=as.character(PL$probe), geneID=as.character(PL$geneID), Z.real=as.numeric(PL$Z.real))
    pvals=mcmapply(getNegEpval, top$probe, top$geneID, top$Z.real, mc.cores=cores)
    PL$pval=pvals
    selected=PL[which(PL$pval<permutation.cutoff),]
    if (dim(selected)[1]<1){
        print("did you run steps correctly? no significant hyper.G- links are found from step3")
    }
    write.table(selected, file=paste("hyper.G-.link.zscore.perm.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.permutation.z.score.for.hyperG_pos.enhancers.bringfiles.realfast.from.step2 <- function(parallel, cores, permutation.cutoff) {
    load("../../settings.rda")
    # library(parallel)
    # let's do for hyper.G+ #
    PL=read.delim("../../step2/hyper.output/hyper.G+.link.zscore.txt", header=F)
    colnames(PL)=c("probe", "geneSymbol", "geneID", "Z.real")
    dir.create("hyper.G+.output")
    setwd("./hyper.G+.output")
    top=list(probe=as.character(PL$probe), geneID=as.character(PL$geneID), Z.real=as.numeric(PL$Z.real))
    pvals=mcmapply(getPosEpval, top$probe, top$geneID, top$Z.real, mc.cores=cores)
    PL$pval=pvals
    selected=PL[which(PL$pval<permutation.cutoff),]
    if (dim(selected)[1]<1){
        print("did you run steps correctly? no significant hyper.G+ links are found from step3")
    }
    write.table(selected, file=paste("hyper.G+.link.zscore.perm.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.permutation.z.score.for.hypermeth.enhancers.bringfiles.realfast.from.step2 <- function(parallel, cores, permutation.cutoff) {
    load("../../settings.rda")
    # library(parallel)
    # let's do for hyper.G+ #
    PL=read.delim("../../step2/hyper.output/hyper.G+.link.zscore.txt", header=F)
    colnames(PL)=c("probe", "geneSymbol", "geneID", "Z.real")
    dir.create("hyper.G+.output")
    setwd("./hyper.G+.output")
    top=list(probe=as.character(PL$probe), geneID=as.character(PL$geneID), Z.real=as.numeric(PL$Z.real))
    pvals=mcmapply(getPosEpval, top$probe, top$geneID, top$Z.real, mc.cores=cores)
    PL$pval=pvals
    selected=PL[which(PL$pval<permutation.cutoff),]
    if (dim(selected)[1]<1){
        print("did you run steps correctly? no significant hyper.G+ links are found from step3")
    }
    write.table(selected, file=paste("hyper.G+.link.zscore.perm.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
    ## let's do for hyper.G- #
    PL=read.delim("../../step2/hyper.output/hyper.G-.link.zscore.txt", header=F)
    colnames(PL)=c("probe", "geneSymbol", "geneID", "Z.real")
    dir.create("hyper.G-.output")
    setwd("./hyper.G-.output")
    top=list(probe=as.character(PL$probe), geneID=as.character(PL$geneID), Z.real=as.numeric(PL$Z.real))
    getNegEpval=function(probe, geneID, Z.real){
        LS=list.files("../", pattern=paste("^", as.character(geneID), ".hyper", sep=""))
        TESTSR=read.delim(paste("../",LS[1], sep=""), header=F)
        colnames(TESTSR)=c("probe", "geneSymbol", "geneID", "Z.real")
        P=TESTSR[which(TESTSR$geneID==as.character(geneID)),]
        P2=P[order(as.numeric(P$Z.real), decreasing=T),]
        P2$rank=c(1:dim(TESTSR)[1])
        epval=as.numeric(P2[which(P2$probe==as.character(probe)),5])/as.numeric(length(as.numeric(Z.real)[which(is.finite(as.numeric(P2$Z.real))==TRUE)]))
        print(epval)
    }
    pvals=mcmapply(getNegEpval, top$probe, top$geneID, top$Z.real, mc.cores=cores)
    PL$pval=pvals
    selected=PL[which(PL$pval<permutation.cutoff),]
    if (dim(selected)[1]<1){
        print("did you run steps correctly? no significant hyper.G- links are found from step3")
    }
    write.table(selected, file=paste("hyper.G-.link.zscore.perm.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.permutation.z.score.for.hypoG_pos.enhancers.bringfiles.realfast.from.step2 <- function(parallel, cores, permutation.cutoff) {
    load("../../settings.rda")
    # library(parallel)
    # let's do for hypo.G+ #
    PL=read.delim("../../step2/hypo.output/hypo.G+.link.zscore.txt", header=F)
    colnames(PL)=c("probe", "geneSymbol", "geneID", "Z.real")
    dir.create("hypo.G+.output")
    setwd("./hypo.G+.output")
    top=list(probe=as.character(PL$probe), geneID=as.character(PL$geneID), Z.real=as.numeric(PL$Z.real))
    getPosEpval=function(probe, geneID, Z.real){
        LS=list.files("../", pattern=paste("^", as.character(geneID), ".hypo", sep=""))
        TESTSR=read.delim(paste("../",LS[1], sep=""), header=F)
        colnames(TESTSR)=c("probe", "geneSymbol", "geneID", "Z.real")
        P=TESTSR[which(TESTSR$geneID==as.character(geneID)),]
        P2=P[order(as.numeric(P$Z.real), decreasing=F),]
        P2$rank=c(1:dim(TESTSR)[1])
        epval=as.numeric(P2[which(P2$probe==as.character(probe)),5])/as.numeric(length(as.numeric(Z.real)[which(is.finite(as.numeric(P2$Z.real))==TRUE)]))
        print(epval)
    }
    pvals=mcmapply(getPosEpval, top$probe, top$geneID, top$Z.real, mc.cores=cores)
    PL$pval=pvals
    selected=PL[which(PL$pval<permutation.cutoff),]
    if (dim(selected)[1]<1){
        print("did you run steps correctly? no significant hypo.G+ links are found from step3")
    }
    write.table(selected, file=paste("hypo.G+.link.zscore.perm.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.permutation.z.score.for.hypoG_neg.enhancers.bringfiles.realfast.from.step2 <- function(parallel, cores, permutation.cutoff) {
    load("../../settings.rda")
    # library(parallel)
    ## let's do for hypo.G- #
    PL=read.delim("../../step2/hypo.output/hypo.G-.link.zscore.txt", header=F)
    colnames(PL)=c("probe", "geneSymbol", "geneID", "Z.real")
    dir.create("hypo.G-.output")
    setwd("./hypo.G-.output")
    top=list(probe=as.character(PL$probe), geneID=as.character(PL$geneID), Z.real=as.numeric(PL$Z.real))
    pvals=mcmapply(getNegEpval, top$probe, top$geneID, top$Z.real, mc.cores=cores)
    PL$pval=pvals
    selected=PL[which(PL$pval<permutation.cutoff),]
    if (dim(selected)[1]<1){
        print("did you run steps correctly? no significant hypo.G- links are found from step3")
    }
    write.table(selected, file=paste("hypo.G-.link.zscore.perm.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

get.permutation.z.score.for.hypometh.enhancers.bringfiles.realfast.from.step2 <- function(parallel, cores, permutation.cutoff) {
    load("../../settings.rda")
    # library(parallel)
    # let's do for hypo.G+ #
    PL=read.delim("../../step2/hypo.output/hypo.G+.link.zscore.txt", header=F)
    colnames(PL)=c("probe", "geneSymbol", "geneID", "Z.real")
    dir.create("hypo.G+.output")
    setwd("./hypo.G+.output")
    top=list(probe=as.character(PL$probe), geneID=as.character(PL$geneID), Z.real=as.numeric(PL$Z.real))
    pvals=mcmapply(getPosEpval, top$probe, top$geneID, top$Z.real, mc.cores=cores)
    PL$pval=pvals
    selected=PL[which(PL$pval<permutation.cutoff),]
    if (dim(selected)[1]<1){
        print("did you run steps correctly? no significant hypo.G+ links are found from step3")
    }
    write.table(selected, file=paste("hypo.G+.link.zscore.perm.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
    ## let's do for hypo.G- #
    PL=read.delim("../../step2/hypo.output/hypo.G-.link.zscore.txt", header=F)
    colnames(PL)=c("probe", "geneSymbol", "geneID", "Z.real")
    dir.create("hypo.G-.output")
    setwd("./hypo.G-.output")
    top=list(probe=as.character(PL$probe), geneID=as.character(PL$geneID), Z.real=as.numeric(PL$Z.real))
    pvals=mcmapply(getNegEpval, top$probe, top$geneID, top$Z.real, mc.cores=cores)
    PL$pval=pvals
    selected=PL[which(PL$pval<permutation.cutoff),]
    if (dim(selected)[1]<1){
        print("did you run steps correctly? no significant hypo.G- links are found from step3")
    }
    write.table(selected, file=paste("hypo.G-.link.zscore.perm.txt", sep="."), row.names=F, col.names=T, quote=F, sep="\t")
    setwd("../")
}

TENET_unique_geneID_subsetting_script_hyper <- function() {
    ## Listing all the previously created *.dataframe.fin.txt files:
    dataframe.fin.txt <- list.files(
        pattern='.dataframe.fin.txt'
    )

    ## Obtaining the list of genes (by ENSG)
    ## which have previously created .dataframe.fin.txt files:
    created_ENSG_dataframes <- substring(
        dataframe.fin.txt,
        1,
        15
    )

    ## Loading the created hyper.link.geneID.uniq.txt file
    ## if .dataframe.fin.txt files have been previously created:
    if(length(dataframe.fin.txt)>0){

        # Loading in the unique geneID file:
        unique_geneID_file <- read.delim(
            file='hyper.link.geneID.uniq.txt',
            header = FALSE,
            stringsAsFactors = FALSE
        )

        # Setting the rownames for this file to also be equal
        # to the gene names (for subsetting purposes):
        rownames(unique_geneID_file) <- unique_geneID_file$V1

        # Obtaining the complete list of genes we expect
        # .dataframe.fin.txt files to be created for:
        total_genes_to_be_created <- unique_geneID_file$V1

        # Getting the genes that haven't been processed yet:
        genes_yet_to_be_created <- setdiff(
            total_genes_to_be_created,
            created_ENSG_dataframes
        )

        # Subsetting unique_geneID_file to the genes that have
        # yet to be created:
        unique_geneID_file_to_be_created <- data.frame(
            'V1'=genes_yet_to_be_created,
            stringsAsFactors = FALSE
        )

        # Write the table of genes that have yet to have dataframe.fin.txt files created
        # as a .txt file similar to the hyper.link.geneID.uniq.txt file:
        write.table(
            x= unique_geneID_file_to_be_created,
            file= 'hyper.link.geneID.uniq.to.create.txt',
            quote= FALSE,
            row.names = FALSE,
            col.names = FALSE
        )
    } else{

        # Loading in the unique geneID file:
        unique_geneID_file <- read.delim(
            file='hyper.link.geneID.uniq.txt',
            header = FALSE,
            stringsAsFactors = FALSE
        )

        # Write the table of genes that have yet to have dataframe.fin.txt files created
        # as a .txt file similar to the hypo.link.geneID.uniq.txt file:
        write.table(
            x= unique_geneID_file,
            file= 'hyper.link.geneID.uniq.to.create.txt',
            quote= FALSE,
            row.names = FALSE,
            col.names = FALSE
        )
    }
}

TENET_unique_geneID_subsetting_script_hypo <- function() {
    ## Listing all the previously created *.dataframe.fin.txt files:
    dataframe.fin.txt <- list.files(
        pattern='.dataframe.fin.txt'
    )

    ## Obtaining the list of genes (by ENSG)
    ## which have previously created .dataframe.fin.txt files:
    created_ENSG_dataframes <- substring(
        dataframe.fin.txt,
        1,
        15
    )

    ## Loading the created hypo.link.geneID.uniq.txt file
    ## if .dataframe.fin.txt files have been previously created:
    if(length(dataframe.fin.txt)>0){

        # Loading in the unique geneID file:
        unique_geneID_file <- read.delim(
            file='hypo.link.geneID.uniq.txt',
            header = FALSE,
            stringsAsFactors = FALSE
        )

        # Setting the rownames for this file to also be equal
        # to the gene names (for subsetting purposes):
        rownames(unique_geneID_file) <- unique_geneID_file$V1

        # Obtaining the complete list of genes we expect
        # .dataframe.fin.txt files to be created for:
        total_genes_to_be_created <- unique_geneID_file$V1

        # Getting the genes that haven't been processed yet:
        genes_yet_to_be_created <- setdiff(
            total_genes_to_be_created,
            created_ENSG_dataframes
        )

        # Subsetting unique_geneID_file to the genes that have
        # yet to be created:
        unique_geneID_file_to_be_created <- data.frame(
            'V1'=genes_yet_to_be_created,
            stringsAsFactors = FALSE
        )

        # Write the table of genes that have yet to have dataframe.fin.txt files created
        # as a .txt file similar to the hypo.link.geneID.uniq.txt file:
        write.table(
            x= unique_geneID_file_to_be_created,
            file= 'hypo.link.geneID.uniq.to.create.txt',
            quote= FALSE,
            row.names = FALSE,
            col.names = FALSE
        )
    } else{

        # Loading in the unique geneID file:
        unique_geneID_file <- read.delim(
            file='hypo.link.geneID.uniq.txt',
            header = FALSE,
            stringsAsFactors = FALSE
        )

        # Write the table of genes that have yet to have dataframe.fin.txt files created
        # as a .txt file similar to the hypo.link.geneID.uniq.txt file:
        write.table(
            x= unique_geneID_file,
            file= 'hypo.link.geneID.uniq.to.create.txt',
            quote= FALSE,
            row.names = FALSE,
            col.names = FALSE
        )
    }
}
