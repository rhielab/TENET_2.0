## Load ggplot2:
library(ggplot2)

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
