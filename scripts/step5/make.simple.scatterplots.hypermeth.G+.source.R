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
          ' expression vs.\n',
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
          'expression',
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