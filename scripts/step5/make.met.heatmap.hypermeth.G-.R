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
  ylab= NULL,
  xlab= NULL
)

## Close the plot:
dev.off()