## Load the matlab library:
library(matlab)

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
jet_color_numeric_color_grad <- matlab::jet.colors(track_top_n_genes)

names(jet_color_numeric_color_grad) <- top_gene_ids

### Now that we have the information, let's assemble the real intersect file:
output_intersect_file <- data.frame(
  'chrom'= intersect_df$TR_chr,
  'chromStart'= intersect_df$TR_start,
  'chromEnd'= intersect_df$TR_end,
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
  'color'= jet_color_numeric_color_grad[
    intersect_df$TR_ID
  ],
  'sourceChrom'= intersect_df$TR_chr,
  'sourceStart'= intersect_df$TR_start,
  'sourceEnd'= intersect_df$TR_end,
  'sourceName'= intersect_df$TR_name,
  'sourceStrand'= rep(
    '.',
    nrow(intersect_df)
  ),
  'targetChrom'= intersect_df$CpG_chr,
  'targetStart'= intersect_df$CpG_start+1,
  'targetEnd'= intersect_df$CpG_end,
  'targetName'= intersect_df$CpG_ID,
  'targetStrand'= rep(
    '.',
    nrow(intersect_df)
  ),
  stringsAsFactors = FALSE
)

## Save the file out as a .bed file:
write.table(
  output_intersect_file,
  file= paste(
    './hyper.G-.output.tracks/',
    'TR.enhancer.probe.links.hg38.bed',
    sep=''
  ),
  row.names= FALSE,
  col.names= FALSE,
  quote= FALSE,
  sep='\t'
)