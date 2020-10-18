## Load settings information from .rda file:
load("../settings.rda")

## Create new directory to contain states file:
dir.create("hypo.G-.output.states")

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

## Read TF link counts from step5:
ordered_TFs_by_link_count <- read.delim(
  file= paste(
    "./hypo.G-.output.histogram/",
    prefix,
    ".hypo.G-.links.all.tf.freq.txt",
    sep=''
  ),
  stringsAsFactors = FALSE
)

## Load the CpGs linked to each gene from step4:
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

## Get the links associated with just those TFs:
CpG_linkage_dataset_of_interest <- CpG_linkage_dataset[
  CpG_linkage_dataset$geneID %in% top_gene_IDs,
]

## Load the original methylation and expression values from out of step1:
load(
  '../step1/output/luad.diff.methylated.datasets.rda'
)

## Write the function to calculate "personal links"
## This determines if a given tumor sample has a given probe-gene link 
## Based on if it has bonferroni 1-sided p-value expression greater 
## Than normal sample mean, and has less than the hypometh cutoff:
LUAD_tumor_sample_link_evaluator <- function(methylation_probe_link, gene_link){
  
  ## Get actual gene name depending on input:
  if(substring(gene_link,1,4)=='ENSG' & nchar(gene_link)==15){
    
    ## Input is in ENSG, get the gene name:
    gene_link_ENSG <- gene_link
    
    gene_name <- gencode_v22_genes_df[
      gene_link_ENSG, 'gene_name'
    ]
    
    ## Get gene ENSG assuming a name is plugged in:
  } else{
    
    ## Assume what was given was the gene name, get the ENSG:
    gene_link_name <- gene_link
    
    gene_link_ENSG <- rownames(
      gencode_v22_genes_df[gencode_v22_genes_df$gene_name==gene_link_name,]
    )
  }
  
  ## extract methylation values from each tumor sample
  tumor_methylation_values_for_probe <- unlist(
    metDataT[
      methylation_probe_link,
    ]
  )
  
  ## extract gene expression values for each tumor sample
  tumor_expression_values_for_gene <- unlist(
    expDataT[
      gene_link_ENSG,
    ]
  )
  
  ## Get expression values from the normal samples:
  normal_expression_values_for_gene <- unlist(
    expDataN[
      gene_link_ENSG,
    ]
  )
  
  ## Use a for loop to determine for which samples the mean in the normal samples
  ## is significantly lower:
  p_values_compared_to_normal <- numeric()
  
  for(i in tumor_expression_values_for_gene){
    
    p_values_compared_to_normal <- c(
      p_values_compared_to_normal,
      t.test(
        normal_expression_values_for_gene,
        mu=i,
        alternative="greater"
      )$p.value
    )
  }
  
  ## Bonferroni correct the p-values:
  p_values_compared_to_normal_bonf <- p.adjust(
    p_values_compared_to_normal,
    method= 'bonferroni',
    n=length(p_values_compared_to_normal)
  )
  
  ## Create a vector that is 1 if the given sample has expression for the gene
  ## significantly different from the normal samples
  significant_from_normal <- ifelse(
    p_values_compared_to_normal_bonf<0.05,
    1,
    0
  )
  
  ## For each methylation sample, check if it is hypomethylated for the given probe:
  cutoff_values <- ifelse(
    tumor_methylation_values_for_probe < hypocutoff,
    1,
    0
  )
  
  ## Create a final vector with 1 for where both conditions are met, and 0 when
  ## they arent both met:
  final_output <- ifelse(
    significant_from_normal==1,
    ifelse(
      cutoff_values==1,
      1,
      0
    ),
    0
  )
  
  ## Assign sample names to the output values:
  names(final_output) <- colnames(expDataT)
  
  ## Return the final output:
  return(final_output)
}

## Calculate if hypometh link is present:
links_dataset <- mapply(
  LUAD_tumor_sample_link_evaluator,
  methylation_probe_link=CpG_linkage_dataset_of_interest$probe,
  gene_link=CpG_linkage_dataset_of_interest$geneID
)

## Change the column names to those of the probe + gene link:
colnames(links_dataset) <- paste(
  CpG_linkage_dataset_of_interest$probe,
  CpG_linkage_dataset_of_interest$geneSymbol,
  sep='_'
)

## Transpose the dataset to make the tumor samples the columns:
links_dataset_transposed <- t(links_dataset)

## Convert it to a dataframe 
##(for potential ease of use in possible downstream analyses)
links_dataset_transposed_df <- as.data.frame(
  links_dataset_transposed,
  stringsAsFactors= FALSE
)

## Write the table out to the states subdirectory:
write.table(
  links_dataset_transposed_df,
  file='./hypo.G-.output.states/links.states.table.tsv',
  quote= FALSE,
  sep='\t'
)