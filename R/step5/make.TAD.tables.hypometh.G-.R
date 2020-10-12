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
