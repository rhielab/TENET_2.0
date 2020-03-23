## Load settings information from .rda file:
load("../settings.rda")

## Load BioCircos:
library('BioCircos')

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