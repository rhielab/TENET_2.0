## Load settings information from .rda file:
load("../settings.rda")

## Load ggplot2:
library('ggplot2')

## Load survival:
library('survival')

## Load parallel:
library('parallel')

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