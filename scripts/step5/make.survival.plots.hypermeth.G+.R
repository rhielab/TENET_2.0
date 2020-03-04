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
rm(gtf, gencode_v22_genes)

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
  
  ## combine expression datasets:
  expData <- cbind(expDataN, expDataT)
  
  ## Get the list of normal and tumor samples:
  normal_samples <- colnames(expDataN)
  tumor_samples <- colnames(expDataT)
  
  ## Get expression values for gene of interest:
  expression_values <- unlist(
    expData[
      gene_ENSG,
    ]
  )
  names(expression_values) <- colnames(expData)
  
  # Setting the normal and tumor datasets:
  if(exists('normal_samples') & exists('tumor_samples')){
    
    # Split gene expression values into normal and tumor samples
    # If both are present in the dataset 
    tumor_expression_values <- expression_values[tumor_samples]
    normal_expression_values <- expression_values[normal_samples]
    
  } else{
    
    # Set the two data groups equal to the whole set (and each other)
    tumor_expression_values <- expression_values
    normal_expression_values <- expression_values
  }
  
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
  
  # Subsetting clinical patient data for individuals with
  # gene expression data in the tumor set:
  relevant_clinical <- clinical[
    clinical[,"bcr_patient_barcode"] %in% names(tumor_expression_values),
  ]
  
  # Keep the patient barcode (saved as rownames and names of the expression values), 
  # days to death, days since last follow up, and vital status columns only
  # other data is unecessary
  relevant_clinical <- relevant_clinical[
    ,c(
      "bcr_patient_barcode",
      "days_to_death",
      "days_to_last_followup",
      "vital_status"
    )
  ]
  
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
    
  ## Sort tumor expression values high to low:
  tumor_expression_decreasing_ordered <- sort(
    tumor_expression_values,
    decreasing=TRUE
  )
    
  # Get expression value in tumor samples that represents the 
  # high threshold specified by the user:
  tumor_expression_high_threshold_value <- quantile(
    tumor_expression_decreasing_ordered,
    high_cutoff
  )[1]
    
  # Get expression value in tumor samples that represents the 
  # low threshold specified by the user:
  tumor_expression_low_threshold_value <- quantile(
    tumor_expression_decreasing_ordered,
    low_cutoff
  )[1]
      
  # Get count of total tumor samples
  tumor_sample_n <- as.numeric(length(tumor_expression_decreasing_ordered))
      
  # Getting sample number which is just higher than the set threshold 
  # for the high expression group
  tumor_expression_high_threshold_sample_number <- max(
    which(
      tumor_expression_decreasing_ordered > tumor_expression_high_threshold_value
    )
  )
      
  # Getting sample number which is just lower than or equal to the set threshold 
  # for the low expression group
  tumor_expression_low_threshold_sample_number <- min(
    which(
      tumor_expression_decreasing_ordered <= tumor_expression_low_threshold_value
    )
  )
      
  # Get the names of the subjects in the high and low expression groups:
  tumor_expression_high_samples <- names(
    tumor_expression_decreasing_ordered[
      1:tumor_expression_high_threshold_sample_number
    ]
  )
  
  tumor_expression_low_samples <- names(
    tumor_expression_decreasing_ordered[
      tumor_expression_low_threshold_sample_number:tumor_sample_n
    ]
  )
      
  # Get the names of the subjects not in either group:
  tumor_expression_intermediate_samples <- names(
    tumor_expression_decreasing_ordered[
      which(
        (tumor_expression_decreasing_ordered) > tumor_expression_low_threshold_value &
        (tumor_expression_decreasing_ordered) <= tumor_expression_high_threshold_value
      )
    ]
  )
  
  # Because some samples were removed from the clinical data due to lack of
  # survival information, we will also have to remove those samples from 
  # the groups listed here:
  tumor_expression_high_samples_clinical_present <- tumor_expression_high_samples[
    tumor_expression_high_samples %in% rownames(relevant_clinical)
  ]
  
  tumor_expression_low_samples_clinical_present <- tumor_expression_low_samples[
    tumor_expression_low_samples %in% rownames(relevant_clinical)
  ]
  
  tumor_expression_intermediate_samples_clinical_present <- tumor_expression_intermediate_samples[
    tumor_expression_intermediate_samples %in% rownames(relevant_clinical)
  ]
      
  # Create individual clinical datasets for the two groups with clinical data
  # and a single combined dataset 
  high_relevant_clinical <- relevant_clinical[
    c(tumor_expression_high_samples_clinical_present),
  ]
  
  low_relevant_clinical <- relevant_clinical[
    c(tumor_expression_low_samples_clinical_present),
  ]
  
  high_low_relevant_clinical <- relevant_clinical[
    c(tumor_expression_high_samples_clinical_present, tumor_expression_low_samples_clinical_present),
  ]
      
  # Obtaining days to death information and death status
  # from high and low expression group and the combined set:
  high_days_to_death <- as.numeric(high_relevant_clinical[, "days_to_death"])

  low_days_to_death <- as.numeric(low_relevant_clinical[, "days_to_death"])
  
  high_low_days_to_death <- as.numeric(high_low_relevant_clinical[, "days_to_death"])
  high_low_death_status <- high_low_days_to_death > 0
      
  # Count the number of patients who had died:
  high_total_deaths <- sum(high_days_to_death >= 0)
  
  low_total_deaths <- sum(low_days_to_death >= 0)
  
  high_low_total_deaths <- sum(high_low_days_to_death >= 0)
      
  # Get expression for the gene of interest in normal samples:
  normal_expression_mean <- mean(normal_expression_values)
  
  # Get expression for the high expression tumor samples in total
  # and only those with clinical data:
  tumor_high_expression_mean <- mean(
    tumor_expression_values[
      tumor_expression_high_samples
    ]
  )
  
  tumor_high_expression_relevant_clinical_mean <- mean(
    tumor_expression_values[
      tumor_expression_high_samples_clinical_present
    ]
  )
  
  # Get expression for the low expression tumor samples in total
  # and only those with clinical data:
  tumor_low_expression_mean <- mean(
    tumor_expression_values[
      tumor_expression_low_samples
      ]
  )
  
  tumor_low_expression_relevant_clinical_mean <- mean(
    tumor_expression_values[
      tumor_expression_low_samples_clinical_present
      ]
  )
 
  # For patients that are not dead replace the days to death with last followup info:
  high_low_days_to_death[high_low_days_to_death<0] <- as.numeric(
    high_low_relevant_clinical[
      high_low_days_to_death<0, 
      "days_to_last_followup"
    ]
  )
  
  # If there are any remaining -Inf values in the 
  # high_low_days_to_death data, replace them with 0
  high_low_days_to_death[high_low_days_to_death==-Inf] <- 0
      
  # Create the survival object for the high + low group
  # Using days to death or days to last follow up data
  # Plus the death status:
  high_low_survival_object <- Surv(
    high_low_days_to_death, 
    high_low_death_status
  )
      
  # Set the rownames of the survival object to be the 
  # patient names for the high/low group smaples:
  rownames(high_low_survival_object) <- rownames(high_low_relevant_clinical)
      
  # Creating legend names for high/low expression groups:
  legend_name_high <- paste(gene_name,"high")

  legend_name_low <- paste(gene_name,"low")
      
  # Calculating the survival p-value with a tryCatch:
  survival_pvalue <- tryCatch(
    {
          
      # Create a vector with "high" repeated once
      # for each sample in the high expression group:
      high_vector <- rep(
        "high", 
        nrow(high_relevant_clinical)
      )
      
      # Create a vector with "low" repeated once
      # for each sample in the low expression group:
      low_vector <- rep(
        "low", 
        nrow(low_relevant_clinical)
      )
          
      # Combine the top and down vectors into a single vector
      # in the order top-down
      expression_group <- c(high_vector, low_vector)
          
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
      
    }, 
    error = function(e){
      
      # if there is an error with the above analysis, return Inf
      # as the p-value
      return(Inf)
      
    }
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
        
    # Create a vector with "high" repeated once
    # for each sample in the high expression group:
    high_vector <- rep(
      "high", 
      nrow(high_relevant_clinical)
    )
    
    # Create a vector with "low" repeated once
    # for each sample in the low expression group:
    low_vector <- rep(
      "low", 
      nrow(low_relevant_clinical)
    )
    
    # Combine the top and down vectors into a single vector
    # in the order top-down
    expression_group <- c(high_vector, low_vector)
    
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
        max(high_low_relevant_clinical$days_to_last_followup) - 2500
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
    
    ## If not graphing, return the p-value instead:
    return(survival_pvalue)
  }
    
}

## Create a function to get the survival p-value or graph for 
##beach probe of interest:
methylation_survival_function_graph <- function(
  probe_of_interest,
  linked_gene_of_interest,
  high_cutoff,
  low_cutoff,
  graph
){
  
  ## Get actual gene name for the linked gene depending on input:
  if(substring(linked_gene_of_interest,1,4)=='ENSG' & nchar(linked_gene_of_interest)==15){
    
    ## Input is in ENSG, get the gene name:
    gene_ENSG <- linked_gene_of_interest
    
    gene_name <- gencode_v22_genes_df[
      gene_ENSG, 'gene_name'
      ]
    
    ## Get gene ENSG assuming a name is plugged in:
  } else{
    
    ## Assume what was given was the gene name, get the ENSG:
    gene_name <- linked_gene_of_interest
    
    gene_ENSG <- rownames(
      gencode_v22_genes_df[gencode_v22_genes_df$gene_name==gene_name,]
    )
  }
  
  ## combine methylation datasets:
  metData <- cbind(metDataN, metDataT)
  
  ## Get the list of normal and tumor samples:
  normal_samples <- colnames(metDataN)
  tumor_samples <- colnames(metDataT)
  
  ## Get methylation values for probe of interest:
  methylation_values <- unlist(
    metData[
      probe_of_interest,
      ]
  )
  names(methylation_values) <- colnames(metData)
  
  # Setting the normal and tumor datasets:
  if(exists('normal_samples') & exists('tumor_samples')){
    
    # Split methylation values into normal and tumor samples
    # If both are present in the dataset 
    tumor_methylation_values <- methylation_values[tumor_samples]
    normal_methylation_values <- methylation_values[normal_samples]
    
  } else{
    
    # Set the two data groups equal to the whole set (and each other)
    tumor_methylation_values <- methylation_values
    normal_methylation_values <- methylation_values
  }
  
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
  
  # Subsetting clinical patient data for individuals with
  # methylation data in the tumor set:
  relevant_clinical <- clinical[
    clinical[,"bcr_patient_barcode"] %in% names(tumor_methylation_values),
    ]
  
  # Keep the patient barcode (saved as rownames and names of the expression values), 
  # days to death, days since last follow up, and vital status columns only
  # other data is unecessary
  relevant_clinical <- relevant_clinical[
    ,c(
      "bcr_patient_barcode",
      "days_to_death",
      "days_to_last_followup",
      "vital_status"
    )
    ]
  
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
  
  ## Sort tumor methylation values high to low:
  tumor_methylation_decreasing_ordered <- sort(
    tumor_methylation_values,
    decreasing=TRUE
  )
  
  # Get expression value in tumor samples that represents the 
  # high threshold specified by the user:
  tumor_methylation_high_threshold_value <- quantile(
    tumor_methylation_decreasing_ordered,
    high_cutoff
  )[1]
  
  # Get expression value in tumor samples that represents the 
  # low threshold specified by the user:
  tumor_methylation_low_threshold_value <- quantile(
    tumor_methylation_decreasing_ordered,
    low_cutoff
  )[1]
  
  # Get count of total tumor samples
  tumor_sample_n <- as.numeric(length(tumor_methylation_decreasing_ordered))
  
  # Getting sample number which is just higher than the set threshold 
  # for the high methylation group
  tumor_methylation_high_threshold_sample_number <- max(
    which(
      tumor_methylation_decreasing_ordered > tumor_methylation_high_threshold_value
    )
  )
  
  # Getting sample number which is just lower than or equal to the set threshold 
  # for the low methylation group
  tumor_methylation_low_threshold_sample_number <- min(
    which(
      tumor_methylation_decreasing_ordered <= tumor_methylation_low_threshold_value
    )
  )
  
  # Get the names of the subjects in the high and low methylation groups:
  tumor_methylation_high_samples <- names(
    tumor_methylation_decreasing_ordered[
      1:tumor_methylation_high_threshold_sample_number
      ]
  )
  
  tumor_methylation_low_samples <- names(
    tumor_methylation_decreasing_ordered[
      tumor_methylation_low_threshold_sample_number:tumor_sample_n
      ]
  )
  
  # Get the names of the subjects not in either group:
  tumor_methylation_intermediate_samples <- names(
    tumor_methylation_decreasing_ordered[
      which(
        (tumor_methylation_decreasing_ordered) > tumor_methylation_low_threshold_value &
          (tumor_methylation_decreasing_ordered) <= tumor_methylation_high_threshold_value
      )
      ]
  )
  
  # Because some samples were removed from the clinical data due to lack of
  # survival information, we will also have to remove those samples from 
  # the groups listed here:
  tumor_methylation_high_samples_clinical_present <- tumor_methylation_high_samples[
    tumor_methylation_high_samples %in% rownames(relevant_clinical)
    ]
  
  tumor_methylation_low_samples_clinical_present <- tumor_methylation_low_samples[
    tumor_methylation_low_samples %in% rownames(relevant_clinical)
    ]
  
  tumor_methylation_intermediate_samples_clinical_present <- tumor_methylation_intermediate_samples[
    tumor_methylation_intermediate_samples %in% rownames(relevant_clinical)
    ]
  
  # Create individual clinical datasets for the two groups with clinical data
  # and a single combined dataset 
  high_relevant_clinical <- relevant_clinical[
    c(tumor_methylation_high_samples_clinical_present),
    ]
  
  low_relevant_clinical <- relevant_clinical[
    c(tumor_methylation_low_samples_clinical_present),
    ]
  
  high_low_relevant_clinical <- relevant_clinical[
    c(tumor_methylation_high_samples_clinical_present, tumor_methylation_low_samples_clinical_present),
    ]
  
  # Obtaining days to death information and death status
  # from high and low methylation group and the combined set:
  high_days_to_death <- as.numeric(high_relevant_clinical[, "days_to_death"])
  
  low_days_to_death <- as.numeric(low_relevant_clinical[, "days_to_death"])
  
  high_low_days_to_death <- as.numeric(high_low_relevant_clinical[, "days_to_death"])
  high_low_death_status <- high_low_days_to_death > 0
  
  # Count the number of patients who had died:
  high_total_deaths <- sum(high_days_to_death >= 0)
  
  low_total_deaths <- sum(low_days_to_death >= 0)
  
  high_low_total_deaths <- sum(high_low_days_to_death >= 0)
  
  # Get methylation for the probe of interest in normal samples:
  normal_methylation_mean <- mean(normal_methylation_values)
  
  # Get methylation for the high methylation tumor samples in total
  # and only those with clinical data:
  tumor_high_methylation_mean <- mean(
    tumor_methylation_values[
      tumor_methylation_high_samples
      ]
  )
  
  tumor_high_methylation_relevant_clinical_mean <- mean(
    tumor_methylation_values[
      tumor_methylation_high_samples_clinical_present
      ]
  )
  
  # Get methylation for the low methylation tumor samples in total
  # and only those with clinical data:
  tumor_low_methylation_mean <- mean(
    tumor_methylation_values[
      tumor_methylation_low_samples
      ]
  )
  
  tumor_low_methylation_relevant_clinical_mean <- mean(
    tumor_methylation_values[
      tumor_methylation_low_samples_clinical_present
      ]
  )
  
  # For patients that are not dead replace the days to death with last followup info:
  high_low_days_to_death[high_low_days_to_death<0] <- as.numeric(
    high_low_relevant_clinical[
      high_low_days_to_death<0, 
      "days_to_last_followup"
      ]
  )
  
  # If there are any remaining -Inf values in the 
  # high_low_days_to_death data, replace them with 0
  high_low_days_to_death[high_low_days_to_death==-Inf] <- 0
  
  # Create the survival object for the high + low group
  # Using days to death or days to last follow up data
  # Plus the death status:
  high_low_survival_object <- Surv(
    high_low_days_to_death, 
    high_low_death_status
  )
  
  # Set the rownames of the survival object to be the 
  # patient names for the high/low group smaples:
  rownames(high_low_survival_object) <- rownames(high_low_relevant_clinical)
  
  # Creating legend names for high/low methylation groups:
  legend_name_high <- paste(probe_of_interest,"high")
  
  legend_name_low <- paste(probe_of_interest,"low")
  
  # Calculating the survival p-value with a tryCatch:
  survival_pvalue <- tryCatch(
    {
      
      # Create a vector with "high" repeated once
      # for each sample in the high methylation group:
      high_vector <- rep(
        "high", 
        nrow(high_relevant_clinical)
      )
      
      # Create a vector with "low" repeated once
      # for each sample in the low methylation group:
      low_vector <- rep(
        "low", 
        nrow(low_relevant_clinical)
      )
      
      # Combine the high and low groups:
      methylation_group <- c(high_vector, low_vector)
      
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
      
    }, 
    error = function(e){
      
      # if there is an error with the above analysis, return Inf
      # as the p-value
      return(Inf)
      
    }
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
    
    # Create a vector with "high" repeated once
    # for each sample in the high expression group:
    high_vector <- rep(
      "high", 
      nrow(high_relevant_clinical)
    )
    
    # Create a vector with "low" repeated once
    # for each sample in the low expression group:
    low_vector <- rep(
      "low", 
      nrow(low_relevant_clinical)
    )
    
    # Combine the high and low groups:
    methylation_group <- c(high_vector, low_vector)
    
    ## Get the path to the step 5 folder to save
    ## the pdf in:
    path_to_folder <- getwd()
    
    ## Create a title for the survival plot pdf:
    ## This is a comination of the probe name with the linked gene:
    survival_plot_pdf_title <- paste(
      path_to_folder,
      '/hypo.G+.output.survival/',
      probe_of_interest,
      '_',
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
        max(high_low_relevant_clinical$days_to_last_followup) - 2500
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
    
    ## If not graphing, return the p-value instead:
    return(survival_pvalue)
  }
  
}

## Create survival plots for each of the genes designated:
mclapply(
  X= top_gene_names,
  FUN= expression_survival_function_graph,
  high_cutoff= high_thresh,
  low_cutoff= low_thresh,
  graph= TRUE,
  mc.cores= cores
)

## Now for each of the genes, we will need to calculate survival-p values
## for each of the probes linked to them, and we will plot the significant ones:
for(gene in top_gene_names){
  
  ## Get the list of probes linked to that gene:
  CpGs_linked <- CpG_linkage_dataset[
    CpG_linkage_dataset$geneID==gene,
    'probe'
  ]
  
  ## Get the p-values for each probe by mclapplying the
  ## methylation survival function to all of the probes
  ## linked to each gene:
  CpGs_linked_p_values <- mclapply(
    X= CpGs_linked,
    FUN= methylation_survival_function_graph,
    linked_gene_of_interest= gene,
    high_cutoff= high_thresh,
    low_cutoff= low_thresh,
    graph= FALSE,
    mc.cores= cores
  )
  
  ## Create a vector of the CpG p-values and add the CpG entries
  ## as names
  CpGs_linked_p_values_vector <- unlist(CpGs_linked_p_values)
  names(CpGs_linked_p_values_vector) <- CpGs_linked
  
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
  
  #### Save plots for the nominally significant CpGs:
  mclapply(
    X= CpGs_linked_nominally_significant_names,
    FUN= methylation_survival_function_graph,
    linked_gene_of_interest= gene,
    high_cutoff= high_thresh,
    low_cutoff= low_thresh,
    graph= TRUE,
    mc.cores= cores
  )
  
}