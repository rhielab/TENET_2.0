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

