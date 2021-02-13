#!/bin/bash
mkdir hyper.output
cd hyper.output/
# this is for all G+ and G-
awk 'BEGIN {OFS="\t"} {print $3}' ../../step2/hyper.output/hyper.G-.link.zscore.txt > hyper.G-.link.geneID.txt
# let's make zscore for unique geneID #
cp hyper.G-.link.geneID.txt hyper.link.geneID.txt
sort hyper.link.geneID.txt | uniq > hyper.link.geneID.uniq.txt
sort -t $'\t' -k 3,3 ../../step2/hyper.output/hyper.zscore.all.genes.rda.txt > hyper.zscore.all.genes.rda.sorted.txt

## Input path to file containing all genes:
GENE_FILE="../../step1/output/TENET_complete_gene_list.txt"

## Input path to file containing Hypo genes of interest:
HYPER_GENES_OF_INTEREST_FILE="../../step3/hyper.output/hyper.link.geneID.uniq.txt"

## Input name of file you want to subset:
HYPER_ORIGINAL_SUBSET_FILE="../../step3/hyper.output/hyper.zscore.all.genes.rda.sorted.txt"

## Path to step1 metadata:
METADATA_PATH="../../step1/output/TENET_step1_overall_metadata.txt"

## Count number of genes present in the data file:
gene_count=$(< ${GENE_FILE} wc -l)

## Count number of hyper_probes present in the data file:
hyper_probes=$(sed -n '3p' ${METADATA_PATH} | awk '{print $2;}')

## Create array with all numbers from 1 to the number of genes:
gene_number_array=($(seq 1 1 $(< ${GENE_FILE} wc -l)))

## Read in the genes as an array:
IFS=$'\n' read -d '' -r -a gene_name_array < ${GENE_FILE}

#### Fix potential issue with "missing" CpGs:
#### Occasionally CpGs will be omitted from the master list in step2, resulting in an incomplete indexing. This is a fix for now and in later TENET versions we will need to create a more formal patch:

## Create dummy placeholder.txt file to properly determine number of cgs present:
head -n ${hyper_probes} ${HYPER_ORIGINAL_SUBSET_FILE} > placeholder.txt

## obtain all cgs present in the placeholder
placeholder_cg_array=($(cut -d$'\t' -f1 placeholder.txt))

## obtain unique cgs present in the placeholder 
## (This is the cgs actually be present in each sliced file)
actual_unique_cgs=($(echo "${placeholder_cg_array[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

## Obtain the unique cg count (this is the number to actually slice on!)
actual_unique_cg_count=(${#actual_unique_cgs[@]})

## Split the .rda.sorted.txt file into individual files starting with 'HYPO_GENE' 
split -d -l ${actual_unique_cg_count} -a ${#gene_count} --additional-suffix=.txt ${HYPER_ORIGINAL_SUBSET_FILE} HYPER_GENE

## Read in an array of all the HYPER_GENE file names:
current_file_names_complete=($(ls -d HYPER_GENE*))

## Get just the split file names: 
# current_file_names=("${current_file_names_complete[@]##*/}")

current_file_names_complete=($(echo HYPO_GENE* | xargs ls))

## Create a vector that contains the names we want to rename the files to:
new_file_names=("${gene_name_array[@]/%/.hyper.zscore.dataframe.fin.txt}")

## Rename the files:
for i in "${!current_file_names[@]}"
do
 	mv "${current_file_names[i]}" "${new_file_names[i]}"
done

## Read in the hyper.G genes of interest:
IFS=$'\n' read -d '' -r -a hyper_genes_of_interest < ${HYPER_GENES_OF_INTEREST_FILE}

## Obtain array of non-hypo genes of interest:
non_genes_of_interest=(`echo ${gene_name_array[@]} ${hyper_genes_of_interest[@]} | tr ' ' '\n' | sort | uniq -u `)

## Create a new vector of gene name files for genes that aren't of interest:
non_genes_of_interest_files=("${non_genes_of_interest[@]/%/.hyper.zscore.dataframe.fin.txt}")

## Remove those files from the folder:
for j in "${non_genes_of_interest_files[@]}"
do
	rm "${j}"
done

## Remove the placeholder file:
# rm "placeholder.txt"

Rscript ../../scripts/step3/get.permutation.z.score.for.hyperG-.enhancers.bringfiles.realfast.from.step2.R
allrda=$(wc -l ./hyper.G-.output/*perm.txt)
allrdac=$(echo $allrda | awk '{print $1}')
if [[ $allrdac > 0 ]];
rm *fin.txt;
then echo "hyper.G- permutaion is done";
fi
