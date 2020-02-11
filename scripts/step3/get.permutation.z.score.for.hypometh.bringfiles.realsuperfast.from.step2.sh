#!/bin/bash
mkdir hypo.output
cd hypo.output
# this is for all G+ and G-
awk 'BEGIN {OFS="\t"} {print $3}' ../../step2/hypo.output/hypo.G+.link.zscore.txt > hypo.G+.link.geneID.txt
awk 'BEGIN {OFS="\t"} {print $3}' ../../step2/hypo.output/hypo.G-.link.zscore.txt > hypo.G-.link.geneID.txt
# let's make zscore for unique geneID #
cat hypo.G+.link.geneID.txt hypo.G-.link.geneID.txt > hypo.link.geneID.txt
sort hypo.link.geneID.txt | uniq > hypo.link.geneID.uniq.txt
Rscript ../../scripts/step3/TENET_unique_geneID_subsetting_script_hypo.R
sort -t $'\t' -k 3,3 ../../step2/hypo.output/hypo.zscore.all.genes.rda.txt > hypo.zscore.all.genes.rda.sorted.txt

## Input path to file containing all genes:
GENE_FILE="../../step1/output/TENET_complete_gene_list.txt"

## Input path to file containing Hypo genes of interest:
HYPO_GENES_OF_INTEREST_FILE="../../step3/hypo.output/hypo.link.geneID.uniq.txt"

## Input name of file you want to subset:
HYPO_ORIGINAL_SUBSET_FILE="../../step3/hypo.output/hypo.zscore.all.genes.rda.sorted.txt"

## Path to step1 metadata:
METADATA_PATH="../../step1/output/TENET_step1_overall_metadata.txt"

## Count number of genes present in the data file:
gene_count=$(< ${GENE_FILE} wc -l)

## Count number of hypo_probes present in the data file:
hypo_probes=$(sed -n '5p' ${METADATA_PATH} | awk '{print $2;}')

## Create array with all numbers from 1 to the number of genes:
gene_number_array=($(seq 1 1 $(< ${GENE_FILE} wc -l)))

## Read in the genes as an array:
IFS=$'\n' read -d '' -r -a gene_name_array < ${GENE_FILE}

#### Fix potential issue with "missing" CpGs:
#### Occasionally CpGs will be omitted from the master list in step2, resulting in an incomplete indexing. This is a fix for now and in later TENET versions we will need to create a more formal patch:

## Create dummy placeholder.txt file to properly determine number of cgs present:
head -n ${hypo_probes} ${HYPO_ORIGINAL_SUBSET_FILE} > placeholder.txt

## obtain all cgs present in the placeholder
placeholder_cg_array=($(cut -d$'\t' -f1 placeholder.txt))

## obtain unique cgs present in the placeholder 
## (This is the cgs actually be present in each sliced file)
actual_unique_cgs=($(echo "${placeholder_cg_array[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

## Obtain the unique cg count (this is the number to actually slice on!)
actual_unique_cg_count=(${#actual_unique_cgs[@]})

## Split the .rda.sorted.txt file into individual files starting with 'HYPO_GENE' 
split -d -l ${hypo_probes} -a ${#gene_count} --additional-suffix=.txt ${HYPO_ORIGINAL_SUBSET_FILE} HYPO_GENE

## Read in an array of all the HYPO_GENE file names:
current_file_names_complete=($(ls -d HYPO_GENE*))

## Get just the split file names: 
current_file_names=("${current_file_names_complete[@]##*/}")

## Create a vector that contains the names we want to rename the files to:
new_file_names=("${gene_name_array[@]/%/.hypo.zscore.dataframe.fin.txt}")

## Rename the files:
for i in "${!current_file_names[@]}"
do
 	mv "${current_file_names[i]}" "${new_file_names[i]}"
done

## Read in the hypo.G genes of interest:
IFS=$'\n' read -d '' -r -a hypo_genes_of_interest < ${HYPO_GENES_OF_INTEREST_FILE}

## Obtain array of non-hypo genes of interest:
non_genes_of_interest=(`echo ${gene_name_array[@]} ${hypo_genes_of_interest[@]} | tr ' ' '\n' | sort | uniq -u `)

## Create a new vector of gene name files for genes that aren't of interest:
non_genes_of_interest_files=("${non_genes_of_interest[@]/%/.hypo.zscore.dataframe.fin.txt}")

## Remove those files from the folder:
for j in "${non_genes_of_interest_files[@]}"
do
	rm "${j}"
done

## Remove the placeholder file:
# rm "placeholder.txt"

Rscript ../../scripts/step3/get.permutation.z.score.for.hypometh.enhancers.bringfiles.realfast.from.step2.R
allrda1=$(wc -l ./hypo.G+.output/*perm.txt)
allrda1c=$(echo $allrda1 | awk '{print $1}')
allrda2=$(wc -l ./hypo.G-.output/*perm.txt)
allrda2c=$(echo $allrda2 | awk '{print $1}')
if [[ $allrda1c > 0 && $allrda2c > 0 ]];
rm *fin.txt;
then echo "hypo permutaion is done";
fi
