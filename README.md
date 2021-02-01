# TENET (Tracing Enhancer Networks using Epigenetic Traits) 2.0

#### Last updated: 1/29/2021

TENET is developed to identify key transcriptional regulators such as transcription factors and enhancers linked to a specific cell type. Please see Rhie et al (PMID: 27833659) and Mullen et al (PMID: 32925947) for its usage. 

## Citation
Mullen DJ et al. TENET 2.0: identification of key transcriptional regulators and enhancers in lung adenocarcinoma. Plos Genetics 2020 16(9): e1009023. https://doi.org/10.1371/journal.pgen.1009023

Rhie SK et al. Identification of activated enhancers and linked transcription factors in breast, prostate, and kidney tumors by tracing enhancer networks using epigenetic traits. Epigenetics Chromatin 2016 Nov 9;9:50. PMID: 27833659

## 1. TENET 2.0 Installation

##### This build of TENET 2.0 depends on the following programs

###### A. bedtools (please install from https://github.com/arq5x/bedtools2/releases)

###### B. R + R packages

TENET 2.0 requires the use of several R packages including: BioCircos, GenomicRanges, ggplot2, htmlwidgets, matlab, parallel, and survival. BiocManager may also be required to install these.

The following code can be used to install these packages in your installation of R

```diff
if (!require("ggplot2", quietly = TRUE)){
  install.packages('ggplot2')
}

if (!require("parallel", quietly = TRUE)){
  install.packages("parallel")
}

if (!require("survival", quietly = TRUE)){
  install.packages("survival")
}

if (!require("htmlwidgets", quietly = TRUE)){
  install.packages("htmlwidgets")
}

if (!require("matlab", quietly = TRUE)){
  install.packages("matlab")
}

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("BioCircos", quietly = TRUE)){
  BiocManager::install("BioCircos")
}

if (!require("GenomicRanges", quietly = TRUE)){
  BiocManager::install("GenomicRanges")
}

```
TENET 2.0 also makes use of heatmap.3 functionality from: 
https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R

TENET 2.0 also makes use of enhancer probes from v1.5.1 of the ELMER package


###### C. Download the TENET 2.0 program from this page into directory of your choosing

## 2. TENET 2.0 Usage

###### A. Get methylation and gene expression data with matched sample IDs

You can either put txt files in ./external.data/data/methylation and ./external.data/data/expression (add case and ctrl in the file name, respectively) or put a rda file in ./external.data/data (for case, name as "metDataT", "expDataT", for ctrl, name as "metDataN", "expDataN"). .rda file can also include a "clinical" object with TCGA clinical data (for use in step5 survival functions). Matching sample names with at least 12 characters (i.e. TCGA-##-####) should be included in the columns of the expData and metData objects, while the GENCODE v22 Ensembl gene IDs (i.e. ENSG###########) or HM450 probe IDs (i.e. cg########) should be contained in the rows of the expData and metData objects, respectively. 

###### B. Fill out parameter settings in settings.txt 

###### C. Run TENET by typing 

```diff
bash ./TENET.sh
``` 
(OPTIONAL 1) put other cells methylation data you want to use for purity adjustment
you can put a txt file in ./external.data/othercells or put a rda file in ./external.data/others (name as "extData")

(OPTIONAL 2) if you have external datasets for enhancer peaks (must be in bed-like format and end with '.bed'), add files in ./external.data/enhancer

(OPTIONAL 3) if you have external datasets for NDR peaks (must be in bed-like format and end with '.bed'), add files in ./external.data/NDR

(OPTIONAL 4) if you have external datasets for other peaks (must be in bed-like format and end with '.bed'), add files in ./external.data/feature

(OPTIONAL 5) if you have Purity, CNV, SM datasets for sample samples, add files in ./external.data/otherinfo (include "CNV", "SM", "Purity" in the file names)

(OPTIONAL 6) if you have TAD file info (must be in bed-like format and can't end with '.rtf'), add files in ./external.data/TAD

## 3. TENET 2.0 Parameters

###### prefix
	Prefix for output files 

###### cores
	Number of cores you want to use in your machine

###### step0
	T: Process data into a .rda file when you input a .txt file for methylation and expression data (when you put other cells' methylation data in a .txt file) 
	
	F: Use when you put a pre-created .rda file containing methylation and expression data

###### step1
	T: Find differentially methylated enhancer regions
	
	F: Do not run step1

###### step2
	T: Select enhancer-gene links by calculating Z scores
	
	F: Do not run step2

###### step3
	T: Select significant enhancer-gene links from step2 by permutation
	
	F: Do not run step3

###### step4
	T: Optimize selection of enhancer-gene lnks
	
	F: Do not run step4

###### step5
	T: Summarize and visualize enhancer-gene links
	
	F: Do not run step5

###### methcutoff
	Set the methylation cutoff of beta-values, range from 0 to 1

###### hypocutoff
	Set the hypomethylation cutoff of beta-values, range from 0 to 1

###### unmethcutoff
	Set the unmethylation cutoff of beta-values, range from 0 to 1

###### hypercutoff
	Set the hypermethylation cutoff of beta-values, range from 0 to 1

###### minTumors
	Set the minimum number of tumors with epigenetic events to determine links

###### leuk
	T: Use leukocyte methylation data to adjust purity
	
	F: Don't use

###### fibro
	T: Use fibroblast methylation data to adjust purity
	
	F: Don't use

###### sm
	T: Use smooth muscle methylation data to adjust purity
	
	F: Don't use

###### othercells
	T: Use other cells methylation data which you loaded in ./external.data/othercells to adjust purity (see Optional 1 above)
	
	F: Don't use

###### purityinfo
	T: Use purity estimates which you loaded in ./external.data/otherinfo to make complex scatterplots
	
	F: Don't use purityinfo

###### udist
	Set minimum upstream distance for enhancers from transcription start sites

###### ddist
	Set minimum downstream distance for enhancers from transcription start sites

###### elmerENH
	T: Use enhancer locations from ELMER R package (T when you want to use ELMER enhancers and external enhancer datasets)
	
	F: Don't use enhancer locations from ELMER R package

###### encodeNDR
	T: Use ENCODE nucleosome depleted regions (obtained from ENCODE Master DNaseI-seq peaks from 125 tissues or cell lines)
	
	F: Don't use ENCODE Master DNaseI-seq peaks

###### extENH
	T: Use external enhancer datasets you loaded in ./external.data/enhancer
	
	F: Don't use external enhancer datasets

###### extNDR
	T: Use external nucleosome depleted regions you loaded in ./external.data/NDR
	
	F: Don't use external NDR datasets

###### onlyextFeature
	T: Use only external datasets for other peaks user has loaded in ./external.data/feature
	
	F: Don't use external datasets for other peaks

###### SMdataset
	T: Use somatic mutation datasets user has loaded in ./external.data/otherinfo to make complex scatterplots
	
	F: Don't use SMdatasets

###### CNVdataset
	T: Use copy number variation datasets user has loaded in ./external.data/otherinfo to make complex scatterplots
	
	F: Don't use CNVdatasets

###### usecaseonly
	T: Use methylation and gene expression data only from cases to select enhancer-gene lnks in step2
	
	F: Use methylation and gene expression data from cases and controls to select enhancer-gene links in step2

###### findhypo
	T: Find hypomethylated enhancer-gene links in step2
	
	F: Don't find hypomethylated enhancer-gene links in step2

###### hypo.strigency
	Set a minimum beta-value for hypomethylated tumors (ranges from 0 to 1)

###### findhypoGpos
	T: Find hypomethylated enhancer to positively linked gene throughout steps
	
	F: Don't find findhypoGpos

###### findhypoGneg
	T: Find hypomethylated enhancer to negatively linked gene throughout steps
	
	F: Don't find findhypoGneg

###### findhyper
	T: Find hypermethylated enhancer-gene links in step2
	
	F: Don't find hypermethylated enhancer-gene links in step2

###### hyper.stringency
	Set a maximum beta-value for hypomethylated tumors (ranges from 0 to 1)

###### findhyperGpos
	T: Find hypermethylated enhancer to positvely linked gene throughout steps
	
	F: Don't find findhyperGpos

###### findhyperGneg
	T: Find hypermethylated enhancer to negatively linked gene throughout steps
	
	F: Don't find findhyperGneg

###### Zcutoff
	Set a p-value cutoff for Z score calculation in step2 (ranges from 0 to 1)

###### permutation.cutoff
	Set an empirical p-value cutoff from permutation tests in step3 (ranges from 0 to 1)

###### adj.pval.cutoff
	Set an adjusted p-value cutoff from Wilcoxon rank sum tests in step4 (ranges from 0 to 1)

###### hypoGposHistogram
	T: Make histograms and tables for hypomethylated enhancers to genes positively linked (hypoGpos)
	
	F: Don't generate histograms and tables for hypoGpos

###### hypoGnegHistogram
	T: Make histograms and tables for hypomethylated enhancers to genes negatively linked (hypoGneg)
	
	F: Don't generate histograms and tables for hypoGneg

###### hyperGposHistogram
	T: Make histograms and tables for hypermethylated enhancers to genes positively linked (hyperGpos)
	
	F: Don't generate histograms and tables for hyperGpos

###### hyperGnegHistogram
	T: Make histograms and tables for hypermethylated enhancers to genes negatively linked (hyperGneg)
	
	F: Don't generate histograms and tables for hyperGneg

###### histcol
	Determine the color of histograms you generate (e.g. for red, histcol=c("red"), for blue, histcol=c("blue"))

###### hypoGposScatter
	T: Make scatterplots displaying the expression of (number specified by the scatterplot_top_n_genes parameter) transcriptional regulators to the methylation level of all positively-linked hypomethylated enhancers (hypoGpos). Allows users to also include transcriptional regulators as listed in a file containing the name 'genes_of_interest' placed in ./external.data/otherinfo/ 
	
	F: top transcriptional regulators and postively-linked hypomethylated enhancers or those specified by the user

###### hypoGnegScatter
	T: Make scatterplots displaying the expression of (number specified by the scatterplot_top_n_genes parameter) transcriptional regulators to the methylation level of all negatively-linked hypomethylated enhancers (hypoGneg). Allows users to also include transcriptional regulators as listed in a file containing the name 'genes_of_interest' placed in ./external.data/otherinfo/ 
	
	F: top transcriptional regulators and negatively-linked hypomethylated enhancers or those specified by the user

###### hyperGposScatter
	T: Make scatterplots displaying the expression of (number specified by the scatterplot_top_n_genes parameter) transcriptional regulators to the methylation level of all positively-linked hypermethylated enhancers (hyperGpos). Allows users to also include transcriptional regulators as listed in a file containing the name 'genes_of_interest' placed in ./external.data/otherinfo/ 
	
	F: top transcriptional regulators and postively-linked hypermethylated enhancers or those specified by the user

###### hyperGnegScatter
	T: Make scatterplots displaying the expression of (number specified by the scatterplot_top_n_genes parameter) transcriptional regulators to the methylation level of all negatively-linked hypermethylated enhancers (hyperGneg). Allows users to also include transcriptional regulators as listed in a file containing the name 'genes_of_interest' placed in ./external.data/otherinfo/ 
	
	F: top transcriptional regulators and negatively-linked hypermethylated enhancers or those specified by the user
	
###### scatterplot_top_n_genes
	Parameter sets the top n transcription factors (by number of linked/negatively-linked enhancer probes of the specified type) to create scatterplots for the top n transcriptional regulators to linked probes. 
	
###### hypoGposCScatter
	T: Make complex scatterplots for hypomethylated enhancers to genes positively linked (hypoGpos) using files in ./external.data/otherinfo (Purity, CNV, SM datasets)

	F: Don't generate complex scatterplots for hypoGpos

###### hypoGnegCScatter
	T: Make complex scatterplots for hypomethylated enhancers to genes negatively linked (hypoGneg) using files in ./external.data/otherinfo (Purity, CNV, SM datasets)
	
	F: Don't generate complex scatterplots for hypoGneg

###### hyperGposCScatter
	T: Make complex scatterplots for hypermethylated enhancers to genes positively linked (hyperGpos) using files in ./external.data/otherinfo (Purity, CNV, SM datasets)
	
	F: Don't generate complex scatterplots for hyperGpos

###### hyperGnegCScatter
	T: Make complex scatterplots for hypermethylated enhancers to genes negatively linked (hyperGneg) using files in ./external.data/otherinfo (Purity, CNV, SM datasets)
	
	F: Don't generate complex scatterplots for hyperGneg

###### complexscatterplot_top_n_genes
	Parameter sets the top n transcription factors (by number of linked/negatively-linked enhancer probes of the specified type) to create complex scatterplots (Purity, CNV, SM datasets) for the top n transcriptional regulators to linked probes. 
	
###### makeScatter4probe
	T: Make simple scatterplots for the enhancer probes users have selected and listed in a file containing the name 'probes_of_interest' placed in ./external.data/otherinfo/ 

	F: Don't generate simple scatterplots for the probes specified by the user in ./external.data/otherinfo

###### hypoGposTracks
	T: Make UCSC genome browser tracks for (number specified by the track_top_n_genes parameter) transcriptional regulators and their linked hypomethylated enhancers (hypoGpos)
	
	F: Don't generate genome browser tracks for top transcriptional regulators and their linked hypomethylated enhancers

###### hypoGnegTracks
	T: Make UCSC genome browser tracks for (number specified by the track_top_n_genes parameter) transcriptional regulators and their negatively-linked hypomethylated enhancers (hypoGneg)
	
	F: Don't generate genome browser tracks for top transcriptional regulators and their negatively-linked hypomethylated enhancers

###### hyperGposTracks
	T: Make UCSC genome browser tracks for (number specified by the track_top_n_genes parameter) transcriptional regulators and their linked hypermethylated enhancers (hyperGneg)
	
	F: Don't generate genome browser tracks for top transcriptional regulators and their linked hypermethylated enhancers

###### hyperGnegTracks
	T: Make UCSC genome browser tracks for (number specified by the track_top_n_genes parameter) transcriptional regulators and their negatively-linked hypermethylated enhancers (hyperGneg)
	
	F: Don't generate genome browser tracks for top transcriptional regulators and their negatively-linked hypermethylated enhancers
	
###### track_top_n_genes
	Parameter sets the top n transcription factors (by number of linked/negatively-linked enhancer probes of the specified type) to create UCSC genome browser tracks for above

###### hypoGposStates
	T: Make a table that lists the activity states of (number specified by the states_top_n_genes parameter) transcriptional regulators and their linked hypomethylated enhancers (hypoGpos)
	
	F: Don't generate the activity status table for top transcriptional regulators and their linked hypomethylated enhancers
	
###### hypoGnegStates
	T: Make a table that lists the activity states of (number specified by the states_top_n_genes parameter) transcriptional regulators and their negatively-linked hypomethylated enhancers (hypoGneg)
	
	F: Don't generate the activity status table for top transcriptional regulators and their negatively-linked hypomethylated enhancers

###### hyperGposStates
	T: Make a table that lists the activity states of (number specified by the states_top_n_genes parameter) transcriptional regulators and their linked hypermethylated enhancers (hyperGpos)
	
	F: Don't generate the activity status table for top transcriptional regulators and their linked hypermethylated enhancers
	
###### hyperGnegStates
	T: Make a table that lists the activity states of (number specified by the states_top_n_genes parameter) transcriptional regulators and their negatively-linked hypermethylated enhancers (hyperGneg)
	
	F: Don't generate the activity status table for top transcriptional regulators and their negatively-linked hypermethylated enhancers

###### states_top_n_genes
	Parameter sets the top n transcription factors (by number of linked/negatively-linked enhancer probes of the specified type) to create the states table for the top n transcriptional regulators to linked probes.

###### hypoGposSurvival
	T: Make Kaplan-Meier survival curves for (number specified by the survival_top_n_genes parameter) transcriptional regulators and their nominally-significant (uncorrected p<0.05) linked hypomethylated enhancers (hypoGpos) using percentile cutoffs determined by high_thresh and low_thresh parameters, respectively. Note: Using this requires that hypoGposHistogram is also set to "T"
	
	F: Don't make Kaplan-Meier survival curves for top transcriptional regulators and linked hypomethylated enhancers.

###### hypoGnegSurvival
	T: Make Kaplan-Meier survival curves for (number specified by the survival_top_n_genes parameter) transcriptional regulators and their nominally-significant (uncorrected p<0.05) negatively-linked hypomethylated enhancers (hypoGneg) using percentile cutoffs determined by high_thresh and low_thresh parameters, respectively. Note: Using this requires that hypoGnegHistogram is also set to "T"
	
	F: Don't make Kaplan-Meier survival curves for top transcriptional regulators and negatively-linked hypomethylated enhancers.

###### hyperGposSurvival
	T: Make Kaplan-Meier survival curves for (number specified by the survival_top_n_genes parameter) transcriptional regulators and their nominally-significant (uncorrected p<0.05) linked hypermethylated enhancers (hyperGpos) using percentile cutoffs determined by high_thresh and low_thresh parameters, respectively. Note: Using this requires that hyperGposHistogram is also set to "T"
	
	F: Don't make Kaplan-Meier survival curves for top transcriptional regulators and linked hypermethylated enhancers.

###### hyperGnegSurvival
	T: Make Kaplan-Meier survival curves for (number specified by the survival_top_n_genes parameter) transcriptional regulators and their nominally-significant (uncorrected p<0.05) negatively-linked hypermethylated enhancers (hyperGneg) using percentile cutoffs determined by high_thresh and low_thresh parameters, respectively. Note: Using this requires that hyperGnegHistogram is also set to "T"
	
	F: Don't make Kaplan-Meier survival curves for top transcriptional regulators and negatively-linked hypermethylated enhancers.

###### survival_top_n_genes
	Parameter sets the top n transcription factors (by number of linked/negatively-linked enhancer probes of the specified type) to create survival plots for above
	
###### visualize_survival_plots_genes
	T: Generate the survival plots for the top n genes specified above (along with a .tsv file containing survival information for the top n genes)
	F: Don't generate the survival plots for the top n genes specified above (still generate the .tsv file containing survival information for the top n genes)
	
###### visualize_survival_plots_probes
	T: Generate the survival plots for the unique probes linked to at least one of the top n genes specified above (along with a .tsv file containing survival information for these probes)
	F: Don't generate the survival plots for the unique probes linked to at least one of the top n genes specified above (still generate the .tsv file containing survival information for these probes)

###### high_thresh
	Set a threshold for proportion of samples to include in the high expression/methylation group (range from 0 to 1, and should be greater than low_thresh to prevent samples from appearing in both groups)

###### low_thresh
	Set a threshold for proportion of samples to include in the low expression group (range from 0 to 1, and should be less than high_thresh to prevent samples from appearing in both groups)

###### hypoGposCircos
	T: Make circos plots for (number specified by the circos_top_n_genes parameter) transcriptional regulators displaying relationship to their linked hypomethylated enhancers (hypoGpos). Note: Using this requires that hypoGposHistogram is also set to "T"
	
	F: Don't make circos plots for transcriptional regulators and their relationship to their linked hypomethylated enhancers

###### hypoGnegCircos
	T: Make circos plots for (number specified by the circos_top_n_genes parameter) transcriptional regulators displaying relationship to their negatively-linked hypomethylated enhancers (hypoGneg). Note: Using this requires that hypoGnegHistogram is also set to "T"
	
	F: Don't make circos plots for transcriptional regulators and their relationship to their negatively-linked hypomethylated enhancers

###### hypoGposCircos
	T: Make circos plots for (number specified by the circos_top_n_genes parameter) transcriptional regulators displaying relationship to their linked hypermethylated enhancers (hyperGpos). Note: Using this requires that hyperGposHistogram is also set to "T"
	
	F: Don't make circos plots for transcriptional regulators and their relationship to their linked hypermethylated enhancers

###### hypoGnegCircos
	T: Make circos plots for (number specified by the circos_top_n_genes parameter) transcriptional regulators displaying relationship to their negatively-linked hypermethylated enhancers (hyperGneg). Note: Using this requires that hyperGnegHistogram is also set to "T"
	
	F: Don't make circos plots for transcriptional regulators and their relationship to their negatively-linked hypermethylated enhancers

###### circos_top_n_genes
	Parameter sets the top n transcription factors (by number of linked/negatively-linked enhancer probes of the specified type) to create circos plots for above
	
###### hypoGposTAD
	T: Make TAD tables for (number specified by the TAD_top_n_genes parameter) top transcriptional regulators in this category, listing probes positively linked to these regulators and the genes that lie within the same TADs as them based on user input files. Note: Using this requires that hypoGposHistogram is also set to "T"
	
	F: Don't make TAD tables for these transcriptional regulators and their listed probes
	
###### hypoGnegTAD
	T: Make TAD tables for (number specified by the TAD_top_n_genes parameter) top transcriptional regulators in this category, listing probes negatively linked to these regulators and the genes that lie within the same TADs as them based on user input files. Note: Using this requires that hypoGnegHistogram is also set to "T"
	
	F: Don't make TAD tables for these transcriptional regulators and their listed probes
	
###### hyperGposTAD
	T: Make TAD tables for (number specified by the TAD_top_n_genes parameter) top transcriptional regulators in this category, listing probes positively linked to these regulators and the genes that lie within the same TADs as them based on user input files. Note: Using this requires that hyperGposHistogram is also set to "T"
	
	F: Don't make TAD tables for these transcriptional regulators and their listed probes
	
###### hyperGnegTAD
	T: Make TAD tables for (number specified by the TAD_top_n_genes parameter) top transcriptional regulators in this category, listing probes negatively linked to these regulators and the genes that lie within the same TADs as them based on user input files. Note: Using this requires that hyperGnegHistogram is also set to "T"
	
	F: Don't make TAD tables for these transcriptional regulators and their listed probes
	
###### TAD_top_n_genes
	Parameter sets the top n transcription factors (by number of linked/negatively-linked enhancer probes of the specified type) to create TAD tables for above
	
###### hypoGposMetHeatmap
	T: Make methylation heatmaps of probes linked to, and gene expression of the top (number specified by the probe_heatmap_top_n_genes parameter) transcriptional regulators in this category, using unique probes positively linked to these regulators based on user input files. Note: Using this requires that hypoGposHistogram is also set to "T"
	
	F: Don't make methylation heatmaps for these transcriptional regulators and their listed probes
	
###### hypoGnegMetHeatmap
	T: Make methylation heatmaps of probes linked to, and gene expression of the top (number specified by the probe_heatmap_top_n_genes parameter) transcriptional regulators in this category, using unique probes negatively linked to these regulators based on user input files. Note: Using this requires that hypoGnegHistogram is also set to "T"
	
	F: Don't make methylation heatmaps for these transcriptional regulators and their listed probes
	
###### hyperGposMetHeatmap
	T: Make methylation heatmaps of probes linked to, and gene expression of the top (number specified by the probe_heatmap_top_n_genes parameter) transcriptional regulators in this category, using unique probes positively linked to these regulators based on user input files. Note: Using this requires that hyperGposHistogram is also set to "T"
	
	F: Don't make methylation heatmaps for these transcriptional regulators and their listed probes
	
###### hyperGnegMetHeatmap
	T: Make methylation heatmaps of probes linked to, and gene expression of the top (number specified by the probe_heatmap_top_n_genes parameter) transcriptional regulators in this category, using unique probes negatively linked to these regulators based on user input files. Note: Using this requires that hyperGnegHistogram is also set to "T"
	
	F: Don't make methylation heatmaps for these transcriptional regulators and their listed probes
	
###### probe_heatmap_top_n_genes
	Parameter sets the top n transcription factors (by number of linked/negatively-linked enhancer probes of the specified type) to create methylation heatmaps for above

##### If you would like to run TENET newly with different settings, we recommend you to move output folders from the previous run (e.g. ./step3/ ./step4/) somewhere else or rename them in order to avoid bringing wrong files in a new run.
