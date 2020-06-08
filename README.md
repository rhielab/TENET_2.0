# TENET (Tracing Enhancer Networks using Epigenetic Traits) 2.0

#### Last updated: 6/7/2020

## 1.TENET 2.0 Installation

##### This build of TENET 2.0 depends on the following programs

###### A. bedtools (please install from https://github.com/arq5x/bedtools2/releases)

###### B. ELMER (R package)

ELMER can be installed from Bioconductor by typing below in R:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ELMER")

###### C. Download the TENET 2.0 program from this page into directory of your choosing

## 2. TENET 2.0 Usage

###### A. Get methylation and gene expression data with matched sample IDs

You can either put txt files in ./external.data/data/methylation and ./external.data/data/expression (add case and ctrl in the file name, respectively) or put a rda file in ./external.data/data (for case, name as "metDataT", "expDataT", for ctrl, name as "metDataN", "expDataN")

###### B. Fill out parameter settings in settings.txt 

###### C. Run TENET by typing bash ./TENET.sh

(OPTIONAL 1) put other cells methylation data you want to use for purity adjustment
you can put a txt file in ./external.data/othercells or put a rda file in ./external.data/others (name as "extData")
see the example othercells data in http://farnhamlab.com/software

(OPTIONAL 2) if you have external datasets for enhancer peaks (bedfile), add files in ./external.data/enhancer
see the example bed file in http://farnhamlab.com/software

(OPTIONAL 3) if you have external datasets for NDR peaks (bedfile), add files in ./external.data/NDR
see the example bed file in http://farnhamlab.com/software

(OPTIONAL 4) if you have external datasets for other peaks (bed file), add files in ./external.data/feature
see the example bed file in http://farnhamlab.com/software

(OPTIONAL 5) if you have Purity, CNV, SM datasets for sample samples, add files in ./external.data/otherinfo (include "CNV", "SM", "Purity" in the file names)
see the example CNV file in http://farnhamlab.com/software

see the example SM file in http://farnhamlab.com/software

see the example Purity file in http://farnhamlab.com/software

## 3. TENET 2.0 Parameters Primer

###### prefix
	Prefix for output files 

###### cores
	Number of cores you want to use in your machine

###### step0
	T: Process data into a .rda file when you input a .txt file for methylation and expression data (when you put other cells' methylation data in a .txt file) 
	
	F: Use when you put a pre-created .rda file containing methylation and expression data

step1
	T: Find differentially methylated enhancer regions
	F: Do not run step1

step2
	T: Select enhancer-gene links by calculating Z scores
	F: Do not run step2

step3
	T: Select significant enhancer-gene links from step2 by permutation
	F: Do not run step3

step4
	T: Optimize selection of enhancer-gene lnks
	F: Do not run step4

step5
	T: Summarize and visualize enhancer-gene links
	F: Do not run step5

methcutoff
	Set the methylation cutoff of beta-values, range from 0 to 1

hypocutoff
	Set the hypomethylation cutoff of beta-values, range from 0 to 1

unmethcutoff
	Set the unmethylation cutoff of beta-values, range from 0 to 1

hypercutoff
	Set the hypermethylation cutoff of beta-values, range from 0 to 1

minTumors
	Set the minimum number of tumors with epigenetic events to determine links

leuk
	T: Use leukocyte methylation data to adjust purity
	F: Don't use

fibro
	T: Use fibroblast methylation data to adjust purity
	F: Don't use

sm
	T: Use smooth muscle methylation data to adjust purity
	F: Don't use

othercells
	T: Use other cells methylation data which you loaded in ./external.data/othercells to adjust purity (see Optional 1 above)
	F: Don't use

purityinfo
	T: Use purity estimates which you loaded in ./external.data/otherinfo to make complex scatterplots
	F: Don't use purityinfo

udist
	Set minimum upstream distance for enhancers from transcription start sites

ddist
	Set minimum downstream distance for enhancers from transcription start sites

elmerENH
	T: Use enhancer locations from ELMER R package (T when you want to use ELMER enhancers and external enhancer datasets)
	F: Don't use enhancer locations from ELMER R package

encodeNDR
	T: Use ENCODE nucleosome depleted regions (obtained from ENCODE Master DNaseI-seq peaks from 125 tissues or cell lines)
	F: Don't use ENCODE Master DNaseI-seq peaks

extENH
	T: Use external enhancer datasets you loaded in ./external.data/enhancer
	F: Don't use external enhancer datasets

extNDR
	T: Use external nucleosome depleted regions you loaded in ./external.data/NDR
	F: Don't use external NDR datasets

onlyextFeature
	T: Use only external datasets for other peaks user has loaded in ./external.data/feature
	F: Don't use external datasets for other peaks

SMdataset
        T: Use somatic mutation datasets purity user has loaded in ./external.data/otherinfo to make complex scatterplots
	F: Don't use SMdatasets

CNVdataset
	T: Use copy number variation datasets purity user has loaded in ./external.data/otherinfo to make complex scatterplots
	F: Don't use CNVdatasets

usecaseonly
	T: Use methylation and gene expression data only from cases to select enhancer-gene lnks in step2
	F: Use methylation and gene expression data from cases and controls to select enhancer-gene links in step2

findhypo
	T: Find hypomethylated enhancer-gene links in step2
	F: Don't find hypomethylated enhancer-gene links in step2

hypo.strigency
	Set a minimum beta-value for hypomethylated tumors (ranges from 0 to 1)

findhypoGpos
	T: Find hypomethylated enhancer to positively linked gene throughout steps
	F: Don't find findhypoGpos

findhypoGneg
	T: Find hypomethylated enhancer to negatively linked gene throughout steps
	F: Don't find findhypoGneg

findhyper
	T: Find hypermethylated enhancer-gene links in step2
	F: Don't find hypermethylated enhancer-gene links in step2

hyper.stringency
	Set a maximum beta-value for hypomethylated tumors (ranges from 0 to 1)

findhyperGpos
	T: Find hypermethylated enhancer to positvely linked gene throughout steps
	F: Don't find findhyperGpos

findhyperGneg
	T: Find hypermethylated enhancer to negatively linked gene throughout steps
	F: Don't find findhyperGneg

Zcutoff
	Set a p-value cutoff for Z score calculation in step2 (ranges from 0 to 1)

permutation.cutoff
	Set an empirical p-value cutoff from permutation tests in step3 (ranges from 0 to 1)

adj.pval.cutoff
	Set an adjusted p-value cutoff from Wilcoxon rank sum tests in step4 (ranges from 0 to 1)

hypoGposHistogram
	T: Make histograms and tables for hypomethylated enhancers to genes positively linked (hypoGpos)
	F: Don't generate histograms and tables for hypoGpos

hypoGnegHistogram
	T: Make histograms and tables for hypomethylated enhancers to genes negatively linked (hypoGneg)
	F: Don't generate histograms and tables for hypoGneg

hyperGposHistogram
	T: Make histograms and tables for hypermethylated enhancers to genes positively linked (hyperGpos)
	F: Don't generate histograms and tables for hyperGpos

hyperGnegHistogram
	T: Make histograms and tables for hypermethylated enhancers to genes negatively linked (hyperGneg)
	F: Don't generate histograms and tables for hyperGneg

histcol
	Determine the color of histograms you generate (e.g. for red, histcol=c("red"), for blue, histcol=c("blue"))

hypoGposScatter
	T: Make simple scatterplots for hypomethylated enhancers to genes positively linked (hypoGpos)
	F: Don't generate simple scatterplots for hypoGpos

hypoGnegScatter
	T: Make simple scatterplots for hypomethylated enhancer to genes negatively linked (hypoGneg)
	F: Don't generate simple scatterplots for hypoGneg

hyperGposScatter
	T: Make simple scatterplots for hypermethylated enhancers to genes positively linked (hyperGpos)
	F: Don't generate simple scatterplots for hyperGpos

hyperGnegScatter
	T: Make simple scatterplots for hypermethylated enhancers to genes negatively linked (hyperGneg)
	F: Don't generate simple scatterplots for hyperGneg

hypoGposTracks
	T: Make genome browser tracks for hypomethylated enhancer to genes positively linked (hypoGpos)
	F: Don't generate genome browser tracks for hypoGpos

hypoGnegStates
	T: Make genome browser tracks for hypomethylated enhancers to genes negatively linked (hypoGneg)
	F: Don't generate genome browser tracks for hypoGneg

hyperGposStates
	T: Make genome browser tracks for hypermethylated enhancers to genes positively linked (hyperGpos)
	F: Don't generate genome browser tracks for hyperGpos

hyperGnegStates
	T: Make genome browser tracks for hypermethylated enhancers to genes negatively linked (hyperGneg)
	F: Don't generate genome browser tracks for hyperGneg

hypoGposCScatter
	T: Make complex scatterplots for hypomethylated enhancers to genes positively linked (hypoGpos) using files in ./external.data/otherinfo (Purity, CNV, SM datasets)
	F: Don't generate complex scatterplots for hypoGpos

hypoGnegCScatter
	T: Make complex scatterplots for hypomethylated enhancers to genes negatively linked (hypoGneg) using files in ./external.data/otherinfo (Purity, CNV, SM datasets)
	F: Don't generate complex scatterplots for hypoGneg

hyperGposCScatter
	T: Make complex scatterplots for hypermethylated enhancers to genes positively linked (hyperGpos) using files in ./external.data/otherinfo (Purity, CNV, SM datasets)
	F: Don't generate complex scatterplots for hyperGpos

hyperGnegCScatter
	T: Make complex scatterplots for hypermethylated enhancers to genes negatively linked (hyperGneg) using files in ./external.data/otherinfo (Purity, CNV, SM datasets)
	F: Don't generate complex scatterplots for hyperGneg

makeScatter4gene
	T: make simple scatterplots for the gene you selected in ./external.data/otherinfo - please make sure you put the list of genes with the correct file names (e.g. for hypomethylated enhancer to genes positively linked, "hypo.G+.gene.txt", for hypermethylated enhancer to genes negatively lnked, "hyper.G-.gene.txt")
	F: don't generate simple scatterplots for the gene you selected in ./external.data/otherinfo

makeScatter4gene
	T: Make simple scatterplots for the enhancer probes you selected in ./external.data/otherinfo - please make sure you put the list of probes with the correct file names (e.g. for hypomethylated enhancer to probes positively linked, "hypo.G+.probe.txt", for hypermethylated enhancer to probes negatively lnked, "hyper.G-.probe.txt")
	F: Don't generate simple scatterplots for the probes you selected in ./external.data/otherinfo

hypoGposSurvival
	T: Make Kaplan-Meier survival curves for (number specified by the survival_top_n_genes parameter) transcriptional regulators and their linked hypomethylated enhancers (hypoGpos) using percentile cutoffs determined by high_thresh and low_thresh parameters, respectively. 
	F: Don't make Kaplan-Meier survival curves for top transcriptional regulators and linked hypomethylated enhancers.

hypoGnegSurvival
	T: Make Kaplan-Meier survival curves for (number specified by the survival_top_n_genes parameter) transcriptional regulators and their negatively-linked hypomethylated enhancers (hypoGneg) using percentile cutoffs determined by high_thresh and low_thresh parameters, respectively. 
	F: Don't make Kaplan-Meier survival curves for top transcriptional regulators and negatively-linked hypomethylated enhancers.

hyperGposSurvival
	T: Make Kaplan-Meier survival curves for (number specified by the survival_top_n_genes parameter) transcriptional regulators and their linked hypermethylated enhancers (hyperGpos) using percentile cutoffs determined by high_thresh and low_thresh parameters, respectively. 
	F: Don't make Kaplan-Meier survival curves for top transcriptional regulators and linked hypermethylated enhancers.

hypoGposSurvival
	T: Make Kaplan-Meier survival curves for (number specified by the survival_top_n_genes parameter) transcriptional regulators and their negatively-linked hypermethylated enhancers (hyperGneg) using percentile cutoffs determined by high_thresh and low_thresh parameters, respectively. 
	F: Don't make Kaplan-Meier survival curves for top transcriptional regulators and negatively-linked hypermethylated enhancers.

survival_top_n_genes
	Parameter sets the top n transcription factors (by number of linked/negatively-linked enhancer probes of the specified type) to create survival plots for above

high_thresh
	Set a threshold for proportion of samples to include in the high expression/methylation group (range from 0 to 1, and should be greater than low_thresh to prevent samples from appearing in both groups)

low_thresh
	Set a threshold for proportion of samples to include in the low expression group (range from 0 to 1, and should be less than high_thresh to prevent samples from appearing in both groups)

hypoGposCircos
	T: Make circos plots for (number specified by the circos_top_n_genes parameter) transcriptional regulators displaying relationship to their linked hypomethylated enhancers (hypoGpos)
	F: Don't make circos plots for transcriptional regulators and their relationship to their linked hypomethylated enhancers

hypoGnegCircos
	T: Make circos plots for (number specified by the circos_top_n_genes parameter) transcriptional regulators displaying relationship to their negatively-linked hypomethylated enhancers (hypoGneg)
	F: Don't make circos plots for transcriptional regulators and their relationship to their negatively-linked hypomethylated enhancers

hypoGposCircos
	T: Make circos plots for (number specified by the circos_top_n_genes parameter) transcriptional regulators displaying relationship to their linked hypermethylated enhancers (hyperGpos)
	F: Don't make circos plots for transcriptional regulators and their relationship to their linked hypermethylated enhancers

hypoGnegCircos
	T: Make circos plots for (number specified by the circos_top_n_genes parameter) transcriptional regulators displaying relationship to their negatively-linked hypermethylated enhancers (hyperGpos)
	F: Don't make circos plots for transcriptional regulators and their relationship to their negatively-linked hypermethylated enhancers

circos_top_n_genes
	Parameter sets the top n transcription factors (by number of linked/negatively-linked enhancer probes of the specified type) to create circos plots for above

*If you would like to run TENET newly with different settings, we recommend you to move output folders from the previous run (e.g. ./step3/ ./step4/) somewhere else or rename them in order to avoid bringing wrong files in a new run.
