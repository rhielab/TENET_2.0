# TENET (Tracing Enhancer Networks using Epigenetic Traits) 2.0 Vignettes

#### Last updated: 9/3/2020

## TENET 2.0 Input Data Structures:

##### Gene Expression Datasets:

Should contain two objects with expression data for normal samples (or control group) in one, and expression data for the tumor samples (or experimental group) in the other. These objects can either be text files or as .rda files containing data frames named expDataN (for the “normal” expression data) and expDataT (for the “tumor” expression data). In either case, both objects should have the sample names listed as the column names in the form of unique ids such as TCGA barcodes (for ease of use we used 19 character sample names, up through the “Portion” element of the barcode - see: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/#:~:text=TCGA%20barcodes%20were%20used%20to,metadata%20values%20for%20a%20sample ). The row names should contain the Ensembl IDs for each of the genes, which should all be included in the Human GENCODE v22 release (see: https://www.gencodegenes.org/human/release_22.html). For recent publication, Upper quartile normalized FPKM values were used for gene expression quantification. Below shows example gene expression datasets (e.g. expDataN, expDataT) 

```diff
> head(expDataN,2)
                TCGA-50-5939-11A-01 TCGA-50-5933-11A-01 TCGA-44-5645-11A-01
ENSG00000000003            16.89178            17.17372           17.098905
ENSG00000000005            11.26400             0.00000            9.837402
                TCGA-44-6778-11A-01 TCGA-44-6144-11A-01 TCGA-49-6745-11A-01
ENSG00000000003            17.05940            16.47293            18.38691
ENSG00000000005            10.12427             0.00000            10.36832
                TCGA-50-5931-11A-01 TCGA-73-4676-11A-01 TCGA-50-5932-11A-01
ENSG00000000003           16.593883           17.360160           16.785994
ENSG00000000005            8.911002            9.170199            9.002069
                TCGA-50-5936-11A-01 TCGA-44-2668-11A-01 TCGA-44-2665-11A-01
ENSG00000000003            16.72215           16.458058           17.108430
ENSG00000000005             0.00000            9.006524            8.828383

> head(expDataT,2)
                TCGA-05-4384-01A-01 TCGA-05-4390-01A-02 TCGA-05-4396-01A-21
ENSG00000000003            17.50848            18.44378            17.59512
ENSG00000000005             0.00000             0.00000             0.00000
                TCGA-05-4405-01A-21 TCGA-05-4410-01A-21 TCGA-05-4415-01A-22
ENSG00000000003            18.50671           18.138609            18.12105
ENSG00000000005             0.00000            9.995623             0.00000
                TCGA-05-4417-01A-22 TCGA-05-4424-01A-22 TCGA-05-4425-01A-01
ENSG00000000003            18.41143            18.51993            17.64216
ENSG00000000005            10.32602             0.00000             0.00000
                TCGA-05-4427-01A-21 TCGA-05-4433-01A-22 TCGA-05-5420-01A-01
ENSG00000000003            17.36598            17.11916            19.54976
ENSG00000000005            14.81771             0.00000             0.00000
                TCGA-05-5423-01A-01 TCGA-05-5425-01A-02 TCGA-05-5428-01A-01
ENSG00000000003            19.47698            18.26736            18.69031
ENSG00000000005            12.42247             0.00000             0.00000
```

##### DNA methylation Datasets:

Similar to expression datasets, should contain two objects with methylation data for normal samples (or control group) in one, and DNA methylation data for the tumor samples (or experimental group) in the other. These objects can either be text files or as .rda files containing data frames named metDataN (for the “normal” DNA methylation data) and metDataT (for the “tumor” DNA methylation data). In either case, both objects should have the sample names listed as the column names in the form of unique ids such as TCGA barcodes (for ease of use we used 19 character sample names, up through the “Portion” element of the barcode - see: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/#:~:text=TCGA%20barcodes%20were%20used%20to,metadata%20values%20for%20a%20sample). The row names should contain the probe IDs (cg########) for each methylation probe. For Mullen et al (Plos Genet 2020), beta values were used for DNA methylation probe quantification. Below shows example DNA methylation datasets (e.g. metDataN, metDataT) 

```diff
> head(metDataN,2)
           TCGA-50-5939-11A-01 TCGA-50-5933-11A-01 TCGA-44-5645-11A-01
cg00000029           0.2434987           0.3032385           0.3042853
cg00000108                  NA                  NA                  NA
           TCGA-44-6778-11A-01 TCGA-44-6144-11A-01 TCGA-49-6745-11A-01
cg00000029           0.3704211           0.3500071           0.3124818
cg00000108                  NA                  NA                  NA
           TCGA-50-5931-11A-01 TCGA-73-4676-11A-01 TCGA-50-5932-11A-01
cg00000029           0.3177163           0.2903119            0.304546
cg00000108                  NA                  NA                  NA
           TCGA-50-5936-11A-01 TCGA-44-2668-11A-01 TCGA-44-2665-11A-01
cg00000029           0.2715483           0.3234379           0.2117347
cg00000108                  NA                  NA                  NA
           TCGA-50-5930-11A-01 TCGA-44-6147-11A-01 TCGA-38-4632-11A-01
cg00000029           0.3781955           0.2741015           0.4017598
cg00000108                  NA                  NA                  NA

> head(metDataT,2)
           TCGA-05-4384-01A-01 TCGA-05-4390-01A-02 TCGA-05-4396-01A-21
cg00000029           0.2714587           0.2299041           0.1785486
cg00000108                  NA                  NA                  NA
           TCGA-05-4405-01A-21 TCGA-05-4410-01A-21 TCGA-05-4415-01A-22
cg00000029            0.421797           0.3563369           0.2667541
cg00000108                  NA                  NA                  NA
           TCGA-05-4417-01A-22 TCGA-05-4424-01A-22 TCGA-05-4425-01A-01
cg00000029           0.3059861           0.3242333           0.3653999
cg00000108                  NA                  NA                  NA
           TCGA-05-4427-01A-21 TCGA-05-4433-01A-22 TCGA-05-5420-01A-01
cg00000029           0.2743896           0.2733237           0.3958237
cg00000108                  NA                  NA                  NA
           TCGA-05-5423-01A-01 TCGA-05-5425-01A-02 TCGA-05-5428-01A-01
cg00000029           0.2139236           0.2635692           0.1836697
cg00000108                  NA                  NA                  NA

```

DNA methylation and expression data should be matched to have same number of samples, and each id should match one sample in both datasets. 

##### Clinical Dataset:

The clinical dataset is used only for the survival functions in step 5 of TENET. Currently, clinical data must be stored as an .rda object in either './external.data/data/clinical' or as as part of a combined .rda file with the DNA methylation and gene expression data. In either case, the clinical data should be contained as a data frame named clinical_data. The clinical data object should contain four specifically named columns: "bcr_patient_barcode", "days_to_death", "days_to_last_followup", and "vital_status". "bcr_patient_barcode" should contain unique ids such as TCGA barcodes matching the DNA methylation and gene expression data names through the first 12 characters (up through the "Participant" element of the barcode - see: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/). "days_to_death" should contain the number of days a given patient survived after the sample was collected (for patients who are still alive, this should be 'NA'). "days_to_last_followup" should contain the number of days a given patient has survived since the data was collected (for patients who have died, this should be 'NA'). Finally "vital_status" should note whether a given patient is "Dead" or "Alive". Below shows an example clinical dataset (e.g. clinical) 

```diff
> head(clinical,4)
             bcr_patient_barcode days_to_death days_to_last_followup
TCGA-05-4244        TCGA-05-4244            NA                     0
TCGA-05-4245        TCGA-05-4245            NA                   395
TCGA-05-4249        TCGA-05-4249            NA                  1158
TCGA-05-4250        TCGA-05-4250           121                    NA
             vital_status
TCGA-05-4244        Alive
TCGA-05-4245        Alive
TCGA-05-4249        Alive
TCGA-05-4250         Dead

```

##### Enhancer and Nucleosome Depleted Regions (NDR) Datasets:

User might also choose to supply information for potential enhancer regions of interest (e.g. H3K27ac or H3K4me1 ChIP-seq peaks), as well as nucleosome depleted regions (NDRs) (e.g. DNaseI hypersensitive sites, ATAC-seq, FAIRE-seq peaks) from relevant cell types to their experiment. Enhancer datasets should be placed in  the './external.data/enhancer/' directory, while NDR datasets should be placed in './external.data/NDR/'. Both enhancer and NDR datasets should be tab-delimited, bed-like files with chromosome information in the first column, listed as 'chr#', the second column containing start positions for the peaks, and the third column containing the end positions. Other columns can be included but will not be considered in the TENET 2.0 analysis. Row and column names should not be included in the files. Enhancer and NDR datasets should be aligned to the human hg38 reference genome, as TENET 2.0 databases are based off that reference genome. Below shows an example enhancer/NDR dataset

```diff
chr1	777850	778178
chr1	778359	779442
chr1	779553	780185
chr1	817095	817495
chr1	826649	826989
chr1	827249	827622

```

## TENET 2.0 Output and recommended functions:

We higly recommended to run TENET 2.0 on a high performance computing cluster as it requires large computer memory. We recommend running all four "quadrants" of analysis (HypoGplus, HypoGminus, HyperGplus, and HyperGminus) though users can choose. HypoGplus and HyperGplus analyses encapsulate the direct effects of potential oncogenes and tumor supressors, respectively. All five steps should be run to identify key transcriptional regulators (e.g. transcription factors) and enhancers; steps 1-5 are vital for generating data, but step 5 includes optional functions (i.e. table and histogram, scatterplot, genome browser track, enhancer probe:gene link state, survival, circos, topologically associating domain (TAD) and heatmap functions). Please see Rhie et al (PMID: 27833659) and Mullen et al (In Press) to find example output figures and tables.

### Recommended functions:

For step 5 functions, it is recommended at the moment to run the histogram, survival, circos, TAD, and heatmap functions. 

##### Histogram output:

Histogram functions output a variety of both .txt and .pdf files containing summarized information about the number of links to relevant transcriptional regulators. The .txt files contain the total number of relevant linked DNA methylation probes to relevant transcription factors organized from the regulators with the largest number of links to the smallest. The .pdf files take the information from the corresponding .txt files and display that information in the form of a histogram, showing a given number of links on the x-axis, and the number of regulators with that number of links on the y-axis. For complete information about all links to identified key transcriptional regulators (e.g. transcription factors), look at the 'links.all.tf.freq' list and historgram. The cis links indicate ones located in the same chromosome ('links.cis.tf.freq') and the trans links indicate ones located in the different chromosomes ('links.trans.tf.freq').

This function is a prerequisite for several other important TENET functions in the step 5 analyses including the survival, circos, TAD, and heatmap functions. 

```diff
gene.nameID	Freq	geneSymbol	geneID
CENPA|ENSG00000115163	875	CENPA	ENSG00000115163
FOXM1|ENSG00000111206	845	FOXM1	ENSG00000111206
TCF24|ENSG00000261787	843	TCF24	ENSG00000261787
MYBL2|ENSG00000101057	840	MYBL2	ENSG00000101057
SOX2|ENSG00000181449	713	SOX2	ENSG00000181449

```
