# TENET (Tracing Enhancer Networks using Epigenetic Traits) 2.0 Vignettes

#### Last updated: 10/18/2020

## TENET 2.0 Input Data Structures:

##### Gene Expression Datasets:

Should contain two objects with expression data for normal samples (or control group) in one, and expression data for the tumor samples (or experimental group) in the other. These objects can either be text files or as .rda files containing data frames named expDataN (for the “normal” expression data) and expDataT (for the “tumor” expression data). In either case, both objects should have the sample names listed as the column names in the form of unique ids such as TCGA barcodes (for ease of use we used 19 character sample names, up through the “Portion” element of the barcode - see: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/#:~:text=TCGA%20barcodes%20were%20used%20to,metadata%20values%20for%20a%20sample). The row names should contain the Ensembl IDs for each of the genes, which should all be included in the Human GENCODE v22 release (see: https://www.gencodegenes.org/human/release_22.html). For recent publication, Upper quartile normalized FPKM values were used for gene expression quantification. Below shows example gene expression datasets (e.g. expDataN, expDataT) 

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

DNA methylation and expression data should be matched to have the same number of samples, and each id should match one sample in both datasets. 

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

##### Genes of interest file:

For simple scatterplot functions, users have the option to include a genes of interest file to specify additional transcriptional regulators they would like to have scatterplots generated for. This file should contain "genes_of_interest" in the filename (such as 'genes_of_interest.txt'), and should be deposited in './external.data/otherinfo/'. This file should list either the gene names or their ensembl IDs out in a single column, with one gene listed per line, and without a header. Note that current functionality of the scatterplot functions will not create plots for genes if they do not have any links to them (scatterplots will not be created if the genes are listed incorrectly, are not accepted human transcriptional factors, or do not have any probes linked to them). 

```diff
ENSG00000197472
ENSG00000179774
ENSG00000152192
ENSG00000129173
ENSG00000204366

```

##### Probes of interest file:

For makeScatter4probe function, users have the option to include a probes of interest file to specify specific probes they would like to have scatterplots generated for. This file should contain "probes_of_interest" in the filename (such as 'probes_of_interest.txt'), and should be deposited in './external.data/otherinfo/'. This file should list the DNA methylation probes on the Illumina HumanMethylation450K beadchip array in a single column, with one probe listed per line, and without a header. Note that current functionality of the scatterplot functions will not create plots for probes if they do not have any links to them (scatterplots will not be created if the probes are listed incorrectly, are not included on the array, or do not have any transcriptional regulators linked to them). 

```diff
cg00068377
cg02202315
cg04855678
cg13873387
cg22953687

```

##### Copy Number Variation file:

For complex scatterplot function, users have the option to include copy number variation (CNV) called file to generated scatterplots with CNV annotation. This file should contain "CNV" in the filename (such as 'CNV.txt'), and should be deposited in './external.data/otherinfo/'. This file should include colnames as sample IDs that matched to above DNA methylation and gene expression datasets for tumor samples (case group) and rownames as gene names with CNV calls. CNV call should have deletion as negative value, amplification as positive value and no change as 0. An example CNV.txt file is shown here

```diff
       TCGA-YL-A8SQ-01 TCGA-YL-A8SR-01 TCGA-YL-A9WH-01 TCGA-YL-A9WI-01
ACAP3               -1               0              -1               0
ACTRT2              -1               0              -1               0

```

##### Somatic Mutation file: 

For complex scatterplot function, users have the option to include a somatic mutations (SM) called file to generate scatterplots with SM annotations. This file should contain "SM" in the filename (such as 'SM.txt'), and should be deposited in './external.data/otherinfo/'. This file should include as colnames the sample IDs that matched to above DNA methylation and gene expression datasets for tumor samples (case group) and as rownames the gene names with SM calls. Samples that have somatic mutation should be indicated as 1, and samples without mutation should be indicated as 0. An example SM.txt file is shown here

```diff
     TCGA-J4-A83N-01 TCGA-J9-A52B-01 TCGA-J9-A8CK-01 TCGA-J9-A8CL-01
TTN                0               0               0               1
TP53               0               1               1               0

```

##### Purity file:

For complex scatterplot function, users have the option to include a sample purity file to generate scatterplots with purity (e.g. tumor purity estimates comparing with normal molecular signatures) annotations. This file should contain "Purity" in the filename (such as 'Purity.txt'), and should be deposited in './external.data/otherinfo/'. This file should include the colnames; ID, and purity, and rownames with sample IDs that matched to above DNA methylation and gene expression datasets for tumor samples (case group). Purity value should be indicated from 0 to 1. An example Purity.txt file is shown here

```diff
             ID  purity
TCGA-2A-A8VL-01  0.510
TCGA-2A-A8VO-01  0.560
TCGA-2A-A8VT-01  0.725
TCGA-2A-A8VV-01  0.835
TCGA-2A-A8W1-01  0.945
TCGA-2A-A8W3-01  0.455

```

## TENET 2.0 Output and recommended functions:

We higly recommended to run TENET 2.0 on a high performance computing cluster as it requires large computer memory. We recommend running all four "quadrants" of analysis (HypoGplus, HypoGminus, HyperGplus, and HyperGminus) though users can choose. HypoGplus and HyperGplus analyses encapsulate the direct effects of potential oncogenes and tumor supressors, respectively. All five steps should be run to identify key transcriptional regulators (e.g. transcription factors) and enhancers; steps 1-5 are vital for generating data, but step 5 includes optional functions (i.e. table and histogram, scatterplot, genome browser track, enhancer probe:gene link state, survival, circos, topologically associating domain (TAD) and heatmap functions). Please see Rhie et al (PMID: 27833659) and Mullen et al (PMID: 32925947) to find example output figures and tables.

### Recommended functions:

For step 5 functions, it is recommended at the moment to run the histogram, genome browser track, states, scatterplot, survival, circos, TAD, and heatmap functions. 

##### Histogram output:

Histogram functions output a variety of both .txt and .pdf files containing summarized information about the number of links to relevant transcriptional regulators. The .txt files contain the total number of relevant linked DNA methylation probes to relevant transcription factors organized from the regulators with the largest number of links to the smallest. The .pdf files take the information from the corresponding .txt files and display that information in the form of a histogram, showing a given number of links on the x-axis, and the number of regulators with that number of links on the y-axis. For complete information about all links to identified key transcriptional regulators (e.g. transcription factors), look at the 'links.all.tf.freq' list and historgram. The cis links indicate ones located in the same chromosome ('links.cis.tf.freq') and the trans links indicate ones located in the different chromosomes ('links.trans.tf.freq').

This function is a prerequisite for several other important TENET functions in the step 5 analyses including the enome browser track, states, scatterplot, survival, circos, TAD, and heatmap functions. 

```diff
gene.nameID	Freq	geneSymbol	geneID
CENPA|ENSG00000115163	875	CENPA	ENSG00000115163
FOXM1|ENSG00000111206	845	FOXM1	ENSG00000111206
TCF24|ENSG00000261787	843	TCF24	ENSG00000261787
MYBL2|ENSG00000101057	840	MYBL2	ENSG00000101057
SOX2|ENSG00000181449	713	SOX2	ENSG00000181449

```

<img src="https://github.com/suhnrhie/TENET_2.0/blob/master/example_images/histogram_example.png" alt="Example histogram output" width="504"/>

##### Genome browser track output:

Genome browser track functions output interact files in .bed format showing interactions between the top n transcriptional regulators (as specified by the user) and their linked DNA methylation probes. These files can be uploaded by the user on the UCSC Genome Browser (https://genome.ucsc.edu/) for visualization on the hg38 human genome using the custom tracks option there. 

```diff
track type=interact name="TENET2.0_Hypo.G+_interactions" description="TENET2.0 top TR to enhancer DNA methylation probe links" 
chr2 26764288 26764290 ENSG00000115163_cg00019511_link 0 0 . #FF0000FF chr2 26764288 26764290 CENPA . chr14 100734950 100734953 cg00019511 .
chr2 26764288 26764290 ENSG00000115163_cg00030432_link 0 0 . #FF0000FF chr2 26764288 26764290 CENPA . chr7 100431350 100431353 cg00030432 .
chr2 26764288 26764290 ENSG00000115163_cg00071210_link 0 0 . #FF0000FF chr2 26764288 26764290 CENPA . chr6 148415805 148415808 cg00071210 .
chr2 26764289 29818472 ENSG00000115163_cg00073794_link 0 0 . #FF0000FF chr2 26764288 26764290 CENPA . chr2 29818471 29818474 cg00073794 .
chr2 26764288 26764290 ENSG00000115163_cg00093522_link 0 0 . #FF0000FF chr2 26764288 26764290 CENPA . chr1 2320579 2320582 cg00093522 .

```

<img src="https://github.com/suhnrhie/TENET_2.0/blob/master/example_images/tracks_example.png" alt="Example UCSC genome browser track output" width="504"/>

##### States output:
States functions output files that list the activation status of each enhancer probe-transcriptional regulator link. For example, below file generated by running TENET using lung adenocarcinoma data indicates that tumor sample, TCGA-05-4415-01A-22 has activated enhancer probe cg00019511 linked to FOXM1, CENPA, and TCF24 while other samples do not. Users can search the activities of top n transcriptional regulators and their linked enhancers by using the parameter, "states_top_n_genes" in the settings.txt file.

```diff
TCGA-05-4384-01A-01	TCGA-05-4390-01A-02	TCGA-05-4396-01A-21	TCGA-05-4405-01A-21	TCGA-05-4410-01A-21	TCGA-05-4415-01A-22
cg00002190_MYBL2	0	0	0	0	0	0
cg00019511_FOXM1	0	0	0	0	0	1
cg00019511_CENPA	0	0	0	0	0	1
cg00019511_TCF24	0	0	0	0	0	1
cg00022235_SOX2	0	0	0	0	0	0

```
##### Simple scatterplot output:

Scatterplot functions output .pdf files with scatterplots displaying the expression of transcriptional regulators on the x-axis, and the methylation of all their linked DNA methylation probes on the y-axis. Individual points are colored red (for tumor samples) or blue (for normal samples). Scatterplots can be generated for the top n transcriptional regulators and their linked enhancers by using the parameter, "scatterplot_top_n_genes" in the settings.txt file. Scatterplots can be generated by inputting genes as specified by the genes of interest file (see above). 

<img src="https://github.com/suhnrhie/TENET_2.0/blob/master/example_images/scatterplot_example.png" alt="Example scatterplot output" width="504"/>

##### Complex scatterplot output:

Complex Scatterplot output .pdf files with scatterplots displaying the expression of transcriptional regulators on the x-axis, and the methylation of all their linked DNA methylation probes on the y-axis. Individual points are colored red (for tumor (case) samples) or blue (for normal (control) samples). Scatterplots can be generated for the top n transcriptional regulators and their linked enhancers by using the parameter, "complexscatterplot_top_n_genes" in the settings.txt file. As the below example complex scatterplot shows, CNV and somatic mutation status is indicated by the shape, and sample purity is indicated by the size.

<img src="https://github.com/suhnrhie/TENET_2.0/blob/master/example_images/complexscatterplot_example.png" alt="Example complex scatterplot output" width="504"/> 

##### MakeScatter4probe:

This function is optional if the user has enhancer probes of particular interest and has supplied them in the Probes of interest (see above). Like the simple scatterplot functions, this function outputs .pdf files with scatterplots displaying the expression of transcriptional regulators on the x-axis, and the methylation of all their linked DNA methylation probes on the y-axis. Individual points are colored red (for tumor (case) samples) or blue (for normal (control) samples). See above for example of output for this function. 

##### Survival output:

Depending on user input, survival functions can output .tsv files containing in-depth survival information about the top n transcriptional regulators (as specified by the user) and the DNA methylation probes linked to any one of of these transcriptional regulators, as well as pdf files displaying the Kaplan-Meier survival curves for the given transcriptional regulators and their linked DNA methylation probes (if the user has selected this option). 

```diff
normal_sample_count	tumor_sample_count	normal_sample_count_missing	tumor_sample_count_missing	mean_normal_expression	mean_tumor_expression	normal_sample_with_clinical_count	tumor_sample_with_clinical_count	tumor_sample_with_clinical_NA_count	tumor_sample_with_clinical_low_count	tumor_sample_with_clinical_intermediate_count	tumor_sample_with_clinical_high_count	mean_tumor_with_clinical_low_expression	mean_tumor_with_clinical_intermediate_expression	mean_tumor_with_clinical_high_expression	proportion_dead_in_low_expression	proportion_dead_in_high_expression	survival_direction_of_effect	survival_p_value
CENPA	21	453	0	0	12.3325086264362	15.5937095161134	21	453	0	109	219	105	13.5368179935483	15.7439186181039	17.3927100762181	0.128440366972477	0.295238095238095	high_expression_low_survival	0.00334844476544105
FOXM1	21	453	0	0	13.9098560452808	16.8571631413882	21	453	0	108	221	104	15.0514790151095	16.9052547479288	18.5686071113993	0.138888888888889	0.336538461538462	high_expression_low_survival	0.000478435670498101
TCF24	21	453	0	0	4.80371882433539	9.08699422170773	21	453	0	108	216	109	2.26963702103498	10.4201103389594	13.3748235506378	0.185185185185185	0.284403669724771	high_expression_low_survival	0.340335620021469
MYBL2	21	453	0	0	14.2238015665364	18.0201935332658	21	453	0	109	215	109	15.7221232722669	18.2055691603018	19.9664164551512	0.110091743119266	0.302752293577982	high_expression_low_survival	0.00209558525594733
SOX2	21	453	0	0	13.3581835433044	15.9366883414569	21	453	0	109	214	110	12.0606586527843	16.3308839057577	19.0555956835861	0.211009174311927	0.181818181818182	low_expression_low_survival	0.0799571894772609

```

<img src="https://github.com/suhnrhie/TENET_2.0/blob/master/example_images/gene_expression_survival_example.png" alt="Example gene expression survival output" width="504"/>

<img src="https://github.com/suhnrhie/TENET_2.0/blob/master/example_images/DNA_methylation_probe_survival_example.png" alt="Example DNA methylation survival output" width="504"/>

##### Circos output:

Circos functions output .html files that contain circos plots demonstrating the links from the top transcriptional regulators (as specified by the user) to their linked DNA methylation probes using the BioCircos package. Currently the links are generated in red, and the plots can be saved as .pdfs by the user by opening the html files and printing the page that opens using the Print -> save as .pdf options

<img src="https://github.com/suhnrhie/TENET_2.0/blob/master/example_images/circos_example.png" alt="Example circos output" width="302"/>

##### TAD output:

TAD functions output .tsv files that contain TAD information for all unique DNA methylation probes linked to the top transcriptional regulators (as specified by the user). These files note which of the top transcriptional regulators each probe is linked to, and lists out all the genes by Ensembl ID and name within the TAD files supplied by the user (in './external.data/TAD/') in a comma-delimited list (these lists can be rather long so a sample is not noted here). 

```diff
probe_ID      seqnames  start     end       CENPA_linked FOXM1_linked  TCF24_linked MYBL2_linked  SOX2_linked  GM12878_Rao_2014.raw_TADs_gene_count_in_TAD  GM12878_Rao_2014.raw_TADs_TAD_gene_ENSGs 
cg00272484     chr8 142319159 142319160         TRUE         TRUE         FALSE       TRUE        FALSE       ENSG00000226807,ENSG00000226490,ENSG00000221123,ENSG00000254183,ENSG00000261710,ENSG00000265247,ENSG00000254008,ENSG00000171045,ENSG00000253602,ENSG00000261693,ENSG00000181790,ENSG00000232722         MROH5,AC138647.1,AC104417.1,RP11-953B20.2,RP11-953B20.1,MIR4472-1,LINC00051,TSNARE1,Metazoa_SRP,RP13-467H17.1,BAI1,MROH4P
cg00389036     chr6 170094374 170094375         TRUE        FALSE         FALSE       FALSE       FALSE       ENSG00000236173,ENSG00000232197,ENSG00000273100,ENSG00000218716,ENSG00000230960,ENSG00000271820,ENSG00000227508,ENSG00000198719,ENSG00000112584,ENSG00000271234,ENSG00000266245,ENSG00000261003,ENSG00000008018,ENSG00000112592,ENSG00000071994,ENSG00000217874         RP1-182D15.2,RP11-302L19.1,RP11-302L19.3,RPL12P23,RP5-1086L22.1,RP5-894D12.5,FLJ38122,DLL1,FAM120B,RP5-894D12.4,MIR4644,RP1-140C12.2,PSMB1,TBP,PDCD2,OR4F7P
```
##### Heatmap output:

Heatmap functions output .pdf files displaying a heatmap with samples in the columns and each unique DNA methylation probes linked to the top transcriptional regulators (as specified by the user). Above the heatmap bars displaying relative expression of each of top transcriptional regulator genes (as specified by the user) are shown (scaled using the function: [X - Xmin] /[Xmax - Xmin], with 0 expression values set to the minimum, non-zero value). Heatmaps are plotted using the heatmap.3 function with euclidean distance and Ward.D2 function for clustering. 

<img src="https://github.com/suhnrhie/TENET_2.0/blob/master/example_images/heatmap_example.png" alt="Example heatmap output" width="504"/>
