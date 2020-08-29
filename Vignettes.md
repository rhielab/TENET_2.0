# TENET (Tracing Enhancer Networks using Epigenetic Traits) 2.0 Vignettes

#### Last updated: 8/28/2020

## TENET 2.0 Data Structures:

##### Expression Datasets:

Should contain two objects with expression data for normal samples (or control group) in one, and expression data for the tumor samples (or experimental group) in the other. These objects can either be text files or as .rda files containing data frames named expDataN (for the “normal” expression data) and expDataT (for the “tumor” expression data). In either case, both objects should have the TCGA sample names listed as the column names in the form of unique TCGA barcodes (for ease of use we used 19 character sample names, up through the “Portion” element of the barcode - see: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/#:~:text=TCGA%20barcodes%20were%20used%20to,metadata%20values%20for%20a%20sample ). The row names should contain the Ensembl IDs for each of the genes, which should all be included in the Human GENCODE v22 release (see: https://www.gencodegenes.org/human/release_22.html ). For recent publication, Upper quartile normalized FPKM values were used for gene expression quantification.

##### Methylation Datasets:

Similar to expression datasets, should contain two objects with methylation data for normal samples (or control group) in one, and expression data for the tumor samples (or experimental group) in the other. These objects can either be text files or as .rda files containing data frames named metDataN (for the “normal” expression data) and metDataT (for the “tumor” expression data). In either case, both objects should have the TCGA sample names listed as the column names in the form of unique TCGA barcodes (for ease of use we used 19 character sample names, up through the “Portion” element of the barcode - see: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/#:~:text=TCGA%20barcodes%20were%20used%20to,metadata%20values%20for%20a%20sample ). The row names should contain the probe IDs (cg########) for each methylation probe. For recent publication, beta values were used for DNA methylation probe quantification. 

Methylation and expression data should be matched, and each barcode should match one sample in both datasets. 

##### Clinical Dataset:

The clinical dataset is used only for the survival functions in step 5 of TENET. Currently, clinical data must be stored as an .rda object in either './external.data/data/clinical' or as as part of a combined .rda file with the methylation and expression data. In either case, the clinical data should be contained as a data frame named clinical_data. The clinical data object should contain four specifically named columns:     "bcr_patient_barcode", "days_to_death", "days_to_last_followup", and "vital_status". "bcr_patient_barcode" should contain TCGA barcodes matching the methylation and expression data names through the first 12 characters (up through the "Participant" element of the barcode - see: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/ ). "days_to_death" should contain the number of days a given patient survived after the sample was collected (for patients who are still alive, this should be 'NA'). "days_to_last_followup" should contain the number of days a given patient has survived since the data was collected (for patients who have died, this should be 'NA'). Finally "vital_status" should note whether a given patient is "Dead" or "Alive". 
