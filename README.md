## Supporting code for "Molecular phenotyping reveals the identity of Barrett’s esophagus and its malignant transition" manuscript

This repository contains all steps required for analysis and production of figures as presented in the "Deep molecular phenotyping reveals the identity of Barrett’s esophagus and its malignant transition" manuscript. Scripts are split into the following folders:

#### Preprocessing

This folder contains low level steps using in the analysis of data. The analysis follows a stepwise process:

##### 1.4_Quality_control.Rmd

Contains script used for the reading of raw count data output from **Cell Ranger v3.0.1 _count_** function. The reads were aligned to GRCh38.92 gencode annotation. Quality control includes assessment of proportion of reads mapping to mitochondrial genes, removal of cells with low read and feature count and high proportion of reads mapped to mitochondrial RNA. The selection criteria are described in **Table S1** of the manuscript.

##### 2.4_Normalization_filtered.Rmd

Read count normalisation within individual samples

##### 3.4_Further_processing_filtered.Rmd

Clustering of cells within individual samples.

##### 4.4_Tissue_correction_filtered.Rmd

Batch correction and clustering of individual cells with individual tissue types

##### 5.4_Reclustering.Rmd

In the case of NSCJ we identified novel cell type. This cell type seem to be comprised of four additional cell types. Code associated with this reclustering is here. 
Further, in the BSCJ NOUROG3 and Goblet cells were not separated into two clusters in step 4.4. Additional reclustering is performed here.

##### 6.4_Full_data_batch_correction_filtered.Rmd

This code includes full data integration across all tissue types and samples. 

##### 7.4_Cluster_annotation.Rmd

This code introduces manually annoted cell types into SingleCellExperiment containing all high quality filtered cell types. 

#### Figures

All code used for generation of figures within the manuscript are located within that folder. The name of hte files corresponds to individual figures. Please note that in additional to the images used within the manuscripts, code also contains information for generation of diagnostic figures within each analysis. 

#### Analysis/Functions

Functions associated with the analysis are located in the _auxiliary.R_ file. 


#### Analysis/Battenberg

Contains steps used in the low level analysis of WGS samples using Battenberg and DPClust
