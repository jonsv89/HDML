# - Expression layers

Expression layers have been generated using RNAseq and microarray data separately. 

- The microarray data has been extracted from:
  https://www.biorxiv.org/content/10.1101/431312v2.abstract

- The RNAseq data has been extracted from:
  https://amp.pharm.mssm.edu/archs4/

# Normalization
Microarray data is generated using the same Affymetrix microarray platform (HG-U133Plus2) to remove differences due to the use of different platforms. Expression data has been normalized using fRMA (PMID: 20097884). 

In the case of RNAseqdata, for each study separately, those genes with 0 expression in the 20% or more of the samples of the study have been removed. Transcript counts have been normalized using DESeq2's median of ratios approach, this is, "counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene".

# Differential expression analysis
- Is done in the same way for microarray and RNAseq data.
- In the cases where a single study analyzes a disease, ind.analysis from the MetaDE package has been used for identifying significantly differentially expressed genes.
 - In the cases where more than one study analyze a disease, ind.cal.Es and Meta.DE functions from the MetaDE package have been used for calculating effect sizes and conducting the meta-analysis. As in our previous meta-analyses (PMID: 28667284), we have used the Random Effects Model (REM). 
 
 # Network generation
 - Using fisher exact tests, we have calculated significant overlaps between significantly differentially up- and down-regulated genes.
   - Two diseases are positively connected if there is a significant overlap between the genes significantly de-regulated in the same direction in both diseases (up-regulated in both diseases and also down-regulated in both diseases).
   - Two diseases are negatively connected if there is a significant overlap between the genes significantly de-regulated in the opposite direction in both diseases (up-regulated in one disease and down-regulated in the other one, and vice versa).
