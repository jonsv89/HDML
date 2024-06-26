# Introduction
Co-occurrence of diseases decreases patient quality of life, complicates treat- ment choices, and increases mortality. Analyses of electronic health records present a complex scenario of comorbidity relationships that vary by age, sex, and cohort under study. The study of similarities between diseases using ‘omics data, such as genes altered in diseases, gene expression, proteome, and microbiome, are fundamental to uncovering the origin of, and potential treat- ment for, comorbidities. Recent studies have produced a first generation of ge- netic interpretations for as much as 46% of the comorbidities described in large cohorts. Integrating different sources of molecular information and using artificial intelligence (AI) methods are promising approaches for the study of co- morbidities. They may help to improve the treatment of comorbidities, including the potential repositioning of drugs (https://www.cell.com/action/showPdf?pii=S0168-9525%2823%2900134-8). The goal of this project is to build a disease-multilayer network integrating different sources of information (genes, miRNAs, drugs, symptoms, metabolome, methylome, microbiome, transcriptome) to better understand comorbidity relationships.

# Dictionary
Disease names, MeSH terms, and DOID have been translated into ICD10 and ICD9 codes (International Classification of Diseases) using the Unified Medical Language System (UMLS), followed by manual curation. The generated dictionary ("DiseaseName_ICD10_dictionary.txt") provides disease name - ICD10 associations indicating the disease category they belong to and provides a palette color.

# Multilayer construction
![Construction of the disease multilayer network](https://github.com/jonsv89/Disease_multilayer/blob/master/Images/LayerConstruction.png)

# Transcriptomic layer
Expression data has been downloaded from GEO, focusing only on those experiments conducted on human samples with the microarray platform HG U133Plus2. After quality filtering, normalization, and filtering of lowly expressed genes differential expression analyses have been conducted comparing expression levels in disease vs. control state. Cosine similarities between diseases were calculated based on the fold changes obtained by differential expression analysis. Using randomizations, the p-value of the similarities was calculated.

# Gene layer
Disease-gene associations have been extracted from DisGeNET (https://www.disgenet.org). We have selected only curated disease-gene associations. After translating disease names into ICD10 codes a weight has been asigned to each gene based on the number of diseases it's associated to. Then, cosine similarities have been calculated between pairs of diseases, calculating a p-value for each pair of diseases through randomizations.

# miRNA layer
Disease-miRNA associations have been extracted from HMDD (https://www.cuilab.cn/hmdd). These associations are extracted from publications and are manually curated. After translating disease names into ICD10 codes a weight has been asigned to each miRNA-ICD10 pair calculating a term-frequency inverse document frequency based on term co-ocurrences. Then, cosine similarities have been calculated between pairs of diseases, calculating a p-value for each pair of diseases through randomizations.

# Drug layer
Disease-drug associations have been extracted from CTD (https://ctdbase.org). We have focused on those drugs described to treat the disease, removing inferred associations. After translating disease names into ICD10 codes a weight has been asigned to each drug based on the number of diseases it's associated to. Then, cosine similarities have been calculated between pairs of diseases, calculating a p-value for each pair of diseases through randomizations.

# Symptom layer
Disease-symptom associations have been extracted from the Human Symptom Disease Network (https://www.nature.com/articles/ncomms5212). The publication by Zhou and colleagues provides in their supplementary material the number of co-ocurrences between diseases and symptoms extracted from pubmed metadata. Through chi-square calculations we have selected significant disease-symptom associations, and after translating disease names into ICD10 codes we have calculatd term-frequency inverse document frequencies based on the number of disease-symptom co-occurences. Then, cosine similarities have been calculated between pairs of diseases, calculating a p-value for each pair of diseases through randomizations.

# Microbiome layer
Disease-microbe associations have been extracted from DISBIOME (https://disbiome.ugent.be/home). We have selected only disease-microbe associations described in the gut. After translating disease names into ICD10 codes a weight has been asigned to each microbe based on the number of diseases it's associated to, giving it a positive or negative value depending on whether the levels of the microorganism are increased or decreased in the disease. Then, cosine similarities have been calculated between pairs of diseases, calculating a p-value for each pair of diseases through randomizations.

# Methylome layer
Disease-cg associations have been extracted from EWASdb (http://www.bioapp.org/ewasdb/). We have selected only disease-cg associations with p-values < 1e-3 and < 1e-7 obtained through Illumina HumanMethylation450 Bead Chip (450K) methylation arrays. disease names have been translated into ICD10 codes and cosine similarities have been calculated between them based on the t-values provided by the database, which are obtained through the use of the EWAS v2.0 software, calculating later a p-value for each pair of diseases through randomizations.

# Metabolome layer
Disease-metabolite associations have been extracted from the Human Metabolome Database (https://hmdb.ca). We have selected disease-metabolite associations, regardless of the concentrations, that were measured in blood. After translating disease names into ICD10 codes a weight has been asigned to each metabolite based on the number of diseases it's associated to. Then, cosine similarities have been calculated between pairs of diseases, calculating a p-value for each pair of diseases through randomizations.

# Diseases by layer
![Number of diseases in each layer](https://github.com/jonsv89/Disease_multilayer/blob/master/Images/DiseaseByLayer.png)
