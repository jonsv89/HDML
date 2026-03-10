#### With this code we are going to generate the multiomic layers (disease-disease similarities) for the BSC-IGC collaboration ####
## Load the needed libraries ##
library("data.table")
library("igraph")
library(lsa)

args = commandArgs(trailingOnly=TRUE)

# print(args[1])
# print(args[2])

#### get the name of the ICD10 three digits codes ####
if(args[1]=="get_icd10_names"){
  tabla<-fread("Data/icd10cm-order-April-2025.txt",stringsAsFactors = F,sep="\t",header=F)
  codigo<-substring(tabla$V1,7,13)
  codigos<-gsub(" .+","",codigo)
  icd10names<-substring(tabla$V1,17,283)
  icd10names<-gsub("  .+","",icd10names)
  cual<-which(nchar(codigos)==3)
  
  newtab<-cbind(codigos[cual],icd10names[cual])
  newtab[which(newtab[,1]=="U07"),2]<-"COVID-19"
  newtab[,2]<-tolower(newtab[,2])
  colnames(newtab)<-c("Code","Name")
  write.table(newtab,"ICD10_three_digits_names.txt",quote=F,sep="\t",row.names=F)
}

#### create a file associating ICD codes to the omic variables under analysis ####
if(args[1]=="get_icd_variable_associations"){
  #### symptoms - 2025 version ####
  ## @@ @@ @@ @@ @ @ @@ @@ @@ @@ ##
  if(args[2]=="get_icd10_symptom_pmid_tables"){
    ## load the dictionary ##
    dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F)
    setkey(dictionary,"disease")
    ## Select the significant associations between diseases and symptoms ##
    ## Calculate chi-squares to select significant disease-symptom associations ##
    ## Read disease-symptom associations ##
    listofdiseasesymptoms<-fread("Data/Symptoms_Zhou2014_20220712/41467_2014_BFncomms5212_MOESM1045_ESM.txt",stringsAsFactors = F,sep="\t")
    
    ## Diseases ##
    ## @@ @@ @@ ##
    ## 16,749,657 manuscripts
    ## Create a list with the list of papers in which each MeSH term has been indexed ##
    diseasefiles<-list.files("Data/Symptoms_MK_31032025/mesh_pmids_043/")
    disl<-list("covid-19"=fread("Data/Symptoms_MK_31032025/mesh_pmids_covid/covid-19_covid.txt",stringsAsFactors = F,sep="\t")$V1)
    disnumbers<-c("covid-19",length(disl$`covid-19`))
    for(a in 1:length(diseasefiles)){
      codes<-unique(fread(paste("Data/Symptoms_MK_31032025/mesh_pmids_043/",diseasefiles[a],sep=""),stringsAsFactors = F,sep="\t")$V1)
      term<-tolower(gsub("_"," ",gsub("_043.txt","",diseasefiles[a])))
      disl[[term]]<-codes
      disnumbers<-rbind(disnumbers,c(term,length(codes)))
      print(paste(round((a/length(diseasefiles))*100,2),"%",sep=""))
    }
    disnumbers<-as.data.table(disnumbers) ; disnumbers$V2<-as.numeric(disnumbers$V2) ; disnumbers<-disnumbers[order(disnumbers$V2,decreasing = T)]
    colnames(disnumbers)<-c("Disease","PMIDs")
    
    ## Symptoms ##
    ## @@ @@ @@ ##
    ## 2,181,093 manuscripts
    symptomfiles<-list.files("Data/Symptoms_MK_31032025/mesh_pmids_044/")
    sympnumbers<-c()
    for(a in 1:length(symptomfiles)){
      codes<-unique(fread(paste("Data/Symptoms_MK_31032025/mesh_pmids_044/",symptomfiles[a],sep=""),stringsAsFactors = F,sep="\t")$V1)
      term<-tolower(gsub("_"," ",gsub("_044.txt","",symptomfiles[a])))
      if(length(codes)>0){
        tabla<-cbind(term,codes)
        write.table(tabla,"Data/Symptoms_MK_31032025/Symptom_PMIDs.txt",quote=F,sep="\t",row.names = F,col.names = F,append = T)
      }
      sympnumbers<-rbind(sympnumbers,c(term,length(codes)))
      print(paste(round((a/length(symptomfiles))*100,2),"%",sep=""))
    }
    sympnumbers<-as.data.table(sympnumbers) ; sympnumbers$V2<-as.numeric(sympnumbers$V2) ; sympnumbers<-sympnumbers[order(sympnumbers$V2,decreasing = T)]
    colnames(sympnumbers)<-c("Symptom","PMIDs")
    
    # codeswitherror<-rbind(cbind("disorders",gsub("_"," ",disnumbers[grep("err.txt",disnumbers$Disease)]$Disease),
    #                             c("Disorder of Sex Development, 46,XY","Alien Limb Phenomenon","Aortic Dissection","Spider Bites","Arthritis, Juvenile","Uterine Cervical Dysplasia","Neurocognitive Disorders","Severe Dengue","","Drug-Related Side Effects and Adverse Reactions",
    #                               "Chemical and Drug Induced Liver Injury","Chemical and Drug Induced Liver Injury, Chronic","Feeding and Eating Disorders","Stevens-Johnson Syndrome","Felty Syndrome","Fetal Alcohol Spectrum Disorders","Gestational Trophoblastic Disease","Castleman Disease","Hyaline Fibromatosis Syndrome","Idiopathic Hypersomnia",
    #                               "Familial Hypophosphatemic Rickets","Disruptive, Impulse Control, and Conduct Disorders","Brief, Resolved, Unexplained Event","Infections","Intervertebral Disc Degeneration","Intervertebral Disc Displacement","Job Syndrome","Klatskin Tumor","Learning Disabilities","Lecithin Cholesterol Acyltransferase Deficiency",
    #                               "Legg-Calve-Perthes Disease","Megalencephaly","Intellectual Disability","Metabolic Syndrome","Mpox, Monkeypox","Dissociative Identity Disorder","Mycobacterium Infections, Nontuberculous","Malformations of Cortical Development, Group II","Drug Overdose","Paraphilic Disorders",
    #                               "Congenital Hyperinsulinism","Diffuse Intrinsic Pontine Glioma","IgA Vasculitis","Laryngopharyngeal Reflux","Chronic Kidney Disease-Mineral and Bone Disorder","Respiratory Distress Syndrome","Retrognathia","Sarcoma, Ewing","Sleep Wake Disorders","Waardenburg Syndrome",
    #                               "Granulomatosis with Polyangiitis")),
    #                       cbind("symptoms",gsub("_"," ",sympnumbers[grep("err.txt",sympnumbers$Symptom)]$Symptom),
    #                             c("Alien Limb Phenomenon","Hypoxia","Brief, Resolved, Unexplained Event","Learning Disabilities","Intellectual Disability","IgA Vasculitis","Sleep Disorders")))
    # colnames(codeswitherror)<-c("class","old","new")
    # write.table(codeswitherror,"Data/Symptoms_MK_31032025/codes_with_errors.txt",quote=F,sep="\t",row.names = F)
    
    if(length(which(sympnumbers$PMIDs==0))>0){sympnumbers<-sympnumbers[-which(sympnumbers$PMIDs==0)]}
    if(length(which(disnumbers$PMIDs==0))>0){disnumbers<-disnumbers[-which(disnumbers$PMIDs==0)]}
    
    ## transform the disease name to ICD10 codes ##
    common_dictio<-dictionary[intersect(dictionary$disease,disnumbers$Disease)]
    icds<-unique(common_dictio$icd10)
    icdnumbers<-c()
    for(a in 1:length(icds)){
      # a<-1
      cuales<-common_dictio$disease[which(common_dictio$icd10==icds[a])]
      terminos<-c() ; for(b in 1:length(cuales)){terminos<-c(terminos,disl[[cuales[b]]])} ; terminos<-unique(terminos)
      tabla<-cbind(icds[a],terminos)
      write.table(tabla,"Data/Symptoms_MK_31032025/ICD10_PMIDs.txt",quote=F,sep="\t",row.names = F,col.names = F,append = T)
      icdnumbers<-rbind(icdnumbers,c(icds[a],length(terminos)))
      print(paste(round((a/length(icds))*100,2),"%",sep=""))
    }
    colnames(icdnumbers)<-c("icd10","PMIDs") ; icdnumbers<-as.data.table(icdnumbers) ; icdnumbers$PMIDs<-as.numeric(icdnumbers$PMIDs)
    ## Write the numbers ##
    write.table(icdnumbers,"Data/Symptoms_MK_31032025/PubMedIdsByICD10.txt",quote=F,sep="\t",row.names=F)
    write.table(sympnumbers,"Data/Symptoms_MK_31032025/PubMedIdsBySymptom.txt",quote=F,sep="\t",row.names=F)
  }
  ## Define the function for calculating the TF-IDF ##
  tfidf<-function(occurrence,diseaseswithsymptom,numberdiseases){
    rtfidf<-occurrence*log(numberdiseases/diseaseswithsymptom)
    return(rtfidf)
  }
  ## Calculate the number of shared symptoms ##
  ini<-Sys.time()
  icdpmid<-fread("Data/Symptoms_MK_31032025/ICD10_PMIDs.txt",stringsAsFactors = F,sep="\t",header=F) ; colnames(icdpmid)<-c("icd10","pubmedid") ; setkey(icdpmid,"icd10")
  sympmid<-fread("Data/Symptoms_MK_31032025/Symptom_PMIDs.txt",stringsAsFactors = F,sep="\t",header=F) ; colnames(sympmid)<-c("symptom","pubmedid") ; setkey(sympmid,"symptom")
  icds<-unique(icdpmid$icd10) ; symptoms<-unique(sympmid$symptom)
  icduniverse<-unique(icdpmid$pubmedid) ; sympuniverse<-unique(sympmid$pubmedid)
  cooccurrences<-matrix(nrow=length(icds)*length(symptoms),ncol=5,NA)
  counter<-0
  for(a in 1:length(icds)){
    # a<-1
    for(b in 1:length(symptoms)){
      counter<-counter+1
      # b<-1
      idicds<-icdpmid[icds[a]]$pubmedid
      idsymp<-sympmid[symptoms[b]]$pubmedid
      cooccurrences[counter,]<-c(icds[a],symptoms[b],length(intersect(idicds,idsymp)),length(idicds),length(idsymp))
    }
    print(paste(round((a/length(icds))*100,2),"%",sep=""))
  }
  cooccurrences<-as.data.table(cooccurrences)
  colnames(cooccurrences)<-c("icd10","symptom","occurrence","n_icd10","n_symptom")
  cooccurrences$occurrence<-as.numeric(cooccurrences$occurrence)
  ## First get the number of ICD10 associated to each symptom ##
  Ndiseasesbysymptom<-table(cooccurrences[which(cooccurrences$occurrence>0)]$symptom)
  ## Now go disease by disease ##
  theicds<-unique(cooccurrences$icd10)
  icdswithtfidf<-c()
  for(a in 1:length(theicds)){
    # a<-1
    thedisease<-cooccurrences[which(cooccurrences$icd10==theicds[a])]
    thetfidfs<-c()
    for(b in 1:length(thedisease$icd10)){
      thevalue<-tfidf(thedisease$occurrence[b],Ndiseasesbysymptom[thedisease$symptom[b]],length(theicds))
      thetfidfs<-c(thetfidfs,as.numeric(thevalue))
    }
    thedisease<-cbind(thedisease,thetfidfs)
    icdswithtfidf<-rbind(icdswithtfidf,thedisease)
    print(paste(round((a/length(theicds))*100,2),"%",sep=""))
  }
  icdswithtfidf<-icdswithtfidf[,c(1,2,6)]
  ## Consider only those real icd-disease associations for the fs and ifs calculations ##
  numsympassocdis<-table(cooccurrences$symptom[cooccurrences$occurrence!=0])
  pesos<-numsympassocdis[icdswithtfidf$symptom]
  ifs_f<-1/pesos ; ifs_f[which(cooccurrences$occurrence==0)]<-0
  # ifs_f<-1/table(icdswithtfidf$symptom)
  fs_f<-pesos ; fs_f[which(cooccurrences$occurrence==0)]<-0
  # fs_f<-table(icdswithtfidf$symptom)/max(table(icdswithtfidf$symptom))
  icdswithtfidf<-cbind(icdswithtfidf,as.numeric(ifs_f),as.numeric(fs_f))
  # icdswithtfidf<-cbind(icdswithtfidf,as.numeric(ifs_f[icdswithtfidf$symptom]),as.numeric(fs_f[icdswithtfidf$symptom]))
  colnames(icdswithtfidf)[3:5]<-c("tfidf","ifs","fs")
  ## Write the table containing the number of ICD10 - symptom occurrences with the associated TF-IDF ##
  icdswithtfidf$symptom<-tolower(icdswithtfidf$symptom)
  write.table(icdswithtfidf,"WorkFiles/ICD10_Omics/ICD10_mksymptom.txt",quote=F,sep="\t",row.names=F)
  fin<-Sys.time()
  
  # #### symptoms ####
  # ##   @@ @@ @@   ##
  # ## load the dictionary ##
  # dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F)
  # setkey(dictionary,"disease")
  # ## Select the significant associations between diseases and symptoms ##
  # ## Calculate chi-squares to select significant disease-symptom associations ##
  # ## Read disease-symptom associations ##
  # listofdiseasesymptoms<-fread("Data/Symptoms_Zhou2014_20220712/41467_2014_BFncomms5212_MOESM1045_ESM.txt",stringsAsFactors = F,sep="\t")
  # ## Get the number of papers associated to a symptom ##
  # symptomfrequency<-fread("Data/Symptoms_Zhou2014_20220712/41467_2014_BFncomms5212_MOESM1044_ESM.txt")
  # setkey(symptomfrequency,"MeSH Symptom Term")
  # ## Get the number of papers associated to a disease ##
  # diseasefrequency<-fread("Data/Symptoms_Zhou2014_20220712/41467_2014_BFncomms5212_MOESM1043_ESM.txt")
  # setkey(diseasefrequency,"MeSH Disease Term")
  # ## Get the p-values and the chi-square estimates ##
  # pvaltabs<-cbind(listofdiseasesymptoms[,1:2],NA,NA)
  # for(a in 1:length(listofdiseasesymptoms$`MeSH Symptom Term`)){
  #   # a<-1
  #   subtab<-rbind(c(listofdiseasesymptoms$`PubMed occurrence`[a],
  #                   symptomfrequency[listofdiseasesymptoms$`MeSH Symptom Term`[a]]$`PubMed occurrence`-listofdiseasesymptoms$`PubMed occurrence`[a]),
  #                 c(diseasefrequency[listofdiseasesymptoms$`MeSH Disease Term`[a]]$`PubMed occurrence`-listofdiseasesymptoms$`PubMed occurrence`[a],
  #                   849103-symptomfrequency[listofdiseasesymptoms$`MeSH Symptom Term`[a]]$`PubMed occurrence`-diseasefrequency[listofdiseasesymptoms$`MeSH Disease Term`[a]]$`PubMed occurrence`+listofdiseasesymptoms$`PubMed occurrence`[a]))
  #   chisq <- chisq.test(subtab)
  #   pvaltabs$V2[a]<-chisq$p.value
  #   pvaltabs$V3[a]<-chisq$statistic
  #   print(paste(round((a/length(listofdiseasesymptoms$`MeSH Symptom Term`))*100,2),"%",sep=""))
  # }
  # pvaltabs<-cbind(pvaltabs,NA)
  # colnames(pvaltabs)[3:5]<-c("pval","chisquare","fdr")
  # pvaltabs$fdr<-pvaltabs$pval*length(pvaltabs$pval)
  # if(length(which(pvaltabs$fdr>1))>0){pvaltabs$fdr[which(pvaltabs$fdr>1)]<-1}
  # listofdiseasesymptoms<-cbind(listofdiseasesymptoms,pvaltabs$pval,pvaltabs$fdr)
  # colnames(listofdiseasesymptoms)[5:6]<-c("pval","fdr")
  # write.table(listofdiseasesymptoms,"Data/Symptoms_Zhou2014_20220712/DiseaseSymptom_fdr.txt",quote=F,sep="\t",row.names=F)
  # ## Load the disease-symptom association file
  # diseasesymptom<-fread("Data/Symptoms_Zhou2014_20220712/DiseaseSymptom_fdr.txt")
  # diseasesymptom2<-diseasesymptom
  # diseasesymptom2$`MeSH Symptom Term`<-tolower(diseasesymptom2$`MeSH Symptom Term`) ; diseasesymptom2$`MeSH Disease Term`<-tolower(diseasesymptom2$`MeSH Disease Term`)
  # colnames(diseasesymptom2)[1:2]<-c("symptom","disease")
  # diseasesymptom2<-diseasesymptom2[,c(2,1,3:6)]
  # write.table(diseasesymptom2,"WorkFiles/DiseaseName_Omics/DiseaseName_symptom.txt",quote=F,sep="\t",row.names=F)
  # ## Select only the significant disease - symptom associations ##
  # diseasesymptom<-diseasesymptom[which(diseasesymptom$fdr<=0.05)]
  # ## Select the diseases for which we have and ICD correspondence and sum their mentions ##
  # common_dictio<-dictionary[intersect(dictionary$disease,tolower(diseasesymptom$`MeSH Disease Term`))]
  # icds<-unique(common_dictio$icd10)
  # icdsymptom<-c()
  # for(a in 1:length(icds)){
  #   # a<-1
  #   diseases<-common_dictio$disease[which(common_dictio$icd10==icds[a])]
  #   icdtab<-c()
  #   for(b in 1:length(diseases)){
  #     # b<-1
  #     icdtab<-rbind(icdtab,diseasesymptom[which(tolower(diseasesymptom$`MeSH Disease Term`)==diseases[b])])
  #   }
  #   icdtab<-icdtab[,lapply(.SD,sum),by=.(`MeSH Symptom Term`), .SDcols=c("PubMed occurrence")]
  #   icdtab<-cbind(icds[a],icdtab)
  #   colnames(icdtab)<-c("icd10","symptom","occurrence")
  #   icdsymptom<-rbind(icdsymptom,icdtab)
  #   print(paste(round((a/length(icds))*100,2),"%",sep=""))
  # }
  # ## Write the table containing the number of ICD10 - symptom occurrences ##
  # write.table(icdsymptom,"Omics/ICD10_symptom_associations.txt",quote=F,sep="\t",row.names=F)
  # ## Calculate the TF-IDF ##
  # if(dim(icdsymptom)[2]==3){
  #   ## First get the number of ICD10 associated to each symptom ##
  #   Ndiseasesbysymptom<-table(icdsymptom$symptom)
  #   ## Now go disease by disease ##
  #   theicds<-unique(icdsymptom$icd10)
  #   icdswithtfidf<-c()
  #   for(a in 1:length(theicds)){
  #     # a<-1
  #     thedisease<-icdsymptom[which(icdsymptom$icd10==theicds[a])]
  #     thetfidfs<-c()
  #     for(b in 1:length(thedisease$icd10)){
  #       thevalue<-tfidf(thedisease$occurrence[b],Ndiseasesbysymptom[thedisease$symptom[b]],length(theicds))
  #       thetfidfs<-c(thetfidfs,as.numeric(thevalue))
  #     }
  #     thedisease<-cbind(thedisease,thetfidfs)
  #     icdswithtfidf<-rbind(icdswithtfidf,thedisease)
  #     print(paste(round((a/length(theicds))*100,2),"%",sep=""))
  #   }
  #   icdswithtfidf<-icdswithtfidf[,c(1,2,4)]
  #   ifs_f<-1/table(icdswithtfidf$symptom)
  #   fs_f<-table(icdswithtfidf$symptom)/max(table(icdswithtfidf$symptom))
  #   icdswithtfidf<-cbind(icdswithtfidf,as.numeric(ifs_f[icdswithtfidf$symptom]),as.numeric(fs_f[icdswithtfidf$symptom]))
  #   colnames(icdswithtfidf)[3:5]<-c("tfidf","ifs","fs")
  #   ## Write the table containing the number of ICD10 - symptom occurrences with the associated TF-IDF ##
  #   write.table(icdswithtfidf,"Omics/ICD10_symptom_associations.txt",quote=F,sep="\t",row.names=F)
  #   icdswithtfidf$symptom<-tolower(icdswithtfidf$symptom)
  #   write.table(icdswithtfidf,"WorkFiles/ICD10_Omics/ICD10_symptom.txt",quote=F,sep="\t",row.names=F)
  # }
  
  #### miRNAs ####
  ##   @@  @@   ##
  ## load the dictionary ##
  dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F)
  setkey(dictionary,"disease")
  ## Load the disease-mirna association file
  diseasemirna<-fread("Data/miRNA_HMDD_20240523/alldata_v4.txt",stringsAsFactors = F,sep="\t")
  diseasemirna<-diseasemirna[,c(4,3,1,2)]
  diseasemirna$disease<-tolower(diseasemirna$disease)
  ## For the future calculation of the TF-IDF remove if a paper associates a disease to a miRNA twice ##
  if(length(which(duplicated(diseasemirna[,c(1,2,4)])))>0){diseasemirna<-diseasemirna[-which(duplicated(diseasemirna[,c(1,2,4)]))]}
  write.table(diseasemirna,"WorkFiles/DiseaseName_Omics/DiseaseName_miRNA.txt",quote=F,sep="\t",row.names=F)
  ## Each paper referencing the association between a disease and a mirna, in this case focused on Alzheimer's disease ##
  # diseasemirna[intersect(which(diseasemirna$disease=="Alzheimer Disease"),which(diseasemirna$miRNA=="hsa-mir-146a")),]
  ## Select common diseases ##
  common_dictio<-dictionary[intersect(dictionary$disease,diseasemirna$disease)]
  icds<-unique(common_dictio$icd10)
  ## Build an ICD10 - mirna association table ##
  icdmirna<-c()
  for(a in 1:length(icds)){
    # a<-1
    diseases<-common_dictio$disease[which(common_dictio$icd10==icds[a])]
    icdtab<-c()
    for(b in 1:length(diseases)){
      # b<-1
      icdtab<-rbind(icdtab,diseasemirna[which(diseasemirna$disease==diseases[b])])
    }
    icdtab<-cbind(icds[a],icdtab)
    colnames(icdtab)<-c("icd10","disease","mirna","code","pmid")
    icdmirna<-rbind(icdmirna,icdtab)
    print(paste(round((a/length(icds))*100,2),"%",sep=""))
  }
  ## Calculate the number of papers associating each icd10 to a miRNA ##
  icdtab<-icdmirna[, .(.N), by = .(icd10,mirna)]
  colnames(icdtab)<-c("icd10","mirna","occurrence")
  ## Write the table containing the number of ICD10 - mirna occurrences ##
  write.table(icdtab,"Omics/ICD10_mirna_associations.txt",quote=F,sep="\t",row.names=F)
  ## Calculate the TF-IDF ##
  ## First get the number of ICD10 associated to each symptom ##
  Ndiseasesbymirna<-table(icdtab$mirna)
  ## Now go disease by disease ##
  theicds<-unique(icdtab$icd10)
  icdswithtfidf<-c()
  for(a in 1:length(theicds)){
    # a<-1
    thedisease<-icdtab[which(icdtab$icd10==theicds[a])]
    thetfidfs<-c()
    for(b in 1:length(thedisease$icd10)){
      thevalue<-tfidf(thedisease$occurrence[b],Ndiseasesbymirna[thedisease$mirna[b]],length(theicds))
      thetfidfs<-c(thetfidfs,as.numeric(thevalue))
    }
    thedisease<-cbind(thedisease,thetfidfs)
    icdswithtfidf<-rbind(icdswithtfidf,thedisease)
    print(paste(round((a/length(theicds))*100,2),"%",sep=""))
  }
  icdswithtfidf<-icdswithtfidf[,c(1,2,4)]
  ifs_f<-1/table(icdswithtfidf$mirna)
  fs_f<-table(icdswithtfidf$mirna)/max(table(icdswithtfidf$mirna))
  icdswithtfidf<-cbind(icdswithtfidf,as.numeric(ifs_f[icdswithtfidf$mirna]),as.numeric(fs_f[icdswithtfidf$mirna]))
  colnames(icdswithtfidf)[3:5]<-c("tfidf","ifs","fs")
  ## Write the table containing the number of ICD10 - mirna occurrences with the associated TF-IDF ##
  write.table(icdswithtfidf,"Omics/ICD10_mirna_associations.txt",quote=F,sep="\t",row.names=F)
  write.table(icdswithtfidf,"WorkFiles/ICD10_Omics/ICD10_miRNA.txt",quote=F,sep="\t",row.names=F)
  
  #### genes ####
  ##   @@ @@   ##
  ## load the dictionary ##
  dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F)
  setkey(dictionary,"disease")
  ## Load the disease-symptom association file
  diseasegene<-fread("Data/Genes_DisGeNET_20230113/curated_gene_disease_associations.tsv",stringsAsFactors = F,sep="\t")
  diseasegene<-diseasegene[,c(2,6)]
  diseasegene$diseaseName<-tolower(diseasegene$diseaseName)
  diseasegene2<-diseasegene[,2:1] ; colnames(diseasegene2)<-c("disease","gene")
  ## We are going to add COVID-19 genes from CTD on 2025/01/30
  covid<-fread("Data/Genes_DisGeNET_20230113/CTD_D000086382_genes_20250130134713.csv",stringsAsFactors = F)
  covid<-covid[grep("marker",covid$`Direct Evidence`),c(3,1)]
  covid$`Disease Name`<-tolower(covid$`Disease Name`) ; colnames(covid)<-c("disease","gene")
  diseasegene2<-rbind(diseasegene2,covid)
  write.table(diseasegene2,"WorkFiles/DiseaseName_Omics/DiseaseName_gene.txt",quote=F,sep="\t",row.names=F)
  ## Select common diseases ##
  covid<-covid[,2:1] ; colnames(covid)<-colnames(diseasegene)
  diseasegene<-rbind(diseasegene,covid)
  common_dictio<-dictionary[intersect(dictionary$disease,diseasegene$diseaseName)]
  icds<-unique(common_dictio$icd10)
  ## Build an ICD10 - gene association table ##
  icdgene<-c()
  for(a in 1:length(icds)){
    # a<-1
    diseases<-dictionary$disease[which(dictionary$icd10==icds[a])]
    icdtab<-c()
    for(b in 1:length(diseases)){
      # b<-1
      icdtab<-rbind(icdtab,diseasegene[which(diseasegene$diseaseName==diseases[b])])
    }
    if(length(which(duplicated(icdtab$geneSymbol)))>0){icdtab<-icdtab[-which(duplicated(icdtab$geneSymbol))]}
    icdtab<-cbind(icds[a],icdtab$geneSymbol)
    colnames(icdtab)<-c("icd10","gene")
    icdgene<-rbind(icdgene,icdtab)
    print(paste(round((a/length(icds))*100,2),"%",sep=""))
  }
  icdgene<-as.data.table(icdgene)
  ## remove miRNAs ##
  icdgene<-icdgene[-grep("mir",icdgene$gene,ignore.case = T)]
  ## calculate ifs and fs ##
  ifs_f<-1/table(icdgene$gene)
  fs_f<-table(icdgene$gene)/max(table(icdgene$gene))
  icdgene<-cbind(icdgene,as.numeric(ifs_f[icdgene$gene]),as.numeric(fs_f[icdgene$gene]))
  colnames(icdgene)[3:4]<-c("ifs","fs")
  ## Write the table containing the number of ICD10 - gene occurrences ##
  write.table(icdgene,"Omics/ICD10_gene_associations.txt",quote=F,sep="\t",row.names=F)
  write.table(icdgene,"WorkFiles/ICD10_Omics/ICD10_gene.txt",quote=F,sep="\t",row.names=F)
  
  #### drugs ####
  ##   @@ @@   ##
  ## load the dictionary ##
  dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F)
  setkey(dictionary,"disease")
  ## Load the disease-drug association file
  ## downloaded on 20250130 ##
  diseasedrug<-fread("Data/Drugs_CTD_20220712/CTD_chemicals_diseases.tsv",stringsAsFactors = F,sep="\t")
  diseasedrug<-diseasedrug[unique(c(which(diseasedrug$DirectEvidence=="marker/mechanism"),which(diseasedrug$DirectEvidence=="therapeutic"))),c(1,2,4:6,10)]
  write.table(diseasedrug,"Data/Drugs_CTD_20220712/Disease_drug_direct_evidence.txt",quote=F,sep="\t",row.names = F)
  ## Select only the drugs with a therapeutic association with the disease ##
  diseasedrug<-diseasedrug[which(diseasedrug$DirectEvidence=="therapeutic")]
  diseasedrug$DiseaseName<-tolower(diseasedrug$DiseaseName)
  diseasedrug2<-diseasedrug[,c(3,1,6,4,2)] ; colnames(diseasedrug2)<-c("disease","drug","PMID","mesh_disease","mesh_drug")
  write.table(diseasedrug2,"WorkFiles/DiseaseName_Omics/DiseaseName_drug.txt",quote=F,sep="\t",row.names=F)
  ## Select common diseases ##
  common_dictio<-dictionary[intersect(dictionary$disease,diseasedrug$DiseaseName)]
  icds<-unique(common_dictio$icd10)
  ## Build an ICD10 - drug association table ##
  icddrug<-c()
  for(a in 1:length(icds)){
    # a<-1
    diseases<-dictionary$disease[which(dictionary$icd10==icds[a])]
    icdtab<-c()
    for(b in 1:length(diseases)){
      # b<-1
      icdtab<-rbind(icdtab,diseasedrug[which(diseasedrug$DiseaseName==diseases[b])][,c(1,3)])
    }
    if(length(which(duplicated(icdtab$drugSymbol)))>0){icdtab<-icdtab[-which(duplicated(icdtab$drugSymbol))]}
    icdtab<-cbind(icds[a],icdtab)
    icddrug<-rbind(icddrug,icdtab)
    print(paste(round((a/length(icds))*100,2),"%",sep=""))
  }
  icddrug<-icddrug[,1:2]
  icddrug<-as.data.table(icddrug)
  colnames(icddrug)<-c("icd10","drug")
  ## calculate ifs and fs ##
  ifs_f<-1/table(icddrug$drug)
  fs_f<-table(icddrug$drug)/max(table(icddrug$drug))
  icddrug<-cbind(icddrug,as.numeric(ifs_f[icddrug$drug]),as.numeric(fs_f[icddrug$drug]))
  colnames(icddrug)[3:4]<-c("ifs","fs")
  ## Write the table containing the number of ICD10 - drug occurrences ##
  write.table(icddrug,"Omics/ICD10_drug_associations.txt",quote=F,sep="\t",row.names=F)
  write.table(icddrug,"WorkFiles/ICD10_Omics/ICD10_drug.txt",quote=F,sep="\t",row.names=F)
  
  #### drugs Chembl ####
  ##   @@ @@   ##
  ## load the dictionary ##
  dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F)
  setkey(dictionary,"disease")
  ## Load the disease-drug association file
  ## downloaded on 20250130 ##
  diseasedrug<-fread("Data/Drugs_ChemBL_20240326/approved_indication.csv",stringsAsFactors = F,sep=";")
  diseasedrug<-diseasedrug[,c(2,7,6,5)]
  colnames(diseasedrug)<-c("drug","disease","mesh_disease","year")
  ## Remove the drugs that have been withdrawn ##
  withdrawn<-fread("Data/Drugs_ChemBL_20240326/ChEMBL_withdrawn drugs.csv",stringsAsFactors = F,sep=";")
  removethem<-unique(withdrawn$"Parent Molecule Name")
  removeind<-c() ; for(a in 1:length(removethem)){removeind<-c(removeind,which(diseasedrug$drug==removethem[a]))}
  diseasedrug<-diseasedrug[-removeind,]
  ## put disease names into lower case
  diseasedrug$disease<-tolower(diseasedrug$disease) ; diseasedrug$drug<-tolower(diseasedrug$drug)
  diseasedrug<-diseasedrug[,c(2,1,3,4)]
  write.table(diseasedrug,"WorkFiles/DiseaseName_Omics/DiseaseName_drugchembl.txt",quote=F,sep="\t",row.names=F)
  ## Select common diseases ##
  common_dictio<-dictionary[intersect(dictionary$disease,diseasedrug$disease)]
  icds<-unique(common_dictio$icd10)
  ## Build an ICD10 - drug association table ##
  icddrug<-c()
  for(a in 1:length(icds)){
    # a<-1
    diseases<-dictionary$disease[which(dictionary$icd10==icds[a])]
    icdtab<-c()
    for(b in 1:length(diseases)){
      # b<-1
      icdtab<-rbind(icdtab,diseasedrug[which(diseasedrug$disease==diseases[b])][,c(1,2)])
    }
    if(length(which(duplicated(icdtab$drug)))>0){icdtab<-icdtab[-which(duplicated(icdtab$drug))]}
    icdtab<-cbind(icds[a],icdtab)
    icddrug<-rbind(icddrug,icdtab)
    print(paste(round((a/length(icds))*100,2),"%",sep=""))
  }
  icddrug<-icddrug[,c(1,3)]
  icddrug<-as.data.table(icddrug)
  colnames(icddrug)<-c("icd10","drug")
  ## calculate ifs and fs ##
  ifs_f<-1/table(icddrug$drug)
  fs_f<-table(icddrug$drug)/max(table(icddrug$drug))
  icddrug<-cbind(icddrug,as.numeric(ifs_f[icddrug$drug]),as.numeric(fs_f[icddrug$drug]))
  colnames(icddrug)[3:4]<-c("ifs","fs")
  ## Write the table containing the number of ICD10 - drug occurrences ##
  write.table(icddrug,"WorkFiles/ICD10_Omics/ICD10_drugchembl.txt",quote=F,sep="\t",row.names=F)
  
  #### drugs Chembl 2 ####
  ##   @@ @@   ##
  ## load the dictionary ##
  dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F)
  setkey(dictionary,"disease")
  ## Load the disease-drug association file
  ## downloaded on 20250130 ##
  diseasedrug<-fread("Data/Drugs_ChemBL_20240326/all_drug_indications_ChEMBL2905.csv",stringsAsFactors = F,sep=";")
  diseasedrug<-diseasedrug[,c(2:7)]
  colnames(diseasedrug)<-c("drug","type","phase","year","mesh_disease","disease")
  ## put disease names into lower case
  diseasedrug$disease<-tolower(diseasedrug$disease) ; diseasedrug$drug<-tolower(diseasedrug$drug)
  diseasedrug<-diseasedrug[,c(6,1,2:5)]
  write.table(diseasedrug,"WorkFiles/DiseaseName_Omics/DiseaseName_drugchembl2.txt",quote=F,sep="\t",row.names=F)
  ## Select common diseases ##
  common_dictio<-dictionary[intersect(dictionary$disease,diseasedrug$disease)]
  icds<-unique(common_dictio$icd10)
  ## Build an ICD10 - drug association table ##
  icddrug<-c()
  for(a in 1:length(icds)){
    # a<-1
    diseases<-dictionary$disease[which(dictionary$icd10==icds[a])]
    icdtab<-c()
    for(b in 1:length(diseases)){
      # b<-1
      icdtab<-rbind(icdtab,diseasedrug[which(diseasedrug$disease==diseases[b])][,c(1,2)])
    }
    if(length(which(duplicated(icdtab$drug)))>0){icdtab<-icdtab[-which(duplicated(icdtab$drug))]}
    icdtab<-cbind(icds[a],icdtab)
    icddrug<-rbind(icddrug,icdtab)
    print(paste(round((a/length(icds))*100,2),"%",sep=""))
  }
  icddrug<-icddrug[,c(1,3)]
  icddrug<-as.data.table(icddrug)
  colnames(icddrug)<-c("icd10","drug")
  ## calculate ifs and fs ##
  ifs_f<-1/table(icddrug$drug)
  fs_f<-table(icddrug$drug)/max(table(icddrug$drug))
  icddrug<-cbind(icddrug,as.numeric(ifs_f[icddrug$drug]),as.numeric(fs_f[icddrug$drug]))
  colnames(icddrug)[3:4]<-c("ifs","fs")
  ## Write the table containing the number of ICD10 - drug occurrences ##
  write.table(icddrug,"WorkFiles/ICD10_Omics/ICD10_drugchembl2.txt",quote=F,sep="\t",row.names=F)
  
  #### metabolome ####
  ##   @ @ @@ @ @   ##
  ## load the dictionary ##
  dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F)
  setkey(dictionary,"disease")
  ## Load the disease-symptom association file
  diseasemetabolite<-fread("Data/Metabolome_HMDB_20220712/disease_metabolite_table.txt",stringsAsFactors = F,sep="\t")
  diseasemetabolite$Disease<-tolower(diseasemetabolite$Disease)
  ## For simplicity we focus only on blood levels ##
  diseasemetabolite<-diseasemetabolite[which(diseasemetabolite$Where=="Blood"),]
  diseasemetabolite$Metabolite<-tolower(diseasemetabolite$Metabolite)
  ## load COVID-19 metabolites ##
  covid<-fread("Data/Metabolome_HMDB_20220712/COVID_datasets/MetaboliteCOVID_Serum Metabolites.tsv",stringsAsFactors = F)
  covid<-covid[which(covid$Indication=="Disease biomarker")] ; colnames(covid)[2]<-"metabolite"
  covid$metabolite<-tolower(covid$metabolite)
  setkey(covid,"metabolite")
  covid<-covid[intersect(unique(diseasemetabolite$Metabolite),unique(covid$metabolite))]
  diseasemetabolite<-diseasemetabolite[,1:2]
  covid<-cbind("covid-19",covid$metabolite) ; colnames(covid)<-c("Disease","Metabolite")
  diseasemetabolite<-rbind(diseasemetabolite,covid)
  colnames(diseasemetabolite)<-c("disease","metabolite")
  write.table(diseasemetabolite,"WorkFiles/DiseaseName_Omics/DiseaseName_metabolome.txt",quote=F,sep="\t",row.names=F)
  ## Select common diseases ##
  common_dictio<-dictionary[intersect(dictionary$disease,diseasemetabolite$disease)]
  icds<-unique(common_dictio$icd10)
  ## Build an ICD10 - metabolite association table ##
  icdmetabolite<-c()
  for(a in 1:length(icds)){
    # a<-1
    diseases<-common_dictio$disease[which(common_dictio$icd10==icds[a])]
    diseases<-intersect(diseasemetabolite$disease,diseases)
    if(length(diseases)>0){
      icdtab<-c()
      for(b in 1:length(diseases)){
        # b<-1
        icdtab<-rbind(icdtab,diseasemetabolite[which(diseasemetabolite$disease==diseases[b])])
      }
      if(length(which(duplicated(icdtab$HMDB_ID)))>0){icdtab<-icdtab[-which(duplicated(icdtab$HMDB_ID))]}
      icdtab<-cbind(icds[a],icdtab)
      icdmetabolite<-rbind(icdmetabolite,icdtab)
      print(paste(round((a/length(icds))*100,2),"%",sep=""))
    }
  }
  icdmetabolite<-icdmetabolite[,-2]
  icdmetabolite<-as.data.table(icdmetabolite)
  colnames(icdmetabolite)<-c("icd10","metabolome")
  ## calculate ifs and fs ##
  ifs_f<-1/table(icdmetabolite$metabolome)
  fs_f<-table(icdmetabolite$metabolome)/max(table(icdmetabolite$metabolome))
  icdmetabolite<-cbind(icdmetabolite,as.numeric(ifs_f[icdmetabolite$metabolome]),as.numeric(fs_f[icdmetabolite$metabolome]))
  colnames(icdmetabolite)[3:4]<-c("ifs","fs")
  ## Write the table containing the number of ICD10 - metabolome occurrences ##
  write.table(icdmetabolite,"Omics/ICD10_metabolome_associations.txt",quote=F,sep="\t",row.names=F)
  write.table(icdmetabolite,"WorkFiles/ICD10_Omics/ICD10_metabolome.txt",quote=F,sep="\t",row.names=F)
  
  #### microbiome ####
  ##   @ @ @@ @ @   ##
  ## load the dictionary ##
  dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F)
  setkey(dictionary,"disease")
  ## Load the disease-symptom association file
  diseasemicrobiome<-fread("Data/microbiome_disbiome_20220714/DiseaseMicrobeTable.txt",stringsAsFactors = F,sep="\t")
  ## For simplicity we focus only on gut microbiome ##
  diseasemicrobiome<-diseasemicrobiome[which(diseasemicrobiome$sample_name=="Faeces")]
  diseasemicrobiome$disease_name<-tolower(diseasemicrobiome$disease_name)
  diseasemicrobiome<-diseasemicrobiome[,c(8,10,3)]
  write.table(diseasemicrobiome,"WorkFiles/DiseaseName_Omics/DiseaseName_microbiome.txt",quote=F,sep="\t",row.names=F)
  ## Select common diseases ##
  common_dictio<-dictionary[intersect(dictionary$disease,diseasemicrobiome$disease_name)]
  icds<-unique(common_dictio$icd10)
  ## Build an ICD10 - microbiome association table ##
  icdmicrobiome<-c()
  for(a in 1:length(icds)){
    # a<-1
    diseases<-common_dictio$disease[which(common_dictio$icd10==icds[a])]
    if(length(diseases)>0){
      icdtab<-c()
      for(b in 1:length(diseases)){
        # b<-1
        icdtab<-rbind(icdtab,diseasemicrobiome[which(diseasemicrobiome$disease_name==diseases[b]),])
      }
      # if(length(which(duplicated(icdtab[,2:3])))>0){icdtab<-icdtab[-which(duplicated(icdtab[,2:3]))]}
      icdtab<-cbind(icds[a],icdtab)
      icdmicrobiome<-rbind(icdmicrobiome,icdtab)
      print(paste(round((a/length(icds))*100,2),"%",sep=""))
    }
  }
  icdmicrobiome<-icdmicrobiome[,-2]
  icdmicrobiome<-as.data.table(icdmicrobiome)
  colnames(icdmicrobiome)<-c("icd10","microbiome","direction")
  ## If it has been described both an increase and a decrease in the disease do the following:
  ## if there are more associations in one direction than the other select it
  ## if there are the same number of elevated and reduced associations remove them ##
  if(length(which(duplicated(icdmicrobiome[,1:2])))>0){
    duplis<-icdmicrobiome[which(duplicated(icdmicrobiome[,1:2])),1:2]
    ## we just need the duplicated pair once, so if they appear more times we remove them ##
    if(length(which(duplicated(duplis)))>0){duplis<-duplis[-which(duplicated(duplis))]}
    removethem<-c()
    for(a in 1:length(duplis$icd10)){
      # a<-1
      theindex<-intersect(which(icdmicrobiome$icd10==duplis$icd10[a]),which(icdmicrobiome$microbiome==duplis$microbiome[a]))
      subtab<-table(icdmicrobiome[theindex,3])
      if(length(subtab)==2){
        if(as.numeric(subtab["Elevated"])==as.numeric(subtab["Reduced"])){
          removethem<-c(removethem,theindex)
        }
        if(as.numeric(subtab["Elevated"])>as.numeric(subtab["Reduced"])){
          # print(a)
          remove_r<-theindex[which(icdmicrobiome$direction[theindex]=="Reduced")]
          remove_e<-setdiff(theindex,remove_r)[-1]
          removethem<-c(removethem,remove_r,remove_e)
        }
        if(as.numeric(subtab["Elevated"])<as.numeric(subtab["Reduced"])){
          # print(a)
          remove_e<-theindex[which(icdmicrobiome$direction[theindex]=="Elevated")]
          remove_r<-setdiff(theindex,remove_e)[-1]
          removethem<-c(removethem,remove_r,remove_e)
        }
      }
      if(length(subtab)==1){
        removethem<-c(removethem,theindex[-1])
      }
    }
    icdmicrobiome<-icdmicrobiome[-removethem]
  }
  ## calculate ifs and fs ##
  ifs_f<-1/table(icdmicrobiome$microbiome)
  fs_f<-table(icdmicrobiome$microbiome)/max(table(icdmicrobiome$microbiome))
  icdmicrobiome<-cbind(icdmicrobiome,as.numeric(ifs_f[icdmicrobiome$microbiome]),as.numeric(fs_f[icdmicrobiome$microbiome]))
  colnames(icdmicrobiome)[4:5]<-c("ifs","fs")
  icdmicrobiome<-icdmicrobiome[,c(1,2,4,5,3)]
  icdmicrobiome<-cbind(icdmicrobiome,icdmicrobiome$ifs,icdmicrobiome$fs)
  ## give a positive or negative value depending on if its elevated or reduced ##
  icdmicrobiome$V2[which(icdmicrobiome$direction=="Reduced")]<-icdmicrobiome$V2[which(icdmicrobiome$direction=="Reduced")]*(-1)
  icdmicrobiome$V3[which(icdmicrobiome$direction=="Reduced")]<-icdmicrobiome$V3[which(icdmicrobiome$direction=="Reduced")]*(-1)
  icdmicrobiome<-icdmicrobiome[,c(1:4,6,7)]
  colnames(icdmicrobiome)[5:6]<-c("ifsd","fsd")
  ## Write the table containing the number of ICD10 - microorganism associations ##
  write.table(icdmicrobiome,"Omics/ICD10_microbiome_associations.txt",quote=F,sep="\t",row.names=F)
  write.table(icdmicrobiome,"WorkFiles/ICD10_Omics/ICD10_microbiome.txt",quote=F,sep="\t",row.names=F)
  
  #### methylome ####
  ##   @@ @ @ @@   ##
  distab<-fread("Data/Methylation_EWASdb_20220701/icd10_phenotype.txt",stringsAsFactors = F,sep="\t")
  fdistab<-distab$class ; names(fdistab)<-distab$ewas_id
  methylation<-fread("Data/Methylation_EWASdb_20220701/ICD10_methylome_associations_duplicated_associations.txt",stringsAsFactors = F)
  ewassize<-fread("Data/Methylation_EWASdb_20220701/EWAS_platform_associations.txt",stringsAsFactors = F)
  ## To avoid problems we remove the EWAS analyzing the big platform ##
  remove_ewas<-ewassize$ewas[which(ewassize$size==850000)]
  for(a in 1:length(remove_ewas)){
    # a<-1
    if(length(which(methylation$EWAS_ID==remove_ewas[a]))>0){
      print(length(which(methylation$EWAS_ID==remove_ewas[a])))
      methylation<-methylation[-which(methylation$EWAS_ID==remove_ewas[a]),]
    }
  }
  if(args[2]=="plot"){
    thresholds<-c(1e-3,1e-4,1e-5,1e-6,1e-7,1e-8)
    thethresholdscgs<-c() ; thethresholdsgenes<-c() ; thethresholdscggenes<-c()
    for(a in 1:length(thresholds)){
      # a<-5
      methyl<-methylation[which(methylation$T_Pvalue<=thresholds[a])]
      ## check the number of cg significantly associated to the diseases with varying thresholds ##
      ## working at the level of cg ##
      methylcgs<-methyl[,lapply(.SD,mean),by=.(icd10,ID_REF), .SDcols=c("T_statistics")]
      thethresholdscgs<-rbind(thethresholdscgs,cbind(thresholds[a],as.data.frame(table(methylcgs$icd10))))
      ## working at the gene level - focused on cg ##
      methylcggenes<-methyl[which(methyl$Gene_symbol!="-")]
      methylcggenes<-methylcggenes[,lapply(.SD,mean),by=.(icd10,ID_REF), .SDcols=c("T_statistics")]
      thethresholdscggenes<-rbind(thethresholdscggenes,cbind(thresholds[a],as.data.frame(table(methylcggenes$icd10))))
      ## working at the gene level - focused on gene ##
      methylgenes<-methyl[which(methyl$Gene_symbol!="-")]
      methylgenes<-methylgenes[,lapply(.SD,mean),by=.(icd10,Gene_symbol), .SDcols=c("T_statistics")]
      thethresholdsgenes<-rbind(thethresholdsgenes,cbind(thresholds[a],as.data.frame(table(methylgenes$icd10))))
    }
    ## cgs level ##
    colnames(thethresholdscgs)<-c("thresholds","icd10","cg")
    thethresholdscgs$thresholds <- as.factor(thethresholdscgs$thresholds)
    dcgs <- ggplot(thethresholdscgs, aes(x=thresholds, y=cg, fill=thresholds)) + 
      geom_violin(trim=FALSE) +
      theme_minimal() +
      geom_boxplot(width=0.1, fill="white") +
      labs(x="p-values", y = "significant cgs") +
      ggtitle("number of significant cgs by disease \n(with different p-value thresholds)") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_fill_brewer(palette="Blues")
    ## gene level - cgs ##
    colnames(thethresholdscggenes)<-c("thresholds","icd10","cg")
    thethresholdscggenes$thresholds <- as.factor(thethresholdscggenes$thresholds)
    dgenescg <- ggplot(thethresholdscggenes, aes(x=thresholds, y=cg, fill=thresholds)) + 
      geom_violin(trim=FALSE) +
      theme_minimal() +
      geom_boxplot(width=0.1, fill="white") +
      labs(x="p-values", y = "significant cgs (genes)") +
      ggtitle("number of significant cgs (genes) by disease \n(with different p-value thresholds)") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_fill_brewer(palette="Blues")
    ## gene level ##
    colnames(thethresholdsgenes)<-c("thresholds","icd10","cg")
    thethresholdsgenes$thresholds <- as.factor(thethresholdsgenes$thresholds)
    dgenes <- ggplot(thethresholdsgenes, aes(x=thresholds, y=cg, fill=thresholds)) + 
      geom_violin(trim=FALSE) +
      theme_minimal() +
      geom_boxplot(width=0.1, fill="white") +
      labs(x="p-values", y = "significant cgs (genes)") +
      ggtitle("number of significant genes by disease \n(with different p-value thresholds)") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_fill_brewer(palette="Blues")
    pdf(file="Plots/Violin_cgs_by_threshold.pdf")
    print(dcgs)
    print(dgenescg)
    print(dgenes)
    dev.off()
  }
  ## Select only those cgs associated to a gene ##
  methylation<-methylation[which(methylation$Gene_symbol!="-")]
  methyldis<-methylation 
  methyldis$icd10<-as.character(fdistab[methylation$EWAS_ID])
  methyldis<-methyldis[,c(1:6,11,12)]
  colnames(methyldis)<-c("disease","ewas_id","cpg","chr","position","gene","t-value","p-value")
  methyldis<-methyldis[which(methyldis$`p-value`<=1e-7)]
  write.table(methyldis,"WorkFiles/DiseaseName_Omics/DiseaseName_methylome.txt",quote=F,sep="\t",row.names=F)
  # ## No filter ##
  # methyl<-methylation[,lapply(.SD,mean),by=.(icd10,ID_REF), .SDcols=c("T_statistics")]
  # ifs_f<-1/table(methyl$ID_REF)
  # fs_f<-table(methyl$ID_REF)/max(table(methyl$ID_REF))
  # methyl<-cbind(methyl,as.numeric(ifs_f[methyl$ID_REF]),as.numeric(fs_f[methyl$ID_REF]))
  # colnames(methyl)<-c("icd10","methylome1e3","tvalue","ifs","fs")
  # methyl<-methyl[,c(1,2,4,5,3)]
  # write.table(methyl,"Omics/ICD10_methylome1e3_associations.txt",quote=F,sep="\t",row.names = F)
  ## Filter disease-cg associations by a p-value 1e-7
  methylation<-methylation[which(methylation$T_Pvalue<=1e-7)]
  methyl<-methylation[,lapply(.SD,mean),by=.(icd10,ID_REF), .SDcols=c("T_statistics")]
  ifs_f<-1/table(methyl$ID_REF)
  fs_f<-table(methyl$ID_REF)/max(table(methyl$ID_REF))
  methyl<-cbind(methyl,as.numeric(ifs_f[methyl$ID_REF]),as.numeric(fs_f[methyl$ID_REF]))
  colnames(methyl)<-c("icd10","methylome","tvalue","ifs","fs")
  methyl<-methyl[,c(1,2,4,5,3)]
  cpgene<-methylation[,c(3,6)] ; if(length(which(duplicated(cpgene)))>0){cpgene<-cpgene[-which(duplicated(cpgene)),]}
  addgenes<-cpgene$Gene_symbol ; names(addgenes)<-cpgene$ID_REF
  methyl<-cbind(methyl,as.character(addgenes[methyl$methylome1e7]))
  colnames(methyl)[6]<-"gene"
  write.table(methyl,"Omics/ICD10_methylome1e7_associations.txt",quote=F,sep="\t",row.names = F)
  write.table(methyl,"WorkFiles/ICD10_Omics/ICD10_methylome.txt",quote=F,sep="\t",row.names = F)
}

#### generate icd - omic variable matrix ####
if(args[1]=="generate_omic_icd_matrixes"){
  # classification<-"ICD10"
  # input<-"miRNA"
  # valor<-"ifs"
  icd_omic_matrix_generator<-function(classification,input,valor){
    ini<-Sys.time()
    associations<-fread(paste("WorkFiles/",classification,"_Omics/",classification,"_",input,".txt",sep=""),stringsAsFactors = F)
    assocs<-as.data.table(cbind(associations[[1]],associations[[2]],associations[[which(colnames(associations)==valor)]])) ; assocs$V3<-as.numeric(assocs$V3)
    colnames(assocs)<-c("icd","omic","value")
    ## First create a matrix with the value of each omic variable in each ICD ##
    diseases<-unique(assocs$icd) ; variables<-unique(assocs$omic)
    icdomic_matrix<-matrix(ncol=length(diseases),nrow=length(variables),0) ; colnames(icdomic_matrix)<-diseases ; rownames(icdomic_matrix)<-variables
    ## If we have less diseases than variables ##
    if(length(diseases)<=length(variables)){
      for(a in 1:length(diseases)){icdomic_matrix[assocs$omic[which(assocs$icd==diseases[a])],diseases[a]]<-assocs$value[which(assocs$icd==diseases[a])]}
    }
    ## If we have more diseases than variables ##
    if(length(diseases)>length(variables)){
      for(a in 1:length(variables)){icdomic_matrix[variables[a],assocs$icd[which(assocs$omic==variables[a])]]<-assocs$value[which(assocs$omic==variables[a])]}
    }
    ## Write the generated table ##
    write.table(icdomic_matrix,paste("WorkFiles/",classification,"_Omics_Matrix/",classification,"_",input,"_",valor,".txt",sep=""),quote=F,sep="\t")
    fin<-Sys.time()
    print(paste("generation of the matrix for:",classification,input,valor,"done in:",fin-ini))
  }
  ## Create the matrixes ##
  ## @@ @@ @@ @ @@ @@ @@ ##
  ## mirna ##
  icd_omic_matrix_generator("ICD10","miRNA","tfidf")
  icd_omic_matrix_generator("ICD10","miRNA","ifs")
  icd_omic_matrix_generator("ICD10","miRNA","fs")
  ## mksymptom ##
  icd_omic_matrix_generator("ICD10","mksymptom","tfidf")
  icd_omic_matrix_generator("ICD10","mksymptom","ifs")
  icd_omic_matrix_generator("ICD10","mksymptom","fs")
  # ## symptom ##
  # icd_omic_matrix_generator("ICD10","symptom","tfidf")
  # icd_omic_matrix_generator("ICD10","symptom","ifs")
  # icd_omic_matrix_generator("ICD10","symptom","fs")
  ## gene ##
  icd_omic_matrix_generator("ICD10","gene","ifs")
  icd_omic_matrix_generator("ICD10","gene","fs")
  ## drug ##
  icd_omic_matrix_generator("ICD10","drug","ifs")
  icd_omic_matrix_generator("ICD10","drug","fs")
  ## drugchembl ##
  icd_omic_matrix_generator("ICD10","drugchembl","ifs")
  icd_omic_matrix_generator("ICD10","drugchembl","fs")
  ## drugchembl2 ##
  icd_omic_matrix_generator("ICD10","drugchembl2","ifs")
  icd_omic_matrix_generator("ICD10","drugchembl2","fs")
  ## metabolome ##
  icd_omic_matrix_generator("ICD10","metabolome","ifs")
  icd_omic_matrix_generator("ICD10","metabolome","fs")
  ## methylome ##
  icd_omic_matrix_generator("ICD10","methylome","tvalue")
  icd_omic_matrix_generator("ICD10","methylome","ifs")
  icd_omic_matrix_generator("ICD10","methylome","fs")
  ## microbiome ##
  icd_omic_matrix_generator("ICD10","microbiome","ifs")
  icd_omic_matrix_generator("ICD10","microbiome","fs")
  icd_omic_matrix_generator("ICD10","microbiome","ifsd")
  icd_omic_matrix_generator("ICD10","microbiome","fsd")
}

#### calculate cosine similarities between icds ####
if(args[1]=="calculate_cosinesimilarities"){
  # classification<-"ICD10"
  # input<-"miRNA"
  # valor<-"ifs"
  cosine_similarities<-function(classification,input,valor){
    ini<-Sys.time()
    icdomic_matrix<-read.csv2(paste("WorkFiles/",classification,"_Omics_Matrix/",classification,"_",input,"_",valor,".txt",sep=""),stringsAsFactors = F,sep="\t")
    icdomic_matrix2<-as.matrix(sapply(icdomic_matrix, as.numeric))
    rownames(icdomic_matrix2)<-rownames(icdomic_matrix)
    remove(icdomic_matrix)
    ## Calculate cosine similarities ##
    cosinesimilarities<-cosine(icdomic_matrix2)
    ind <- which(upper.tri(cosinesimilarities, diag = TRUE), arr.ind = TRUE)
    netmode2<-cbind(colnames(cosinesimilarities)[ind[,1]],colnames(cosinesimilarities)[ind[,2]], cosinesimilarities[ind])
    netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
    write.table(netmode2,paste("Similarities/OriginalAssociations/",classification,"_",input,"_",valor,"_cosinesimilarities.txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
    fin<-Sys.time()
    print(paste("Calculating similarities based on ",classification," - ",input," - ",valor," done, it took: ",fin-ini,sep=""))
  }
  ## Calculate cosine similarities with the ifs or the tf-idf ##
  ## mirna ##
  cosine_similarities("ICD10","miRNA","tfidf")
  cosine_similarities("ICD10","miRNA","ifs")
  cosine_similarities("ICD10","miRNA","fs")
  ## symptom ##
  cosine_similarities("ICD10","symptom","tfidf")
  cosine_similarities("ICD10","symptom","ifs")
  cosine_similarities("ICD10","symptom","fs")
  ## gene ##
  cosine_similarities("ICD10","gene","ifs")
  cosine_similarities("ICD10","gene","fs")
  ## drug ##
  cosine_similarities("ICD10","drug","ifs")
  cosine_similarities("ICD10","drug","fs")
  ## drugchembl ##
  cosine_similarities("ICD10","drugchembl","ifs")
  cosine_similarities("ICD10","drugchembl","fs")
  ## drugchembl2 ##
  cosine_similarities("ICD10","drugchembl2","ifs")
  cosine_similarities("ICD10","drugchembl2","fs")
  ## metabolome ##
  cosine_similarities("ICD10","metabolome","ifs")
  cosine_similarities("ICD10","metabolome","fs")
  ## methylome1e7 ##
  cosine_similarities("ICD10","methylome","tvalue")
  cosine_similarities("ICD10","methylome","ifs")
  cosine_similarities("ICD10","methylome","fs")
  ## microbiome ##
  cosine_similarities("ICD10","microbiome","ifs")
  cosine_similarities("ICD10","microbiome","fs")
  cosine_similarities("ICD10","microbiome","ifsd")
  cosine_similarities("ICD10","microbiome","fsd")
}

#### calculate jaccard indexes ####
if(args[1]=="calculate_fuzzyjaccards"){
  # classification<-"ICD10"
  # input<-"mirna"
  # valor<-"ifs"
  calculatejaccards<-function(classification,input,valor){
    ini<-Sys.time()
    icdomic_matrix<-read.csv2(paste("WorkFiles/",classification,"_Omics_Matrix/",classification,"_",input,"_",valor,".txt",sep=""),stringsAsFactors = F,sep="\t")
    icdomic_matrix2<-as.matrix(sapply(icdomic_matrix, as.numeric))
    rownames(icdomic_matrix2)<-rownames(icdomic_matrix)
    remove(icdomic_matrix)
    ## Get the number of comparisons that will be made ##
    ncompas<-((length(icdomic_matrix2[1,])*length(icdomic_matrix2[1,]))-length(icdomic_matrix2[1,]))/2
    normaljaccard<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
    ## Start a counter to register the loops ##
    contador<-c()
    for(a in 1:(length(icdomic_matrix2[1,])-1)){
      # a<-1
      for(b in (a+1):length(icdomic_matrix2[1,])){
        # b<-3
        contador<-sum(contador,1)
        submat<-icdomic_matrix2[,c(a,b)]
        theintersection<-intersect(which(submat[,1]>0),which(submat[,2]>0))
        theunion<-unique(c(which(submat[,1]>0),which(submat[,2]>0)))
        normaljaccard[contador,]<-c(colnames(submat),length(theintersection)/length(theunion))
        if(length(theintersection)==0){
          fuzzyjaccard[contador,]<-c(colnames(submat),0)
        }
        fuzzyjaccard[contador,]<-c(colnames(submat),sum(apply(submat,1,min))/sum(apply(submat,1,max)))
      }
      print(paste(round((a/(length(icdomic_matrix2[1,])-1))*100,2),"%",sep=""))
    }
    write.table(normaljaccard,paste("Similarities/OriginalAssociations/",classification,"_",input,"_",valor,"_jaccardindex.txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
    write.table(fuzzyjaccard,paste("Similarities/OriginalAssociations/",classification,"_",input,"_",valor,"_weightedjaccardindex.txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
    fin<-Sys.time()
    print(paste("Calculating similarities based on ",input," done, it took: ",fin-ini,sep=""))
  }
  ## Calculate weighted and normal jaccard indexes ##
  ## mirna ##
  calculatejaccards("ICD10","miRNA","tfidf")
  calculatejaccards("ICD10","miRNA","ifs")
  calculatejaccards("ICD10","miRNA","fs")
  ## symptom ##
  calculatejaccards("ICD10","symptom","tfidf")
  calculatejaccards("ICD10","symptom","ifs")
  calculatejaccards("ICD10","symptom","fs")
  ## gene ##
  calculatejaccards("ICD10","gene","ifs")
  calculatejaccards("ICD10","gene","fs")
  ## drug ##
  calculatejaccards("ICD10","drug","ifs")
  calculatejaccards("ICD10","drug","fs")
  ## drugchembl ##
  calculatejaccards("ICD10","drugchembl","ifs")
  calculatejaccards("ICD10","drugchembl","fs")
  ## drugchembl2 ##
  calculatejaccards("ICD10","drugchembl2","ifs")
  calculatejaccards("ICD10","drugchembl2","fs")
  ## metabolome ##
  calculatejaccards("ICD10","metabolome","ifs")
  calculatejaccards("ICD10","metabolome","fs")
  ## methylome1e7 ##
  calculatejaccards("ICD10","methylome","ifs")
  calculatejaccards("ICD10","methylome","fs")
  ## microbiome ##
  calculatejaccards("ICD10","microbiome","ifs")
  calculatejaccards("ICD10","microbiome","fs")
}

#### calculate "random" cosine similarities shuffling the disease-variable associations ####
if(args[1]=="random_cosinesimilarities"){
  # classification<-"ICD10"
  # input<-"mksymptom"
  # valor<-"ifs"
  # semilla<-1
  # randomnum<-200
  random_distances<-function(classification,input,valor,semilla,randomnum){
    inicio<-Sys.time()
    set.seed(semilla)
    icdomic_matrix<-read.csv2(paste("WorkFiles/",classification,"_Omics_Matrix/",classification,"_",input,"_",valor,".txt",sep=""),stringsAsFactors = F,sep="\t")
    icdomic_matrix2<-as.matrix(sapply(icdomic_matrix, as.numeric))
    rownames(icdomic_matrix2)<-rownames(icdomic_matrix)
    remove(icdomic_matrix)
    for(a in 1:randomnum){
      ini<-Sys.time()
      ricdomic_matrix<-apply(icdomic_matrix2,2,sample)
      rcosinesimilarities<-cosine(ricdomic_matrix)
      ind <- which(upper.tri(rcosinesimilarities, diag = TRUE), arr.ind = TRUE)
      netmode2<-cbind(colnames(rcosinesimilarities)[ind[,1]],colnames(rcosinesimilarities)[ind[,2]], rcosinesimilarities[ind])
      netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
      write.table(netmode2,paste("Similarities/RandomAssociations/",classification,"_",input,"_",valor,"_cosinesimilarities/Seed_",semilla,"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
      fin<-Sys.time()
      print(paste("Seed ",semilla," round ",a," time: ",fin-ini,sep=""))
    }
    final<-Sys.time()
    print(final-inicio)
  }
  random_distances(args[2],args[3],args[4],args[5],args[6])
}

#### calculate "random" jaccard indexes shuffling the disease-variable associations ####
if(args[1]=="random_jaccardindexes"){
  # classification<-"ICD10"
  # input<-"mirna"
  # valor<-"ifs"
  # semilla<-1
  # randomnum<-10000
  random_distances<-function(classification,input,valor,semilla,randomnum){
    inicio<-Sys.time()
    set.seed(semilla)
    icdomic_matrix<-read.csv2(paste("WorkFiles/",classification,"_Omics_Matrix/",classification,"_",input,"_",valor,".txt",sep=""),stringsAsFactors = F,sep="\t")
    icdomic_matrix2<-as.matrix(sapply(icdomic_matrix, as.numeric))
    rownames(icdomic_matrix2)<-rownames(icdomic_matrix)
    remove(icdomic_matrix)
    for(a in 1:randomnum){
      ini<-Sys.time()
      ricdomic_matrix<-apply(icdomic_matrix2,2,sample)
      
      ## Get the number of comparisons that will be made ##
      ncompas<-((length(ricdomic_matrix[1,])*length(ricdomic_matrix[1,]))-length(ricdomic_matrix[1,]))/2
      normaljaccard<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
      ## Start a counter to register the loops ##
      contador<-c()
      for(d in 1:(length(ricdomic_matrix[1,])-1)){
        # d<-1
        for(e in (d+1):length(ricdomic_matrix[1,])){
          # e<-3
          contador<-sum(contador,1)
          submat<-ricdomic_matrix[,c(d,e)]
          theintersection<-intersect(which(submat[,1]>0),which(submat[,2]>0))
          theunion<-unique(c(which(submat[,1]>0),which(submat[,2]>0)))
          normaljaccard[contador,]<-c(colnames(submat),length(theintersection)/length(theunion))
          if(length(theintersection)==0){
            fuzzyjaccard[contador,]<-c(colnames(submat),0)
          }
          fuzzyjaccard[contador,]<-c(colnames(submat),sum(apply(submat,1,min))/sum(apply(submat,1,max)))
        }
        print(paste(round((d/(length(ricdomic_matrix[1,])-1))*100,2),"%",sep=""))
      }
      write.table(normaljaccard,paste("Similarities/RandomAssociations/",classification,"_",input,"_",valor,"_jaccardindex/NSeed_",semilla,"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
      write.table(fuzzyjaccard,paste("Similarities/RandomAssociations/",classification,"_",input,"_",valor,"_weightedjaccardindex/NSeed_",semilla,"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
      fin<-Sys.time()
      print(paste("Seed ",semilla," round ",a," time: ",fin-ini,sep=""))
    }
    final<-Sys.time()
    print(final-inicio)
  }
  random_distances(args[2],args[3],args[4],args[5],args[6])
}

#### New Symptoms layer - MK - calculate random TFIDF for disease-symptom associations by shuffling disease - paper associations ####
if(args[1]=="generate_random_tfidf"){
  if("Random_TFIDF"%in%list.files("Data/Symptoms_MK_31032025/")==FALSE){dir.create("Data/Symptoms_MK_31032025/Random_TFIDF")}
  if("Random_Omics_Matrix"%in%list.files("Data/Symptoms_MK_31032025/")==FALSE){dir.create("Data/Symptoms_MK_31032025/Random_Omics_Matrix")}
  ## Define the function for calculating the TF-IDF ##
  tfidf<-function(occurrence,diseaseswithsymptom,numberdiseases){
    rtfidf<-occurrence*log(numberdiseases/diseaseswithsymptom)
    return(rtfidf)
  }
  ## Calculate the number of shared symptoms ##
  ini<-Sys.time()
  sympmid<-fread("Data/Symptoms_MK_31032025/Symptom_PMIDs.txt",stringsAsFactors = F,sep="\t",header=F) ; colnames(sympmid)<-c("symptom","pubmedid") ; setkey(sympmid,"symptom")
  icdpmid<-fread("Data/Symptoms_MK_31032025/ICD10_PMIDs.txt",stringsAsFactors = F,sep="\t",header=F) ; colnames(icdpmid)<-c("icd10","pubmedid") ; setkey(icdpmid,"icd10")
  icduniverse<-unique(icdpmid$pubmedid)
  icds<-table(icdpmid$icd10) ; symptoms<-unique(sympmid$symptom)
  set.seed(as.numeric(args[2]))
  for(z in 1:20){
    cooccurrences<-matrix(nrow=length(icds)*length(symptoms),ncol=5,NA)
    counter<-0
    for(a in 1:length(icds)){
      # a<-1
      idicds<-sample(icduniverse,as.numeric(icds[a]))
      for(b in 1:length(symptoms)){
        counter<-counter+1
        # b<-1
        idsymp<-sympmid[symptoms[b]]$pubmedid
        cooccurrences[counter,]<-c(names(icds)[a],symptoms[b],length(intersect(idicds,idsymp)),length(idicds),length(idsymp))
      }
      print(paste(round((a/length(icds))*100,2),"%",sep=""))
    }
    cooccurrences<-as.data.table(cooccurrences)
    colnames(cooccurrences)<-c("icd10","symptom","occurrence","n_icd10","n_symptom")
    cooccurrences$occurrence<-as.numeric(cooccurrences$occurrence)
    ## First get the number of ICD10 associated to each symptom ##
    Ndiseasesbysymptom<-table(cooccurrences[which(cooccurrences$occurrence>0)]$symptom)
    ## Now go disease by disease ##
    theicds<-unique(cooccurrences$icd10)
    icdswithtfidf<-c()
    for(a in 1:length(theicds)){
      # a<-1
      thedisease<-cooccurrences[which(cooccurrences$icd10==theicds[a])]
      thetfidfs<-c()
      for(b in 1:length(thedisease$icd10)){
        thevalue<-tfidf(thedisease$occurrence[b],Ndiseasesbysymptom[thedisease$symptom[b]],length(theicds))
        thetfidfs<-c(thetfidfs,as.numeric(thevalue))
      }
      thedisease<-cbind(thedisease,thetfidfs)
      icdswithtfidf<-rbind(icdswithtfidf,thedisease)
      print(paste(round((a/length(theicds))*100,2),"%",sep=""))
    }
    icdswithtfidf<-icdswithtfidf[,c(1,2,6)]
    ## Consider only those real icd-disease associations for the fs and ifs calculations ##
    numsympassocdis<-table(cooccurrences$symptom[cooccurrences$occurrence!=0])
    pesos<-numsympassocdis[icdswithtfidf$symptom]
    ifs_f<-1/pesos ; ifs_f[which(cooccurrences$occurrence==0)]<-0
    # ifs_f<-1/table(icdswithtfidf$symptom)
    fs_f<-pesos ; fs_f[which(cooccurrences$occurrence==0)]<-0
    # fs_f<-table(icdswithtfidf$symptom)/max(table(icdswithtfidf$symptom))
    icdswithtfidf<-cbind(icdswithtfidf,as.numeric(ifs_f),as.numeric(fs_f))
    # icdswithtfidf<-cbind(icdswithtfidf,as.numeric(ifs_f[icdswithtfidf$symptom]),as.numeric(fs_f[icdswithtfidf$symptom]))
    colnames(icdswithtfidf)[3:5]<-c("tfidf","ifs","fs")
    ## Write the table containing the number of ICD10 - symptom occurrences with the associated TF-IDF ##
    icdswithtfidf$symptom<-tolower(icdswithtfidf$symptom)
    write.table(icdswithtfidf,paste("Data/Symptoms_MK_31032025/Random_TFIDF/Random_TFIDF_",as.numeric(args[2]),"_",z,".txt",sep=""),quote=F,sep="\t",row.names=F)
    ## convert them into a matrix ##
    valores<-c("tfidf","ifs","fs")
    for(y in 1:length(valores)){
      # y<-1
      assocs<-as.data.table(cbind(icdswithtfidf[[1]],icdswithtfidf[[2]],icdswithtfidf[[which(colnames(icdswithtfidf)==valores[y])]])) ; assocs$V3<-as.numeric(assocs$V3)
      colnames(assocs)<-c("icd","omic","value")
      ## First create a matrix with the value of each omic variable in each ICD ##
      diseases<-unique(assocs$icd) ; variables<-unique(assocs$omic)
      icdomic_matrix<-matrix(ncol=length(diseases),nrow=length(variables),0) ; colnames(icdomic_matrix)<-diseases ; rownames(icdomic_matrix)<-variables
      ## If we have less diseases than variables ##
      if(length(diseases)<=length(variables)){
        for(a in 1:length(diseases)){icdomic_matrix[assocs$omic[which(assocs$icd==diseases[a])],diseases[a]]<-assocs$value[which(assocs$icd==diseases[a])]}
      }
      ## If we have more diseases than variables ##
      if(length(diseases)>length(variables)){
        for(a in 1:length(variables)){icdomic_matrix[variables[a],assocs$icd[which(assocs$omic==variables[a])]]<-assocs$value[which(assocs$omic==variables[a])]}
      }
      ## Write the generated table ##
      write.table(icdomic_matrix,paste("Data/Symptoms_MK_31032025/Random_Omics_Matrix/Random_TFIDF_",as.numeric(args[2]),"_",z,"_",valores[y],".txt",sep=""),quote=F,sep="\t")
    }
    fin<-Sys.time()
    print(paste("generation of the random associations done in:",fin-ini))
  }
}

#### New Symptoms layer - MK - calculate cosine similarities between icds ####
if(args[1]=="symptom_mk_calculate_cosinesimilarities"){
  ## Calculate cosine similarities with the ifs or the tf-idf ##
  if(args[2]=="original"){
    # classification<-"ICD10"
    # input<-"mksymptom"
    # valor<-"fs"
    cosine_similarities<-function(classification,input,valor){
      ini<-Sys.time()
      icdomic_matrix<-read.csv2(paste("WorkFiles/",classification,"_Omics_Matrix/",classification,"_",input,"_",valor,".txt",sep=""),stringsAsFactors = F,sep="\t")
      icdomic_matrix2<-as.matrix(sapply(icdomic_matrix, as.numeric))
      rownames(icdomic_matrix2)<-rownames(icdomic_matrix)
      remove(icdomic_matrix)
      ## Calculate cosine similarities ##
      cosinesimilarities<-cosine(icdomic_matrix2)
      ind <- which(upper.tri(cosinesimilarities, diag = TRUE), arr.ind = TRUE)
      netmode2<-cbind(colnames(cosinesimilarities)[ind[,1]],colnames(cosinesimilarities)[ind[,2]], cosinesimilarities[ind])
      netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
      write.table(netmode2,paste("Similarities/OriginalAssociations/",classification,"_",input,"_",valor,"_cosinesimilarities.txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
      fin<-Sys.time()
      print(paste("Calculating similarities based on ",classification," - ",input," - ",valor," done, it took: ",fin-ini,sep=""))
    }
    cosine_similarities("ICD10","mksymptom","tfidf")
    cosine_similarities("ICD10","mksymptom","ifs")
    cosine_similarities("ICD10","mksymptom","fs")
  }
  if(args[2]=="random"){
    if("Random_similarities"%in%list.files("Data/Symptoms_MK_31032025")==FALSE){dir.create("Data/Symptoms_MK_31032025/Random_similarities")}
    # input<-"Random_TFIDF_9_78_ifs.txt"
    rcosine_similarities<-function(input){
      ini<-Sys.time()
      icdomic_matrix<-read.csv2(paste("Data/Symptoms_MK_31032025/Random_Omics_Matrix/",input,sep=""),stringsAsFactors = F,sep="\t")
      icdomic_matrix2<-as.matrix(sapply(icdomic_matrix, as.numeric))
      rownames(icdomic_matrix2)<-rownames(icdomic_matrix)
      remove(icdomic_matrix)
      ## Calculate cosine similarities ##
      cosinesimilarities<-cosine(icdomic_matrix2)
      ind <- which(upper.tri(cosinesimilarities, diag = TRUE), arr.ind = TRUE)
      netmode2<-cbind(colnames(cosinesimilarities)[ind[,1]],colnames(cosinesimilarities)[ind[,2]], cosinesimilarities[ind])
      netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
      write.table(netmode2,paste("Data/Symptoms_MK_31032025/Random_similarities/",gsub(".txt","_cosinesimilarities.txt",input),sep=""),quote=F,sep="\t",row.names = F,col.names = F)
      fin<-Sys.time()
      print(paste("Calculating similarities based on ",gsub(".txt","",input)," done, it took: ",fin-ini,sep=""))
    }
    files<-list.files("Data/Symptoms_MK_31032025/Random_Omics_Matrix/")
    files<-files[as.numeric(args[3]):as.numeric(args[4])]
    for(a in 1:length(files)){
      rcosine_similarities(files[a])
    }
  }
}

#### New Symptoms layer - MK - calculate jaccard indexes between icds ####
if(args[1]=="symptom_mk_calculate_fuzzyjaccards"){
  if(args[2]=="original"){
    # classification<-"ICD10"
    # input<-"mirna"
    # valor<-"ifs"
    calculatejaccards<-function(classification,input,valor){
      ini<-Sys.time()
      icdomic_matrix<-read.csv2(paste("WorkFiles/",classification,"_Omics_Matrix/",classification,"_",input,"_",valor,".txt",sep=""),stringsAsFactors = F,sep="\t")
      icdomic_matrix2<-as.matrix(sapply(icdomic_matrix, as.numeric))
      rownames(icdomic_matrix2)<-rownames(icdomic_matrix)
      remove(icdomic_matrix)
      ## Get the number of comparisons that will be made ##
      ncompas<-((length(icdomic_matrix2[1,])*length(icdomic_matrix2[1,]))-length(icdomic_matrix2[1,]))/2
      normaljaccard<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
      ## Start a counter to register the loops ##
      contador<-c()
      for(a in 1:(length(icdomic_matrix2[1,])-1)){
        # a<-1
        for(b in (a+1):length(icdomic_matrix2[1,])){
          # b<-3
          contador<-sum(contador,1)
          submat<-icdomic_matrix2[,c(a,b)]
          theintersection<-intersect(which(submat[,1]>0),which(submat[,2]>0))
          theunion<-unique(c(which(submat[,1]>0),which(submat[,2]>0)))
          normaljaccard[contador,]<-c(colnames(submat),length(theintersection)/length(theunion))
          if(length(theintersection)==0){
            fuzzyjaccard[contador,]<-c(colnames(submat),0)
          }
          fuzzyjaccard[contador,]<-c(colnames(submat),sum(apply(submat,1,min))/sum(apply(submat,1,max)))
        }
        print(paste(round((a/(length(icdomic_matrix2[1,])-1))*100,2),"%",sep=""))
      }
      write.table(normaljaccard,paste("Similarities/OriginalAssociations/",classification,"_",input,"_",valor,"_jaccardindex.txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
      write.table(fuzzyjaccard,paste("Similarities/OriginalAssociations/",classification,"_",input,"_",valor,"_weightedjaccardindex.txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
      fin<-Sys.time()
      print(paste("Calculating similarities based on ",input," done, it took: ",fin-ini,sep=""))
    }
    ## Calculate weighted and normal jaccard indexes ##
    ## symptom ##
    calculatejaccards("ICD10","mksymptom","tfidf")
    calculatejaccards("ICD10","mksymptom","ifs")
    calculatejaccards("ICD10","mksymptom","fs")
  }
  if(args[2]=="random"){
    if("Random_similarities"%in%list.files("Data/Symptoms_MK_31032025")==FALSE){dir.create("Data/Symptoms_MK_31032025/Random_similarities")}
    # input<-"Random_TFIDF_9_78_ifs.txt"
    rcalculatejaccards<-function(input){
      ini<-Sys.time()
      icdomic_matrix<-read.csv2(paste("Data/Symptoms_MK_31032025/Random_Omics_Matrix/",input,sep=""),stringsAsFactors = F,sep="\t")
      icdomic_matrix2<-as.matrix(sapply(icdomic_matrix, as.numeric))
      rownames(icdomic_matrix2)<-rownames(icdomic_matrix)
      remove(icdomic_matrix)
      ## Get the number of comparisons that will be made ##
      ncompas<-((length(icdomic_matrix2[1,])*length(icdomic_matrix2[1,]))-length(icdomic_matrix2[1,]))/2
      normaljaccard<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
      ## Start a counter to register the loops ##
      contador<-c()
      for(a in 1:(length(icdomic_matrix2[1,])-1)){
        # a<-1
        for(b in (a+1):length(icdomic_matrix2[1,])){
          # b<-3
          contador<-sum(contador,1)
          submat<-icdomic_matrix2[,c(a,b)]
          theintersection<-intersect(which(submat[,1]>0),which(submat[,2]>0))
          theunion<-unique(c(which(submat[,1]>0),which(submat[,2]>0)))
          normaljaccard[contador,]<-c(colnames(submat),length(theintersection)/length(theunion))
          if(length(theintersection)==0){
            fuzzyjaccard[contador,]<-c(colnames(submat),0)
          }
          fuzzyjaccard[contador,]<-c(colnames(submat),sum(apply(submat,1,min))/sum(apply(submat,1,max)))
        }
        print(paste(round((a/(length(icdomic_matrix2[1,])-1))*100,2),"%",sep=""))
      }
      write.table(normaljaccard,paste("Data/Symptoms_MK_31032025/Random_similarities/",gsub(".txt","_jaccardindex.txt",input),sep=""),quote=F,sep="\t",row.names = F,col.names = F)
      write.table(fuzzyjaccard,paste("Data/Symptoms_MK_31032025/Random_similarities/",gsub(".txt","_weightedjaccardindex.txt",input),sep=""),quote=F,sep="\t",row.names = F,col.names = F)
      fin<-Sys.time()
      print(paste("Calculating similarities based on ",gsub(".txt","",input)," done, it took: ",fin-ini,sep=""))
    }
    ## Calculate weighted and normal jaccard indexes ##
    ## symptom ##
    files<-list.files("Data/Symptoms_MK_31032025/Random_Omics_Matrix/")
    files<-files[as.numeric(args[3]):as.numeric(args[4])]
    for(a in 1:length(files)){
      rcalculatejaccards(files[a])
    }
  }
}

#### get p-values for the symptoms - 1 ####
if(args[1]=="get_pvals_for_symptoms"){
  coords<-cbind(seq(1,10000,1000),seq(0,10000,1000)[-1])
  archivos <- list.files(path = "Data/Symptoms_MK_31032025/Random_similarities", pattern = "\\.txt$", full.names = TRUE)
  archivos <- archivos[grep(gsub(".txt","",gsub("ICD10_mksymptom","",args[2])),archivos)]
  for(a in 1:length(coords[,1])){
    subarch<-archivos[coords[a,1]:coords[a,2]]
    # Leer solo la tercera columna de cada archivo
    columnas <- lapply(subarch, function(f) fread(f, select = 3)[[1]])
    # Combinar todas las columnas en un data.table o data.frame
    resultado <- as.data.table(columnas)
    write.table(resultado,paste("Data/Symptoms_MK_31032025/Random_",gsub(".txt","",gsub("ICD10_mksymptom_","",args[2])),"_",a,".txt",sep=""),quote = F,sep="\t",row.names=F,col.names=F)
    print(paste(round((a/length(coords[,1]))*100,2),"%",sep=""))
  }
}

#### get p-values for the symptoms - 2 ####
if(args[1]=="get_pvals_for_symptoms_art_2"){
  ## using bash put together all the tables for the same weight and similarity metric ##
  # paste Random_fs_cosinesimilarities_1.txt Random_fs_cosinesimilarities_2.txt Random_fs_cosinesimilarities_3.txt Random_fs_cosinesimilarities_4.txt Random_fs_cosinesimilarities_5.txt Random_fs_cosinesimilarities_6.txt Random_fs_cosinesimilarities_7.txt Random_fs_cosinesimilarities_8.txt Random_fs_cosinesimilarities_9.txt Random_fs_cosinesimilarities_10.txt > Random_fs_cosinesimilarities.txt
  # paste Random_fs_jaccardindex_1.txt Random_fs_jaccardindex_2.txt Random_fs_jaccardindex_3.txt Random_fs_jaccardindex_4.txt Random_fs_jaccardindex_5.txt Random_fs_jaccardindex_6.txt Random_fs_jaccardindex_7.txt Random_fs_jaccardindex_8.txt Random_fs_jaccardindex_9.txt Random_fs_jaccardindex_10.txt > Random_fs_jaccardindex.txt
  # paste Random_fs_weightedjaccardindex_1.txt Random_fs_weightedjaccardindex_2.txt Random_fs_weightedjaccardindex_3.txt Random_fs_weightedjaccardindex_4.txt Random_fs_weightedjaccardindex_5.txt Random_fs_weightedjaccardindex_6.txt Random_fs_weightedjaccardindex_7.txt Random_fs_weightedjaccardindex_8.txt Random_fs_weightedjaccardindex_9.txt Random_fs_weightedjaccardindex_10.txt > Random_fs_weightedjaccardindex.txt
  
  # paste Random_ifs_cosinesimilarities_1.txt Random_ifs_cosinesimilarities_2.txt Random_ifs_cosinesimilarities_3.txt Random_ifs_cosinesimilarities_4.txt Random_ifs_cosinesimilarities_5.txt Random_ifs_cosinesimilarities_6.txt Random_ifs_cosinesimilarities_7.txt Random_ifs_cosinesimilarities_8.txt Random_ifs_cosinesimilarities_9.txt Random_ifs_cosinesimilarities_10.txt > Random_ifs_cosinesimilarities.txt
  # paste Random_ifs_jaccardindex_1.txt Random_ifs_jaccardindex_2.txt Random_ifs_jaccardindex_3.txt Random_ifs_jaccardindex_4.txt Random_ifs_jaccardindex_5.txt Random_ifs_jaccardindex_6.txt Random_ifs_jaccardindex_7.txt Random_ifs_jaccardindex_8.txt Random_ifs_jaccardindex_9.txt Random_ifs_jaccardindex_10.txt > Random_ifs_jaccardindex.txt
  # paste Random_ifs_weightedjaccardindex_1.txt Random_ifs_weightedjaccardindex_2.txt Random_ifs_weightedjaccardindex_3.txt Random_ifs_weightedjaccardindex_4.txt Random_ifs_weightedjaccardindex_5.txt Random_ifs_weightedjaccardindex_6.txt Random_ifs_weightedjaccardindex_7.txt Random_ifs_weightedjaccardindex_8.txt Random_ifs_weightedjaccardindex_9.txt Random_ifs_weightedjaccardindex_10.txt > Random_ifs_weightedjaccardindex.txt
  
  # paste Random_tfidf_cosinesimilarities_1.txt Random_tfidf_cosinesimilarities_2.txt Random_tfidf_cosinesimilarities_3.txt Random_tfidf_cosinesimilarities_4.txt Random_tfidf_cosinesimilarities_5.txt Random_tfidf_cosinesimilarities_6.txt Random_tfidf_cosinesimilarities_7.txt Random_tfidf_cosinesimilarities_8.txt Random_tfidf_cosinesimilarities_9.txt Random_tfidf_cosinesimilarities_10.txt > Random_tfidf_cosinesimilarities.txt
  # paste Random_tfidf_jaccardindex_1.txt Random_tfidf_jaccardindex_2.txt Random_tfidf_jaccardindex_3.txt Random_tfidf_jaccardindex_4.txt Random_tfidf_jaccardindex_5.txt Random_tfidf_jaccardindex_6.txt Random_tfidf_jaccardindex_7.txt Random_tfidf_jaccardindex_8.txt Random_tfidf_jaccardindex_9.txt Random_tfidf_jaccardindex_10.txt > Random_tfidf_jaccardindex.txt
  # paste Random_tfidf_weightedjaccardindex_1.txt Random_tfidf_weightedjaccardindex_2.txt Random_tfidf_weightedjaccardindex_3.txt Random_tfidf_weightedjaccardindex_4.txt Random_tfidf_weightedjaccardindex_5.txt Random_tfidf_weightedjaccardindex_6.txt Random_tfidf_weightedjaccardindex_7.txt Random_tfidf_weightedjaccardindex_8.txt Random_tfidf_weightedjaccardindex_9.txt Random_tfidf_weightedjaccardindex_10.txt > Random_tfidf_weightedjaccardindex.txt
  
  ## Create a table with one column showing the real values ##
  # tt<-fread("Similarities/OriginalAssociations/ICD10_mksymptom_fs_cosinesimilarities.txt",stringsAsFactors = F,sep="\t")$V3
  # write.table(tt,"Data/Symptoms_MK_31032025/Real_fs_cosinesimilarities.txt",quote=F,sep="\t",row.names=F,col.names=F)
  # tt<-fread("Similarities/OriginalAssociations/ICD10_mksymptom_fs_jaccardindex.txt",stringsAsFactors = F,sep="\t")$V3
  # write.table(tt,"Data/Symptoms_MK_31032025/Real_fs_jaccardindex.txt",quote=F,sep="\t",row.names=F,col.names=F)
  # tt<-fread("Similarities/OriginalAssociations/ICD10_mksymptom_fs_weightedjaccardindex.txt",stringsAsFactors = F,sep="\t")$V3
  # write.table(tt,"Data/Symptoms_MK_31032025/Real_fs_weightedjaccardindex.txt",quote=F,sep="\t",row.names=F,col.names=F)
  
  # tt<-fread("Similarities/OriginalAssociations/ICD10_mksymptom_ifs_cosinesimilarities.txt",stringsAsFactors = F,sep="\t")$V3
  # write.table(tt,"Data/Symptoms_MK_31032025/Real_ifs_cosinesimilarities.txt",quote=F,sep="\t",row.names=F,col.names=F)
  # tt<-fread("Similarities/OriginalAssociations/ICD10_mksymptom_ifs_jaccardindex.txt",stringsAsFactors = F,sep="\t")$V3
  # write.table(tt,"Data/Symptoms_MK_31032025/Real_ifs_jaccardindex.txt",quote=F,sep="\t",row.names=F,col.names=F)
  # tt<-fread("Similarities/OriginalAssociations/ICD10_mksymptom_ifs_weightedjaccardindex.txt",stringsAsFactors = F,sep="\t")$V3
  # write.table(tt,"Data/Symptoms_MK_31032025/Real_ifs_weightedjaccardindex.txt",quote=F,sep="\t",row.names=F,col.names=F)
  
  # tt<-fread("Similarities/OriginalAssociations/ICD10_mksymptom_tfidf_cosinesimilarities.txt",stringsAsFactors = F,sep="\t")$V3
  # write.table(tt,"Data/Symptoms_MK_31032025/Real_tfidf_cosinesimilarities.txt",quote=F,sep="\t",row.names=F,col.names=F)
  # tt<-fread("Similarities/OriginalAssociations/ICD10_mksymptom_tfidf_jaccardindex.txt",stringsAsFactors = F,sep="\t")$V3
  # write.table(tt,"Data/Symptoms_MK_31032025/Real_tfidf_jaccardindex.txt",quote=F,sep="\t",row.names=F,col.names=F)
  # tt<-fread("Similarities/OriginalAssociations/ICD10_mksymptom_tfidf_weightedjaccardindex.txt",stringsAsFactors = F,sep="\t")$V3
  # write.table(tt,"Data/Symptoms_MK_31032025/Real_tfidf_weightedjaccardindex.txt",quote=F,sep="\t",row.names=F,col.names=F)
  
  ## Ver el numero de filas ##
  # wc -l Real_fs_cosinesimilarities.txt 
  # wc -l Random_fs_cosinesimilarities.txt 
  
  ## Ver el numero de columnas ##
  # head -n 1 Real_fs_cosinesimilarities.txt | awk -F'\t' '{print NF}'
  # head -n 1 Random_fs_cosinesimilarities.txt | awk -F'\t' '{print NF}'
  
  ## Cosine similarities ##
  ## @@ @@ @@ @ @@ @@ @@ ##
  ## frequency score ##
#   awk -F'\t' '
#   NR == FNR {real[FNR] = $1; next} 
#   {
#     count = 0
#     for(i=1; i<=NF; i++) if ($i >= real[FNR]) count++
#     p_value = (count) / (NF)  # Método Monte Carlo (+1 incluye el valor real)
#     print p_value
#   }
# ' Real_fs_cosinesimilarities.txt Random_fs_cosinesimilarities.txt > pvalues_fs_cosinesimilarities.tsv
  ## inverse frequency score ##
#   awk -F'\t' '
#   NR == FNR {real[FNR] = $1; next} 
#   {
#     count = 0
#     for(i=1; i<=NF; i++) if ($i >= real[FNR]) count++
#     p_value = (count) / (NF)  # Método Monte Carlo (+1 incluye el valor real)
#     print p_value
#   }
# ' Real_ifs_cosinesimilarities.txt Random_ifs_cosinesimilarities.txt > pvalues_ifs_cosinesimilarities.tsv
#   ## term frequency-inverse document frequency
#   awk -F'\t' '
#   NR == FNR {real[FNR] = $1; next} 
#   {
#     count = 0
#     for(i=1; i<=NF; i++) if ($i >= real[FNR]) count++
#     p_value = (count) / (NF)  # Método Monte Carlo (+1 incluye el valor real)
#     print p_value
#   }
# ' Real_tfidf_cosinesimilarities.txt Random_tfidf_cosinesimilarities.txt > pvalues_tfidf_cosinesimilarities.tsv
  
  ## jaccard index ##
  ## @@ @@ @ @@ @@ ##
  ## frequency score ##
#   awk -F'\t' '
#   NR == FNR {real[FNR] = $1; next} 
#   {
#     count = 0
#     for(i=1; i<=NF; i++) if ($i >= real[FNR]) count++
#     p_value = (count) / (NF)  # Método Monte Carlo (+1 incluye el valor real)
#     print p_value
#   }
# ' Real_fs_jaccardindex.txt Random_fs_jaccardindex.txt > pvalues_fs_jaccardindex.tsv
  ## inverse frequency score ##
#   awk -F'\t' '
#   NR == FNR {real[FNR] = $1; next} 
#   {
#     count = 0
#     for(i=1; i<=NF; i++) if ($i >= real[FNR]) count++
#     p_value = (count) / (NF)  # Método Monte Carlo (+1 incluye el valor real)
#     print p_value
#   }
# ' Real_ifs_jaccardindex.txt Random_ifs_jaccardindex.txt > pvalues_ifs_jaccardindex.tsv
  ## term frequency-inverse document frequency
#   awk -F'\t' '
#   NR == FNR {real[FNR] = $1; next} 
#   {
#     count = 0
#     for(i=1; i<=NF; i++) if ($i >= real[FNR]) count++
#     p_value = (count) / (NF)  # Método Monte Carlo (+1 incluye el valor real)
#     print p_value
#   }
# ' Real_tfidf_jaccardindex.txt Random_tfidf_jaccardindex.txt > pvalues_tfidf_jaccardindex.tsv
  
  ## weighted jaccard index ##
  ## @@ @@ @@ @@@@ @@ @@ @@ ##
  ## frequency score ##
#   awk -F'\t' '
#   NR == FNR {real[FNR] = $1; next} 
#   {
#     count = 0
#     for(i=1; i<=NF; i++) if ($i >= real[FNR]) count++
#     p_value = (count) / (NF)  # Método Monte Carlo (+1 incluye el valor real)
#     print p_value
#   }
# ' Real_fs_weightedjaccardindex.txt Random_fs_weightedjaccardindex.txt > pvalues_fs_weightedjaccardindex.tsv
  ## inverse frequency score ##
#   awk -F'\t' '
#   NR == FNR {real[FNR] = $1; next} 
#   {
#     count = 0
#     for(i=1; i<=NF; i++) if ($i >= real[FNR]) count++
#     p_value = (count) / (NF)  # Método Monte Carlo (+1 incluye el valor real)
#     print p_value
#   }
# ' Real_ifs_weightedjaccardindex.txt Random_ifs_weightedjaccardindex.txt > pvalues_ifs_weightedjaccardindex.tsv
  ## term frequency-inverse document frequency
#   awk -F'\t' '
#   NR == FNR {real[FNR] = $1; next} 
#   {
#     count = 0
#     for(i=1; i<=NF; i++) if ($i >= real[FNR]) count++
#     p_value = (count) / (NF)  # Método Monte Carlo (+1 incluye el valor real)
#     print p_value
#   }
# ' Real_tfidf_weightedjaccardindex.txt Random_tfidf_weightedjaccardindex.txt > pvalues_tfidf_weightedjaccardindex.tsv
}

#### calculate the p-value and remove all the individual files ####
if(args[1]=="calculate_pval"){
  if("Similarities_with_pval"%in%list.files("Similarities")==FALSE){dir.create("Similarities/Similarities_with_pval")}
  # classification<-"ICD10"
  # input<-"mksymptom"
  # valor<-"ifs"
  # metric<-"weightedjaccardindex"
  getpval<-function(classification,input,valor,metric){
    thefile<-paste(classification,"_",input,"_",valor,"_",metric,sep="")
    originalnetwork<-fread(paste("Similarities/OriginalAssociations/",thefile,".txt",sep=""),stringsAsFactors=F,sep="\t")
    ## continue only if we haven't created the pvalue file yet ## 
    if(paste(thefile,".txt",sep="")%in%list.files("Similarities/Similarities_with_pval")==F){
      randomfiles<-list.files(paste("Similarities/RandomAssociations/",thefile,sep=""))
      ## continue only if we have 10.000 random associations ##
      if(length(randomfiles)==10000){
        randommatrix<-matrix(ncol=10000,nrow=length(originalnetwork$V1),NA)
        for(b in 1:length(randomfiles)){
          # b<-1
          randomnetwork<-fread(paste("Similarities/RandomAssociations/",thefile,"/",randomfiles[b],sep=""),stringsAsFactors = F,sep="\t")
          randommatrix[,b]<-randomnetwork$V3
          print(paste(round((b/length(randomfiles))*100,2),"%",sep=""))
        }
        resultmatrix<-matrix(ncol=4,nrow=length(originalnetwork$V1),NA)
        for(a in 1:length(originalnetwork$V1)){
          # a<-1
          thepvalup<-length(which(randommatrix[a,]>=originalnetwork$V3[a]))/length(randommatrix[a,])
          thepvaldown<-length(which(randommatrix[a,]<=originalnetwork$V3[a]))/length(randommatrix[a,])
          meandistance<-mean(randommatrix[a,])
          sddistance<-sd(randommatrix[a,])
          resultmatrix[a,]<-c(thepvalup,thepvaldown,meandistance,sddistance)
          print(paste(round((a/length(originalnetwork$V1))*100,2),"%",sep=""))
        }
        resultmatrix<-cbind(originalnetwork,resultmatrix)
        colnames(resultmatrix)<-c("disease1","disease2",metric,"pval_up","pval_down","mean","sd")
        write.table(resultmatrix,paste("Similarities/Similarities_with_pval/",thefile,".txt",sep=""),quote=F,sep="\t",row.names=F)
      }
      if(length(randomfiles)!=10000){print(paste(thefile,"has only",length(randomfiles),"random files"))}
    }
    if(paste(thefile,".txt",sep="")%in%list.files("Similarities/Similarities_with_pval")==T){print("The file has been already created")}
    print("Finished!")
  }
  ## example: getpval("ICD10","drug","ifs","weightedjaccardindex")
  ## get the p-values ##
  ## @@ @@ @@@@ @@ @@ ##
  getpval(args[2],args[3],args[4],args[5])
  # getpval(classification,input,valor,metric)
}

#### calculate FDR (Bonferroni correction) and sort associations ####
if(args[1]=="calculate_fdr"){
  if("Similarities_with_fdr"%in%list.files("Similarities")==FALSE){dir.create("Similarities/Similarities_with_fdr")}
  # classification<-"ICD10"
  # input<-"microbiome"
  # valor<-"ifs"
  # metric<-"weightedjaccardindex"
  getfdr<-function(classification,input,valor,metric){
    thefile<-paste(classification,"_",input,"_",valor,"_",metric,sep="")
    originalnetwork<-fread(paste("Similarities/Similarities_with_pval/",thefile,".txt",sep=""),stringsAsFactors = F,sep="\t")
    sortednet<-t(apply(originalnetwork[,1:2], 1, sort))
    originalnetwork$disease1<-sortednet[,1] ; originalnetwork$disease2<-sortednet[,2]
    fdr_up<-originalnetwork$pval_up*length(originalnetwork$pval_up)
    if(length(which(fdr_up>1))>0){fdr_up[which(fdr_up>1)]<-1}
    fdr_down<-originalnetwork$pval_down*length(originalnetwork$pval_down)
    if(length(which(fdr_down>1))>0){fdr_down[which(fdr_down>1)]<-1}
    originalnetwork<-cbind(originalnetwork,fdr_up,fdr_down)
    colnames(originalnetwork)[6:7]<-c("random_mean","random_sd")
    originalnetwork<-originalnetwork[,c(1:4,8,5,9,6,7)]
    write.table(originalnetwork,paste("Similarities/Similarities_with_fdr/",thefile,".txt",sep=""),quote=F,sep="\t",row.names=F)
    print(paste(input,"finished!"))
  }
  ## get the FDR ##
  ## @@ @@ @@ @@ ##
  files<-list.files("Similarities/Similarities_with_pval/")
  dones<-list.files("Similarities/Similarities_with_fdr/")
  ## run it only for those networks that have not been generated yet ##
  files<-setdiff(files,dones)
  for(a in 1:length(files)){
    # a<-1
    elements<-strsplit(files[a],"_")[[1]]
    getfdr(elements[1],elements[2],elements[3],gsub(".txt","",elements[4]))
    print(paste(round((a/length(files))*100,2),"%",sep=""))
  }
}

#### put, for each omic layer, all disease-disease similarities in the same order #### 
# if(args[1]=="sameorderofdiseasesineachlayer"){
#   if("HDMLN"%in%list.files()==F){dir.create("HDMLN")}
#   if("Intralayers"%in%list.files("HDMLN")==FALSE){dir.create("HDMLN/Intralayers")}
#   if("Interlayers"%in%list.files("HDMLN")==FALSE){dir.create("HDMLN/Interlayers")}
#   files<-list.files("Similarities/Similarities_with_fdr/")
#   omicas<-gsub("_.+","",gsub("ICD10_","",files))
#   uomic<-unique(omicas)
#   tabla<-c()
#   for(a in 1:length(uomic)){
#     # a<-1
#     sonomica<-files[which(omicas==uomic[a])]
#     tablal<-list()
#     pastel<-list()
#     for(b in 1:length(sonomica)){
#       # b<-1
#       tt<-fread(paste("Similarities/Similarities_with_fdr/",sonomica[b],sep=""),stringsAsFactors = F,sep="\t")
#       tt<-cbind(paste(tt$disease1,tt$disease2,sep="_"),tt)
#       setkey(tt,"V1")
#       tablal[[sonomica[b]]]<-tt
#       pastel[[sonomica[b]]]<-tt$V1
#       print(length(unique(tt$V1)))
#     }
#     print(paste("Shared interactions:",length(Reduce(intersect,pastel))))
#     # if(uomic[a]!="epidemiology" && uomic[a]!="transcriptome"){
#     #   theorder<-sort(tablal[[1]]$V1)
#     #   latabla<-c()
#     #   for(b in 1:length(names(tablal))){
#     #     # b<-1
#     #     tt<-tablal[[b]][theorder]
#     #     tt<-tt[,-1]
#     #     write.table(tt,paste("Similarities/Similarities_with_fdr/",sonomica[b],sep=""),quote=F,sep="\t",row.names=F)
#     #     if(b==1){latabla<-tt[,1:5] ; colnames(latabla)[3:5]<-c(gsub(paste(uomic[a],"_",sep=""),"",gsub(".txt","",gsub("ICD10_","",sonomica[b]))),paste(c("pval","fdr"),gsub(paste(uomic[a],"_",sep=""),"",gsub(".txt","",gsub("ICD10_","",sonomica[b]))),sep="_"))}
#     #     if(b>1){
#     #       newtab<-tt[,3:5]
#     #       colnames(newtab)<-c(gsub(paste(uomic[a],"_",sep=""),"",gsub(".txt","",gsub("ICD10_","",sonomica[b]))),paste(c("pval","fdr"),gsub(paste(uomic[a],"_",sep=""),"",gsub(".txt","",gsub("ICD10_","",sonomica[b]))),sep="_"))
#     #       latabla<-cbind(latabla,newtab)
#     #     }
#     #   }
#     #   write.table(latabla,paste("HDMLN/",uomic[a],".txt",sep=""),quote=F,sep="\t",row.names=F)
#     # }
#   }
# }


#### Add the two extra layers previously generated: transcriptomics and epidemiology ####
if(args[1]=="transcriptomicsandepidemiology"){
  ## transcriptomics ##
  ## @@ @@ @ @ @@ @@ ##
  transcrip<-fread("Data/Transcriptomics_JSV/ICD10_Adjusted_pvals_network.txt",stringsAsFactors = F,sep="\t")
  transcrip<-transcrip[,c(1,2,16,17)]
  transcriptomics<-fread("Data/Transcriptomics_JSV/ICD10_Adjusted_fdr_network.txt",stringsAsFactors = F,sep="\t")
  transcriptomics<-transcriptomics[,c(1,2,16,17)]
  transcriptomics<-cbind(transcriptomics$Disease1,transcriptomics$Disease2,transcriptomics$AllCosineSimilarity,transcrip$AllPvalCosineSimilarity,transcriptomics$AllFDRCosineSimilarity)
  colnames(transcriptomics)<-c("disease1","disease2","cosinesimilarities","pval","fdr")
  write.table(transcriptomics,"Similarities/Similarities_with_fdr/ICD10_transcriptome_fc_cosinesimilarities.txt",quote=F,sep="\t",row.names=F)
  ## epidemiology ##
  ## @@ @@  @@ @@ ##
  epidemiology<-fread("Data/Epidemiology_westergaard/41467_2019_8475_MOESM6_ESM.txt",stringsAsFactors = F,sep="\t")
  epidemiology<-epidemiology[,c(1,2,7)]
  epidemiology<-cbind(t(apply(epidemiology[,1:2], 1, sort)),epidemiology$adjustedRR2.5)
  colnames(epidemiology)<-c("disease1","disease2","RR")
  ## we don't take into consideration the NAs ##
  if(length(which(is.na(epidemiology[,3])))>0){epidemiology<-epidemiology[-which(is.na(epidemiology[,3])),]}
  epidemiology<-as.data.table(epidemiology) ; epidemiology$RR<-as.numeric(epidemiology$RR)
  ## Remove duplicated associations selecting always the layers RR
  if(length(which(duplicated(epidemiology)))>0){
    duplis<-epidemiology[which(duplicated(epidemiology))]
    for(a in 1:length(duplis$disease1)){
      # a<-1
      indices<-intersect(which(epidemiology$disease1==duplis$disease1[a]),which(epidemiology$disease2==duplis$disease2[a]))
      epidemiology<-epidemiology[-setdiff(indices,indices[which(epidemiology$RR[indices]==max(epidemiology$RR[indices]))[1]])]
    }
  }
  pasted<-paste(epidemiology$disease1,epidemiology$disease2)
  epidemiology<-epidemiology[order(pasted)]
  write.table(epidemiology,"Similarities/Similarities_with_fdr/ICD10_epidemiology_relativerisk.txt",quote=F,sep="\t",row.names=F)
}
















