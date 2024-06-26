#### With this code we are going to generate the multiomic layers (disease-disease similarities) for the BSC-IGC collaboration ####
## Load the needed libraries ##
library("data.table")
library(VennDiagram)
library(RColorBrewer)
library("igraph")
library("UpSetR")
library("dendextend")
library(lsa)
library("geometry")
library("ggplot2")
library("gplots")

args = commandArgs(trailingOnly=TRUE)

# print(args[1])
# print(args[2])

if("Summaries"%in%list.files()==FALSE){dir.create("Summaries")}
if("Final"%in%list.files("Interlayers")==FALSE){dir.create("Interlayers/Final")}

#### Build all the intralayers ####
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##

#### Symptoms #### 
## Cosine Similarity done ##
## Jaccard index pending ##
if(args[1]=="symptoms"){
  if(args[2]=="getthenewvectors"){
    ## Sum the disease - symptom co-occurrences ##
    if("ICD10_symptom_associations.txt"%in%list.files("Symptoms_Zhou2014_20220712/")==FALSE){
      ## Load the disease name - icd10 association file
      symptomnametoicd<-fread("Symptoms_Zhou2014_20220712/Disease_ICD10_translation.txt",stringsAsFactors = F,sep="\t")
      symptomnametoicd<-symptomnametoicd[which(nchar(symptomnametoicd$ICD10)==3)]
      ## Load the disease-symptom association file
      diseasesymptom<-fread("Symptoms_Zhou2014_20220712/41467_2014_BFncomms5212_MOESM1045_ESM.txt")
      ## Select the diseases for which we have and ICD correspondence and sum their mentions ##
      icds<-unique(symptomnametoicd$ICD10)
      icdsymptom<-c()
      for(a in 1:length(icds)){
        # a<-1
        diseases<-symptomnametoicd$Name[which(symptomnametoicd$ICD10==icds[a])]
        icdtab<-c()
        for(b in 1:length(diseases)){
          # b<-1
          icdtab<-rbind(icdtab,diseasesymptom[which(diseasesymptom$`MeSH Disease Term`==diseases[b])])
        }
        icdtab<-icdtab[,lapply(.SD,sum),by=.(`MeSH Symptom Term`), .SDcols=c("PubMed occurrence")]
        icdtab<-cbind(icds[a],icdtab)
        colnames(icdtab)<-c("icd10","symptom","occurrence")
        icdsymptom<-rbind(icdsymptom,icdtab)
        print(paste(round((a/length(icds))*100,2),"%",sep=""))
      }
      ## Write the table containing the number of ICD10 - symptom occurrences ##
      write.table(icdsymptom,"Symptoms_Zhou2014_20220712/ICD10_symptom_associations.txt",quote=F,sep="\t",row.names=F)
    }
    ## Define the function for calculating the TF-IDF ##
    tfidf<-function(occurrence,diseaseswithsymptom,numberdiseases){
      rtfidf<-occurrence*log(numberdiseases/diseaseswithsymptom)
      return(rtfidf)
    }
    ## Calculate the TF-IDF ##
    if(dim(icdsymptom)[2]==3){
      ## First get the number of ICD10 associated to each symptom ##
      Ndiseasesbysymptom<-table(icdsymptom$symptom)
      ## Now go disease by disease ##
      theicds<-unique(icdsymptom$icd10)
      icdswithtfidf<-c()
      for(a in 1:length(theicds)){
        # a<-1
        thedisease<-icdsymptom[which(icdsymptom$icd10==theicds[a])]
        thetfidfs<-c()
        for(b in 1:length(thedisease$icd10)){
          thevalue<-tfidf(thedisease$occurrence[b],Ndiseasesbysymptom[thedisease$symptom[b]],length(theicds))
          thetfidfs<-c(thetfidfs,as.numeric(thevalue))
        }
        thedisease<-cbind(thedisease,thetfidfs)
        icdswithtfidf<-rbind(icdswithtfidf,thedisease)
        print(paste(round((a/length(theicds))*100,2),"%",sep=""))
      }
      ## Write the table containing the number of ICD10 - symptom occurrences with the associated TF-IDF ##
      write.table(icdswithtfidf,"Symptoms_Zhou2014_20220712/ICD10_symptom_associations.txt",quote=F,sep="\t",row.names=F)
    }
  }
  if(args[2]=="calculatesimilarities"){
    if("RandomAssociations"%in%list.files("Symptoms_Zhou2014_20220712/")==F){dir.create("Symptoms_Zhou2014_20220712/RandomAssociations")}
    associations<-fread("Symptoms_Zhou2014_20220712/ICD10_symptom_associations.txt",stringsAsFactors = F,sep="\t")
    ## First create a matrix with the TF-IDF of each symptom in each disease ##
    diseases<-unique(associations$icd10) ; symptoms<-unique(associations$symptom)
    icdsymptom_matrix<-matrix(ncol=length(diseases),nrow=length(symptoms),0) ; colnames(icdsymptom_matrix)<-diseases ; rownames(icdsymptom_matrix)<-symptoms
    # icdsymptom_number<-icdsymptom_matrix
    for(a in 1:length(symptoms)){
      icdsymptom_matrix[symptoms[a],associations$icd10[which(associations$symptom==symptoms[a])]]<-associations$thetfidfs[which(associations$symptom==symptoms[a])]
      # icdsymptom_number[symptoms[a],associations$icd10[which(associations$symptom==symptoms[a])]]<-associations$occurrence[which(associations$symptom==symptoms[a])]
    }
    if(args[3]=="CreateTFIDFmatrix"){
      write.table(icdsymptom_matrix,"Symptoms_Zhou2014_20220712/ICD10_symptom_tfidf_matrix.txt",quote=F,sep="\t")
      ## Calculate cosine similarities ##
      cosinesimilarities<-cosine(icdsymptom_matrix)
      ind <- which(upper.tri(cosinesimilarities, diag = TRUE), arr.ind = TRUE)
      netmode2<-cbind(colnames(cosinesimilarities)[ind[,1]],colnames(cosinesimilarities)[ind[,2]], cosinesimilarities[ind])
      netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
      write.table(netmode2,"Symptoms_Zhou2014_20220712/ICD10_CosineSimilarities.txt",quote=F,sep="\t",row.names = F,col.names = F)
    }
    if(args[3]=="random_distances"){
      ## Calculate random cosine similarities ##
      set.seed(as.numeric(args[4]))
      for(a in 1:2000){
        ini<-Sys.time()
        ricdsymptom_matrix<-apply(icdsymptom_matrix,2,sample)
        rcosinesimilarities<-cosine(ricdsymptom_matrix)
        ind <- which(upper.tri(rcosinesimilarities, diag = TRUE), arr.ind = TRUE)
        netmode2<-cbind(colnames(rcosinesimilarities)[ind[,1]],colnames(rcosinesimilarities)[ind[,2]], rcosinesimilarities[ind])
        netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
        write.table(netmode2,paste("Symptoms_Zhou2014_20220712/RandomAssociations/Seed_",args[4],"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
        fin<-Sys.time()
        print(paste("Seed ",args[4]," round ",a," time: ",fin-ini,sep=""))
      }
    }
  }
  if(args[2]=="calculatesignificance"){
    originalnetwork<-fread("Symptoms_Zhou2014_20220712/ICD10_CosineSimilarities.txt",stringsAsFactors = F,sep="\t")
    files<-list.files("Symptoms_Zhou2014_20220712/RandomAssociations/")
    randommatrix<-matrix(ncol=10000,nrow=length(originalnetwork$V1),NA)
    for(b in 1:length(files)){
      # b<-1
      randomnetwork<-fread(paste("Symptoms_Zhou2014_20220712/RandomAssociations/",files[b],sep=""),stringsAsFactors = F,sep="\t")
      randommatrix[,b]<-randomnetwork$V3
      print(paste(round((b/length(files))*100,2),"%",sep=""))
    }
    resultmatrix<-matrix(ncol=6,nrow=length(originalnetwork$V1),NA)
    for(a in 1:length(originalnetwork$V1)){
      # a<-1
      thepvalup<-length(which(randommatrix[a,]>=originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvaldown<-length(which(randommatrix[a,]<=originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvalupni<-length(which(randommatrix[a,]>originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvaldownni<-length(which(randommatrix[a,]<originalnetwork$V3[a]))/length(randommatrix[a,])
      meandistance<-mean(randommatrix[a,])
      sddistance<-sd(randommatrix[a,])
      resultmatrix[a,]<-c(thepvalup,thepvaldown,thepvalupni,thepvaldownni,meandistance,sddistance)
      print(paste(round((a/length(originalnetwork$V1))*100,2),"%",sep=""))
    }
    resultmatrix<-cbind(originalnetwork,resultmatrix)
    colnames(resultmatrix)<-c("disease1","disease2","cosine_similarity","pval_up","pval_down","pval_up_not_equal","pval_down_not_equal","mean","sd")
    write.table(resultmatrix,"Symptoms_Zhou2014_20220712/ICD10_CosineSimilarities_pvals.txt",quote=F,sep="\t",row.names=F)
  }
  ## Work only with the significant disease-symptom associations ##
  if(args[2]=="significantdiseasesymptomassociations"){
    ## Select the significant associations between diseases and symptoms ##
    ## Calculate chi-squares to select significant disease-symptom associations ##
    ## Read disease-symptom associations ##
    listofdiseasesymptoms<-fread("Symptoms_Zhou2014_20220712/41467_2014_BFncomms5212_MOESM1045_ESM.txt",stringsAsFactors = F,sep="\t")
    ## Get the number of papers associated to a symptom ##
    symptomfrequency<-fread("Symptoms_Zhou2014_20220712/41467_2014_BFncomms5212_MOESM1044_ESM.txt")
    setkey(symptomfrequency,"MeSH Symptom Term")
    ## Get the number of papers associated to a disease ## 
    diseasefrequency<-fread("Symptoms_Zhou2014_20220712/41467_2014_BFncomms5212_MOESM1043_ESM.txt")
    setkey(diseasefrequency,"MeSH Disease Term")
    ## Get the p-values and the chi-square estimates ##
    pvaltabs<-cbind(listofdiseasesymptoms[,1:2],NA,NA)
    for(a in 1:length(listofdiseasesymptoms$`MeSH Symptom Term`)){
      # a<-1
      subtab<-rbind(c(listofdiseasesymptoms$`PubMed occurrence`[a],
                      symptomfrequency[listofdiseasesymptoms$`MeSH Symptom Term`[a]]$`PubMed occurrence`-listofdiseasesymptoms$`PubMed occurrence`[a]),
                    c(diseasefrequency[listofdiseasesymptoms$`MeSH Disease Term`[a]]$`PubMed occurrence`-listofdiseasesymptoms$`PubMed occurrence`[a],
                      849103-symptomfrequency[listofdiseasesymptoms$`MeSH Symptom Term`[a]]$`PubMed occurrence`-diseasefrequency[listofdiseasesymptoms$`MeSH Disease Term`[a]]$`PubMed occurrence`+listofdiseasesymptoms$`PubMed occurrence`[a]))
      chisq <- chisq.test(subtab)
      pvaltabs$V2[a]<-chisq$p.value
      pvaltabs$V3[a]<-chisq$statistic
      print(paste(round((a/length(listofdiseasesymptoms$`MeSH Symptom Term`))*100,2),"%",sep=""))
    }
    pvaltabs<-cbind(pvaltabs,NA)
    colnames(pvaltabs)[3:5]<-c("pval","chisquare","fdr")
    pvaltabs$fdr<-pvaltabs$pval*length(pvaltabs$pval)
    if(length(which(pvaltabs$fdr>1))>0){
      pvaltabs$fdr[which(pvaltabs$fdr>1)]<-1
    }
    listofdiseasesymptoms<-cbind(listofdiseasesymptoms,pvaltabs$pval,pvaltabs$fdr)
    colnames(listofdiseasesymptoms)[5:6]<-c("pval","fdr")
    write.table(listofdiseasesymptoms,"Symptoms_Zhou2014_20220712/DiseaseSymptom_fdr.txt",quote=F,sep="\t",row.names=F)
  }
  if(args[2]=="sign_getthenewvectors"){
    ## Sum the disease - symptom co-occurrences ##
    if("ICD10_symptom_significant_associations.txt"%in%list.files("Symptoms_Zhou2014_20220712/")==FALSE){
      ## Load the disease name - icd10 association file
      symptomnametoicd<-fread("Symptoms_Zhou2014_20220712/Disease_ICD10_translation.txt",stringsAsFactors = F,sep="\t")
      symptomnametoicd<-symptomnametoicd[which(nchar(symptomnametoicd$ICD10)==3)]
      ## Load the disease-symptom association file
      diseasesymptom<-fread("Symptoms_Zhou2014_20220712/DiseaseSymptom_fdr.txt")
      diseasesymptom<-diseasesymptom[which(diseasesymptom$fdr<=0.05)]
      ## Select the diseases for which there is still information available ##
      commondis<-intersect(symptomnametoicd$Name,diseasesymptom$`MeSH Disease Term`)
      setkey(symptomnametoicd,"Name")
      symptomnametoicd<-symptomnametoicd[commondis]
      ## Select the diseases for which we have and ICD correspondence and sum their mentions ##
      icds<-unique(symptomnametoicd$ICD10)
      icdsymptom<-c()
      for(a in 1:length(icds)){
        # a<-1
        diseases<-symptomnametoicd$Name[which(symptomnametoicd$ICD10==icds[a])]
        icdtab<-c()
        for(b in 1:length(diseases)){
          # b<-1
          icdtab<-rbind(icdtab,diseasesymptom[which(diseasesymptom$`MeSH Disease Term`==diseases[b])])
        }
        icdtab<-icdtab[,lapply(.SD,sum),by=.(`MeSH Symptom Term`), .SDcols=c("PubMed occurrence")]
        icdtab<-cbind(icds[a],icdtab)
        colnames(icdtab)<-c("icd10","symptom","occurrence")
        icdsymptom<-rbind(icdsymptom,icdtab)
        print(paste(round((a/length(icds))*100,2),"%",sep=""))
      }
      ## Write the table containing the number of ICD10 - symptom occurrences ##
      write.table(icdsymptom,"Symptoms_Zhou2014_20220712/ICD10_symptom_significant_associations.txt",quote=F,sep="\t",row.names=F)
    }
    ## Define the function for calculating the TF-IDF ##
    tfidf<-function(occurrence,diseaseswithsymptom,numberdiseases){
      rtfidf<-occurrence*log(numberdiseases/diseaseswithsymptom)
      return(rtfidf)
    }
    ## Calculate the TF-IDF ##
    if(dim(icdsymptom)[2]==3){
      ## First get the number of ICD10 associated to each symptom ##
      Ndiseasesbysymptom<-table(icdsymptom$symptom)
      ## Now go disease by disease ##
      theicds<-unique(icdsymptom$icd10)
      icdswithtfidf<-c()
      for(a in 1:length(theicds)){
        # a<-1
        thedisease<-icdsymptom[which(icdsymptom$icd10==theicds[a])]
        thetfidfs<-c()
        for(b in 1:length(thedisease$icd10)){
          thevalue<-tfidf(thedisease$occurrence[b],Ndiseasesbysymptom[thedisease$symptom[b]],length(theicds))
          thetfidfs<-c(thetfidfs,as.numeric(thevalue))
        }
        thedisease<-cbind(thedisease,thetfidfs)
        icdswithtfidf<-rbind(icdswithtfidf,thedisease)
        print(paste(round((a/length(theicds))*100,2),"%",sep=""))
      }
      ## Write the table containing the number of ICD10 - symptom occurrences with the associated TF-IDF ##
      write.table(icdswithtfidf,"Symptoms_Zhou2014_20220712/ICD10_symptom_significant_associations.txt",quote=F,sep="\t",row.names=F)
    }
  }
  if(args[2]=="sign_calculatesimilarities"){
    if("RandomAssociations"%in%list.files("Symptoms_Zhou2014_20220712/")==F){dir.create("Symptoms_Zhou2014_20220712/RandomAssociations")}
    associations<-fread("Symptoms_Zhou2014_20220712/ICD10_symptom_significant_associations.txt",stringsAsFactors = F,sep="\t")
    ## First create a matrix with the TF-IDF of each symptom in each disease ##
    diseases<-unique(associations$icd10) ; symptoms<-unique(associations$symptom)
    icdsymptom_matrix<-matrix(ncol=length(diseases),nrow=length(symptoms),0) ; colnames(icdsymptom_matrix)<-diseases ; rownames(icdsymptom_matrix)<-symptoms
    # icdsymptom_number<-icdsymptom_matrix
    for(a in 1:length(symptoms)){
      icdsymptom_matrix[symptoms[a],associations$icd10[which(associations$symptom==symptoms[a])]]<-associations$thetfidfs[which(associations$symptom==symptoms[a])]
      # icdsymptom_number[symptoms[a],associations$icd10[which(associations$symptom==symptoms[a])]]<-associations$occurrence[which(associations$symptom==symptoms[a])]
    }
    if(args[3]=="CreateTFIDFmatrix"){
      write.table(icdsymptom_matrix,"Symptoms_Zhou2014_20220712/ICD10_symptom_significant_tfidf_matrix.txt",quote=F,sep="\t")
      ## Calculate cosine similarities ##
      cosinesimilarities<-cosine(icdsymptom_matrix)
      ind <- which(upper.tri(cosinesimilarities, diag = TRUE), arr.ind = TRUE)
      netmode2<-cbind(colnames(cosinesimilarities)[ind[,1]],colnames(cosinesimilarities)[ind[,2]], cosinesimilarities[ind])
      netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
      write.table(netmode2,"Symptoms_Zhou2014_20220712/ICD10_significant_CosineSimilarities.txt",quote=F,sep="\t",row.names = F,col.names = F)
    }
    if(args[3]=="random_distances"){
      ## Calculate random cosine similarities ##
      set.seed(as.numeric(args[4]))
      for(a in 1:2000){
        ini<-Sys.time()
        ricdsymptom_matrix<-apply(icdsymptom_matrix,2,sample)
        rcosinesimilarities<-cosine(ricdsymptom_matrix)
        ind <- which(upper.tri(rcosinesimilarities, diag = TRUE), arr.ind = TRUE)
        netmode2<-cbind(colnames(rcosinesimilarities)[ind[,1]],colnames(rcosinesimilarities)[ind[,2]], rcosinesimilarities[ind])
        netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
        write.table(netmode2,paste("Symptoms_Zhou2014_20220712/RandomAssociations/Seed_",args[4],"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
        fin<-Sys.time()
        print(paste("Seed ",args[4]," round ",a," time: ",fin-ini,sep=""))
      }
    }
  }
  if(args[2]=="sign_calculatesignificance"){
    originalnetwork<-fread("Symptoms_Zhou2014_20220712/ICD10_significant_CosineSimilarities.txt",stringsAsFactors = F,sep="\t")
    files<-list.files("Symptoms_Zhou2014_20220712/RandomAssociations/")
    randommatrix<-matrix(ncol=10000,nrow=length(originalnetwork$V1),NA)
    for(b in 1:length(files)){
      # b<-1
      randomnetwork<-fread(paste("Symptoms_Zhou2014_20220712/RandomAssociations/",files[b],sep=""),stringsAsFactors = F,sep="\t")
      randommatrix[,b]<-randomnetwork$V3
      print(paste(round((b/length(files))*100,2),"%",sep=""))
    }
    resultmatrix<-matrix(ncol=6,nrow=length(originalnetwork$V1),NA)
    for(a in 1:length(originalnetwork$V1)){
      # a<-1
      thepvalup<-length(which(randommatrix[a,]>=originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvaldown<-length(which(randommatrix[a,]<=originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvalupni<-length(which(randommatrix[a,]>originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvaldownni<-length(which(randommatrix[a,]<originalnetwork$V3[a]))/length(randommatrix[a,])
      meandistance<-mean(randommatrix[a,])
      sddistance<-sd(randommatrix[a,])
      resultmatrix[a,]<-c(thepvalup,thepvaldown,thepvalupni,thepvaldownni,meandistance,sddistance)
      print(paste(round((a/length(originalnetwork$V1))*100,2),"%",sep=""))
    }
    resultmatrix<-cbind(originalnetwork,resultmatrix)
    colnames(resultmatrix)<-c("disease1","disease2","cosine_similarity","pval_up","pval_down","pval_up_not_equal","pval_down_not_equal","mean","sd")
    write.table(resultmatrix,"Symptoms_Zhou2014_20220712/ICD10_significant_CosineSimilarities_pvals.txt",quote=F,sep="\t",row.names=F)
  }
  if(args[2]=="compare_sign_nonsign_assoc"){
    allass<-fread("Symptoms_Zhou2014_20220712/ICD10_symptom_associations.txt",stringsAsFactors = F)
    sigass<-fread("Symptoms_Zhou2014_20220712/ICD10_symptom_significant_associations.txt",stringsAsFactors = F)
    
    allass[which(allass$symptom=="Pain")]
    sigass[which(sigass$symptom=="Pain")]
  }
}

#### Genes ####
## Cosine Similarity done ##
## Jaccard index pending ##
if(args[1]=="genes"){
  if(args[2]=="disgenet"){
    if(args[3]=="getthenewvectors"){
      ## Sum the disease - symptom co-occurrences ##
      if("ICD10_genes_associations_DisGeNET.txt"%in%list.files("Symptoms_Zhou2014_20220712/")==FALSE){
        ## Load the disease name - icd10 association file
        genenametoicd<-fread("Genes_DisGeNET_20230113/Disease_ICD10_translation.txt",stringsAsFactors = F,sep="\t")
        genenametoicd<-genenametoicd[which(nchar(genenametoicd$threedigits)==3)]
        ## Load the disease-symptom association file
        diseasegene<-fread("Genes_DisGeNET_20230113/curated_gene_disease_associations.tsv",stringsAsFactors = F,sep="\t")
        diseasegene<-diseasegene[,c(2:5,10,14:16)]
        ## Histogram of the scores ##
        hist(diseasegene$score,breaks = 100,main="score")
        ## The Disease Specificity Index (DSI) ranges from 0 to 1, and is inversely proportional to the number of diseases 
        ## associated to a particular gene. A gene associated to a large number of diseases (e.g. TNF, associated to more 
        ## than 1500 diseases) will have a DSI close to zero, while a gene associated to only one disease, is more ‘specific’
        ## for that disease and has DSI of 1
        hist(diseasegene$DSI,breaks = 100,main="DSI")
        ## The Disease Pleiotropy Index (DPI) ranges from 0 to 1 and is proportional to the number of different (MeSH) disease 
        ## classes a gene is associated to. Thus, a gene associated to diseases of diverse classes (such as APOE, associated to 
        ## Cardiovascular Diseases, Mental Disorders, Neoplasms, Respiratory Tract Diseases, etc), will have a DPI close to 1. 
        ## Conversely, the PSCA, associated to 58 diseases, most of which are neoplasms has a relatively low DPI.
        # hist(diseasegene$DPI,breaks = 20,main="DPI")
        head(sort(table(diseasegene$source),decreasing = T))
        ## Select the diseases for which we have and ICD correspondence and do the union of their genes ##
        icds<-unique(genenametoicd$threedigits)
        
        
        ## Build an ICD10 - gene association table ##
        icdgene<-c()
        for(a in 1:length(icds)){
          # a<-1
          diseases<-genenametoicd$CUI[which(genenametoicd$threedigits==icds[a])]
          icdtab<-c()
          for(b in 1:length(diseases)){
            # b<-1
            icdtab<-rbind(icdtab,diseasegene[which(diseasegene$diseaseId==diseases[b])][,c(1,2,5)])
          }
          if(length(which(duplicated(icdtab$geneSymbol)))>0){icdtab<-icdtab[-which(duplicated(icdtab$geneSymbol))]}
          icdtab<-cbind(icds[a],icdtab)
          colnames(icdtab)<-c("icd10","gene","dsi","score")
          icdgene<-rbind(icdgene,icdtab)
          print(paste(round((a/length(icds))*100,2),"%",sep=""))
        }
        ## Write the table containing the number of ICD10 - symptom occurrences ##
        write.table(icdgene,"Genes_DisGeNET_20230113/ICD10_gene_associations.txt",quote=F,sep="\t",row.names=F)
      }
    }
    if(args[3]=="calculatesimilarities"){
      if("RandomAssociations"%in%list.files("Genes_DisGeNET_20230113/")==F){dir.create("Genes_DisGeNET_20230113//RandomAssociations")}
      associations<-fread("Genes_DisGeNET_20230113/ICD10_gene_associations.txt",stringsAsFactors = F,sep="\t")
      ## First create a matrix with the DSI of each gene in each disease ##
      diseases<-unique(associations$icd10) ; genes<-unique(associations$gene)
      icdgene_matrix<-matrix(ncol=length(diseases),nrow=length(genes),0) ; colnames(icdgene_matrix)<-diseases ; rownames(icdgene_matrix)<-genes
      for(a in 1:length(genes)){
        icdgene_matrix[genes[a],associations$icd10[which(associations$gene==genes[a])]]<-associations$dsi[which(associations$gene==genes[a])]
      }
      if(args[4]=="CreateDSImatrix"){
        write.table(icdgene_matrix,"Genes_DisGeNET_20230113/ICD10_gene_dsi_matrix.txt",quote=F,sep="\t")
        ## Calculate cosine similarities ##
        cosinesimilarities<-cosine(icdgene_matrix)
        ind <- which(upper.tri(cosinesimilarities, diag = TRUE), arr.ind = TRUE)
        netmode2<-cbind(colnames(cosinesimilarities)[ind[,1]],colnames(cosinesimilarities)[ind[,2]], cosinesimilarities[ind])
        netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
        write.table(netmode2,"Genes_DisGeNET_20230113/ICD10_CosineSimilarities.txt",quote=F,sep="\t",row.names = F,col.names = F)
      }
      if(args[4]=="random_distances"){
        ## Calculate random cosine similarities ##
        set.seed(as.numeric(args[5]))
        for(a in 1:1000){
          ini<-Sys.time()
          ricdgene_matrix<-apply(icdgene_matrix,2,sample)
          rcosinesimilarities<-cosine(ricdgene_matrix)
          ind <- which(upper.tri(rcosinesimilarities, diag = TRUE), arr.ind = TRUE)
          netmode2<-cbind(colnames(rcosinesimilarities)[ind[,1]],colnames(rcosinesimilarities)[ind[,2]], rcosinesimilarities[ind])
          netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
          write.table(netmode2,paste("Genes_DisGeNET_20230113/RandomAssociations/Seed_",args[5],"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
          fin<-Sys.time()
          print(paste("Seed ",args[5]," round ",a," time: ",fin-ini,sep=""))
        }
      }
    }
    if(args[3]=="calculatesignificance"){
      originalnetwork<-fread("Genes_DisGeNET_20230113/ICD10_CosineSimilarities.txt",stringsAsFactors = F,sep="\t")
      files<-list.files("Genes_DisGeNET_20230113/RandomAssociations/")
      randommatrix<-matrix(ncol=10000,nrow=length(originalnetwork$V1),NA)
      for(b in 1:length(files)){
        # b<-1
        randomnetwork<-fread(paste("Genes_DisGeNET_20230113/RandomAssociations/",files[b],sep=""),stringsAsFactors = F,sep="\t")
        randommatrix[,b]<-randomnetwork$V3
        print(paste(round((b/length(files))*100,2),"%",sep=""))
      }
      resultmatrix<-matrix(ncol=6,nrow=length(originalnetwork$V1),NA)
      for(a in 1:length(originalnetwork$V1)){
        # a<-1
        thepvalup<-length(which(randommatrix[a,]>=originalnetwork$V3[a]))/length(randommatrix[a,])
        thepvaldown<-length(which(randommatrix[a,]<=originalnetwork$V3[a]))/length(randommatrix[a,])
        thepvalupni<-length(which(randommatrix[a,]>originalnetwork$V3[a]))/length(randommatrix[a,])
        thepvaldownni<-length(which(randommatrix[a,]<originalnetwork$V3[a]))/length(randommatrix[a,])
        meandistance<-mean(randommatrix[a,])
        sddistance<-sd(randommatrix[a,])
        resultmatrix[a,]<-c(thepvalup,thepvaldown,thepvalupni,thepvaldownni,meandistance,sddistance)
        print(paste(round((a/length(originalnetwork$V1))*100,2),"%",sep=""))
      }
      resultmatrix<-cbind(originalnetwork,resultmatrix)
      colnames(resultmatrix)<-c("disease1","disease2","cosine_similarity","pval_up","pval_down","pval_up_not_equal","pval_down_not_equal","mean","sd")
      write.table(resultmatrix,"Genes_DisGeNET_20230113/ICD10_CosineSimilarities_pvals.txt",quote=F,sep="\t",row.names=F)
    }
  }
  if(args[2]=="disgenet_sinmirna"){
    if(args[3]=="calculatesimilarities"){
      if("RandomAssociations"%in%list.files("Genes_DisGeNET_20230113/")==F){dir.create("Genes_DisGeNET_20230113//RandomAssociations")}
      associations<-fread("Genes_DisGeNET_20230113/ICD10_gene_associations.txt",stringsAsFactors = F,sep="\t")
      if(length(grep("mir",associations$gene,ignore.case = T))>0){associations<-associations[-grep("mir",associations$gene,ignore.case = T)]}
      ## First create a matrix with the DSI of each gene in each disease ##
      diseases<-unique(associations$icd10) ; genes<-unique(associations$gene)
      icdgene_matrix<-matrix(ncol=length(diseases),nrow=length(genes),0) ; colnames(icdgene_matrix)<-diseases ; rownames(icdgene_matrix)<-genes
      for(a in 1:length(genes)){
        icdgene_matrix[genes[a],associations$icd10[which(associations$gene==genes[a])]]<-associations$dsi[which(associations$gene==genes[a])]
      }
      if(args[4]=="CreateDSImatrix"){
        write.table(icdgene_matrix,"Genes_DisGeNET_20230113/ICD10_gene_dsi_matrix_without_mirna.txt",quote=F,sep="\t")
        ## Calculate cosine similarities ##
        cosinesimilarities<-cosine(icdgene_matrix)
        ind <- which(upper.tri(cosinesimilarities, diag = TRUE), arr.ind = TRUE)
        netmode2<-cbind(colnames(cosinesimilarities)[ind[,1]],colnames(cosinesimilarities)[ind[,2]], cosinesimilarities[ind])
        netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
        write.table(netmode2,"Genes_DisGeNET_20230113/ICD10_CosineSimilarities_without_mirna.txt",quote=F,sep="\t",row.names = F,col.names = F)
      }
      if(args[4]=="random_distances"){
        ## Calculate random cosine similarities ##
        set.seed(as.numeric(args[5]))
        for(a in 1:1000){
          ini<-Sys.time()
          ricdgene_matrix<-apply(icdgene_matrix,2,sample)
          rcosinesimilarities<-cosine(ricdgene_matrix)
          ind <- which(upper.tri(rcosinesimilarities, diag = TRUE), arr.ind = TRUE)
          netmode2<-cbind(colnames(rcosinesimilarities)[ind[,1]],colnames(rcosinesimilarities)[ind[,2]], rcosinesimilarities[ind])
          netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
          write.table(netmode2,paste("Genes_DisGeNET_20230113/RandomAssociations/Seed_",args[5],"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
          fin<-Sys.time()
          print(paste("Seed ",args[5]," round ",a," time: ",fin-ini,sep=""))
        }
      }
    }
    if(args[3]=="calculatesignificance"){
      originalnetwork<-fread("Genes_DisGeNET_20230113/ICD10_CosineSimilarities_without_mirna.txt",stringsAsFactors = F,sep="\t")
      files<-list.files("Genes_DisGeNET_20230113/RandomAssociations/")
      randommatrix<-matrix(ncol=10000,nrow=length(originalnetwork$V1),NA)
      for(b in 1:length(files)){
        # b<-1
        randomnetwork<-fread(paste("Genes_DisGeNET_20230113/RandomAssociations/",files[b],sep=""),stringsAsFactors = F,sep="\t")
        randommatrix[,b]<-randomnetwork$V3
        print(paste(round((b/length(files))*100,2),"%",sep=""))
      }
      resultmatrix<-matrix(ncol=6,nrow=length(originalnetwork$V1),NA)
      for(a in 1:length(originalnetwork$V1)){
        # a<-1
        thepvalup<-length(which(randommatrix[a,]>=originalnetwork$V3[a]))/length(randommatrix[a,])
        thepvaldown<-length(which(randommatrix[a,]<=originalnetwork$V3[a]))/length(randommatrix[a,])
        thepvalupni<-length(which(randommatrix[a,]>originalnetwork$V3[a]))/length(randommatrix[a,])
        thepvaldownni<-length(which(randommatrix[a,]<originalnetwork$V3[a]))/length(randommatrix[a,])
        meandistance<-mean(randommatrix[a,])
        sddistance<-sd(randommatrix[a,])
        resultmatrix[a,]<-c(thepvalup,thepvaldown,thepvalupni,thepvaldownni,meandistance,sddistance)
        print(paste(round((a/length(originalnetwork$V1))*100,2),"%",sep=""))
      }
      resultmatrix<-cbind(originalnetwork,resultmatrix)
      colnames(resultmatrix)<-c("disease1","disease2","cosine_similarity","pval_up","pval_down","pval_up_not_equal","pval_down_not_equal","mean","sd")
      write.table(resultmatrix,"Genes_DisGeNET_20230113/ICD10_CosineSimilarities_without_mirna_pvals.txt",quote=F,sep="\t",row.names=F)
    }
  }
}

#### Drugs ####
## Cosine Similarity done ##
## Jaccard index pending ##
if(args[1]=="drugs"){
  if(args[2]=="getthenewvectors"){
    ## Sum the disease - symptom co-occurrences ##
    if("ICD10_drug_associations.txt"%in%list.files("Drugs_CTD_20220712/")==FALSE){
      ## Load the disease name - icd10 association file
      drugnametoicd<-fread("Drugs_CTD_20220712/Disease_ICD10_translation.txt",stringsAsFactors = F,sep="\t")
      drugnametoicd<-drugnametoicd[which(nchar(drugnametoicd$threedigits)==3)]
      ## Load the disease-drug association file
      diseasedrug<-fread("Drugs_CTD_20220712/CTD_chemicals_diseases.csv",stringsAsFactors = F,sep=",")
      diseasedrug<-diseasedrug[unique(c(which(diseasedrug$DirectEvidence=="marker/mechanism"),which(diseasedrug$DirectEvidence=="therapeutic"))),c(1,2,4:6)]
      write.table(diseasedrug,"Drugs_CTD_20220712/Disease_drug_direct_evidence.txt",quote=F,sep="\t",row.names = F)
      ## Select only the drugs with a therapeutic association with the disease ##
      diseasedrug<-diseasedrug[which(diseasedrug$DirectEvidence=="therapeutic")]
      ## Calculate a score for each drug based on the number of diseases it's associated to ##
      frequencyscore<-1/table(diseasedrug$ChemicalName)
      diseasedrug<-cbind(diseasedrug,as.character(frequencyscore[diseasedrug$ChemicalName]))
      colnames(diseasedrug)[6]<-"IFS" ## IFS states for Inverse Frequency Score ##
      write.table(diseasedrug,"Drugs_CTD_20220712/Disease_drug_therapeutic_use.txt",quote=F,sep="\t",row.names = F)
      
      ## Histogram of the scores ##
      hist(as.numeric(diseasedrug$IFS),breaks = 20,main="Inverse Frequency Score\n(drugs)",ylim=c(0,20000),xlab="IFS")
      ## Select the diseases for which we have and ICD correspondence and do the union of their drugs ##
      icds<-unique(drugnametoicd$threedigits)
      
      ## Build an ICD10 - drug association table ##
      icddrug<-c()
      for(a in 1:length(icds)){
        # a<-1
        diseases<-drugnametoicd$disease[which(drugnametoicd$threedigits==icds[a])]
        icdtab<-c()
        for(b in 1:length(diseases)){
          # b<-1
          icdtab<-rbind(icdtab,diseasedrug[which(diseasedrug$DiseaseName==diseases[b])][,c(1,6)])
        }
        if(length(which(duplicated(icdtab$drugSymbol)))>0){icdtab<-icdtab[-which(duplicated(icdtab$drugSymbol))]}
        icdtab<-cbind(icds[a],icdtab)
        colnames(icdtab)<-c("icd10","drug","ifs")
        icddrug<-rbind(icddrug,icdtab)
        print(paste(round((a/length(icds))*100,2),"%",sep=""))
      }
      ## Write the table containing the number of ICD10 - symptom occurrences ##
      write.table(icddrug,"Drugs_CTD_20220712/ICD10_drug_associations.txt",quote=F,sep="\t",row.names=F)
    }
  }
  if(args[2]=="calculatesimilarities"){
    if("RandomAssociations"%in%list.files("Drugs_CTD_20220712/")==F){dir.create("Drugs_CTD_20220712//RandomAssociations")}
    associations<-fread("Drugs_CTD_20220712/ICD10_drug_associations.txt",stringsAsFactors = F,sep="\t")
    ## First create a matrix with the DSI of each drug in each disease ##
    diseases<-unique(associations$icd10) ; drugs<-unique(associations$drug)
    icddrug_matrix<-matrix(ncol=length(diseases),nrow=length(drugs),0) ; colnames(icddrug_matrix)<-diseases ; rownames(icddrug_matrix)<-drugs
    for(a in 1:length(drugs)){
      icddrug_matrix[drugs[a],associations$icd10[which(associations$drug==drugs[a])]]<-associations$ifs[which(associations$drug==drugs[a])]
    }
    if(args[3]=="CreateIFSmatrix"){
      write.table(icddrug_matrix,"Drugs_CTD_20220712/ICD10_drug_ifs_matrix.txt",quote=F,sep="\t")
      ## Calculate cosine similarities ##
      cosinesimilarities<-cosine(icddrug_matrix)
      ind <- which(upper.tri(cosinesimilarities, diag = TRUE), arr.ind = TRUE)
      netmode2<-cbind(colnames(cosinesimilarities)[ind[,1]],colnames(cosinesimilarities)[ind[,2]], cosinesimilarities[ind])
      netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
      write.table(netmode2,"Drugs_CTD_20220712/ICD10_CosineSimilarities.txt",quote=F,sep="\t",row.names = F,col.names = F)
    }
    if(args[3]=="random_distances"){
      ## Calculate random cosine similarities ##
      set.seed(as.numeric(args[4]))
      for(a in 1:2000){
        ini<-Sys.time()
        ricddrug_matrix<-apply(icddrug_matrix,2,sample)
        rcosinesimilarities<-cosine(ricddrug_matrix)
        ind <- which(upper.tri(rcosinesimilarities, diag = TRUE), arr.ind = TRUE)
        netmode2<-cbind(colnames(rcosinesimilarities)[ind[,1]],colnames(rcosinesimilarities)[ind[,2]], rcosinesimilarities[ind])
        netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
        write.table(netmode2,paste("Drugs_CTD_20220712/RandomAssociations/Seed_",args[4],"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
        fin<-Sys.time()
        print(paste("Seed ",args[4]," round ",a," time: ",fin-ini,sep=""))
      }
    }
  }
  if(args[2]=="calculatesignificance"){
    originalnetwork<-fread("Drugs_CTD_20220712/ICD10_CosineSimilarities.txt",stringsAsFactors = F,sep="\t")
    files<-list.files("Drugs_CTD_20220712/RandomAssociations/")
    randommatrix<-matrix(ncol=10000,nrow=length(originalnetwork$V1),NA)
    for(b in 1:length(files)){
      # b<-1
      randomnetwork<-fread(paste("Drugs_CTD_20220712/RandomAssociations/",files[b],sep=""),stringsAsFactors = F,sep="\t")
      randommatrix[,b]<-randomnetwork$V3
      print(paste(round((b/length(files))*100,2),"%",sep=""))
    }
    resultmatrix<-matrix(ncol=6,nrow=length(originalnetwork$V1),NA)
    for(a in 1:length(originalnetwork$V1)){
      # a<-1
      thepvalup<-length(which(randommatrix[a,]>=originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvaldown<-length(which(randommatrix[a,]<=originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvalupni<-length(which(randommatrix[a,]>originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvaldownni<-length(which(randommatrix[a,]<originalnetwork$V3[a]))/length(randommatrix[a,])
      meandistance<-mean(randommatrix[a,])
      sddistance<-sd(randommatrix[a,])
      resultmatrix[a,]<-c(thepvalup,thepvaldown,thepvalupni,thepvaldownni,meandistance,sddistance)
      print(paste(round((a/length(originalnetwork$V1))*100,2),"%",sep=""))
    }
    resultmatrix<-cbind(originalnetwork,resultmatrix)
    colnames(resultmatrix)<-c("disease1","disease2","cosine_similarity","pval_up","pval_down","pval_up_not_equal","pval_down_not_equal","mean","sd")
    write.table(resultmatrix,"Drugs_CTD_20220712/ICD10_CosineSimilarities_pvals.txt",quote=F,sep="\t",row.names=F)
  }
}

#### miRNAs ####
## Cosine Similarity done ##
## Jaccard index pending ##
if(args[1]=="mirnas"){
  ## Define the function for calculating the TF-IDF ##
  tfidf<-function(occurrence,diseaseswithmirna,numberdiseases){
    rtfidf<-occurrence*log(numberdiseases/diseaseswithmirna)
    return(rtfidf)
  }
  if(args[2]=="getthenewvectors"){
    ## Sum the disease - symptom co-occurrences ##
    if("ICD10_mirna_associations.txt"%in%list.files("miRNA_HMDD_20240523/")==FALSE){
      ## Load the disease name - icd10 association file
      mirnanametoicd<-fread("miRNA_HMDD_20240523/Disease_ICD10_translation.txt",stringsAsFactors = F,sep="\t")
      mirnanametoicd<-mirnanametoicd[which(nchar(mirnanametoicd$threedigits)==3)]
      ## Load the disease-mirna association file
      diseasemirna<-fread("miRNA_HMDD_20240523/alldata_v4.txt",stringsAsFactors = F,sep="\t")
      diseasemirnainfo<-diseasemirna$description
      diseasemirna<-diseasemirna[,c(4,3,1,2)]
      ## For the future calculation of the TF-IDF remove if a paper associates a disease to a miRNA twice ##
      if(length(which(duplicated(diseasemirna[,c(1,2,4)])))>0){
        diseasemirna<-diseasemirna[-which(duplicated(diseasemirna[,c(1,2,4)]))]
      }
      ## Each paper referencing the association between a disease and a mirna, in this case focused on Alzheimer's disease ##
      # diseasemirna[intersect(which(diseasemirna$disease=="Alzheimer Disease"),which(diseasemirna$miRNA=="hsa-mir-146a")),]
      
      ## Build an ICD10 - mirna association table ##
      icds<-unique(mirnanametoicd$threedigits)
      icdmirna<-c()
      for(a in 1:length(icds)){
        # a<-1
        diseases<-mirnanametoicd$disease[which(mirnanametoicd$threedigits==icds[a])]
        icdtab<-c()
        for(b in 1:length(diseases)){
          # b<-1
          icdtab<-rbind(icdtab,diseasemirna[which(tolower(diseasemirna$disease)==diseases[b])])
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
      write.table(icdtab,"miRNA_HMDD_20240523/ICD10_mirna_associations.txt",quote=F,sep="\t",row.names=F)
    }
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
    ## Write the table containing the number of ICD10 - mirna occurrences with the associated TF-IDF ##
    write.table(icdswithtfidf,"miRNA_HMDD_20240523/ICD10_mirna_associations.txt",quote=F,sep="\t",row.names=F)
  }
  if(args[2]=="calculatesimilarities"){
    if("RandomAssociations"%in%list.files("miRNA_HMDD_20240523/")==F){dir.create("miRNA_HMDD_20240523//RandomAssociations")}
    associations<-fread("miRNA_HMDD_20240523/ICD10_mirna_associations.txt",stringsAsFactors = F,sep="\t")
    ## First create a matrix with the DSI of each mirna in each disease ##
    diseases<-unique(associations$icd10) ; mirnas<-unique(associations$mirna)
    icdmirna_matrix<-matrix(ncol=length(diseases),nrow=length(mirnas),0) ; colnames(icdmirna_matrix)<-diseases ; rownames(icdmirna_matrix)<-mirnas
    for(a in 1:length(mirnas)){
      icdmirna_matrix[mirnas[a],associations$icd10[which(associations$mirna==mirnas[a])]]<-associations$thetfidfs[which(associations$mirna==mirnas[a])]
    }
    if(args[3]=="CreateTFIDFmatrix"){
      write.table(icdmirna_matrix,"miRNA_HMDD_20240523/ICD10_mirna_tfidfs_matrix.txt",quote=F,sep="\t")
      ## Calculate cosine similarities ##
      cosinesimilarities<-cosine(icdmirna_matrix)
      ind <- which(upper.tri(cosinesimilarities, diag = TRUE), arr.ind = TRUE)
      netmode2<-cbind(colnames(cosinesimilarities)[ind[,1]],colnames(cosinesimilarities)[ind[,2]], cosinesimilarities[ind])
      netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
      write.table(netmode2,"miRNA_HMDD_20240523/ICD10_CosineSimilarities.txt",quote=F,sep="\t",row.names = F,col.names = F)
    }
    if(args[3]=="random_distances"){
      ## Calculate random cosine similarities ##
      set.seed(as.numeric(args[4]))
      for(a in 1:2000){
        ini<-Sys.time()
        ricdmirna_matrix<-apply(icdmirna_matrix,2,sample)
        rcosinesimilarities<-cosine(ricdmirna_matrix)
        ind <- which(upper.tri(rcosinesimilarities, diag = TRUE), arr.ind = TRUE)
        netmode2<-cbind(colnames(rcosinesimilarities)[ind[,1]],colnames(rcosinesimilarities)[ind[,2]], rcosinesimilarities[ind])
        netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
        write.table(netmode2,paste("miRNA_HMDD_20240523/RandomAssociations/Seed_",args[4],"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
        fin<-Sys.time()
        print(paste("Seed ",args[4]," round ",a," time: ",fin-ini,sep=""))
      }
    }
  }
  if(args[2]=="calculatesignificance"){
    originalnetwork<-fread("miRNA_HMDD_20240523/ICD10_CosineSimilarities.txt",stringsAsFactors = F,sep="\t")
    files<-list.files("miRNA_HMDD_20240523/RandomAssociations/")
    randommatrix<-matrix(ncol=10000,nrow=length(originalnetwork$V1),NA)
    for(b in 1:length(files)){
      # b<-1
      randomnetwork<-fread(paste("miRNA_HMDD_20240523/RandomAssociations/",files[b],sep=""),stringsAsFactors = F,sep="\t")
      randommatrix[,b]<-randomnetwork$V3
      print(paste(round((b/length(files))*100,2),"%",sep=""))
    }
    resultmatrix<-matrix(ncol=6,nrow=length(originalnetwork$V1),NA)
    for(a in 1:length(originalnetwork$V1)){
      # a<-1
      thepvalup<-length(which(randommatrix[a,]>=originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvaldown<-length(which(randommatrix[a,]<=originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvalupni<-length(which(randommatrix[a,]>originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvaldownni<-length(which(randommatrix[a,]<originalnetwork$V3[a]))/length(randommatrix[a,])
      meandistance<-mean(randommatrix[a,])
      sddistance<-sd(randommatrix[a,])
      resultmatrix[a,]<-c(thepvalup,thepvaldown,thepvalupni,thepvaldownni,meandistance,sddistance)
      print(paste(round((a/length(originalnetwork$V1))*100,2),"%",sep=""))
    }
    resultmatrix<-cbind(originalnetwork,resultmatrix)
    colnames(resultmatrix)<-c("disease1","disease2","cosine_similarity","pval_up","pval_down","pval_up_not_equal","pval_down_not_equal","mean","sd")
    write.table(resultmatrix,"miRNA_HMDD_20240523/ICD10_CosineSimilarities_pvals.txt",quote=F,sep="\t",row.names=F)
  }
}


#### Metabolome ####
## Cosine Similarity done ##
## Jaccard index pending ##
if(args[1]=="metabolome"){
  if(args[2]=="getthenewvectors"){
    ## Sum the disease - symptom co-occurrences ##
    if("ICD10_metabolite_associations.txt"%in%list.files("Metabolome_HMDB_20220712/")==FALSE){
      ## Load the disease name - icd10 association file
      metabolitenametoicd<-fread("Metabolome_HMDB_20220712/Disease_ICD10_translation.txt",stringsAsFactors = F,sep="\t")
      metabolitenametoicd<-metabolitenametoicd[which(nchar(metabolitenametoicd$threedigits)==3)]
      ## Load the disease-symptom association file
      diseasemetabolite<-fread("Metabolome_HMDB_20220712/disease_metabolite_table.txt",stringsAsFactors = F,sep="\t")
      diseasemetabolite$Disease<-tolower(diseasemetabolite$Disease)
      ## For simplicity we focus only on blood levels ##
      diseasemetabolite<-diseasemetabolite[which(diseasemetabolite$Where=="Blood"),]
      ## Calculate for each metabolite an inverse score
      metabolites<-table(diseasemetabolite$HMDB_ID)
      metabolites<-1/metabolites
      setkey(diseasemetabolite,"HMDB_ID")
      diseasemetabolite<-cbind(diseasemetabolite,as.numeric(metabolites[diseasemetabolite$HMDB_ID]))
      colnames(diseasemetabolite)[9]<-"IFS"
      ## Histogram of the scores ##
      hist(diseasemetabolite$IFS,breaks = 100,main="score")
      ## Select the diseases for which we have and ICD correspondence and do the union of their metabolites ##
      icds<-unique(metabolitenametoicd$threedigits)
      ## Build an ICD10 - metabolite association table ##
      icdmetabolite<-c()
      for(a in 1:length(icds)){
        # a<-1
        diseases<-metabolitenametoicd$Name[which(metabolitenametoicd$threedigits==icds[a])]
        diseases<-intersect(diseasemetabolite$Disease,diseases)
        if(length(diseases)>0){
          icdtab<-c()
          for(b in 1:length(diseases)){
            # b<-1
            icdtab<-rbind(icdtab,diseasemetabolite[which(diseasemetabolite$Disease==diseases[b])][,c(1,3,9)])
          }
          if(length(which(duplicated(icdtab$HMDB_ID)))>0){icdtab<-icdtab[-which(duplicated(icdtab$HMDB_ID))]}
          icdtab<-cbind(icds[a],icdtab)
          colnames(icdtab)<-c("icd10","disease","metabolite","ifs")
          icdmetabolite<-rbind(icdmetabolite,icdtab)
          print(paste(round((a/length(icds))*100,2),"%",sep=""))
        }
      }
      icdmetabolite<-icdmetabolite[,-2]
      ## Write the table containing the number of ICD10 - symptom occurrences ##
      write.table(icdmetabolite,"Metabolome_HMDB_20220712/ICD10_metabolite_associations.txt",quote=F,sep="\t",row.names=F)
    }
  }
  if(args[2]=="calculatesimilarities"){
    if("RandomAssociations"%in%list.files("Metabolome_HMDB_20220712/")==F){dir.create("Metabolome_HMDB_20220712//RandomAssociations")}
    associations<-fread("Metabolome_HMDB_20220712/ICD10_metabolite_associations.txt",stringsAsFactors = F,sep="\t")
    ## First create a matrix with the ifs of each metabolite in each disease ##
    diseases<-unique(associations$icd10) ; metabolites<-unique(associations$metabolite)
    icdmetabolite_matrix<-matrix(ncol=length(diseases),nrow=length(metabolites),0) ; colnames(icdmetabolite_matrix)<-diseases ; rownames(icdmetabolite_matrix)<-metabolites
    for(a in 1:length(metabolites)){
      icdmetabolite_matrix[metabolites[a],associations$icd10[which(associations$metabolite==metabolites[a])]]<-associations$ifs[which(associations$metabolite==metabolites[a])]
    }
    if(args[3]=="CreateIFSmatrix"){
      write.table(icdmetabolite_matrix,"Metabolome_HMDB_20220712/ICD10_metabolite_ifs_matrix.txt",quote=F,sep="\t")
      ## Calculate cosine similarities ##
      cosinesimilarities<-cosine(icdmetabolite_matrix)
      ind <- which(upper.tri(cosinesimilarities, diag = TRUE), arr.ind = TRUE)
      netmode2<-cbind(colnames(cosinesimilarities)[ind[,1]],colnames(cosinesimilarities)[ind[,2]], cosinesimilarities[ind])
      netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
      write.table(netmode2,"Metabolome_HMDB_20220712/ICD10_CosineSimilarities.txt",quote=F,sep="\t",row.names = F,col.names = F)
    }
    if(args[3]=="random_distances"){
      ## Calculate random cosine similarities ##
      set.seed(as.numeric(args[4]))
      for(a in 1:2000){
        ini<-Sys.time()
        ricdmetabolite_matrix<-apply(icdmetabolite_matrix,2,sample)
        rcosinesimilarities<-cosine(ricdmetabolite_matrix)
        ind <- which(upper.tri(rcosinesimilarities, diag = TRUE), arr.ind = TRUE)
        netmode2<-cbind(colnames(rcosinesimilarities)[ind[,1]],colnames(rcosinesimilarities)[ind[,2]], rcosinesimilarities[ind])
        netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
        write.table(netmode2,paste("Metabolome_HMDB_20220712/RandomAssociations/Seed_",args[4],"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
        fin<-Sys.time()
        print(paste("Seed ",args[4]," round ",a," time: ",fin-ini,sep=""))
      }
    }
  }
  if(args[2]=="calculatesignificance"){
    originalnetwork<-fread("Metabolome_HMDB_20220712/ICD10_CosineSimilarities.txt",stringsAsFactors = F,sep="\t")
    files<-list.files("Metabolome_HMDB_20220712/RandomAssociations/")
    randommatrix<-matrix(ncol=10000,nrow=length(originalnetwork$V1),NA)
    for(b in 1:length(files)){
      # b<-1
      randomnetwork<-fread(paste("Metabolome_HMDB_20220712/RandomAssociations/",files[b],sep=""),stringsAsFactors = F,sep="\t")
      randommatrix[,b]<-randomnetwork$V3
      print(paste(round((b/length(files))*100,2),"%",sep=""))
    }
    resultmatrix<-matrix(ncol=6,nrow=length(originalnetwork$V1),NA)
    for(a in 1:length(originalnetwork$V1)){
      # a<-1
      thepvalup<-length(which(randommatrix[a,]>=originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvaldown<-length(which(randommatrix[a,]<=originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvalupni<-length(which(randommatrix[a,]>originalnetwork$V3[a]))/length(randommatrix[a,])
      thepvaldownni<-length(which(randommatrix[a,]<originalnetwork$V3[a]))/length(randommatrix[a,])
      meandistance<-mean(randommatrix[a,])
      sddistance<-sd(randommatrix[a,])
      resultmatrix[a,]<-c(thepvalup,thepvaldown,thepvalupni,thepvaldownni,meandistance,sddistance)
      print(paste(round((a/length(originalnetwork$V1))*100,2),"%",sep=""))
    }
    resultmatrix<-cbind(originalnetwork,resultmatrix)
    colnames(resultmatrix)<-c("disease1","disease2","cosine_similarity","pval_up","pval_down","pval_up_not_equal","pval_down_not_equal","mean","sd")
    write.table(resultmatrix,"Metabolome_HMDB_20220712/ICD10_CosineSimilarities_pvals.txt",quote=F,sep="\t",row.names=F)
  }
}

#### Microbiome ####
## Cosine Similarity done ##
## Jaccard index pending ##
if(args[1]=="microbiome"){
  if(args[2]=="disbiome"){
    if(args[3]=="getthenewvectors"){
      ## Sum the disease - symptom co-occurrences ##
      if("ICD10_microbiome_associations.txt"%in%list.files("microbiome_disbiome_20220714/")==FALSE){
        ## Load the disease name - icd10 association file
        microbiomenametoicd<-fread("microbiome_disbiome_20220714/Disease_ICD10_translation.txt",stringsAsFactors = F,sep="\t")
        microbiomenametoicd<-microbiomenametoicd[which(nchar(microbiomenametoicd$ICD10)==3)]
        ## Load the disease-symptom association file
        diseasemicrobiome<-fread("microbiome_disbiome_20220714/DiseaseMicrobeTable.txt",stringsAsFactors = F,sep="\t")
        ## For simplicity we focus only on blood levels ##
        diseasemicrobiome<-diseasemicrobiome[which(diseasemicrobiome$sample_name=="Faeces")]
        diseasemicrobiome$disease_name<-tolower(diseasemicrobiome$disease_name)
        diseasemicrobiome<-diseasemicrobiome[,c(8,10,3)]
        
        ## Calculate for each microbiome an inverse score
        microbiomes<-table(diseasemicrobiome$organism_name)
        microbiomes<-1/microbiomes
        setkey(diseasemicrobiome,"organism_name")
        diseasemicrobiome<-cbind(diseasemicrobiome,as.numeric(microbiomes[diseasemicrobiome$organism_name]))
        colnames(diseasemicrobiome)[4]<-"IFS"
        ## Select the diseases for which we have and ICD correspondence and do the union of their microbiomes ##
        icds<-unique(microbiomenametoicd$ICD10)
        ## Build an ICD10 - microbiome association table ##
        icdmicrobiome<-c()
        for(a in 1:length(icds)){
          # a<-1
          diseases<-microbiomenametoicd$Name[which(microbiomenametoicd$ICD10==icds[a])]
          diseases<-intersect(diseasemicrobiome$disease_name,diseases)
          if(length(diseases)>0){
            icdtab<-c()
            for(b in 1:length(diseases)){
              # b<-1
              icdtab<-rbind(icdtab,diseasemicrobiome[which(diseasemicrobiome$disease_name==diseases[b]),])
            }
            if(length(which(duplicated(icdtab[,1:3])))>0){icdtab<-icdtab[-which(duplicated(icdtab[,]))]}
            icdtab<-cbind(icds[a],icdtab)
            colnames(icdtab)<-c("icd10","disease","microbiome","direction","ifs")
            icdmicrobiome<-rbind(icdmicrobiome,icdtab)
            print(paste(round((a/length(icds))*100,2),"%",sep=""))
          }
        }
        icdmicrobiome<-icdmicrobiome[,-2]
        ## Write the table containing the number of ICD10 - symptom occurrences ##
        write.table(icdmicrobiome,"microbiome_disbiome_20220714/ICD10_microbiome_associations.txt",quote=F,sep="\t",row.names=F)
      }
    }
    if(args[3]=="calculatesimilarities"){
      if("RandomAssociations"%in%list.files("microbiome_disbiome_20220714/")==F){dir.create("microbiome_disbiome_20220714/RandomAssociations")}
      associations<-fread("microbiome_disbiome_20220714/ICD10_microbiome_associations.txt",stringsAsFactors = F,sep="\t")
      ## Work with ifs and forget about direction ##
      if(length(which(duplicated(associations[,c(1,2,4)])))>0){associations<-associations[-which(duplicated(associations[,c(1,2,4)])),c(1,2,4)]}
      ## First create a matrix with the ifs of each microbiome in each disease ##
      diseases<-unique(associations$icd10) ; microbiomes<-unique(associations$microbiome)
      icdmicrobiome_matrix<-matrix(ncol=length(diseases),nrow=length(microbiomes),0) ; colnames(icdmicrobiome_matrix)<-diseases ; rownames(icdmicrobiome_matrix)<-microbiomes
      for(a in 1:length(microbiomes)){
        icdmicrobiome_matrix[microbiomes[a],associations$icd10[which(associations$microbiome==microbiomes[a])]]<-associations$ifs[which(associations$microbiome==microbiomes[a])]
      }
      if(args[4]=="CreateIFSmatrix"){
        write.table(icdmicrobiome_matrix,"microbiome_disbiome_20220714/ICD10_microbiome_ifs_matrix.txt",quote=F,sep="\t")
        ## Calculate cosine similarities ##
        cosinesimilarities<-cosine(icdmicrobiome_matrix)
        ind <- which(upper.tri(cosinesimilarities, diag = TRUE), arr.ind = TRUE)
        netmode2<-cbind(colnames(cosinesimilarities)[ind[,1]],colnames(cosinesimilarities)[ind[,2]], cosinesimilarities[ind])
        netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
        write.table(netmode2,"microbiome_disbiome_20220714/ICD10_CosineSimilarities.txt",quote=F,sep="\t",row.names = F,col.names = F)
      }
      if(args[4]=="random_distances"){
        ## Calculate random cosine similarities ##
        set.seed(as.numeric(args[5]))
        for(a in 1:2000){
          ini<-Sys.time()
          ricdmicrobiome_matrix<-apply(icdmicrobiome_matrix,2,sample)
          rcosinesimilarities<-cosine(ricdmicrobiome_matrix)
          ind <- which(upper.tri(rcosinesimilarities, diag = TRUE), arr.ind = TRUE)
          netmode2<-cbind(colnames(rcosinesimilarities)[ind[,1]],colnames(rcosinesimilarities)[ind[,2]], rcosinesimilarities[ind])
          netmode2<-netmode2[-which(netmode2[,1]==netmode2[,2]),]
          write.table(netmode2,paste("microbiome_disbiome_20220714/RandomAssociations/Seed_",args[5],"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
          fin<-Sys.time()
          print(paste("Seed ",args[5]," round ",a," time: ",fin-ini,sep=""))
        }
      }
    }
    if(args[3]=="calculatesignificance"){
      originalnetwork<-fread("microbiome_disbiome_20220714/ICD10_CosineSimilarities.txt",stringsAsFactors = F,sep="\t")
      files<-list.files("microbiome_disbiome_20220714/RandomAssociations/")
      randommatrix<-matrix(ncol=10000,nrow=length(originalnetwork$V1),NA)
      for(b in 1:length(files)){
        # b<-1
        randomnetwork<-fread(paste("microbiome_disbiome_20220714/RandomAssociations/",files[b],sep=""),stringsAsFactors = F,sep="\t")
        randommatrix[,b]<-randomnetwork$V3
        print(paste(round((b/length(files))*100,2),"%",sep=""))
      }
      resultmatrix<-matrix(ncol=6,nrow=length(originalnetwork$V1),NA)
      for(a in 1:length(originalnetwork$V1)){
        # a<-1
        thepvalup<-length(which(randommatrix[a,]>=originalnetwork$V3[a]))/length(randommatrix[a,])
        thepvaldown<-length(which(randommatrix[a,]<=originalnetwork$V3[a]))/length(randommatrix[a,])
        thepvalupni<-length(which(randommatrix[a,]>originalnetwork$V3[a]))/length(randommatrix[a,])
        thepvaldownni<-length(which(randommatrix[a,]<originalnetwork$V3[a]))/length(randommatrix[a,])
        meandistance<-mean(randommatrix[a,])
        sddistance<-sd(randommatrix[a,])
        resultmatrix[a,]<-c(thepvalup,thepvaldown,thepvalupni,thepvaldownni,meandistance,sddistance)
        print(paste(round((a/length(originalnetwork$V1))*100,2),"%",sep=""))
      }
      resultmatrix<-cbind(originalnetwork,resultmatrix)
      colnames(resultmatrix)<-c("disease1","disease2","cosine_similarity","pval_up","pval_down","pval_up_not_equal","pval_down_not_equal","mean","sd")
      write.table(resultmatrix,"microbiome_disbiome_20220714/ICD10_CosineSimilarities_pvals.txt",quote=F,sep="\t",row.names=F)
    }
  }
  if(args[2]=="GMMAD"){
    if(args[3]=="getthenewvectors"){
      ## Sum the disease - symptom co-occurrences ##
      if("ICD10_microbiome_associations.txt"%in%list.files("Microbiome_GMMAD_20240613/")==FALSE){
        ## Load the disease name - icd10 association file
        microbiomenametoicd<-fread("Microbiome_GMMAD_20240613/Disease_ICD10_translation.txt",stringsAsFactors = F,sep="\t")
        microbiomenametoicd<-microbiomenametoicd[which(nchar(microbiomenametoicd$icd10)==3)] ; setkey(microbiomenametoicd,"icd10")
        ## Load the disease-symptom association file
        diseasemicrobiome<-fread("Microbiome_GMMAD_20240613/disease_species.csv",stringsAsFactors = F,sep=",")
        diseasemicrobiome<-diseasemicrobiome[,c(2,3,16)]
        ## Select one representative for ICDs with more than 1 associated diseases ##
        duplis<-c("M45","P07","K90","K76","K05","F84","E66","B20")

        repes<-microbiomenametoicd[duplis]
        setkey(repes,"disease")
        diseases<-c("Spondylarthritis","Infant, Premature","Celiac Disease","Non-alcoholic Fatty Liver Disease","Gingivitis","Autism Spectrum Disorder",
                    "Obesity","HIV")
        repeselection<-repes[diseases]
        for(a in 1:length(duplis)){
          microbiomenametoicd<-microbiomenametoicd[-which(microbiomenametoicd$icd10==duplis[a])]
        }
        microbiomenametoicd<-rbind(microbiomenametoicd,repeselection)
        icdmicrobtab<-c()
        for(a in 1:length(microbiomenametoicd$disease)){
          # a<-1
          icdmicrobtab<-rbind(icdmicrobtab,cbind(microbiomenametoicd$icd10[a],diseasemicrobiome[which(diseasemicrobiome$disease==microbiomenametoicd$disease[a]),c(2:3)]))
        }
        colnames(icdmicrobtab)[1]<-"icd10"
        ## Write the table containing the ICD10 - microbe associations ##
        write.table(icdmicrobtab,"Microbiome_GMMAD_20240613/ICD10_microbiome_associations.txt",quote=F,sep="\t",row.names=F)
      }
    }
  }
}




#### Layer overlaps ####
## @@ @@ @@  @@ @@ @@ ##
if(args[1]=="layer_overlaps"){
  if(args[2]=="preparetheintralayers"){
    ## mirna ##
    ## @@ @@ ##
    mirna<-fread("miRNA_HMDD_20240523/ICD10_CosineSimilarities_pvals.txt",stringsAsFactors = F,sep="\t")
    mirna<-cbind(mirna,as.numeric(mirna$pval_up)*length(mirna$pval_up),as.numeric(mirna$pval_down)*length(mirna$pval_down))
    colnames(mirna)[c(10,11)]<-c("FDR_up","FDR_down")
    if(length(which(mirna$FDR_up>1))>0){mirna$FDR_up[which(mirna$FDR_up>1)]<-1}
    if(length(which(mirna$FDR_down>1))>0){mirna$FDR_down[which(mirna$FDR_down>1)]<-1}
    mirnasorted<-t(apply(mirna[,1:2], 1, sort))
    mirna$disease1<-mirnasorted[,1] ; mirna$disease2<-mirnasorted[,2]
    mirna<-mirna[,c(1:4,10)] ; colnames(mirna)[4:5]<-c("pval","fdr")
    write.table(mirna,"IntraLayers/miRNA_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names=F)
    
    ## genes ##
    ## @@ @@ ##
    gene<-fread("Genes_DisGeNET_20230113/ICD10_CosineSimilarities_pvals.txt",stringsAsFactors = F,sep="\t")
    gene<-cbind(gene,as.numeric(gene$pval_up)*length(gene$pval_up),as.numeric(gene$pval_down)*length(gene$pval_down))
    colnames(gene)[c(10,11)]<-c("FDR_up","FDR_down")
    if(length(which(gene$FDR_up>1))>0){gene$FDR_up[which(gene$FDR_up>1)]<-1}
    if(length(which(gene$FDR_down>1))>0){gene$FDR_down[which(gene$FDR_down>1)]<-1}
    genesorted<-t(apply(gene[,1:2], 1, sort))
    gene$disease1<-genesorted[,1] ; gene$disease2<-genesorted[,2]
    gene<-gene[,c(1:4,10)] ; colnames(gene)[4:5]<-c("pval","fdr")
    write.table(gene,"IntraLayers/Genes_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names=F)
    
    ## drugs ##
    ## @@ @@ ##
    drug<-fread("Drugs_CTD_20220712/ICD10_CosineSimilarities_pvals.txt",stringsAsFactors = F,sep="\t")
    drug<-cbind(drug,as.numeric(drug$pval_up)*length(drug$pval_up),as.numeric(drug$pval_down)*length(drug$pval_down))
    colnames(drug)[c(10,11)]<-c("FDR_up","FDR_down")
    if(length(which(drug$FDR_up>1))>0){drug$FDR_up[which(drug$FDR_up>1)]<-1}
    if(length(which(drug$FDR_down>1))>0){drug$FDR_down[which(drug$FDR_down>1)]<-1}
    drugsorted<-t(apply(drug[,1:2], 1, sort))
    drug$disease1<-drugsorted[,1] ; drug$disease2<-drugsorted[,2]
    drug<-drug[,c(1:4,10)] ; colnames(drug)[4:5]<-c("pval","fdr")
    write.table(drug,"IntraLayers/Drugs_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names=F)
    
    # ## symptoms ##
    # ## @@ @@ @@ ##
    # symptom<-fread("Symptoms_Zhou2014_20220712/ICD10_CosineSimilarities_pvals.txt",stringsAsFactors = F,sep="\t")
    # symptom<-cbind(symptom,as.numeric(symptom$pval_up)*length(symptom$pval_up),as.numeric(symptom$pval_down)*length(symptom$pval_down))
    # colnames(symptom)[c(10,11)]<-c("FDR_up","FDR_down")
    # if(length(which(symptom$FDR_up>1))>0){symptom$FDR_up[which(symptom$FDR_up>1)]<-1}
    # if(length(which(symptom$FDR_down>1))>0){symptom$FDR_down[which(symptom$FDR_down>1)]<-1}
    # symptomsorted<-t(apply(symptom[,1:2], 1, sort))
    # symptom$disease1<-symptomsorted[,1] ; symptom$disease2<-symptomsorted[,2]
    # symptom<-symptom[,c(1:4,10)] ; colnames(symptom)[4:5]<-c("pval","fdr")
    # write.table(symptom,"IntraLayers/Symptoms_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names=F)
    
    ## symptoms significantly associated to diseases ##
    ## @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ ##
    symptom<-fread("Symptoms_Zhou2014_20220712/ICD10_significant_CosineSimilarities_pvals.txt",stringsAsFactors = F,sep="\t")
    symptom<-cbind(symptom,as.numeric(symptom$pval_up)*length(symptom$pval_up),as.numeric(symptom$pval_down)*length(symptom$pval_down))
    colnames(symptom)[c(10,11)]<-c("FDR_up","FDR_down")
    if(length(which(symptom$FDR_up>1))>0){symptom$FDR_up[which(symptom$FDR_up>1)]<-1}
    if(length(which(symptom$FDR_down>1))>0){symptom$FDR_down[which(symptom$FDR_down>1)]<-1}
    symptomsorted<-t(apply(symptom[,1:2], 1, sort))
    symptom$disease1<-symptomsorted[,1] ; symptom$disease2<-symptomsorted[,2]
    symptom<-symptom[,c(1:4,10)] ; colnames(symptom)[4:5]<-c("pval","fdr")
    write.table(symptom,"IntraLayers/SymptomsSignificant_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names=F)
    
    ## transcriptomics ##
    ## @@ @@ @ @ @@ @@ ##
    transcriptomics<-fread("Transcriptomics/Microarrays_JSV/ICD10_Adjusted_fdr_network.txt",stringsAsFactors = F,sep="\t")
    transcriptomics<-transcriptomics[,c(1,2,16,17)]
    colnames(transcriptomics)<-c("disease1","disease2","cosine_similarity","fdr")
    write.table(transcriptomics,"IntraLayers/Transcriptomics_ICD10_Adjusting_fdr.txt",quote=F,sep="\t",row.names=F)
    
    ## epidemiology ##
    ## @@ @@  @@ @@ ##
    epidemiology<-fread("Epidemiology_westergaard/41467_2019_8475_MOESM6_ESM.txt",stringsAsFactors = F,sep="\t")
    epidemiology<-epidemiology[,c(1,2,7)]
    epidemiology<-cbind(t(apply(epidemiology[,1:2], 1, sort)),epidemiology$adjustedRR2.5)
    colnames(epidemiology)<-c("disease1","disease2","RR")
    ## we don't take into consideration the NAs ##
    if(length(which(is.na(epidemiology[,3])))>0){epidemiology<-epidemiology[-which(is.na(epidemiology[,3])),]}
    epidemiology<-as.data.table(epidemiology) ; epidemiology$RR<-as.numeric(epidemiology$RR)
    ## Remove duplicated associations
    if(length(which(duplicated(epidemiology)))>0){
      duplis<-epidemiology[which(duplicated(epidemiology))]
      for(a in 1:length(duplis$disease1)){
        # a<-1
        indices<-intersect(which(epidemiology$disease1==duplis$disease1[a]),which(epidemiology$disease2==duplis$disease2[a]))
        epidemiology<-epidemiology[-setdiff(indices,indices[which(epidemiology$RR[indices]==max(epidemiology$RR[indices]))[1]])]
      }
    }
    write.table(epidemiology,"IntraLayers/Epidemiology_ICD10_RelativeRisk.txt",quote=F,sep="\t",row.names=F)
    
    ## metabolome ##
    ## @@ @@@@ @@ ##
    metabolome<-fread("Metabolome_HMDB_20220712/ICD10_CosineSimilarities_pvals.txt",stringsAsFactors = F,sep="\t")
    metabolome<-cbind(metabolome,as.numeric(metabolome$pval_up)*length(metabolome$pval_up),as.numeric(metabolome$pval_down)*length(metabolome$pval_down))
    colnames(metabolome)[c(10,11)]<-c("FDR_up","FDR_down")
    if(length(which(metabolome$FDR_up>1))>0){metabolome$FDR_up[which(metabolome$FDR_up>1)]<-1}
    if(length(which(metabolome$FDR_down>1))>0){metabolome$FDR_down[which(metabolome$FDR_down>1)]<-1}
    metabolomesorted<-t(apply(metabolome[,1:2], 1, sort))
    metabolome$disease1<-metabolomesorted[,1] ; metabolome$disease2<-metabolomesorted[,2]
    metabolome<-metabolome[,c(1:4,10)] ; colnames(metabolome)[4:5]<-c("pval","fdr")
    write.table(metabolome,"IntraLayers/Metabolome_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names=F)
    
    ## microbiome ##
    ## @ @ @@ @ @ ##
    microbiome<-fread("microbiome_disbiome_20220714/ICD10_CosineSimilarities_pvals.txt",stringsAsFactors = F,sep="\t")
    microbiome<-cbind(microbiome,as.numeric(microbiome$pval_up)*length(microbiome$pval_up),as.numeric(microbiome$pval_down)*length(microbiome$pval_down))
    colnames(microbiome)[c(10,11)]<-c("FDR_up","FDR_down")
    if(length(which(microbiome$FDR_up>1))>0){microbiome$FDR_up[which(microbiome$FDR_up>1)]<-1}
    if(length(which(microbiome$FDR_down>1))>0){microbiome$FDR_down[which(microbiome$FDR_down>1)]<-1}
    microbiomesorted<-t(apply(microbiome[,1:2], 1, sort))
    microbiome$disease1<-microbiomesorted[,1] ; microbiome$disease2<-microbiomesorted[,2]
    microbiome<-microbiome[,c(1:4,10)] ; colnames(microbiome)[4:5]<-c("pval","fdr")
    write.table(microbiome,"IntraLayers/Microbiome_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names=F)
    
  }
  if(args[2]=="plotoverlaps"){
    #### Load the layers ####
    ## mirna ##
    mirna<-fread("IntraLayers/miRNA_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    mirnadiseases<-unique(c(mirna$disease1,mirna$disease2))
    mirnasig<-mirna[which(mirna$fdr<=0.05)]
    ## genes ##
    gene<-fread("IntraLayers/Genes_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    genediseases<-unique(c(gene$disease1,gene$disease2))
    genesig<-gene[which(gene$fdr<=0.05)]
    ## drugs ##
    drug<-fread("IntraLayers/Drugs_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    drugdiseases<-unique(c(drug$disease1,drug$disease2))
    drugsig<-drug[which(drug$fdr<=0.05)]
    ## symptoms ##
    symptom<-fread("IntraLayers/SymptomsSignificant_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    symptomdiseases<-unique(c(symptom$disease1,symptom$disease2))
    symptomsig<-symptom[which(symptom$fdr<=0.05)]
    ## transcriptomics ##
    transcriptomics<-fread("IntraLayers/Transcriptomics_ICD10_Adjusting_fdr.txt",stringsAsFactors = F,sep="\t")
    transcriptomicspos<-transcriptomics[which(transcriptomics$cosine_similarity>0),]
    transcriptomicsdiseases<-unique(c(transcriptomics$disease1,transcriptomics$disease2))
    transcriptomicspossig<-transcriptomicspos[intersect(which(transcriptomicspos$cosine_similarity>0),which(transcriptomicspos$fdr<=0.05)),]
    ## epidemiology ##
    epidemiology<-fread("IntraLayers/Epidemiology_ICD10_RelativeRisk.txt",stringsAsFactors = F,sep="\t")
    epidemiologydiseases<-unique(c(epidemiology$disease1,epidemiology$disease2))
    epidemiologysig<-epidemiology[which(epidemiology$RR>1.01)]
    ## metabolome ##
    metabolome<-fread("IntraLayers/Metabolome_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    metabolomediseases<-unique(c(metabolome$disease1,metabolome$disease2))
    metabolomesig<-metabolome[which(metabolome$fdr<=0.05)]
    ## microbiome ##
    microbiome<-fread("IntraLayers/Microbiome_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    microbiomediseases<-unique(c(microbiome$disease1,microbiome$disease2))
    microbiomesig<-microbiome[which(microbiome$fdr<=0.05)]
    
    
    ## Write a summary of the edges in each network ##
    edgesummary<-rbind(c(length(mirnasig$cosine_similarity),length(which(mirna$cosine_similarity!=0)),length(mirna$cosine_similarity)),
                       c(length(drugsig$cosine_similarity),length(which(drug$cosine_similarity!=0)),length(drug$cosine_similarity)),
                       c(length(genesig$cosine_similarity),length(which(gene$cosine_similarity!=0)),length(gene$cosine_similarity)),
                       c(length(symptomsig$cosine_similarity),length(which(symptom$cosine_similarity!=0)),length(symptom$cosine_similarity)),
                       c(length(which(transcriptomicspossig$cosine_similarity!=0)),length(which(transcriptomicspos$cosine_similarity!=0)),length(transcriptomics$cosine_similarity)),
                       c(length(metabolomesig$cosine_similarity),length(which(metabolome$cosine_similarity!=0)),length(metabolome$cosine_similarity)),
                       c(length(microbiomesig$cosine_similarity),length(which(microbiome$cosine_similarity!=0)),length(microbiome$cosine_similarity)),
                       c(length(epidemiologysig$disease1),length(epidemiologysig$disease1),length(epidemiology$disease1)))
    edgesummary<-cbind(edgesummary,(edgesummary[,1]/edgesummary[,3])*100,(edgesummary[,2]/edgesummary[,3])*100)
    colnames(edgesummary)<-c("significant edges","edges","potential edges","% significant","% edges")
    rownames(edgesummary)<-c("miRNA","drug","gene","symptom","transcriptomic","epidemiology")
    edgesummary<-cbind(c(length(mirnadiseases),length(drugdiseases),length(genediseases),length(symptomdiseases),
                         length(transcriptomicsdiseases),length(metabolomediseases),length(microbiomediseases),length(epidemiologydiseases)),
                       edgesummary)
    colnames(edgesummary)[1]<-"nodes"
    edgesummary<-cbind(c("miRNAs","drugs","genes","symptoms","transcriptomics","metabolome","microbiome","EHRs"),edgesummary)
    edgesummary<-edgesummary[order(as.numeric(edgesummary[,2]),decreasing = T),]
    edgesummary<-as.data.frame(edgesummary)
    colnames(edgesummary)[1]<-"layer"
    write.table(edgesummary,"Number_of_nodes_and_edges.txt",quote=F,sep="\t")
    
    #### Excluding transcriptomic data ####
    #### @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ####
    if(args[3]=="exclude_transcriptomics"){
      notranscriptomicdiseases<-list("miRNA"=mirnadiseases,"genes"=genediseases,"drugs"=drugdiseases,"symptoms"=symptomdiseases,"epidemiology"=epidemiologydiseases)
      notranscriptomicshareddiseases<-Reduce(intersect,notranscriptomicdiseases)
      ## Get the interactions between shared diseases ##
      indmirna1<-c() ; indmirna2<-c() ; indgene1<-c() ; indgene2<-c() ; inddrug1<-c() ; inddrug2<-c() ; indsymptom1<-c() ; indsymptom2<-c() ; indepidemiology1<-c() ; indepidemiology2<-c()
      for(a in 1:length(notranscriptomicshareddiseases)){
        # a<-1
        indmirna1<-c(indmirna1,which(mirna$disease1==notranscriptomicshareddiseases[a]))
        indmirna2<-c(indmirna2,which(mirna$disease2==notranscriptomicshareddiseases[a]))
        indgene1<-c(indgene1,which(gene$disease1==notranscriptomicshareddiseases[a]))
        indgene2<-c(indgene2,which(gene$disease2==notranscriptomicshareddiseases[a]))
        inddrug1<-c(inddrug1,which(drug$disease1==notranscriptomicshareddiseases[a]))
        inddrug2<-c(inddrug2,which(drug$disease2==notranscriptomicshareddiseases[a]))
        indsymptom1<-c(indsymptom1,which(symptom$disease1==notranscriptomicshareddiseases[a]))
        indsymptom2<-c(indsymptom2,which(symptom$disease2==notranscriptomicshareddiseases[a]))
        indepidemiology1<-c(indepidemiology1,which(epidemiology$disease1==notranscriptomicshareddiseases[a]))
        indepidemiology2<-c(indepidemiology2,which(epidemiology$disease2==notranscriptomicshareddiseases[a]))
        print(paste(round((a/length(notranscriptomicshareddiseases))*100,2),"%",sep=""))
      }
      ## All interactions ##
      notrans_mirna<-mirna[intersect(indmirna1,indmirna2)]
      notrans_gene<-gene[intersect(indgene1,indgene2)]
      notrans_drug<-drug[intersect(inddrug1,inddrug2)]
      notrans_symptom<-symptom[intersect(indsymptom1,indsymptom2)]
      notrans_epidemiology<-epidemiology[intersect(indepidemiology1,indepidemiology2)]
      ## Positive interactions ##
      notrans_mirna_pos<-notrans_mirna[which(notrans_mirna$cosine_similarity>0)]
      notrans_gene_pos<-notrans_gene[which(notrans_gene$cosine_similarity>0)]
      notrans_drug_pos<-notrans_drug[which(notrans_drug$cosine_similarity>0)]
      notrans_symptom_pos<-notrans_symptom[which(notrans_symptom$cosine_similarity>0)]
      notrans_epidemiology_pos<-notrans_epidemiology[which(notrans_epidemiology$RR>1.01)]
      ## Significant interactions ##
      notrans_mirna_pos_sig<-notrans_mirna_pos[which(notrans_mirna_pos$fdr<=0.05)]
      notrans_gene_pos_sig<-notrans_gene_pos[which(notrans_gene_pos$fdr<=0.05)]
      notrans_drug_pos_sig<-notrans_drug_pos[which(notrans_drug_pos$fdr<=0.05)]
      notrans_symptom_pos_sig<-notrans_symptom_pos[which(notrans_symptom_pos$fdr<=0.05)]
      ## List with the positive links ##
      notrans_posint<-list("miRNA"=paste(notrans_mirna_pos$disease1,notrans_mirna_pos$disease2,sep="_"),
                           "gene"=paste(notrans_gene_pos$disease1,notrans_gene_pos$disease2,sep="_"),
                           "drug"=paste(notrans_drug_pos$disease1,notrans_drug_pos$disease2,sep="_"),
                           "symptom"=paste(notrans_symptom_pos$disease1,notrans_symptom_pos$disease2,sep="_"),
                           "epidemiology"=paste(notrans_epidemiology_pos$disease1,notrans_epidemiology_pos$disease2,sep="_"))
      ## List with significant links ##
      notrans_posint_sig<-list("miRNA"=paste(notrans_mirna_pos_sig$disease1,notrans_mirna_pos_sig$disease2,sep="_"),
                               "gene"=paste(notrans_gene_pos_sig$disease1,notrans_gene_pos_sig$disease2,sep="_"),
                               "drug"=paste(notrans_drug_pos_sig$disease1,notrans_drug_pos_sig$disease2,sep="_"),
                               "symptom"=paste(notrans_symptom_pos_sig$disease1,notrans_symptom_pos_sig$disease2,sep="_"),
                               "epidemiology"=paste(notrans_epidemiology_pos$disease1,notrans_epidemiology_pos$disease2,sep="_"))
      
      ## Upset plot excluding the transcriptomics ##
      ## @@ @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ @@ ##
      ## Diseases ##
      input<-list("miRNA"=mirnadiseases,"genes"=genediseases,"drugs"=drugdiseases,"symptoms"=symptomdiseases,"epidemiology"=epidemiologydiseases)
      tabla<-fromList(input)
      pdf(file = "VennDiagrams/Upsetplot_no_transcriptomic_layer_diseases.pdf")
        upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
      dev.off()
      ## Interactions ## 
      tabla<-fromList(notrans_posint)
      pdf(file = "VennDiagrams/Upsetplot_no_transcriptomic_layer_interactions.pdf")
        upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
      dev.off()
      ## Significant interactions ## 
      tabla<-fromList(notrans_posint_sig)
      pdf(file = "VennDiagrams/Upsetplot_no_transcriptomic_layer_significant_interactions.pdf")
        upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
      dev.off()
      
      ## Write a summary of the edges in each network ##
      notrans_edgesummary<-cbind(c(length(unique(c(notrans_mirna_pos_sig$disease1,notrans_mirna_pos_sig$disease2))),length(unique(c(notrans_drug_pos_sig$disease1,notrans_drug_pos_sig$disease2))),length(unique(c(notrans_gene_pos_sig$disease1,notrans_drug_pos_sig$disease2))),length(unique(c(notrans_symptom_pos_sig$disease1,notrans_symptom_pos_sig$disease2))),length(unique(c(notrans_epidemiology_pos$disease1,notrans_epidemiology_pos$disease2)))),
                                 c(length(unique(c(notrans_mirna_pos$disease1,notrans_mirna_pos$disease2))),length(unique(c(notrans_drug_pos$disease1,notrans_drug_pos$disease2))),length(unique(c(notrans_gene_pos$disease1,notrans_drug_pos$disease2))),length(unique(c(notrans_symptom_pos$disease1,notrans_symptom_pos$disease2))),length(unique(c(notrans_epidemiology_pos$disease1,notrans_epidemiology_pos$disease2)))),
                                 c(length(unique(c(notrans_mirna$disease1,notrans_mirna$disease2))),length(unique(c(notrans_drug$disease1,notrans_drug$disease2))),length(unique(c(notrans_gene$disease1,notrans_drug$disease2))),length(unique(c(notrans_symptom$disease1,notrans_symptom$disease2))),length(unique(c(notrans_epidemiology$disease1,notrans_epidemiology$disease2)))),
                                 c(length(notrans_mirna_pos_sig$disease1),length(notrans_drug_pos_sig$disease1),length(notrans_gene_pos_sig$disease1),length(notrans_symptom_pos_sig$disease1),length(notrans_epidemiology_pos$disease1)),
                                 c(length(notrans_mirna_pos$disease1),length(notrans_drug_pos$disease1),length(notrans_gene_pos$disease1),length(notrans_symptom_pos$disease1),length(notrans_epidemiology_pos$disease1)),
                                 c(length(notrans_mirna$disease1),length(notrans_drug$disease1),length(notrans_gene$disease1),length(notrans_symptom$disease1),length(notrans_epidemiology$disease1)))
      colnames(notrans_edgesummary)<-c("nodes_se","nodes_e","nodes_pe","edges_se","edges_e","edges_pe")
      rownames(notrans_edgesummary)<-c("miRNA","drug","gene","symptom","epidemiology")
      notrans_edgesummary<-cbind(notrans_edgesummary,(notrans_edgesummary[,4]/notrans_edgesummary[,6])*100,
                                 (notrans_edgesummary[,5]/notrans_edgesummary[,6])*100)
      colnames(notrans_edgesummary)[c(7,8)]<-c("%_se","%_e")
      write.table(notrans_edgesummary,"Number_of_nodes_and_edges_excluding_transcriptomics.txt",quote=F,sep="\t")
    }
    #### Including transcriptomic data ####
    #### @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ####
    if(args[3]=="include_transcriptomics"){
      alldiseases<-list("transcriptomics"=transcriptomicsdiseases,"miRNA"=mirnadiseases,"genes"=genediseases,"drugs"=drugdiseases,"symptoms"=symptomdiseases,"epidemiology"=epidemiologydiseases,"microbiome"=microbiomediseases,"metabolome"=metabolomediseases)
      allshareddiseases<-Reduce(intersect,alldiseases)
      ## Get the interactions between shared diseases ##
      indmirna1<-c() ; indmirna2<-c() ; indgene1<-c() ; indgene2<-c() ; inddrug1<-c() ; inddrug2<-c()
      indsymptom1<-c() ; indsymptom2<-c() ; indtranscriptomics1<-c() ; indtranscriptomics2<-c()
      indepidemiology1<-c() ; indepidemiology2<-c() ; indmetabolome1<-c() ; indmetabolome2<-c() ; indmicrobiome1<-c() ; indmicrobiome2<-c()
      for(a in 1:length(allshareddiseases)){
        # a<-1
        indmirna1<-c(indmirna1,which(mirna$disease1==allshareddiseases[a]))
        indmirna2<-c(indmirna2,which(mirna$disease2==allshareddiseases[a]))
        indgene1<-c(indgene1,which(gene$disease1==allshareddiseases[a]))
        indgene2<-c(indgene2,which(gene$disease2==allshareddiseases[a]))
        inddrug1<-c(inddrug1,which(drug$disease1==allshareddiseases[a]))
        inddrug2<-c(inddrug2,which(drug$disease2==allshareddiseases[a]))
        indsymptom1<-c(indsymptom1,which(symptom$disease1==allshareddiseases[a]))
        indsymptom2<-c(indsymptom2,which(symptom$disease2==allshareddiseases[a]))
        indtranscriptomics1<-c(indtranscriptomics1,which(transcriptomics$disease1==allshareddiseases[a]))
        indtranscriptomics2<-c(indtranscriptomics2,which(transcriptomics$disease2==allshareddiseases[a]))
        indepidemiology1<-c(indepidemiology1,which(epidemiology$disease1==allshareddiseases[a]))
        indepidemiology2<-c(indepidemiology2,which(epidemiology$disease2==allshareddiseases[a]))
        indmetabolome1<-c(indmetabolome1,which(metabolome$disease1==allshareddiseases[a]))
        indmetabolome2<-c(indmetabolome2,which(metabolome$disease2==allshareddiseases[a]))
        indmicrobiome1<-c(indmicrobiome1,which(microbiome$disease1==allshareddiseases[a]))
        indmicrobiome2<-c(indmicrobiome2,which(microbiome$disease2==allshareddiseases[a]))
        print(paste(round((a/length(allshareddiseases))*100,2),"%",sep=""))
      }
      ## All interactions ##
      alldis_mirna<-mirna[intersect(indmirna1,indmirna2)]
      alldis_gene<-gene[intersect(indgene1,indgene2)]
      alldis_drug<-drug[intersect(inddrug1,inddrug2)]
      alldis_symptom<-symptom[intersect(indsymptom1,indsymptom2)]
      alldis_transcriptomics<-transcriptomics[intersect(indtranscriptomics1,indtranscriptomics2)]
      alldis_epidemiology<-epidemiology[intersect(indepidemiology1,indepidemiology2)]
      alldis_microbiome<-microbiome[intersect(indmicrobiome1,indmicrobiome2)]
      alldis_metabolome<-metabolome[intersect(indmetabolome1,indmetabolome2)]
      ## Positive interactions ##
      alldis_mirna_pos<-alldis_mirna[which(alldis_mirna$cosine_similarity>0)]
      alldis_gene_pos<-alldis_gene[which(alldis_gene$cosine_similarity>0)]
      alldis_drug_pos<-alldis_drug[which(alldis_drug$cosine_similarity>0)]
      alldis_symptom_pos<-alldis_symptom[which(alldis_symptom$cosine_similarity>0)]
      alldis_transcriptomics_pos<-alldis_transcriptomics[which(alldis_transcriptomics$cosine_similarity>0)]
      alldis_epidemiology_pos<-alldis_epidemiology[which(alldis_epidemiology$RR>1.01)]
      alldis_metabolome_pos<-alldis_metabolome[which(alldis_metabolome$cosine_similarity>0)]
      alldis_microbiome_pos<-alldis_microbiome[which(alldis_microbiome$cosine_similarity>0)]
      ## Significant interactions ##
      alldis_mirna_pos_sig<-alldis_mirna_pos[which(alldis_mirna_pos$fdr<=0.05)]
      alldis_gene_pos_sig<-alldis_gene_pos[which(alldis_gene_pos$fdr<=0.05)]
      alldis_drug_pos_sig<-alldis_drug_pos[which(alldis_drug_pos$fdr<=0.05)]
      alldis_symptom_pos_sig<-alldis_symptom_pos[which(alldis_symptom_pos$fdr<=0.05)]
      alldis_transcriptomics_pos_sig<-alldis_transcriptomics_pos[which(alldis_transcriptomics_pos$fdr<=0.05)]
      alldis_microbiome_pos_sig<-alldis_microbiome_pos[which(alldis_microbiome_pos$fdr<=0.05)]
      alldis_metabolome_pos_sig<-alldis_metabolome_pos[which(alldis_metabolome_pos$fdr<=0.05)]
      ## List with the positive links ##
      alldis_posint<-list("miRNA"=paste(alldis_mirna_pos$disease1,alldis_mirna_pos$disease2,sep="_"),
                          "gene"=paste(alldis_gene_pos$disease1,alldis_gene_pos$disease2,sep="_"),
                          "drug"=paste(alldis_drug_pos$disease1,alldis_drug_pos$disease2,sep="_"),
                          "symptom"=paste(alldis_symptom_pos$disease1,alldis_symptom_pos$disease2,sep="_"),
                          "transcriptomics"=paste(alldis_transcriptomics_pos$disease1,alldis_transcriptomics_pos$disease2,sep="_"),
                          "metabolome"=paste(alldis_metabolome_pos$disease1,alldis_metabolome_pos$disease2,sep="_"),
                          "microbiome"=paste(alldis_microbiome_pos$disease1,alldis_microbiome_pos$disease2,sep="_"),
                          "epidemiology"=paste(alldis_epidemiology_pos$disease1,alldis_epidemiology_pos$disease2,sep="_"))
      ## List with significant links ##
      alldis_posint_sig<-list("miRNA"=paste(alldis_mirna_pos_sig$disease1,alldis_mirna_pos_sig$disease2,sep="_"),
                              "gene"=paste(alldis_gene_pos_sig$disease1,alldis_gene_pos_sig$disease2,sep="_"),
                              "drug"=paste(alldis_drug_pos_sig$disease1,alldis_drug_pos_sig$disease2,sep="_"),
                              "symptom"=paste(alldis_symptom_pos_sig$disease1,alldis_symptom_pos_sig$disease2,sep="_"),
                              "transcriptomics"=paste(alldis_transcriptomics_pos_sig$disease1,alldis_transcriptomics_pos_sig$disease2,sep="_"),
                              "metabolome"=paste(alldis_metabolome_pos_sig$disease1,alldis_metabolome_pos_sig$disease2,sep="_"),
                              "microbiome"=paste(alldis_microbiome_pos_sig$disease1,alldis_microbiome_pos_sig$disease2,sep="_"),
                              "epidemiology"=paste(alldis_epidemiology_pos$disease1,alldis_epidemiology_pos$disease2,sep="_"))
      
      ## Upset plot with all the networks ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
      ## Diseases ##
      input<-list("transcriptomics"=transcriptomicsdiseases,"microbiome"=microbiomediseases,"metabolome"=metabolomediseases,"miRNA"=mirnadiseases,"genes"=genediseases,"drugs"=drugdiseases,"symptoms"=symptomdiseases,"epidemiology"=epidemiologydiseases)
      tabla<-fromList(input)
      pdf(file = "VennDiagrams/Upsetplot_all_layers_diseases.pdf")
        upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
      dev.off()
      ## Interactions ## 
      tabla<-fromList(alldis_posint)
      pdf(file = "VennDiagrams/Upsetplot_all_layers_interactions.pdf")
        upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
      dev.off()
      ## Significant interactions ## 
      tabla<-fromList(alldis_posint_sig)
      pdf(file = "VennDiagrams/Upsetplot_all_layers_significant_interactions.pdf")
        upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
      dev.off()
      
      ## Write a summary of the edges in each network ##
      alldis_edgesummary<-cbind(c(length(unique(c(alldis_mirna_pos_sig$disease1,alldis_mirna_pos_sig$disease2))),length(unique(c(alldis_drug_pos_sig$disease1,alldis_drug_pos_sig$disease2))),length(unique(c(alldis_gene_pos_sig$disease1,alldis_drug_pos_sig$disease2))),length(unique(c(alldis_symptom_pos_sig$disease1,alldis_symptom_pos_sig$disease2))),length(unique(c(alldis_epidemiology_pos$disease1,alldis_epidemiology_pos$disease2))),length(unique(c(alldis_transcriptomics_pos_sig$disease1,alldis_transcriptomics_pos_sig$disease2)))),
                                c(length(unique(c(alldis_mirna_pos$disease1,alldis_mirna_pos$disease2))),length(unique(c(alldis_drug_pos$disease1,alldis_drug_pos$disease2))),length(unique(c(alldis_gene_pos$disease1,alldis_drug_pos$disease2))),length(unique(c(alldis_symptom_pos$disease1,alldis_symptom_pos$disease2))),length(unique(c(alldis_epidemiology_pos$disease1,alldis_epidemiology_pos$disease2))),length(unique(c(alldis_transcriptomics_pos$disease1,alldis_transcriptomics_pos$disease2)))),
                                c(length(unique(c(alldis_mirna$disease1,alldis_mirna$disease2))),length(unique(c(alldis_drug$disease1,alldis_drug$disease2))),length(unique(c(alldis_gene$disease1,alldis_drug$disease2))),length(unique(c(alldis_symptom$disease1,alldis_symptom$disease2))),length(unique(c(alldis_epidemiology$disease1,alldis_epidemiology$disease2))),length(unique(c(alldis_transcriptomics$disease1,alldis_transcriptomics$disease2)))),
                                c(length(alldis_mirna_pos_sig$disease1),length(alldis_drug_pos_sig$disease1),length(alldis_gene_pos_sig$disease1),length(alldis_symptom_pos_sig$disease1),length(alldis_epidemiology_pos$disease1),length(alldis_transcriptomics_pos_sig$disease1)),
                                c(length(alldis_mirna_pos$disease1),length(alldis_drug_pos$disease1),length(alldis_gene_pos$disease1),length(alldis_symptom_pos$disease1),length(alldis_epidemiology_pos$disease1),length(alldis_transcriptomics_pos$disease1)),
                                c(length(alldis_mirna$disease1),length(alldis_drug$disease1),length(alldis_gene$disease1),length(alldis_symptom$disease1),length(alldis_epidemiology$disease1),length(alldis_transcriptomics$disease1)))
      colnames(alldis_edgesummary)<-c("nodes_se","nodes_e","nodes_pe","edges_se","edges_e","edges_pe")
      rownames(alldis_edgesummary)<-c("miRNA","drug","gene","symptom","epidemiology","transcriptomics")
      alldis_edgesummary<-cbind(alldis_edgesummary,(alldis_edgesummary[,4]/alldis_edgesummary[,6])*100,
                                (alldis_edgesummary[,5]/alldis_edgesummary[,6])*100)
      colnames(alldis_edgesummary)[c(7,8)]<-c("%_se","%_e")
      write.table(alldis_edgesummary,"Number_of_nodes_and_edges_including_transcriptomics.txt",quote=F,sep="\t")
    }
  }
  if(args[2]=="safefinalnetworks"){
    if("Final"%in%list.files("IntraLayers/")==F){dir.create("IntraLayers/Final")}
    ## mirna ##
    mirna<-fread("IntraLayers/miRNA_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    mirna<-mirna[which(mirna$cosine_similarity!=0)]
    write.table(mirna,"IntraLayers/Final/miRNA_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names = F)
    ## genes ##
    gene<-fread("IntraLayers/Genes_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    gene<-gene[which(gene$cosine_similarity!=0)]
    write.table(gene,"IntraLayers/Final/Genes_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names = F)
    ## drugs ##
    drug<-fread("IntraLayers/Drugs_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    drug<-drug[which(drug$cosine_similarity!=0)]
    write.table(drug,"IntraLayers/Final/Drugs_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names = F)
    # ## symptoms ##
    # symptom<-fread("IntraLayers/Symptoms_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    # symptom<-symptom[which(symptom$cosine_similarity!=0)]
    # write.table(symptom,"IntraLayers/Final/Symptoms_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names = F)
    ## symptoms significantly associated to diseases ##
    symptom<-fread("IntraLayers/SymptomsSignificant_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    symptom<-symptom[which(symptom$cosine_similarity!=0)]
    write.table(symptom,"IntraLayers/Final/Symptoms_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names = F)
    ## transcriptomics ##
    transcriptomics<-fread("IntraLayers/Transcriptomics_ICD10_Adjusting_fdr.txt",stringsAsFactors = F,sep="\t")
    transcriptomics<-transcriptomics[which(transcriptomics$cosine_similarity>0),]
    write.table(transcriptomics,"IntraLayers/Final/Transcriptomics_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names = F)
    ## epidemiology ##
    epidemiology<-fread("IntraLayers/Epidemiology_ICD10_RelativeRisk.txt",stringsAsFactors = F,sep="\t")
    epidemiology<-epidemiology[which(epidemiology$RR>1.01)]
    write.table(epidemiology,"IntraLayers/Final/Epidemiology_ICD10_RelativeRisk.txt",quote=F,sep="\t",row.names = F)
    ## metabolome ##
    metabolome<-fread("IntraLayers/Metabolome_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    metabolome<-metabolome[which(metabolome$cosine_similarity!=0)]
    write.table(metabolome,"IntraLayers/Final/Metabolome_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names = F)
    ## microbiome ##
    microbiome<-fread("IntraLayers/Microbiome_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
    microbiome<-microbiome[which(microbiome$cosine_similarity!=0)]
    write.table(microbiome,"IntraLayers/Final/Microbiome_ICD10_CosineSimilarities_fdrs.txt",quote=F,sep="\t",row.names = F)
  }
}

#### Build all the interlayers ####
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
if(args[1]=="buildinterlayers"){
  #### Connect miRNA to gene ####
  ## Read files ##
  disease_mirna<-fread("mirna_HMDD_20240523/ICD10_mirna_associations.txt",stringsAsFactors = F,sep="\t") ; disease_mirna$mirna<-tolower(disease_mirna$mirna)
  disease_gene<-fread("Genes_DisGeNET_20230113/ICD10_gene_associations.txt",stringsAsFactors = F,sep="\t") ; disease_gene$gene<-toupper(disease_gene$gene)
  mirna_gene_associations<-fread("Interlayers/miRTarBase_MTI.csv",stringsAsFactors = F,sep=";") ; mirna_gene_associations$miRNA<-tolower(mirna_gene_associations$miRNA) ; mirna_gene_associations$`Target Gene`<-toupper(mirna_gene_associations$`Target Gene`)
  ## Select only interactions between human miRNAs and genes ##
  mirna_gene_associations<-mirna_gene_associations[which(mirna_gene_associations$`Species (miRNA)`=="Homo sapiens")] ; mirna_gene_associations<-mirna_gene_associations[which(mirna_gene_associations$`Species (Target Gene)`=="Homo sapiens")]
  mirna_gene_associations<-mirna_gene_associations[,c(2,4)]
  colnames(mirna_gene_associations)<-c("mirna","gene")
  beggining<-c(length(unique(disease_mirna$icd10)),length(unique(disease_mirna$mirna)),length(disease_mirna$icd10),
               length(unique(disease_gene$icd10)),length(unique(disease_gene$gene)),length(disease_gene$icd10),
               length(mirna_gene_associations$mirna))
  ## Get the miRNAs common to both files
  mirnas_comunes<-intersect(disease_mirna$mirna,mirna_gene_associations$mirna)
  ## Only 432 out of 1,465 miRNAs in the disease network have an interaction with at least one gene ## 29.48%
  length(intersect(disease_mirna$mirna,mirna_gene_associations$mirna))/length(unique(disease_mirna$mirna))
  setkey(mirna_gene_associations,"miRNA")
  mirna_gene_associations<-mirna_gene_associations[mirnas_comunes]
  ## Get the genes common to both files
  genes_comunes<-intersect(disease_gene$gene,mirna_gene_associations$gene)
  ## Only 3,818 out of 7,360 genes in the disease network have an interaction with at least one miRNA ## 51.87%
  length(intersect(disease_gene$gene,mirna_gene_associations$gene))/length(unique(disease_gene$gene))
  setkey(mirna_gene_associations,"Target Gene")
  mirna_gene_associations<-mirna_gene_associations[genes_comunes]
  ## Get the diseases associated to genes and mirnas that have at least one association to each other (gene - miRNA association)
  commongenes<-unique(mirna_gene_associations$gene) ; commonmirnas<-unique(mirna_gene_associations$mirna)
  setkey(disease_mirna,"mirna") ; setkey(disease_gene,"gene")
  disease_gene<-disease_gene[commongenes] ; disease_mirna<-disease_mirna[commonmirnas]
  ## At the end we have interactions between 430 miRNAs and 3,818 genes ##
  ending<-c(length(unique(disease_mirna$icd10)),length(unique(disease_mirna$mirna)),length(disease_mirna$icd10),
            length(unique(disease_gene$icd10)),length(unique(disease_gene$gene)),length(disease_gene$icd10),
            length(mirna_gene_associations$mirna))
  association_summary<-cbind(beggining,ending)
  rownames(association_summary)<-c("ICD10(mirna)","miRNAs","ICD10-mirna associations","ICD10(gene)","Genes","ICD10-gene associations","miRNA-gene associations")
  association_summary<-cbind(association_summary,(association_summary[,2]/association_summary[,1])*100)
  colnames(association_summary)<-c("original","post_filtering","percentage_kept")
  write.table(association_summary,"Summaries/Summary_mirna_gene_associations.txt",quote=F,sep="\t")
  ## 430 miRNAs are connected to 319 ICD10 - 1,847 associations ##
  ## 3,818 genes are connected to 631 ICD10 - 12,796 associations ##
  ## Build the network connecting diseases in one layer to diseases in the other_one ##
  vmirna_gene_associations<-paste(mirna_gene_associations$mirna,mirna_gene_associations$gene,sep="__")
  genediseases<-unique(disease_gene$icd10) ; mirnadiseases<-unique(disease_mirna$icd10)
  mirna_gene_layer<-c()
  for(a in 1:length(genediseases)){
    # a<-1
    genes<-unique(disease_gene$gene[which(disease_gene$icd10==genediseases[a])])
    for(b in 1:length(mirnadiseases)){
      # b<-1
      mirnas<-unique(disease_mirna$mirna[which(disease_mirna$icd10==mirnadiseases[b])])
      mirna_gene_layer<-rbind(mirna_gene_layer,c(genediseases[a],mirnadiseases[b],length(intersect(apply(expand.grid(mirnas, genes), 1, paste, collapse="__"),vmirna_gene_associations))))
    }
    print(paste(round((a/length(genediseases))*100,2),"%",sep=""))
  }
  mirna_gene_layer<-cbind(mirna_gene_layer,as.numeric(mirna_gene_layer[,3])/max(as.numeric(mirna_gene_layer[,3])))
  colnames(mirna_gene_layer)<-c("ICD10_gene","ICD10_mirna","interactions","normalized") ## 201,289 combinations
  mirna_gene_interlayer<-mirna_gene_layer[which(as.numeric(mirna_gene_layer[,3])!=0),]
  mirna_gene_interlayer<-as.data.table(mirna_gene_interlayer)
  mirna_gene_interlayer$interactions<-as.numeric(mirna_gene_interlayer$interactions)
  mirna_gene_interlayer$normalized<-as.numeric(mirna_gene_interlayer$normalized)
  write.table(mirna_gene_interlayer,"Interlayers/Final/Genes_miRNA_ICD10_interactions.txt",row.names = F,sep="\t",quote=F)
  
  
  #### Connect genes to drug ####
  disease_gene<-fread("Genes_DisGeNET_20230113/ICD10_gene_associations.txt",stringsAsFactors = F,sep="\t")
  disease_gene$gene<-toupper(disease_gene$gene)
  disease_drug<-fread("Drugs_CTD_20220712/ICD10_drug_associations.txt",stringsAsFactors = F,sep="\t")
  disease_drug$drug<-tolower(disease_drug$drug)
  drug_gene_associations<-fread("Interlayers/CTD_chem_gene_ixns.csv",stringsAsFactors = F,sep=",")
  ## Select only interactions between drugs and human genes ##
  drug_gene_associations<-drug_gene_associations[which(drug_gene_associations$Organism=="Homo sapiens")]
  drug_gene_associations<-drug_gene_associations[,c(1,4,10,11)]
  colnames(drug_gene_associations)<-c("drug","gene","interaction","pmid")
  beggining<-c(length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
               length(unique(disease_gene$icd10)),length(unique(disease_gene$gene)),length(disease_gene$icd10),
               length(drug_gene_associations$drug))
  ## Get the drugs common to both files
  drugs_comunes<-intersect(disease_drug$drug,drug_gene_associations$drug)
  ## Only 1,120 out of 5,351 drugs in the disease network have an interaction with at least one gene ## 20.93%
  length(intersect(disease_drug$drug,drug_gene_associations$drug))/length(unique(disease_drug$drug))
  setkey(drug_gene_associations,"drug")
  drug_gene_associations<-drug_gene_associations[drugs_comunes]
  ## Get the genes common to both files
  genes_comunes<-intersect(disease_gene$gene,drug_gene_associations$gene)
  ## Only 6,449 out of 7,360 genes in the disease network have an interaction with at least one drug ## 87.62%
  length(intersect(disease_gene$gene,drug_gene_associations$gene))/length(unique(disease_gene$gene))
  setkey(drug_gene_associations,"Target Gene")
  drug_gene_associations<-drug_gene_associations[genes_comunes]
  ## Get the diseases associated to genes and drugs that have at least one association to each other (gene - drug association)
  commongenes<-unique(drug_gene_associations$gene) ; commondrugs<-unique(drug_gene_associations$drug)
  setkey(disease_drug,"drug") ; setkey(disease_gene,"gene")
  disease_gene<-disease_gene[commongenes] ; disease_drug<-disease_drug[commondrugs]
  ## At the end we have interactions between 430 drugs and 3,818 genes ##
  ending<-c(length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
            length(unique(disease_gene$icd10)),length(unique(disease_gene$gene)),length(disease_gene$icd10),
            length(drug_gene_associations$drug))
  association_summary<-cbind(beggining,ending)
  rownames(association_summary)<-c("ICD10(drug)","drugs","ICD10-drug associations","ICD10(gene)","Genes","ICD10-gene associations","drug-gene associations")
  association_summary<-cbind(association_summary,(association_summary[,2]/association_summary[,1])*100)
  colnames(association_summary)<-c("original","post_filtering","percentage_kept")
  write.table(association_summary,"Summaries/Summary_drug_gene_associations.txt",quote=F,sep="\t")
  ## 430 drugs are connected to 319 ICD10 - 1,847 associations ##
  ## 3,818 genes are connected to 631 ICD10 - 12,796 associations ##
  ## Build the network connecting diseases in one layer to diseases in the other_one ##
  tt<-sort(table(drug_gene_associations$drug),decreasing = T)
  vdrug_gene_associations<-unique(paste(drug_gene_associations$drug,drug_gene_associations$gene,sep="__"))
  genediseases<-unique(disease_gene$icd10) ; drugdiseases<-unique(disease_drug$icd10)
  for(a in 1:length(genediseases)){
    # a<-1
    genes<-unique(disease_gene$gene[which(disease_gene$icd10==genediseases[a])])
    drug_gene_layer<-c()
    for(b in 1:length(drugdiseases)){
      # b<-1
      drugs<-unique(disease_drug$drug[which(disease_drug$icd10==drugdiseases[b])])
      drug_gene_layer<-rbind(drug_gene_layer,c(genediseases[a],drugdiseases[b],length(intersect(apply(expand.grid(drugs, genes), 1, paste, collapse="__"),vdrug_gene_associations))))
    }
    write.table(drug_gene_layer,"Interlayers/Final/Genes_drug_ICD10_interactions.txt",quote=F,sep="\t",row.names=F,col.names = F,append = T)
    print(paste(round((a/length(genediseases))*100,2),"%",sep=""))
  }
  drug_gene_layer<-fread("Interlayers/Final/Genes_drug_ICD10_interactions.txt",stringsAsFactors = F,sep="\t",header=F)
  drug_gene_layer<-cbind(drug_gene_layer,as.numeric(drug_gene_layer$V3)/max(as.numeric(drug_gene_layer$V3)))
  colnames(drug_gene_layer)<-c("ICD10_gene","ICD10_drug","interactions","normalized") ## 273,309 combinations
  drug_gene_interlayer<-drug_gene_layer[which(as.numeric(drug_gene_layer$interactions)!=0),]
  drug_gene_interlayer$interactions<-as.numeric(drug_gene_interlayer$interactions)
  write.table(drug_gene_interlayer,"Interlayers/Final/Genes_drug_ICD10_interactions.txt",row.names = F,sep="\t",quote=F)
  #### Connect mirnas to drugs ####
  disease_mirna<-fread("miRNA_HMDD_20240523/ICD10_mirna_associations.txt",stringsAsFactors = F,sep="\t")
  disease_mirna$mirna<-toupper(disease_mirna$mirna)
  disease_drug<-fread("Drugs_CTD_20220712/ICD10_drug_associations.txt",stringsAsFactors = F,sep="\t")
  disease_drug$drug<-tolower(disease_drug$drug)
  ## CTD ##
  ## @@@ ##
  drug_mirna_associations<-fread("Interlayers/CTD_chem_gene_ixns.csv",stringsAsFactors = F,sep=",")
  ## Select only interactions between drugs and human mirnas ##
  drug_mirna_associations<-drug_mirna_associations[which(drug_mirna_associations$Organism=="Homo sapiens")]
  drug_mirna_associations<-drug_mirna_associations[,c(1,4,10,11)]
  colnames(drug_mirna_associations)<-c("drug","mirna","interaction","pmid")
  ## Select the drug-mirna associations ##
  drug_mirna_ctd_associations<-drug_mirna_associations[grep("mir",drug_mirna_associations$mirna,ignore.case = T)]
  drug_mirna_ctd_associations$mirna<-paste("HSA-MIR-",gsub("MIR","",drug_mirna_ctd_associations$mirna),sep="")
  drug_mirna_ctd_associations$drug<-tolower(drug_mirna_ctd_associations$drug)
  ## SM2miR3 ##
  ## @ @ @ @ ##
  drug_mirna_associations<-fread("Interlayers/SM2miR3.csv",stringsAsFactors = F,sep=";")
  ## Select only interactions between drugs and human mirnas ##
  drug_mirna_associations<-drug_mirna_associations[which(drug_mirna_associations$Species=="Homo sapiens")]
  drug_mirna_associations<-drug_mirna_associations[,c(1,3,4:6,9:10,14)]
  colnames(drug_mirna_associations)<-c("mirna","drug","fdr","db","cid","disease","pmid","expression")
  ## Select the drug-mirna associations ##
  drug_mirna_sm2mir3_associations<-drug_mirna_associations
  drug_mirna_sm2mir3_associations$mirna<-paste("HSA-MIR",gsub("miR","",drug_mirna_sm2mir3_associations$mirna),sep="")
  drug_mirna_sm2mir3_associations$drug<-tolower(drug_mirna_sm2mir3_associations$drug)
  ## Put both databases together ##
  ## @@ @@ @@ @@ @ @ @@ @@ @@ @@ ##
  drug_mirna_associations<-rbind(drug_mirna_sm2mir3_associations[,c(1,2)],drug_mirna_ctd_associations[,c(2,1)])
  if(length(which(duplicated(drug_mirna_associations)))>0){drug_mirna_associations<-drug_mirna_associations[-which(duplicated(drug_mirna_associations)),]}
  beggining<-c(length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
               length(unique(disease_mirna$icd10)),length(unique(disease_mirna$mirna)),length(disease_mirna$icd10),
               length(drug_mirna_associations$drug))
  
  ## Get the drugs common to both files
  drugs_comunes<-intersect(disease_drug$drug,drug_mirna_associations$drug)
  ## Only 255 out of 5,351 drugs in the disease network have an interaction with at least one mirna ## 4.76%
  length(intersect(disease_drug$drug,drug_mirna_associations$drug))/length(unique(disease_drug$drug))
  setkey(drug_mirna_associations,"drug")
  drug_mirna_associations<-drug_mirna_associations[drugs_comunes]
  ## Get the mirnas common to both files
  mirnas_comunes<-intersect(disease_mirna$mirna,drug_mirna_associations$mirna)
  ## Only 748 out of 1,465 mirnas in the disease network have an interaction with at least one drug ## 51.05%
  length(intersect(disease_mirna$mirna,drug_mirna_associations$mirna))/length(unique(disease_mirna$mirna))
  setkey(drug_mirna_associations,"mirna")
  drug_mirna_associations<-drug_mirna_associations[mirnas_comunes]
  ## Get the diseases associated to mirnas and drugs that have at least one association to each other (mirna - drug association)
  commonmirnas<-unique(drug_mirna_associations$mirna) ; commondrugs<-unique(drug_mirna_associations$drug)
  setkey(disease_drug,"drug") ; setkey(disease_mirna,"mirna")
  disease_mirna<-disease_mirna[commonmirnas] ; disease_drug<-disease_drug[commondrugs]
  ## At the end we have interactions between 430 drugs and 3,818 mirnas ##
  ending<-c(length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
            length(unique(disease_mirna$icd10)),length(unique(disease_mirna$mirna)),length(disease_mirna$icd10),
            length(drug_mirna_associations$drug))
  association_summary<-cbind(beggining,ending)
  rownames(association_summary)<-c("ICD10(drug)","drugs","ICD10-drug associations","ICD10(mirna)","mirnas","ICD10-mirna associations","drug-mirna associations")
  association_summary<-cbind(association_summary,(association_summary[,2]/association_summary[,1])*100)
  colnames(association_summary)<-c("original","post_filtering","percentage_kept")
  write.table(association_summary,"Summaries/Summary_drug_mirna_associations.txt",quote=F,sep="\t")
  ## 211 drugs are connected to 491 ICD10 - 3,087 associations ##
  ## 748 mirnas are connected to 589 ICD10 - 12,147 associations ##
  ## Build the network connecting diseases in one layer to diseases in the other_one ##
  tt<-sort(table(drug_mirna_associations$drug),decreasing = T)
  vdrug_mirna_associations<-unique(paste(drug_mirna_associations$drug,drug_mirna_associations$mirna,sep="__"))
  mirnadiseases<-unique(disease_mirna$icd10) ; drugdiseases<-unique(disease_drug$icd10)
  for(a in 1:length(mirnadiseases)){
    # a<-1
    mirnas<-unique(disease_mirna$mirna[which(disease_mirna$icd10==mirnadiseases[a])])
    drug_mirna_layer<-c()
    for(b in 1:length(drugdiseases)){
      # b<-1
      drugs<-unique(disease_drug$drug[which(disease_drug$icd10==drugdiseases[b])])
      drug_mirna_layer<-rbind(drug_mirna_layer,c(mirnadiseases[a],drugdiseases[b],length(intersect(apply(expand.grid(drugs, mirnas), 1, paste, collapse="__"),vdrug_mirna_associations))))
    }
    write.table(drug_mirna_layer,"Interlayers/Final/mirnas_drug_ICD10_interactions.txt",quote=F,sep="\t",row.names=F,col.names = F,append = T)
    print(paste(round((a/length(mirnadiseases))*100,2),"%",sep=""))
  }
  drug_mirna_layer<-fread("Interlayers/Final/mirnas_drug_ICD10_interactions.txt",stringsAsFactors = F,sep="\t",header=F)
  drug_mirna_layer<-cbind(drug_mirna_layer,as.numeric(drug_mirna_layer$V3)/max(as.numeric(drug_mirna_layer$V3)))
  colnames(drug_mirna_layer)<-c("ICD10_mirna","ICD10_drug","interactions","normalized") ## 289,199 combinations
  drug_mirna_interlayer<-drug_mirna_layer[which(as.numeric(drug_mirna_layer$interactions)!=0),]
  drug_mirna_interlayer$interactions<-as.numeric(drug_mirna_interlayer$interactions) ## 191,300 interactions
  write.table(drug_mirna_interlayer,"Interlayers/Final/miRNAs_drug_ICD10_interactions.txt",row.names = F,sep="\t",quote=F)
}


if(args[1]=="presentation"){
  ## Disease - variable associations ##
  tt_mirna<-fread("Associations/ICD10_mirna_associations.txt",stringsAsFactors = F,sep="\t")
  tt_genes<-fread("Associations/ICD10_gene_associations.txt",stringsAsFactors = F,sep="\t")
  tt_drugs<-fread("Associations/ICD10_drug_associations.txt",stringsAsFactors = F,sep="\t")
  tt_symptoms<-fread("Associations/ICD10_symptom_associations.txt",stringsAsFactors = F,sep="\t")
  ## Intra-layer ##
  mirnalayer<-fread("IntraLayers/Final/miRNA_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
  druglayer<-fread("IntraLayers/Final/Drugs_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
  genelayer<-fread("IntraLayers/Final/Genes_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
  ## Inter-layer ##
  gene_drug_layer<-fread("Interlayers/Final/Genes_drug_ICD10_interactions.txt",stringsAsFactors = F,sep="\t")
  gene_mirna_layer<-fread("Interlayers/Final/Genes_miRNA_ICD10_interactions.txt",stringsAsFactors = F,sep="\t")
  mirna_drug_layer<-fread("Interlayers/Final/mirnas_drug_ICD10_interactions.txt",stringsAsFactors = F,sep="\t")
  
  intralayerconnections<-rbind(c(length(unique(c(mirnalayer$disease1,mirnalayer$disease2))),length(mirnalayer$disease1),length(which(mirnalayer$fdr<=0.05))),
                               c(length(unique(c(druglayer$disease1,druglayer$disease2))),length(druglayer$disease1),length(which(druglayer$fdr<=0.05))),
                               c(length(unique(c(genelayer$disease1,genelayer$disease2))),length(genelayer$disease1),length(which(genelayer$fdr<=0.05))))
  colnames(intralayerconnections)<-c("diseases","edges","significant_edges")
  rownames(intralayerconnections)<-c("mirna","drug","gene")
  
  
  interlayerconnections<-rbind(c("mirna","drug",length(unique(mirna_drug_layer$ICD10_mirna)),length(unique(mirna_drug_layer$ICD10_drug)),length(mirna_drug_layer$normalized)),
                               c("mirna","gene",length(unique(gene_mirna_layer$ICD10_mirna)),length(unique(gene_mirna_layer$ICD10_gene)),length(gene_mirna_layer$normalized)),
                               c("drug","gene",length(unique(gene_drug_layer$ICD10_drug)),length(unique(gene_drug_layer$ICD10_gene)),length(gene_drug_layer$normalized)))
  colnames(interlayerconnections)<-c("layer1","layer2","disease1","diseases2","edges")
  interlayerconnections<-as.data.frame(interlayerconnections)
  
  
  if("Plots"%in%list.files()==F){dir.create("Plots")}
  pdf(file="Plots/Edge_values_mirna_gene_drug.pdf")
    par(mfrow=c(2,3))
    ## Intra-layer edges ##
    hist(mirnalayer$cosine_similarity,main="miRNA layer",xlab="cosine similarity",breaks = 100)
    hist(druglayer$cosine_similarity,main="drug layer",xlab="cosine similarity",breaks = 100)
    hist(genelayer$cosine_similarity,main="gene layer",xlab="cosine similarity",breaks = 100)
    ## Inter-layer edges ##
    hist(mirna_drug_layer$normalized,main="miRNA-drug\ninter-layer",xlab="IFS",breaks = 100)
    hist(gene_mirna_layer$normalized,main="miRNA-gene\ninter-layer",xlab="IFS",breaks = 100)
    hist(gene_drug_layer$normalized,main="gene-drug\ninter-layer",xlab="IFS",breaks = 100)
  dev.off()
  
  
  #### Look for specific diseases ####
  ## mirna ##
  mirna<-fread("IntraLayers/Final/miRNA_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
  ## genes ##
  gene<-fread("IntraLayers/Final/Genes_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
  ## drugs ##
  drug<-fread("IntraLayers/Final/Drugs_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
  ## symptoms significantly associated to diseases ##
  symptom<-fread("IntraLayers/Final/Symptoms_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
  ## transcriptomics ##
  transcriptomics<-fread("IntraLayers/Final/Transcriptomics_ICD10_CosineSimilarities_fdrs.txt",stringsAsFactors = F,sep="\t")
  ## epidemiology ##
  epidemiology<-fread("IntraLayers/Final/Epidemiology_ICD10_RelativeRisk.txt",stringsAsFactors = F,sep="\t")
  ## metabolome ##
  metabolome<-fread("Metabolome_HMDB_20220712/ICD10_CosineSimilarities.txt",stringsAsFactors = F,sep="\t")
  colnames(metabolome)<-c("disease1","disease2","cosine_similarity")
  metabolomesorted<-t(apply(metabolome[,1:2], 1, sort))
  metabolome$disease1<-metabolomesorted[,1] ; metabolome$disease2<-metabolomesorted[,2]
  metabolome<-metabolome[which(metabolome$cosine_similarity!=0)]
  ## microbiome ##
  microbiome<-fread("microbiome_disbiome_20220714/ICD10_CosineSimilarities.txt",stringsAsFactors = F,sep="\t")
  colnames(microbiome)<-c("disease1","disease2","cosine_similarity")
  microbiomesorted<-t(apply(microbiome[,1:2], 1, sort))
  microbiome$disease1<-microbiomesorted[,1] ; microbiome$disease2<-microbiomesorted[,2]
  microbiome<-microbiome[which(microbiome$cosine_similarity!=0)]
  
  ## List with all the sources ##
  lnet<-list("mirna"=mirna,"gene"=gene,"drug"=drug,"symptom"=symptom,"metabolome"=metabolome,"microbiome"=microbiome,"transcriptomics"=transcriptomics)
  
  
  #### Associations ####
  asomirna<-fread("miRNA_HMDD_20240523/ICD10_mirna_associations.txt") ; colnames(asomirna)[2]<-"var"
  asogene<-fread("Genes_DisGeNET_20230113/ICD10_gene_associations.txt") ; colnames(asogene)[2]<-"var"
  asodrug<-fread("Drugs_CTD_20220712/ICD10_drug_associations.txt") ; colnames(asodrug)[2]<-"var"
  asosymp<-fread("Symptoms_Zhou2014_20220712/ICD10_symptom_significant_associations.txt") ; colnames(asosymp)[2]<-"var"
  asomicr<-fread("microbiome_disbiome_20220714/ICD10_microbiome_associations.txt") ; colnames(asomicr)[2]<-"var"
  asometa<-fread("Metabolome_HMDB_20220712/ICD10_metabolite_associations.txt") ; colnames(asometa)[2]<-"var"
  
  ## List with all the associations ##
  lasso<-list("mirna"=asomirna,"gene"=asogene,"drug"=asodrug,"symptom"=asosymp,"metabolome"=asometa,"microbiome"=asomicr)
  
  
  
  #### Function to get the involved variables ####
  retrieveinteractions<-function(d1,d2,lnet,lasso){
    valores<-c() ; fuentes<-c()
    for(a in 1:length(names(lnet))){
      # a<-1
      tt<-lnet[[a]][intersect(which(lnet[[a]]$disease1==d1),which(lnet[[a]]$disease2==d2))]
      if(dim(tt)[1]>0){
        valores<-rbind(valores,tt[,1:3])
        fuentes<-c(fuentes,names(lnet)[a])
      }
    }
    if(length(fuentes)>0){
      valores<-cbind(fuentes,valores)
      newval<-c()
      for(a in 1:length(fuentes)){
        if(fuentes[a]!="transcriptomics"){
          # a<-1
          ad1<-lasso[[fuentes[a]]][which(lasso[[fuentes[a]]]$icd10==d1)]$var
          ad2<-lasso[[fuentes[a]]][which(lasso[[fuentes[a]]]$icd10==d2)]$var
          ad12<-intersect(ad1,ad2)
          newval<-rbind(newval,cbind(valores[a],ad12))
        }
      }
    }
    resultados<-list("scores"=valores,"variables"=newval)
    return(resultados)
  }
  net<-retrieveinteractions("E11","G40",lnet,lasso)
  
  
  ## Network mode ##
  netmod<-c()
  for(a in 1:length(net$variables$fuentes)){
    # a<-1
    netmod<-rbind(netmod,cbind(net$variables[a]$fuentes,net$variables[a]$disease1,net$variables[a]$ad12),cbind(net$variables[a]$fuentes,net$variables[a]$disease2,net$variables[a]$ad12))
  }
  type<-netmod[,c(1,3)]
  if(length(which(duplicated(type)))>0){type<-type[-which(duplicated(type)),]}
  type<-rbind(type,c("disease","E11"),c("disease","G40"))
  colnames(type)<-c("type","node")
  write.table(type,"test_node.txt",quote=F,sep="\t",row.names=F)
  netmod<-rbind(netmod,c("comorbidity","E11","G40"))
  colnames(netmod)<-c("layer","d1","d2")
  write.table(netmod,"test_edges.txt",quote=F,sep="\t",row.names=F)
}

metaasso<-fread("Metabolome_HMDB_20220712/ICD10_metabolite_associations.txt")



























