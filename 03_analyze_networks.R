#### With this code we are going to generate the multiomic layers (disease-disease similarities) for the BSC-IGC collaboration ####
## Load the needed libraries ##
library("data.table")
library(VennDiagram)
# library(RColorBrewer)
library("igraph")
library("UpSetR")
# library("dendextend")
library(lsa)
library("geometry")
library("ggplot2")
library("plyr")
library("rgexf")
# library("gplots")
# library(GEOquery)
# library(metapro)

args = commandArgs(trailingOnly=TRUE)

## Select the final networks to be represented 
if(args[1]=="Final_networks"){
  #### Intralayers ####
  ##  @@ @@ @ @@ @@  ##
  summary<-c()
  ## drug ##
  red<-fread("Similarities/Similarities_with_fdr/ICD10_drug_fs_weightedjaccardindex.txt",stringsAsFactors = F,sep="\t")
  diseases<-unique(c(red[[1]],red[[2]]))
  net<-red[unique(c(which(red$fdr_down<=0.05),which(red$fdr_up<=0.05)))]
  if(length(which(net[[3]]==0))>0){net<-net[-which(net[[3]]==0),]}
  if(length(which(net$fdr_down<=0.05))>0){net[which(net$fdr_down<=0.05)][[3]]<--1*net[which(net$fdr_down<=0.05)][[3]]}
  net<-net[,1:3] ; colnames(net)<-c("icd10_1","icd10_2","weight")
  write.table(net,"HDMLN/Intralayers/drug.txt",quote=F,sep="\t",row.names=F)
  summary<-rbind(summary,c(length(which(red$fdr_down<=0.05)),length(which(red$fdr_up<=0.05)),length(unique(c(red[[1]],red[[2]]))),length(red[[1]])))
  
  ## gene ##
  red<-fread("Similarities/Similarities_with_fdr/ICD10_gene_fs_weightedjaccardindex.txt",stringsAsFactors = F,sep="\t")
  diseases<-c(diseases,unique(c(red[[1]],red[[2]])))
  net<-red[unique(c(which(red$fdr_down<=0.05),which(red$fdr_up<=0.05)))]
  if(length(which(net[[3]]==0))>0){net<-net[-which(net[[3]]==0),]}
  if(length(which(net$fdr_down<=0.05))>0){net[which(net$fdr_down<=0.05)][[3]]<--1*net[which(net$fdr_down<=0.05)][[3]]}
  net<-net[,1:3] ; colnames(net)<-c("icd10_1","icd10_2","weight")
  write.table(net,"HDMLN/Intralayers/gene.txt",quote=F,sep="\t",row.names=F)
  summary<-rbind(summary,c(length(which(red$fdr_down<=0.05)),length(which(red$fdr_up<=0.05)),length(unique(c(red[[1]],red[[2]]))),length(red[[1]])))
  
  ## metabolome ##
  red<-fread("Similarities/Similarities_with_fdr/ICD10_metabolome_fs_cosinesimilarities.txt",stringsAsFactors = F,sep="\t")
  diseases<-c(diseases,unique(c(red[[1]],red[[2]])))
  net<-red[unique(c(which(red$fdr_down<=0.05),which(red$fdr_up<=0.05)))]
  if(length(which(net[[3]]==0))>0){net<-net[-which(net[[3]]==0),]}
  if(length(which(net$fdr_down<=0.05))>0){net[which(net$fdr_down<=0.05)][[3]]<--1*net[which(net$fdr_down<=0.05)][[3]]}
  net<-net[,1:3] ; colnames(net)<-c("icd10_1","icd10_2","weight")
  write.table(net,"HDMLN/Intralayers/metabolome.txt",quote=F,sep="\t",row.names=F)
  summary<-rbind(summary,c(length(which(red$fdr_down<=0.05)),length(which(red$fdr_up<=0.05)),length(unique(c(red[[1]],red[[2]]))),length(red[[1]])))
  
  ## methylome ##
  red<-fread("Similarities/Similarities_with_fdr/ICD10_methylome_ifs_cosinesimilarities.txt",stringsAsFactors = F,sep="\t")
  diseases<-c(diseases,unique(c(red[[1]],red[[2]])))
  net<-red[unique(c(which(red$fdr_down<=0.05),which(red$fdr_up<=0.05)))]
  if(length(which(net[[3]]==0))>0){net<-net[-which(net[[3]]==0),]}
  if(length(which(net$fdr_down<=0.05))>0){net[which(net$fdr_down<=0.05)][[3]]<--1*net[which(net$fdr_down<=0.05)][[3]]}
  net<-net[,1:3] ; colnames(net)<-c("icd10_1","icd10_2","weight")
  write.table(net,"HDMLN/Intralayers/methylome.txt",quote=F,sep="\t",row.names=F)
  summary<-rbind(summary,c(length(which(red$fdr_down<=0.05)),length(which(red$fdr_up<=0.05)),length(unique(c(red[[1]],red[[2]]))),length(red[[1]])))
  
  ## microbiome ##
  red<-fread("Similarities/Similarities_with_fdr/ICD10_microbiome_fsd_cosinesimilarities.txt",stringsAsFactors = F,sep="\t")
  diseases<-c(diseases,unique(c(red[[1]],red[[2]])))
  net<-red[unique(c(which(red$fdr_down<=0.05),which(red$fdr_up<=0.05)))]
  if(length(which(net[[3]]==0))>0){net<-net[-which(net[[3]]==0),]}
  if(length(which(net$fdr_down<=0.05))>0){net[which(net$fdr_down<=0.05)][[3]]<--1*net[which(net$fdr_down<=0.05)][[3]]}
  net<-net[,1:3] ; colnames(net)<-c("icd10_1","icd10_2","weight")
  write.table(net,"HDMLN/Intralayers/microbiome.txt",quote=F,sep="\t",row.names=F)
  summary<-rbind(summary,c(length(which(red$fdr_down<=0.05)),length(which(red$fdr_up<=0.05)),length(unique(c(red[[1]],red[[2]]))),length(red[[1]])))
  
  ## mirna ##
  red<-fread("Similarities/Similarities_with_fdr/ICD10_miRNA_fs_cosinesimilarities.txt",stringsAsFactors = F,sep="\t")
  diseases<-c(diseases,unique(c(red[[1]],red[[2]])))
  net<-red[unique(c(which(red$fdr_down<=0.05),which(red$fdr_up<=0.05)))]
  if(length(which(net[[3]]==0))>0){net<-net[-which(net[[3]]==0),]}
  if(length(which(net$fdr_down<=0.05))>0){net[which(net$fdr_down<=0.05)][[3]]<--1*net[which(net$fdr_down<=0.05)][[3]]}
  net<-net[,1:3] ; colnames(net)<-c("icd10_1","icd10_2","weight")
  write.table(net,"HDMLN/Intralayers/miRNA.txt",quote=F,sep="\t",row.names=F)
  summary<-rbind(summary,c(length(which(red$fdr_down<=0.05)),length(which(red$fdr_up<=0.05)),length(unique(c(red[[1]],red[[2]]))),length(red[[1]])))
  
  ## symptom ##
  red<-fread("Similarities/Similarities_with_fdr/ICD10_symptom_tfidf_cosinesimilarities.txt",stringsAsFactors = F,sep="\t")
  diseases<-c(diseases,unique(c(red[[1]],red[[2]])))
  net<-red[unique(c(which(red$fdr_down<=0.05),which(red$fdr_up<=0.05)))]
  if(length(which(net[[3]]==0))>0){net<-net[-which(net[[3]]==0),]}
  if(length(which(net$fdr_down<=0.05))>0){net[which(net$fdr_down<=0.05)][[3]]<--1*net[which(net$fdr_down<=0.05)][[3]]}
  net<-net[,1:3] ; colnames(net)<-c("icd10_1","icd10_2","weight")
  write.table(net,"HDMLN/Intralayers/symptom.txt",quote=F,sep="\t",row.names=F)
  summary<-rbind(summary,c(length(which(red$fdr_down<=0.05)),length(which(red$fdr_up<=0.05)),length(unique(c(red[[1]],red[[2]]))),length(red[[1]])))
  
  ## transcriptome ##
  red<-fread("Similarities/Similarities_with_fdr/ICD10_transcriptome_fc_cosinesimilarities.txt",stringsAsFactors = F,sep="\t")
  diseases<-c(diseases,unique(c(red[[1]],red[[2]])))
  net<-red[which(red$fdr<=0.05)]
  net<-net[,1:3] ; colnames(net)<-c("icd10_1","icd10_2","weight")
  write.table(net,"HDMLN/Intralayers/transcriptome.txt",quote=F,sep="\t",row.names=F)
  summary<-rbind(summary,c(length(which(net$weight<0)),length(which(net$weight>0)),length(unique(c(red[[1]],red[[2]]))),length(red[[1]])))
  ## add row and column names 
  colnames(summary)<-c("negative","positive","diseases","interactions")
  rownames(summary)<-c("drug","gene","metabolome","methylome","microbiome","mirna","symptom","transcriptome")
  ## get index for row names
  indfila<-dim(summary)[1]+1
  
  #### Interlayers ####
  ##  @@ @@ @ @@ @@  ##
  ficheros<-list.files("Similarities/Similarities_Interlayer_with_fdr/")
  nombres<-c()
  for(a in 1:length(ficheros)){
    # a<-1
    red<-fread(paste("Similarities/Similarities_Interlayer_with_fdr/",ficheros[a],sep=""),quote=F,sep="\t")
    net<-red[unique(c(which(red$fdr_down<=0.05),which(red$fdr_up<=0.05)))]
    if(length(which(net[[3]]==0))>0){net<-net[-which(net[[3]]==0),]}
    if(length(which(net$fdr_down<=0.05))>0){net[which(net$fdr_down<=0.05)][[3]]<--1*net[which(net$fdr_down<=0.05)][[3]]}
    net<-net[,1:3] ; colnames(net)[3]<-c("weight")
    write.table(net,paste("HDMLN/Interlayers/",gsub("_associations","",gsub("ICD10_","",ficheros[a])),sep=""),quote=F,sep="\t",row.names=F)
    nombres<-c(nombres,paste(gsub("icd10_","",colnames(net)[1:2]),collapse = " - "))
    summary<-rbind(summary,c(length(which(net$weight<0)),length(which(net$weight>0)),paste(length(unique(red[[1]])),"-",length(unique(red[[2]]))),length(red[[1]])))
  }
  rownames(summary)[9:dim(summary)[1]]<-nombres
  write.table(summary,"Summaries/HDMLN_interactions.txt",quote=F,sep="\t")
}


## Calculate the overlap of the different layers with epidemiology
if(args[1]=="calculate_overlaps"){
  # comparison<-gsub("_.+","",gsub(".+ICD10_","","Similarities/Similarities_with_fdr/ICD10_drug_fs_weightedjaccardindex.txt"))
  # thenet<-fread("Similarities/Similarities_with_fdr/ICD10_drug_fs_weightedjaccardindex.txt",stringsAsFactors = F,sep="\t")
  # numberrandom<-10000
  # agerange<-"all"
  # epinet<-"Denmark"
  ## The function ##
  calculateoverlaps<-function(comparison,thenet,numberrandom,agerange,epinet){
    set.seed(1)
    ## Get the icds ##
    icds2<-unique(c(thenet[[1]],thenet[[2]]))
    if(agerange=="all"){
      if(epinet=="Denmark"){
        epidemiology2<-fread("Epidemiology/ICD10_epidemiology_denmark.txt",stringsAsFactors = F,sep="\t")
        ## Diseases in both networks ##
        icds<-intersect(icds2,unique(c(epidemiology2$disease1,epidemiology2$disease2)))
        ## Select the comorbidities between shared diseases ##
        ind1<-c() ; ind2<-c()
        for(a in 1:length(icds)){
          ind1<-c(ind1,which(epidemiology2$disease1==icds[a]))
          ind2<-c(ind2,which(epidemiology2$disease2==icds[a]))
        }
        ## Select comorbidities between shared diseases ##
        epidemiology<-epidemiology2[intersect(ind1,ind2),]
        ## Select the co-occurrences that are significant ##
        aresignificant<-which(as.numeric(epidemiology$RR)>1.01)
        epi<-epidemiology[aresignificant,] ### Changed from 1 to 1.01
        if(length(which(duplicated(epi[,1:2])))>0){epi<-epi[-which(duplicated(epi[,1:2])),]}
        ## Which are the diseases that should be included as isolated nodes in the epidemiology? ##
        nodesthatshouldbethere_epidemiology<-setdiff(icds,unique(c(epi$disease1,epi$disease2)))
      }
      if(epinet=="Catalonia"){
        epidemiology2<-fread("Epidemiology/Catalonia/N0_5_Adjusted.txt",stringsAsFactors = F,sep="\t")
        colnames(epidemiology2)[1:2]<-c("disease1","disease2")
        ## Diseases in both networks ##
        icds<-intersect(icds2,unique(c(epidemiology2$disease1,epidemiology2$disease2)))
        ## Select the comorbidities between shared diseases ##
        ind1<-c() ; ind2<-c()
        for(a in 1:length(icds)){
          ind1<-c(ind1,which(epidemiology2$disease1==icds[a]))
          ind2<-c(ind2,which(epidemiology2$disease2==icds[a]))
        }
        ## Select comorbidities between shared diseases ##
        epidemiology<-epidemiology2[intersect(ind1,ind2),]
        ## Select the co-occurrences that are significant ##
        epi<-epidemiology[intersect(intersect(which(as.numeric(epidemiology$OR)>1),which(as.numeric(epidemiology$FDR)<=0.05)),which(as.numeric(epidemiology$CasesWith)>=100)),]
        sortednet<-t(apply(epi[,1:2], 1, sort))
        epi$disease1<-sortednet[,1] ; epi$disease2<-sortednet[,2]
        if(length(which(duplicated(epi[,1:2])))>0){epi<-epi[-which(duplicated(epi[,1:2])),]}
        ## Which are the diseases that should be included as isolated nodes in the epidemiology? ##
        nodesthatshouldbethere_epidemiology<-setdiff(icds,unique(c(epi$disease1,epi$disease2)))
      }
      ## All the interactions or just the significant?? FDR ##
      if(comparison=="transcriptome"){thenet<-thenet[intersect(which(as.numeric(thenet[[3]])>0),which(as.numeric(thenet$fdr)<=0.05))]}
      if(comparison!="transcriptome"){thenet<-thenet[intersect(which(as.numeric(thenet[[3]])>0),which(as.numeric(thenet$fdr_up)<=0.05))]}
      thenet<-thenet[,1:2]
      ## Select common diseases ##
      colnames(thenet)<-c("disease1","disease2")
      ## Select the similarities between shared diseases ##
      ind1<-c() ; ind2<-c()
      for(a in 1:length(icds)){
        ind1<-c(ind1,which(thenet$disease1==icds[a]))
        ind2<-c(ind2,which(thenet$disease2==icds[a]))
      }
      thenet<-thenet[intersect(ind1,ind2)]

      ## Which are the diseases that should be included as isolated nodes in the omics layer? ##
      nodesthatshouldbethere_omics<-setdiff(icds,unique(c(thenet$disease1,thenet$disease2)))
    }
    if(agerange!="all"){
      gsub("N0_5_Adjusted_","",gsub(".txt","",list.files("Epidemiology/Catalonia/")[grep("-",list.files("Epidemiology/Catalonia/"))]))
      if(epinet=="Catalonia"){
        epidemiology2<-fread(paste("Epidemiology/Catalonia/N0_5_Adjusted_",agerange,".txt",sep=""),stringsAsFactors = F,sep="\t")
        colnames(epidemiology2)[1:2]<-c("disease1","disease2")
        ## Diseases in both networks ##
        icds<-intersect(icds2,unique(c(epidemiology2$disease1,epidemiology2$disease2)))
        ## Select the comorbidities between shared diseases ##
        ind1<-c() ; ind2<-c()
        for(a in 1:length(icds)){
          ind1<-c(ind1,which(epidemiology2$disease1==icds[a]))
          ind2<-c(ind2,which(epidemiology2$disease2==icds[a]))
        }
        ## Select comorbidities between shared diseases ##
        epidemiology<-epidemiology2[intersect(ind1,ind2),]
        ## Select the co-occurrences that are significant ##
        epi<-epidemiology[intersect(intersect(which(as.numeric(epidemiology$OR)>1),which(as.numeric(epidemiology$FDR)<=0.05)),which(as.numeric(epidemiology$CasesWith)>=5)),]
        sortednet<-t(apply(epi[,1:2], 1, sort))
        epi$disease1<-sortednet[,1] ; epi$disease2<-sortednet[,2]
        if(length(which(duplicated(epi[,1:2])))>0){epi<-epi[-which(duplicated(epi[,1:2])),]}
        ## Which are the diseases that should be included as isolated nodes in the epidemiology? ##
        nodesthatshouldbethere_epidemiology<-setdiff(icds,unique(c(epi$disease1,epi$disease2)))
      }
      ## All the interactions or just the significant?? FDR ##
      if(comparison=="transcriptome"){thenet<-thenet[intersect(which(as.numeric(thenet[[3]])>0),which(as.numeric(thenet$fdr)<=0.05))]}
      if(comparison!="transcriptome"){thenet<-thenet[intersect(which(as.numeric(thenet[[3]])>0),which(as.numeric(thenet$fdr_up)<=0.05))]}
      thenet<-thenet[,1:2]
      ## Select common diseases ##
      colnames(thenet)<-c("disease1","disease2")
      ## Select the similarities between shared diseases ##
      ind1<-c() ; ind2<-c()
      for(a in 1:length(icds)){
        ind1<-c(ind1,which(thenet$disease1==icds[a]))
        ind2<-c(ind2,which(thenet$disease2==icds[a]))
      }
      thenet<-thenet[intersect(ind1,ind2)]
      
      ## Which are the diseases that should be included as isolated nodes in the omics layer? ##
      nodesthatshouldbethere_omics<-setdiff(icds,unique(c(thenet$disease1,thenet$disease2)))
    }
    ## Calculate overlaps ##
    ## @@ @@ @@  @@ @@ @@ ##
    print(paste("Calculating significance of the overlap in ",comparison,":",sep=""))
    ## Transform the table into an igraph object - Epidemiology ##
    if(length(nodesthatshouldbethere_epidemiology)>0){allepi<-graph_from_data_frame(epi,directed = FALSE,vertices = unique(c(epi$disease1,epi$disease2,nodesthatshouldbethere_epidemiology)))}
    if(length(nodesthatshouldbethere_epidemiology)==0){allepi<-graph_from_data_frame(epi,directed = FALSE)}
    ## Transform the table into an igraph object - Omics ##
    if(length(nodesthatshouldbethere_omics)>0){omics<-graph_from_data_frame(thenet,directed = FALSE,vertices = unique(c(thenet$disease1,thenet$disease2,nodesthatshouldbethere_omics)))}
    if(length(nodesthatshouldbethere_omics)==0){omics<-graph_from_data_frame(thenet,directed = FALSE)}
    ## Overlapping edges ##
    overlapingedges<-E(intersection(allepi,omics))
    theoverlappingedges<-as_ids(overlapingedges)
    ## Number of overlapping edges ##
    noverlapingedges<-length(overlapingedges)
    ## Random similarities ##
    roverlap<-c()
    for(b in 1:numberrandom){
      omics2<-rewire(omics, with=keeping_degseq(loops=FALSE, niter=ecount(omics)*10))
      roverlap<-c(roverlap,length(E(intersection(allepi,omics2))))
      print(paste(round((b/numberrandom)*100,2),"%",sep=""))
    }
    ## Obtain the pvalues of the overlap ##
    pvalue<-length(which(roverlap>=noverlapingedges))/length(roverlap)
    resultados<-list("Tabla"=c(comparison,length(E(omics)),length(E(allepi)),noverlapingedges,pvalue,mean(roverlap),length(V(omics)),length(V(allepi)),length(V(graph.intersection(omics,allepi))),
                               paste(sort(components(omics)$csize,decreasing = T),collapse = "-"),paste(sort(components(allepi)$csize,decreasing = T),collapse = "-"),
                               paste(sort(components(graph.intersection(omics,allepi))$csize,decreasing = T),collapse = "-")),
                     "omics"=E(omics),"epidemiology"=E(allepi),"overlapping"=theoverlappingedges)
    return(resultados)
  }
  
  # args<-c("calculate_overlaps","inter","Similarities/Similarities_Interlayer_with_fdr/ICD10_drug_gene_associations.txt","Denmark","all")
  ## Calculate the overlaps ##
  if(args[2]=="intra"){
    thecomparison<-gsub("_.+","",gsub(".+ICD10_","",args[3]))
    thenet<-fread(args[3],stringsAsFactors = F,sep="\t")
    if(args[4]=="Denmark"){
      resultados<-calculateoverlaps(thecomparison,thenet,10000,args[5],args[4])
      save(resultados,file=paste("EpidemiologyOverlap/",gsub(".txt","",gsub(".+ICD10_","",args[3])),"_all_Denmark.RData",sep=""))
    }
    if(args[4]=="Catalonia"){
      resultados<-calculateoverlaps(thecomparison,thenet,10000,args[5],args[4])
      save(resultados,file=paste("EpidemiologyOverlap/",gsub(".txt","",gsub(".+ICD10_","",args[3])),"_",args[5],"_Catalonia.RData",sep=""))
    }
  }
  if(args[2]=="inter"){
    thecomparison<-gsub("_associations.+","",gsub(".+ICD10_","",args[3]))
    thenet<-fread(args[3],stringsAsFactors = F,sep="\t") ; colnames(thenet)[1:2]<-c("disease1","disease2")
    sortednet<-t(apply(thenet[,1:2], 1, sort))
    thenet$disease1<-sortednet[,1] ; thenet$disease2<-sortednet[,2]
    
    if(args[4]=="Denmark"){
      resultados<-calculateoverlaps(thecomparison,thenet,10000,args[5],args[4])
      save(resultados,file=paste("EpidemiologyOverlap/",thecomparison,"_all_Denmark.RData",sep=""))
    }
    if(args[4]=="Catalonia"){
      resultados<-calculateoverlaps(thecomparison,thenet,10000,args[5],args[4])
      save(resultados,file=paste("EpidemiologyOverlap/",thecomparison,"_",args[5],"_Catalonia.RData",sep=""))
    }
  }
}

## Calculate the overlap of the different layers with epidemiology
if(args[1]=="calculate_overlaps_shared_diseases"){
  # ## Save the comorbidity networks for the shared diseases ##
  # moesme<-fread("Epidemiology/41467_2019_8475_MOESM6_ESM.txt",stringsAsFactors = F)
  # dendis<-unique(c(moesme$A,moesme$B))
  # moesme<-cbind(moesme[,c(1,2,7)],moesme$nWomen+moesme$nMen)
  # colnames(moesme)[4]<-"caseswith"
  # catalonia<-fread("Epidemiology/Catalonia/N0_5_Adjusted.txt",stringsAsFactors = F,sep="\t")
  # catdis<-unique(c(catalonia$Disease1,catalonia$Disease2))
  # 
  # comunes<-intersect(dendis,catdis)
  # 
  # ind1c<-c() ; ind2c<-c() ; ind1d<-c() ; ind2d<-c()
  # for(a in 1:length(comunes)){
  #   # a<-1
  #   ind1c<-unique(c(ind1c,which(catalonia$Disease1==comunes[a])))
  #   ind2c<-unique(c(ind2c,which(catalonia$Disease2==comunes[a])))
  #   ind1d<-unique(c(ind1d,which(moesme$A==comunes[a])))
  #   ind2d<-unique(c(ind2d,which(moesme$B==comunes[a])))
  # }
  # comcat<-catalonia[intersect(ind1c,ind2c)]
  # comden<-moesme[intersect(ind1d,ind2d)] ## it doesn't change
  # write.table(comcat,"Epidemiology/Catalonia_shared_diseases.txt",quote=F,sep="\t",row.names = F)
  
  ## The function ##
  calculateoverlaps<-function(comparison,thenet,numberrandom){
    set.seed(1)
    ## Get the icds ##
    icds2<-unique(c(thenet[[1]],thenet[[2]]))
    epidemiology2<-fread("Epidemiology/Catalonia_shared_diseases.txt",stringsAsFactors = F,sep="\t")
    colnames(epidemiology2)[1:2]<-c("disease1","disease2")
    ## Diseases in both networks ##
    icds<-intersect(icds2,unique(c(epidemiology2$disease1,epidemiology2$disease2)))
    ## Select the comorbidities between shared diseases ##
    ind1<-c() ; ind2<-c()
    for(a in 1:length(icds)){
      ind1<-c(ind1,which(epidemiology2$disease1==icds[a]))
      ind2<-c(ind2,which(epidemiology2$disease2==icds[a]))
    }
    ## Select comorbidities between shared diseases ##
    epidemiology<-epidemiology2[intersect(ind1,ind2),]
    ## Select the co-occurrences that are significant ##
    epi<-epidemiology[intersect(intersect(which(as.numeric(epidemiology$OR)>1),which(as.numeric(epidemiology$FDR)<=0.05)),which(as.numeric(epidemiology$CasesWith)>=100)),]
    sortednet<-t(apply(epi[,1:2], 1, sort))
    epi$disease1<-sortednet[,1] ; epi$disease2<-sortednet[,2]
    if(length(which(duplicated(epi[,1:2])))>0){epi<-epi[-which(duplicated(epi[,1:2])),]}
    ## Which are the diseases that should be included as isolated nodes in the epidemiology? ##
    nodesthatshouldbethere_epidemiology<-setdiff(icds,unique(c(epi$disease1,epi$disease2)))
    ## All the interactions or just the significant?? FDR ##
    if(comparison=="transcriptome"){thenet<-thenet[intersect(which(as.numeric(thenet[[3]])>0),which(as.numeric(thenet$fdr)<=0.05))]}
    if(comparison!="transcriptome"){thenet<-thenet[intersect(which(as.numeric(thenet[[3]])>0),which(as.numeric(thenet$fdr_up)<=0.05))]}
    thenet<-thenet[,1:2]
    ## Select common diseases ##
    colnames(thenet)<-c("disease1","disease2")
    ## Select the similarities between shared diseases ##
    ind1<-c() ; ind2<-c()
    for(a in 1:length(icds)){
      ind1<-c(ind1,which(thenet$disease1==icds[a]))
      ind2<-c(ind2,which(thenet$disease2==icds[a]))
    }
    thenet<-thenet[intersect(ind1,ind2)]
    
    ## Which are the diseases that should be included as isolated nodes in the omics layer? ##
    nodesthatshouldbethere_omics<-setdiff(icds,unique(c(thenet$disease1,thenet$disease2)))
    ## Calculate overlaps ##
    ## @@ @@ @@  @@ @@ @@ ##
    print(paste("Calculating significance of the overlap in ",comparison,":",sep=""))
    ## Transform the table into an igraph object - Epidemiology ##
    if(length(nodesthatshouldbethere_epidemiology)>0){allepi<-graph_from_data_frame(epi,directed = FALSE,vertices = unique(c(epi$disease1,epi$disease2,nodesthatshouldbethere_epidemiology)))}
    if(length(nodesthatshouldbethere_epidemiology)==0){allepi<-graph_from_data_frame(epi,directed = FALSE)}
    ## Transform the table into an igraph object - Omics ##
    if(length(nodesthatshouldbethere_omics)>0){omics<-graph_from_data_frame(thenet,directed = FALSE,vertices = unique(c(thenet$disease1,thenet$disease2,nodesthatshouldbethere_omics)))}
    if(length(nodesthatshouldbethere_omics)==0){omics<-graph_from_data_frame(thenet,directed = FALSE)}
    ## Overlapping edges ##
    overlapingedges<-E(intersection(allepi,omics))
    theoverlappingedges<-as_ids(overlapingedges)
    ## Number of overlapping edges ##
    noverlapingedges<-length(overlapingedges)
    ## Random similarities ##
    roverlap<-c()
    for(b in 1:numberrandom){
      omics2<-rewire(omics, with=keeping_degseq(loops=FALSE, niter=ecount(omics)*10))
      roverlap<-c(roverlap,length(E(intersection(allepi,omics2))))
      print(paste(round((b/numberrandom)*100,2),"%",sep=""))
    }
    ## Obtain the pvalues of the overlap ##
    pvalue<-length(which(roverlap>=noverlapingedges))/length(roverlap)
    resultados<-list("Tabla"=c(comparison,length(E(omics)),length(E(allepi)),noverlapingedges,pvalue,mean(roverlap),length(V(omics)),length(V(allepi)),length(V(graph.intersection(omics,allepi))),
                               paste(sort(components(omics)$csize,decreasing = T),collapse = "-"),paste(sort(components(allepi)$csize,decreasing = T),collapse = "-"),
                               paste(sort(components(graph.intersection(omics,allepi))$csize,decreasing = T),collapse = "-")),
                     "omics"=E(omics),"epidemiology"=E(allepi),"overlapping"=theoverlappingedges)
    return(resultados)
  }
  intranets<-c("Similarities/Similarities_with_fdr/ICD10_drug_fs_weightedjaccardindex.txt",
               "Similarities/Similarities_with_fdr/ICD10_gene_fs_weightedjaccardindex.txt",
               "Similarities/Similarities_with_fdr/ICD10_metabolome_fs_cosinesimilarities.txt",
               "Similarities/Similarities_with_fdr/ICD10_methylome_ifs_cosinesimilarities.txt",
               "Similarities/Similarities_with_fdr/ICD10_microbiome_fsd_cosinesimilarities.txt",
               "Similarities/Similarities_with_fdr/ICD10_miRNA_fs_cosinesimilarities.txt",
               "Similarities/Similarities_with_fdr/ICD10_symptom_fs_weightedjaccardindex.txt",
               "Similarities/Similarities_with_fdr/ICD10_transcriptome_fc_cosinesimilarities.txt")
  for(qq in 1:length(intranets)){
    # qq<-1
    thecomparison<-gsub("_.+","",gsub(".+ICD10_","",intranets[qq]))
    thenet<-fread(intranets[qq],stringsAsFactors = F,sep="\t")
    resultados<-calculateoverlaps(thecomparison,thenet,10000)
    save(resultados,file=paste("EpidemiologyOverlap/",thecomparison,"_Catalonia_shared.RData",sep=""))
  }
}

if(args[1]=="overlap_epidemiological_networks"){
  epi<-fread("Epidemiology/ICD10_epidemiology_denmark.txt",stringsAsFactors = F)
  moesme<-fread("Epidemiology/41467_2019_8475_MOESM6_ESM.txt",stringsAsFactors = F)
  dendis<-unique(c(moesme$A,moesme$B))
  moesme<-cbind(moesme[,c(1,2,7)],moesme$nWomen+moesme$nMen)
  colnames(moesme)[4]<-"caseswith"
  denmark<-moesme[which(moesme$adjustedRR2.5>1.01)]
  catalonia<-fread("Epidemiology/Catalonia/N0_5_Adjusted.txt",stringsAsFactors = F,sep="\t")
  catdis<-unique(c(catalonia$Disease1,catalonia$Disease2))
  catalonia<-catalonia[intersect(intersect(which(catalonia$OR>1),which(catalonia$FDR<=0.05)),which(catalonia$CasesWith>=100))]
  
  comunes<-intersect(dendis,catdis)
  
  ind1c<-c() ; ind2c<-c() ; ind1d<-c() ; ind2d<-c()
  for(a in 1:length(comunes)){
    # a<-1
    ind1c<-unique(c(ind1c,which(catalonia$Disease1==comunes[a])))
    ind2c<-unique(c(ind2c,which(catalonia$Disease2==comunes[a])))
    ind1d<-unique(c(ind1d,which(denmark$A==comunes[a])))
    ind2d<-unique(c(ind2d,which(denmark$B==comunes[a])))
  }
  comcat<-catalonia[intersect(ind1c,ind2c)]
  comden<-denmark[intersect(ind1d,ind2d)]
  
  pcomcat<-paste(comcat$Disease1,comcat$Disease2,sep="_")
  pcomden<-paste(comden$A,comden$B,sep="_")
  
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  venn.diagram(list("Denmark"=pcomden,"Catalonia"=pcomcat),filename = 'Plots/Venn_diagram_comorbidities_populations.png',category.names = c("Denmark","Catalonia"),
               main="Comorbidities shared between\nDenmark and Catalonia (adjusting for gender)",output=TRUE,
               lwd = 2,lty = 'blank',fill = c("#3A5382","#19787F"),cex=2)
}

## Summarize results ##
if(args[1]=="summaryze_overlaps"){
  overlaps<-list.files("EpidemiologyOverlap/")
  entire_population<-overlaps[grep("all",overlaps)]
  catalonia<-entire_population[grep("Catalonia",entire_population)]
  denmark<-entire_population[grep("Denmark",entire_population)]
  theoverlapsd<-c() ; theoverlapsc<-c() ; loverlaps<-list() ; value<-c() ; population<-c() ; layer<-c() ; colors<-c()
  files<-gsub("_Denmark.RData","",denmark)
  files<-setdiff(files,"epidemiology_relativerisk_all")
  for(a in 1:length(files)){
    # a<-1
    load(paste("EpidemiologyOverlap/",files[a],"_Denmark.RData",sep=""))
    theoverlapsd<-rbind(theoverlapsd,as.numeric(resultados$Tabla[c(2:5,7)]))
    loverlaps[[gsub("_all.+","",files[a])]]$Denmark<-intersect(as_ids(resultados$epidemiology),as_ids(resultados$omics))
    load(paste("EpidemiologyOverlap/",files[a],"_Catalonia.RData",sep=""))
    theoverlapsc<-rbind(theoverlapsc,as.numeric(resultados$Tabla[c(2:5,7)]))
    loverlaps[[gsub("_all.+","",files[a])]]$Catalonia<-intersect(as_ids(resultados$epidemiology),as_ids(resultados$omics))
    population<-c(population,"Denmark","Catalonia","Shared")
    value<-c(value,length(setdiff(loverlaps[[a]]$Denmark,loverlaps[[a]]$Catalonia)),length(setdiff(loverlaps[[a]]$Catalonia,loverlaps[[a]]$Denmark)),length(intersect(loverlaps[[a]]$Catalonia,loverlaps[[a]]$Denmark)))
    layer<-c(layer,rep(gsub("_all.+","",files[a]),3))
    colors<-c(colors,"#9AA2B8","#94B4B7","#CC9E9A")
  }
  rownames(theoverlapsc)<-gsub("_all","",files)
  rownames(theoverlapsd)<-gsub("_all","",files)
  theoverlapsc<-theoverlapsc[c(setdiff(1:length(rownames(theoverlapsc)),grep("_",rownames(theoverlapsc))),grep("_",rownames(theoverlapsc))),]
  theoverlapsd<-theoverlapsd[c(setdiff(1:length(rownames(theoverlapsd)),grep("_",rownames(theoverlapsd))),grep("_",rownames(theoverlapsd))),]
  theoverlapsc<-cbind(theoverlapsc,round((theoverlapsc[,3]/theoverlapsc[,2])*100,2),round((theoverlapsc[,3]/theoverlapsc[,1])*100,2))
  colnames(theoverlapsc)<-c("similarities (omics)","comorbidities","overlap","p-value","shared diseases","recall","precision")
  theoverlapsd<-cbind(theoverlapsd,round((theoverlapsd[,3]/theoverlapsd[,2])*100,2),round((theoverlapsd[,3]/theoverlapsd[,1])*100,2))
  colnames(theoverlapsd)<-c("similarities (omics)","comorbidities","overlap","p-value","shared diseases","recall","precision")
  if("Results"%in%list.files()==F){dir.create("Results")}
  write.table(theoverlapsd,"Results/denmark_recall_and_precision_all.txt",quote=F,sep="\t")
  write.table(theoverlapsc,"Results/catalonia_recall_and_precision_all.txt",quote=F,sep="\t")
  
  ## Select the best weight and metric - Denmark ##
  ## @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ ##
  theoverlapsd<-read.csv2("Results/denmark_recall_and_precision_all.txt",stringsAsFactors = F,sep="\t")
  
  ## layer by layer select the best method combining pvalues with recall and precission ##
  layers<-gsub("_.+","",rownames(theoverlapsd)) ; ulayers<-unique(layers)
  ## calculate a score for selecting the most relevant layers ##
  sumtabd<-c()
  for(a in 1:(length(ulayers)-1)){
    # a<-1
    submat<-theoverlapsd[which(layers==ulayers[a]),]
    ## select only those significantly recovering comorbidity relationships ##
    submat<-submat[which(submat$p.value<=0.05),]
    submat2<-cbind(submat,as.numeric(submat$recall)*as.numeric(submat$precision))
    colnames(submat2)[8]<-"score"
    print(submat2)
    sumtabd<-rbind(sumtabd,submat2[which(submat2$score==max(submat2$score)),])
  }
  transcriptomics<-theoverlapsd[grep("transcriptome",rownames(theoverlapsd)),]
  transcriptomics<-cbind(transcriptomics,as.numeric(transcriptomics[,6])*as.numeric(transcriptomics[,7]))
  colnames(transcriptomics)[8]<-"score"
  sumtabd<-rbind(sumtabd,transcriptomics)
  write.table(sumtabd,"Data/Intermediate/Best_layers_for_comorbidites.txt",quote=F,sep="\t")
  splited<-strsplit(rownames(sumtabd),"_")
  caracteristicas<-c() ; for(a in 1:length(splited)){caracteristicas<-rbind(caracteristicas,splited[[a]])}
  resultadod<-cbind(caracteristicas,sumtabd)
  rownames(resultadod)<-NULL
  colnames(resultadod)<-c("layer","weight","metric","similarities","comorbidities","overlap","pvalue","shared_diseases","recall","precision","score")
  write.table(resultadod,"Results/Multiplex_overlap_comorbidities_denmark.txt",quote=F,sep="\t",row.names=F)
  
  ## Select the best weight and metric - Catalonia ##
  ## @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ ##
  theoverlapsc<-read.csv2("Results/catalonia_recall_and_precision_all.txt",stringsAsFactors = F,sep="\t")
  
  ## layer by layer select the best method combining pvalues with recall and precission ##
  layers<-gsub("_.+","",rownames(theoverlapsc)) ; ulayers<-unique(layers)
  ## calculate a score for selecting the most relevant layers ##
  sumtabc<-c()
  for(a in 1:(length(ulayers)-1)){
    # a<-1
    submat<-theoverlapsc[which(layers==ulayers[a]),]
    ## select only those significantly recovering comorbidity relationships ##
    submat<-submat[which(submat$p.value<=0.05),]
    submat2<-cbind(submat,as.numeric(submat$recall)*as.numeric(submat$precision))
    colnames(submat2)[8]<-"score"
    print(submat2)
    sumtabc<-rbind(sumtabc,submat2[which(submat2$score==max(submat2$score)),])
  }
  transcriptomics<-theoverlapsc[grep("transcriptome",rownames(theoverlapsc)),]
  transcriptomics<-cbind(transcriptomics,as.numeric(transcriptomics[,6])*as.numeric(transcriptomics[,7]))
  colnames(transcriptomics)[8]<-"score"
  sumtabc<-rbind(sumtabc,transcriptomics)
  write.table(sumtabc,"Data/Intermediate/Best_layers_for_comorbidites_catalonia.txt",quote=F,sep="\t")
  splited<-strsplit(rownames(sumtabc),"_")
  caracteristicas<-c() ; for(a in 1:length(splited)){caracteristicas<-rbind(caracteristicas,splited[[a]])}
  resultadoc<-cbind(caracteristicas,sumtabc)
  rownames(resultadoc)<-NULL
  colnames(resultadoc)<-c("layer","weight","metric","similarities","comorbidities","overlap","pvalue","shared_diseases","recall","precision","score")
  write.table(resultadoc,"Results/Multiplex_overlap_comorbidities_catalonia.txt",quote=F,sep="\t",row.names=F)
  
  ## Plot recalls vs. precisions for all weights and similarity metrics ##
  ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
  ## put in red non-significant overlaps ##
  pdf(file="Plots/precision_vs_recall_all_combinations.pdf")
    colores<-rep("black",length(theoverlapsc$similarities..omics.))
    colores[which(theoverlapsc$p.value>0.05)]<-"red"
    plot(theoverlapsc$recall,theoverlapsc$precision,xlim=c(0,100),ylim=c(0,100),xlab="recall",ylab="precision",col=colores,main="Catalonia")
    colores<-rep("black",length(theoverlapsd$similarities..omics.))
    colores[which(theoverlapsd$p.value>0.05)]<-"red"
    plot(theoverlapsd$recall,theoverlapsd$precision,xlim=c(0,100),ylim=c(0,100),xlab="recall",ylab="precision",col=colores,main="Denmark")
  dev.off()
  
  ## Barplot Recalls Denmark and Catalonia ##
  ## @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ ##
  resultadod$layer[which(resultadod$layer=="symptom")]<-"symptom\n(2014)"
  resultadod$layer[which(resultadod$layer=="mksymptom")]<-"symptom\n(2025)"
  resultadod$layer[which(resultadod$layer=="drug")]<-"drug\n(CTD)"
  resultadod$layer[which(resultadod$layer=="drugchembl")]<-"drug\n(ChEMBL)"
  resultadod$layer[which(resultadod$layer=="drugchembl2")]<-"drug\n(ChEMBL-all)"
  resultadoc$layer[which(resultadoc$layer=="symptom")]<-"symptom\n(2014)"
  resultadoc$layer[which(resultadoc$layer=="mksymptom")]<-"symptom\n(2025)"
  resultadoc$layer[which(resultadoc$layer=="drug")]<-"drug\n(CTD)"
  resultadoc$layer[which(resultadoc$layer=="drugchembl")]<-"drug\n(ChEMBL)"
  
  ## Prepare the data.frame for the barplot ##
  layer<-rep(resultadod$layer,2)
  population<-c(rep("Catalonia",11),rep("Denmark",11))
  colors<-c(rep("#94B4B7",11),rep("#9AA2B8",11))
  recall<-as.numeric(c(resultadoc$recall,resultadod$recall))
  
  df2<-data.frame(layer,population,recall,colors)
  pdf(file="Plots/Barplot_recalls_catalonia_denmark.pdf",width = 14)
  ggplot(data=df2, aes(x=layer, y=recall, fill=population)) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=recall), vjust=1.6, color="black",
              position = position_dodge(0.9), size=5)+
    scale_fill_manual(values=c("#94B4B7","#9AA2B8"))+
    theme_minimal() +
    ylim(c(0,100)) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(size=14),
          legend.title=element_text(size=18))
  dev.off()
  
  ## Barplot Recalls only Denmark ##
  ## @@ @@ @@ @@ @@@@ @@ @@ @@ @@ ##
  solodenmark<-resultadod[-10,]
  layer<-solodenmark$layer
  population<-rep("Denmark",10)
  colors<-rep("#9AA2B8",10)
  recall<-as.numeric(solodenmark$recall)
  
  df2<-data.frame(layer,population,recall,colors)
  pdf(file="Plots/Barplot_recalls_denmark.pdf",width = 14)
    ggplot(data=df2, aes(x=layer, y=recall, fill=population)) +
      geom_bar(stat="identity", position=position_dodge())+
      geom_text(aes(label=recall), vjust=1.6, color="black",
                position = position_dodge(0.9), size=5)+
      scale_fill_manual(values="#9AA2B8")+
      theme_minimal() +
      ylim(c(0,100)) +
      theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
            axis.title.x = element_text(color = "grey20", size = 18, angle = 0, hjust = .5, vjust = 0, face = "plain"),
            axis.title.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            legend.text = element_text(size=14),
            legend.title=element_text(size=18))
  dev.off()
  
  ## Plot recalls by each layer in each population (and shared) ##
  ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@@@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
  denmarks<-paste(rownames(sumtabd),"_all_Denmark.RData",sep="")
  catalonias<-paste(rownames(sumtabd),"_all_Catalonia.RData",sep="")
  theoverlapsd<-c() ; theoverlapsc<-c() ; loverlaps<-list() ; value<-c() ; population<-c() ; layer<-c() ; colors<-c()
  for(a in 1:length(denmarks)){
    load(paste("EpidemiologyOverlap/",denmarks[a],sep=""))
    theoverlapsd<-rbind(theoverlapsd,as.numeric(resultados$Tabla[c(2:5,7)]))
    loverlaps[[gsub("_all.+","",denmarks[a])]]$Denmark<-intersect(as_ids(resultados$epidemiology),as_ids(resultados$omics))
    load(paste("EpidemiologyOverlap/",catalonias[a],sep=""))
    theoverlapsc<-rbind(theoverlapsc,as.numeric(resultados$Tabla[c(2:5,7)]))
    loverlaps[[gsub("_all.+","",catalonias[a])]]$Catalonia<-intersect(as_ids(resultados$epidemiology),as_ids(resultados$omics))
    population<-c(population,"Denmark","Catalonia","Shared")
    value<-c(value,length(setdiff(loverlaps[[a]]$Denmark,loverlaps[[a]]$Catalonia)),length(setdiff(loverlaps[[a]]$Catalonia,loverlaps[[a]]$Denmark)),length(intersect(loverlaps[[a]]$Catalonia,loverlaps[[a]]$Denmark)))
    layer<-c(layer,rep(gsub("_.+","",denmarks[a]),3))
    colors<-c(colors,"#9AA2B8","#94B4B7","#CC9E9A")
    print(a)
  }
  
  ## Prepare the data.frame for the barplot ##
  layer[which(layer=="drug")]<-"drug\n(CTD)"
  layer[which(layer=="drugchembl")]<-"drug\n(ChEMBL)"
  layer[which(layer=="drugchembl2")]<-"drug\n(ChEMBL-all)"
  layer[which(layer=="symptom")]<-"symptom\n(2014)"
  layer[which(layer=="mksymptom")]<-"symptom\n(2025)"
  data <- data.frame(layer,population,value,colors)
  data$population<-factor(data$population,levels=c("Denmark","Shared","Catalonia"))
  pdf(file = "Plots/Barplot_overlap_epidemiology.pdf")
    # Stacked + percent
    ggplot(data, aes(fill=population, y=value, x=layer)) + 
      geom_bar(position="fill", stat="identity") + ## position=fill es lo que te dice que agrupe por categoria (layer) y te muestre en porcentajes
      scale_fill_manual(values=c("#9AA2B8","#CC9E9A","#94B4B7")) +
      theme(axis.text.x = element_text(angle = 90, hjust = .5),
            panel.background = element_blank(),
            legend.title=element_blank(),
            axis.ticks.y=element_blank()) +
      geom_text(aes(label = value), 
                position = position_fill(vjust = 0.5), size = 3) +
      scale_y_continuous(labels = scales::percent_format(scale = 100)) ## scale hace que en lugar de mostrar valores del 0 al 1 los muestre del 0 a 100
  dev.off()
}
if(args[1]=="Create_final_layers"){
  ## Create a directory with the networks ##
  sumtab<-read.csv2("Data/Intermediate/Best_layers_for_comorbidites.txt",stringsAsFactors = F,sep="\t")
  if("Networks"%in%list.files()==F){dir.create("Networks")}
  if("Intralayer"%in%list.files("Networks/")==F){dir.create("Networks/Intralayer")}
  if("Interlayer"%in%list.files("Networks/")==F){dir.create("Networks/Interlayer")}
  densidades<-c()
  edges<-c()
  for(a in 1:length(sumtab$comorbidities)){
    # a<-1
    if(a!=11){
      tt<-fread(paste("Similarities/Similarities_with_fdr/ICD10_",rownames(sumtab)[a],".txt",sep=""),stringsAsFactors = F,sep="\t")
      diseases<-unique(c(tt$disease1,tt$disease2))
      ## Calculate the maximum number of potential interactions ##
      potential_interactions<-(length(diseases)*(length(diseases)-1))/2
      if(length(which(tt$fdr_up<=0.05))>0){
        tt2<-tt[which(tt$fdr_up<=0.05)]
      }
      tt2<-tt2[,1:3]
      colnames(tt2)<-c("disease1","disease2","weight")
      densidades<-c(densidades,length(tt2$disease1)/potential_interactions)
      diffes<-setdiff(diseases,c(tt2$disease1,tt2$disease2))
      if(length(diffes)>0){
        diftab<-cbind(diffes,NA,NA) ; colnames(diftab)<-colnames(tt2)
        tt2<-rbind(tt2,diftab)
      }
      edges<-rbind(edges,c(length(unique(c(tt2$disease1,tt2$disease2))),length(which(is.na(tt2$disease2)==FALSE))))
      if(gsub("_.+",".txt",rownames(sumtab)[a])=="drug.txt"){
        write.table(tt2,paste("Networks/Intralayer/drugs_ctd.txt",sep=""),quote=F,sep="\t",row.names=F,col.names = F)
      }
      if(gsub("_.+",".txt",rownames(sumtab)[a])=="drugchembl.txt"){
        write.table(tt2,paste("Networks/Intralayer/drugs_chembl.txt",sep=""),quote=F,sep="\t",row.names=F,col.names = F)
      }
      if(gsub("_.+",".txt",rownames(sumtab)[a])=="drugchembl2.txt"){
        write.table(tt2,paste("Networks/Intralayer/drugs_chemblall.txt",sep=""),quote=F,sep="\t",row.names=F,col.names = F)
      }
      if(gsub("_.+",".txt",rownames(sumtab)[a])=="symptom.txt"){
        write.table(tt2,paste("Networks/Intralayer/symptom_2014.txt",sep=""),quote=F,sep="\t",row.names=F,col.names = F)
      }
      if(gsub("_.+",".txt",rownames(sumtab)[a])=="mksymptom.txt"){
        write.table(tt2,paste("Networks/Intralayer/symptom_2025.txt",sep=""),quote=F,sep="\t",row.names=F,col.names = F)
      }
      if(gsub("_.+",".txt",rownames(sumtab)[a])!="drug.txt" && gsub("_.+",".txt",rownames(sumtab)[a])!="drugchembl.txt" && gsub("_.+",".txt",rownames(sumtab)[a])!="drugchembl2.txt" && gsub("_.+",".txt",rownames(sumtab)[a])!="symptom.txt" && gsub("_.+",".txt",rownames(sumtab)[a])!="mksymptom.txt"){
        write.table(tt2,paste("Networks/Intralayer/",gsub("_.+",".txt",rownames(sumtab)[a]),sep=""),quote=F,sep="\t",row.names=F,col.names = F)
      }
    }
    if(a==11){
      tt<-fread(paste("Similarities/Similarities_with_fdr/ICD10_",rownames(sumtab)[a],".txt",sep=""),stringsAsFactors = F,sep="\t")
      if(length(which(tt$fdr<=0.05))>0){
        tt2<-tt[which(tt$fdr<=0.05)]
        tt2<-tt2[which(tt2$cosinesimilarities>0)]
      }
      tt2<-tt2[,1:3]
      colnames(tt2)<-c("disease1","disease2","weight")
      densidades<-c(densidades,length(tt2$disease1)/potential_interactions)
      diffes<-setdiff(c(tt$disease1,tt$disease2),c(tt2$disease1,tt2$disease2))
      if(length(diffes)>0){
        diftab<-cbind(diffes,NA,NA) ; colnames(diftab)<-colnames(tt2)
        tt2<-rbind(tt2,diftab)
      }
      edges<-rbind(edges,c(length(unique(c(tt2$disease1,tt2$disease2))),length(which(is.na(tt2$disease2)==FALSE))))
      write.table(tt2,paste("Networks/Intralayer/",gsub("_.+",".txt",rownames(sumtab)[a]),sep=""),quote=F,sep="\t",row.names=F,col.names = F)
    }
  }
  splited<-strsplit(rownames(sumtab),"_")
  tabla<-c()
  for(a in 1:length(splited)){tabla<-rbind(tabla,c(splited[[a]][1],splited[[a]][2],splited[[a]][3]))}
  tabla<-cbind(tabla,sumtab,densidades)
  rownames(tabla)<-NULL
  tabla<-cbind(tabla[,c(1:3)],edges,tabla[,12],tabla[,c(4:11)])
  tabla[which(tabla[,1]=="drug"),1]<-"drug (CTD)" ; tabla[which(tabla[,1]=="drugchembl"),1]<-"drug (ChEMBL)" ; tabla[which(tabla[,1]=="drugchembl2"),1]<-"drug (ChEMBL-all)"
  tabla[which(tabla[,1]=="symptom"),1]<-"symptom (2014)" ; tabla[which(tabla[,1]=="mksymptom"),1]<-"symptom (2025)"
  colnames(tabla)<-c("layer","weight","metric","nodes","edges","density","similarities_shared_diseases","comorbidities_shared_diseases","overlap","pvalue","shared_diseases","recall","precision","score")
  write.table(tabla,"Results/Tabla_nodos_y_recall_densidad.txt",quote=F,sep="\t")
}

if(args[1]=="ChEMBL_vs_CTD"){
  drug_chembl<-read.csv2("Networks/Intralayer/drugs_chembl.txt",stringsAsFactors = F,sep="\t",header=F)
  drug_ctd<-read.csv2("Networks/Intralayer/drugs_ctd.txt",stringsAsFactors = F,sep="\t",header=F)
  
  gchembl<-graph_from_data_frame(drug_chembl,directed = F)
  gctd<-graph_from_data_frame(drug_ctd,directed = F)
  
  # Specify nodes of interest by name
  nodes_of_interest <- intersect(c(drug_chembl$V1,drug_chembl$V2),c(drug_ctd$V1,drug_ctd$V2))
  if(length(which(is.na(nodes_of_interest)))>0){nodes_of_interest<-nodes_of_interest[-which(is.na(nodes_of_interest))]}
  
  # Select vertices by name
  selnodeschembl <- V(gchembl)[name %in% nodes_of_interest]
  selnodesctd <- V(gctd)[name %in% nodes_of_interest]
  
  # Create the induced subgraph
  gschembl <- induced_subgraph(gchembl, selnodeschembl)
  gsctd <- induced_subgraph(gctd, selnodesctd)
  
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



if(args[1]=="distances"){
  ## Drug ##
  ## @@@@ ##
  red<-fread("HDMLN/Intralayers/drug.txt",stringsAsFactors = F)
  rednet<-graph_from_data_frame(red, directed = FALSE)
  ## it is considering the weight of the edges
  shortest_paths(
    rednet,
    from = "U07",
    to = "G30",
    mode ="all",
    output = "epath"
    # algorithm = c("automatic", "unweighted", "dijkstra", "bellman-ford", "johnson")
  )
  ## Get all paths of length X
  buf2<-all_simple_paths(rednet, "U07", "G30", cutoff=2)
  buf3<-all_simple_paths(rednet, "U07", "G30", cutoff=3)
  buf4<-all_simple_paths(rednet, "U07", "G30", cutoff=4)
  
  enfermedades<-c() ; for(a in 1:length(buf2)){enfermedades<-c(enfermedades,as_ids(buf2[[a]][2]))}
  
  ## Gene ##
  ## @@@@ ##
  red<-fread("HDMLN/Intralayers/gene.txt",stringsAsFactors = F)
  rednet<-graph_from_data_frame(red, directed = FALSE)
  ## it is considering the weight of the edges
  shortest_paths(
    rednet,
    from = "U07",
    to = "G30",
    mode ="all",
    output = "epath"
    # algorithm = c("automatic", "unweighted", "dijkstra", "bellman-ford", "johnson")
  )
  ## Get all paths of length X
  buf2<-all_simple_paths(rednet, "U07", "G30", cutoff=2)
  buf3<-all_simple_paths(rednet, "U07", "G30", cutoff=3)
  buf4<-all_simple_paths(rednet, "U07", "G30", cutoff=4)
  
  ## miRNA ##
  ## @@ @@ ##
  red<-fread("HDMLN/Intralayers/mirna.txt",stringsAsFactors = F)
  rednet<-graph_from_data_frame(red, directed = FALSE)
  ## it is considering the weight of the edges
  shortest_paths(
    rednet,
    from = "U07",
    to = "G30",
    mode ="all",
    output = "epath"
    # algorithm = c("automatic", "unweighted", "dijkstra", "bellman-ford", "johnson")
  )
  ## Get all paths of length X
  buf2<-all_simple_paths(rednet, "U07", "G30", cutoff=2)
  buf3<-all_simple_paths(rednet, "U07", "G30", cutoff=3)
  buf4<-all_simple_paths(rednet, "U07", "G30", cutoff=4)
  
  red<-fread("HDMLN/Intralayers/microbiome.txt",stringsAsFactors = F)
  rednet<-graph_from_data_frame(red, directed = FALSE)
  ## it is considering the weight of the edges
  shortest_paths(
    rednet,
    from = "U07",
    to = "G30",
    mode ="all",
    output = "epath"
    # algorithm = c("automatic", "unweighted", "dijkstra", "bellman-ford", "johnson")
  )
  
}


if(args[1]=="Shared_diseases"){
  il<-list.files("HDMLN/Intralayers/")
  diseases<-c()
  for(a in 1:length(il)){
    # a<-1
    tt<-fread(paste("HDMLN/Intralayers/",il[a],sep=""),stringsAsFactors = F,sep="\t")
    diseases<-c(diseases,unique(c(tt$icd10_1,tt$icd10_2)))
  }
  diseases<-sort(table(diseases),decreasing = T)
}


if(args[1]=="Upsets"){
  files<-c("Similarities/Similarities_with_fdr/ICD10_symptom_fs_weightedjaccardindex.txt",
           "Similarities/Similarities_with_fdr/ICD10_drug_fs_weightedjaccardindex.txt",
           "Similarities/Similarities_with_fdr/ICD10_gene_fs_weightedjaccardindex.txt",
           "Similarities/Similarities_with_fdr/ICD10_miRNA_fs_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_metabolome_fs_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_microbiome_fsd_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_transcriptome_fc_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_methylome_ifs_cosinesimilarities.txt")
  
  omics<-gsub("_.+","",gsub("Similarities/Similarities_with_fdr/ICD10_","",files))
  nodel<-list()
  edgel<-list()
  diseases<-c()
  for(a in 1:length(files)){
    # a<-1
    tt<-fread(files[a],stringsAsFactors = F,sep="\t")
    nodel[[omics[a]]]<-unique(c(tt[[1]],tt[[2]]))
    diseases<-c(diseases,unique(c(tt[[1]],tt[[2]])))
    print(paste(omics[a],length(unique(c(tt[[1]],tt[[2]])))))
    if(length(which(tt[[5]]<0))>0){tt<-tt[-which(tt[[3]]==0)]}
    tt<-tt[which(tt[[5]]<=0.05)]
    tt<-tt[which(tt[[3]]>0)]
    edgel[[omics[a]]]<-paste(tt[[1]],tt[[2]],sep="_")
  }
  diseases<-sort(table(diseases),decreasing = T)
  as.data.frame(diseases[which(diseases==8)])
  input<-nodel
  tabla<-fromList(input)
  tabla<-tabla[8:1]
  pdf(file = "Plots/Upsetplot_shared_nodes.pdf",width = 12,height = 8)
    upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  dev.off()
  ## Overlaps in terms of similarities between the 5 main layers ##
  ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
  ## Add epidemiology ##
  tt<-fread("Epidemiology/ICD10_epidemiology_denmark.txt",stringsAsFactors = F,sep="\t")
  nodel$epidemiology<-unique(c(tt$disease1,tt$disease2))
  tt<-tt[which(tt$RR>1.01)]
  edgel$epidemiology<-paste(tt[[1]],tt[[2]],sep="_")
  nodel5big<-list("epidemiology"=nodel$epidemiology,"symptom"=nodel$symptom,"drug"=nodel$drug,"gene"=nodel$gene,"miRNA"=nodel$miRNA)
  ## Upset nodes 5big ##
  input<-nodel5big
  tabla<-fromList(input)
  tabla<-tabla[,5:1]
  pdf(file = "Plots/Upsetplot_shared_nodes_big5.pdf",width = 12,height = 8)
    upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  dev.off()
  commonnodes5<-Reduce(intersect,nodel5big)
  edgel5big<-list()
  for(a in 1:length(names(nodel5big))){
    # a<-1
    red<-edgel[[names(nodel5big)[a]]]
    ind1<-c() ; ind2<-c()
    for(b in 1:length(commonnodes5)){
      # b<-1
      ind1<-c(ind1,grep(paste(commonnodes5[b],"_",sep=""),red)) ; ind2<-c(ind2,grep(paste("_",commonnodes5[b],sep=""),red))
    }
    edgel5big[[names(nodel5big)[a]]]<-red[intersect(ind1,ind2)]
  }
  ## Upset edges 5big ##
  input<-edgel5big
  tabla<-fromList(input)
  tabla<-tabla[c(1,4,3,5,2)]
  pdf(file = "Plots/Upsetplot_shared_edges_big5.pdf",width = 12,height = 8)
    upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  dev.off()
  
  ## Overlaps in terms of similarities between the 5 main layers + transcriptomics ##
  ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
  # nodel6big<-list("epidemiology"=nodel$epidemiology,"symptom"=nodel$symptom,"drug"=nodel$drug,"gene"=nodel$gene,"mirna"=nodel$mirna,"transcriptome"=nodel$transcriptome)
  # ## Upset nodes 6big ##
  # input<-nodel6big
  # tabla<-fromList(input)
  # tabla<-tabla[,6:1]
  # pdf(file = "Plots/Upsetplot_shared_nodes_big5_transcriptomics.pdf",width = 12,height = 8)
  # upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  # dev.off()
  # commonnodes5<-Reduce(intersect,nodel6big)
  # edgel6big<-list()
  # for(a in 1:length(names(nodel6big))){
  #   # a<-1
  #   red<-edgel[[names(nodel6big)[a]]]
  #   ind1<-c() ; ind2<-c()
  #   for(b in 1:length(commonnodes5)){
  #     # b<-1
  #     ind1<-c(ind1,grep(paste(commonnodes5[b],"_",sep=""),red)) ; ind2<-c(ind2,grep(paste("_",commonnodes5[b],sep=""),red))
  #   }
  #   edgel6big[[names(nodel6big)[a]]]<-red[intersect(ind1,ind2)]
  # }
  # ## Upset edges 6big ##
  # input<-edgel6big
  # tabla<-fromList(input)
  # tabla<-tabla[c(1,3,4,6,5,2)]
  # pdf(file = "Plots/Upsetplot_shared_edges_big5_transcriptomics.pdf",width = 12,height = 8)
  # upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  # dev.off()
  # tabla<-tabla[which(tabla$epidemiology==1),c(4,3,2,5,6,1)]
  # pdf(file = "Plots/Upsetplot_epidemiology_edges_big5_transcriptomics.pdf",width = 12,height = 8)
  # upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  # dev.off()
  
  ## Overlaps in terms of similarities between the epidemiology & transcriptomics ##
  ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
  nodel2<-list("epidemiology"=nodel$epidemiology,"transcriptome"=nodel$transcriptome)
  commonnodes2<-Reduce(intersect,nodel2)
  edgel2<-list()
  for(a in 1:length(names(edgel))){
    # a<-1
    red<-edgel[[names(edgel)[a]]]
    ind1<-c() ; ind2<-c()
    for(b in 1:length(commonnodes2)){
      # b<-1
      ind1<-c(ind1,grep(paste(commonnodes2[b],"_",sep=""),red)) ; ind2<-c(ind2,grep(paste("_",commonnodes2[b],sep=""),red))
    }
    edgel2[[names(edgel)[a]]]<-red[intersect(ind1,ind2)]
  }
  ## Upset edges 5big ##
  input<-edgel2
  tabla<-fromList(input)
  tabla<-tabla[which(tabla$epidemiology==1),]
  tabla<-tabla[c(5,8,6,2,3,4,1,7)]
  pdf(file = "Plots/Upsetplot_epidemiology_edges_transcriptomics_all.pdf",width = 12,height = 8)
    upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  dev.off()
  ## barplot with the number of nodes in the different layers ##
  nodenumber<-c()
  for(a in 1:length((names(nodel)))){
    # a<-1
    nodenumber<-c(nodenumber,length(intersect(nodel[[a]],commonnodes2)))
  }
  percentage<-(nodenumber/length(commonnodes2))*100
  data<-as.data.frame(cbind(percentage,nodenumber))
  data<-cbind(names(nodel),data)
  data<-data[c(1,2,8,5:3,7,9,6),]
  tsimply<-data$nodenumber ; names(tsimply)<-data$`names(nodel)`
  pdf(file="Plots/Barplot_nodes_in_transcriptomic_epidemiology.pdf",width = 12,height = 12)
  p <- barplot(tsimply,col="black",ylim=c(0,100),main="",las=2)
  text(x = p, y = tsimply + 7, labels = data$nodenumber,srt=90,cex = 4)
  dev.off()
  
  
  
  ## Networks for common diseases ##
  commondis<-names(diseases[which(diseases==8)])
  files<-c("Similarities/Similarities_with_fdr/ICD10_symptom_fs_weightedjaccardindex.txt",
           "Similarities/Similarities_with_fdr/ICD10_drug_fs_weightedjaccardindex.txt",
           "Similarities/Similarities_with_fdr/ICD10_gene_fs_weightedjaccardindex.txt",
           "Similarities/Similarities_with_fdr/ICD10_miRNA_fs_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_metabolome_fs_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_microbiome_fsd_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_transcriptome_fc_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_methylome_ifs_cosinesimilarities.txt")
  
  sharedmultilayer<-c()
  for(a in 1:length(files)){
    # a<-1
    tt<-fread(files[a],stringsAsFactors = F,sep="\t")
    tt<-tt[which(tt[[5]]<=0.05)]
    tt<-tt[which(tt[[3]]>0)]
    ind1<-c() ; ind2<-c()
    for(b in 1:length(commondis)){
      # b<-1
      ind1<-c(ind1,which(tt[[1]]==commondis[b])) ; ind2<-c(ind2,which(tt[[2]]==commondis[b]))
    }
    comunes<-intersect(ind1,ind2)
    if(length(comunes)>0){
      subnet<-cbind(tt[comunes,1:3],gsub("_.+","",gsub("Similarities/Similarities_with_fdr/ICD10_","",files[a])))
      colnames(subnet)<-c("ICD10_1","ICD10_2","weight","layer")
      sharedmultilayer<-rbind(sharedmultilayer,subnet)
    }
  }
  write.table(sharedmultilayer,"Plots/multilayer_shared_nodes.txt",quote=F,sep="\t",row.names=F)
  
  sharedmultilayer2<-cbind(sharedmultilayer,sharedmultilayer$layer)
  colnames(sharedmultilayer2)<-c("SourceNode","TargetNode","weight","SourceLayer","TargetLayer")
  write.table(sharedmultilayer2,"Plots/multilayer_shared_nodes2.txt",quote=F,sep="\t",row.names=F)
  
  
  # epi<-fread("Epidemiology/ICD10_epidemiology_denmark.txt",stringsAsFactors = F)
  # epi<-epi[which(epi$RR>1.01)]
  # ind1<-c() ; ind2<-c()
  # for(b in 1:length(commondis)){
  #   # b<-1
  #   ind1<-c(ind1,which(epi[[1]]==commondis[b])) ; ind2<-c(ind2,which(epi[[2]]==commondis[b]))
  # }
  # comunes<-intersect(ind1,ind2)
  # epi<-epi[comunes]
  # for(a in 1:length(sharedmultilayer$ICD10_1)){
  #   a<-1
  #   
  # }
  
  
}


if(args[1]=="Upsets_20250531"){
  files<-paste("Networks/Intralayer/",list.files("Networks/Intralayer/"),sep="")
  files<-files[c(1:7,9,10)]
  omics<-gsub(".txt","",gsub(".+er/","",files))
  nodel<-list()
  edgel<-list()
  diseases<-c()
  for(a in 1:length(files)){
    # a<-1
    tt<-fread(files[a],stringsAsFactors = F,sep="\t")
    nodes<-unique(c(tt[[1]],tt[[2]])) ; if(length(which(is.na(nodes)))>0){nodes<-nodes[-which(is.na(nodes))]}
    nodel[[omics[a]]]<-nodes
    diseases<-c(diseases,nodes)
    print(paste(omics[a],length(nodes)))
    edges<-paste(tt[[1]],tt[[2]],sep="_") ; if(length(which(is.na(tt[[2]])))>0){edges<-edges[-which(is.na(tt[[2]]))]}
    edgel[[omics[a]]]<-edges
  }
  diseases<-sort(table(diseases),decreasing = T)
  as.data.frame(diseases[which(diseases==9)])
  input<-nodel
  tabla<-fromList(input)
  tabla<-tabla[c(8,2,3,7,1,4,6,9,5)][9:1]
  pdf(file = "Plots/Upsetplot_shared_nodes.pdf",width = 12,height = 8)
  upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  dev.off()
  ## Overlaps in terms of similarities between the 5 main layers ##
  ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
  ## Add epidemiology ##
  tt<-fread("Epidemiology/ICD10_epidemiology_denmark.txt",stringsAsFactors = F,sep="\t")
  nodel$epidemiology<-unique(c(tt$disease1,tt$disease2))
  tt<-tt[which(tt$RR>1.01)]
  edgel$epidemiology<-paste(tt[[1]],tt[[2]],sep="_")
  nodel5big<-list("epidemiology"=nodel$epidemiology,"symptom"=nodel$symptom,"drug"=nodel$drug,"gene"=nodel$gene,"miRNA"=nodel$miRNA)
  ## Upset nodes 5big ##
  input<-nodel5big
  tabla<-fromList(input)
  tabla<-tabla[,5:1]
  pdf(file = "Plots/Upsetplot_shared_nodes_big5.pdf",width = 12,height = 8)
  upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  dev.off()
  commonnodes5<-Reduce(intersect,nodel5big)
  edgel5big<-list()
  for(a in 1:length(names(nodel5big))){
    # a<-1
    red<-edgel[[names(nodel5big)[a]]]
    ind1<-c() ; ind2<-c()
    for(b in 1:length(commonnodes5)){
      # b<-1
      ind1<-c(ind1,grep(paste(commonnodes5[b],"_",sep=""),red)) ; ind2<-c(ind2,grep(paste("_",commonnodes5[b],sep=""),red))
    }
    edgel5big[[names(nodel5big)[a]]]<-red[intersect(ind1,ind2)]
  }
  ## Upset edges 5big ##
  input<-edgel5big
  tabla<-fromList(input)
  tabla<-tabla[c(1,4,3,5,2)]
  pdf(file = "Plots/Upsetplot_shared_edges_big5.pdf",width = 12,height = 8)
  upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  dev.off()
  
  ## Overlaps in terms of similarities between the 5 main layers + transcriptomics ##
  ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
  # nodel6big<-list("epidemiology"=nodel$epidemiology,"symptom"=nodel$symptom,"drug"=nodel$drug,"gene"=nodel$gene,"mirna"=nodel$mirna,"transcriptome"=nodel$transcriptome)
  # ## Upset nodes 6big ##
  # input<-nodel6big
  # tabla<-fromList(input)
  # tabla<-tabla[,6:1]
  # pdf(file = "Plots/Upsetplot_shared_nodes_big5_transcriptomics.pdf",width = 12,height = 8)
  # upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  # dev.off()
  # commonnodes5<-Reduce(intersect,nodel6big)
  # edgel6big<-list()
  # for(a in 1:length(names(nodel6big))){
  #   # a<-1
  #   red<-edgel[[names(nodel6big)[a]]]
  #   ind1<-c() ; ind2<-c()
  #   for(b in 1:length(commonnodes5)){
  #     # b<-1
  #     ind1<-c(ind1,grep(paste(commonnodes5[b],"_",sep=""),red)) ; ind2<-c(ind2,grep(paste("_",commonnodes5[b],sep=""),red))
  #   }
  #   edgel6big[[names(nodel6big)[a]]]<-red[intersect(ind1,ind2)]
  # }
  # ## Upset edges 6big ##
  # input<-edgel6big
  # tabla<-fromList(input)
  # tabla<-tabla[c(1,3,4,6,5,2)]
  # pdf(file = "Plots/Upsetplot_shared_edges_big5_transcriptomics.pdf",width = 12,height = 8)
  # upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  # dev.off()
  # tabla<-tabla[which(tabla$epidemiology==1),c(4,3,2,5,6,1)]
  # pdf(file = "Plots/Upsetplot_epidemiology_edges_big5_transcriptomics.pdf",width = 12,height = 8)
  # upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  # dev.off()
  
  ## Overlaps in terms of similarities between the epidemiology & transcriptomics ##
  ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
  nodel2<-list("epidemiology"=nodel$epidemiology,"transcriptome"=nodel$transcriptome)
  commonnodes2<-Reduce(intersect,nodel2)
  edgel2<-list()
  for(a in 1:length(names(edgel))){
    # a<-1
    red<-edgel[[names(edgel)[a]]]
    ind1<-c() ; ind2<-c()
    for(b in 1:length(commonnodes2)){
      # b<-1
      ind1<-c(ind1,grep(paste(commonnodes2[b],"_",sep=""),red)) ; ind2<-c(ind2,grep(paste("_",commonnodes2[b],sep=""),red))
    }
    edgel2[[names(edgel)[a]]]<-red[intersect(ind1,ind2)]
  }
  ## Upset edges 5big ##
  input<-edgel2
  tabla<-fromList(input)
  tabla<-tabla[which(tabla$epidemiology==1),]
  tabla<-tabla[c(5,8,6,2,3,4,1,7)]
  pdf(file = "Plots/Upsetplot_epidemiology_edges_transcriptomics_all.pdf",width = 12,height = 8)
  upset(tabla,sets=colnames(tabla),order.by = "freq",nsets = length(colnames(tabla)),keep.order=T,set_size.show = TRUE)
  dev.off()
  ## barplot with the number of nodes in the different layers ##
  nodenumber<-c()
  for(a in 1:length((names(nodel)))){
    # a<-1
    nodenumber<-c(nodenumber,length(intersect(nodel[[a]],commonnodes2)))
  }
  percentage<-(nodenumber/length(commonnodes2))*100
  data<-as.data.frame(cbind(percentage,nodenumber))
  data<-cbind(names(nodel),data)
  data<-data[c(1,2,8,5:3,7,9,6),]
  tsimply<-data$nodenumber ; names(tsimply)<-data$`names(nodel)`
  pdf(file="Plots/Barplot_nodes_in_transcriptomic_epidemiology.pdf",width = 12,height = 12)
  p <- barplot(tsimply,col="black",ylim=c(0,100),main="",las=2)
  text(x = p, y = tsimply + 7, labels = data$nodenumber,srt=90,cex = 4)
  dev.off()
  
  
  
  ## Networks for common diseases ##
  commondis<-names(diseases[which(diseases==8)])
  files<-c("Similarities/Similarities_with_fdr/ICD10_symptom_fs_weightedjaccardindex.txt",
           "Similarities/Similarities_with_fdr/ICD10_drug_fs_weightedjaccardindex.txt",
           "Similarities/Similarities_with_fdr/ICD10_gene_fs_weightedjaccardindex.txt",
           "Similarities/Similarities_with_fdr/ICD10_miRNA_fs_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_metabolome_fs_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_microbiome_fsd_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_transcriptome_fc_cosinesimilarities.txt",
           "Similarities/Similarities_with_fdr/ICD10_methylome_ifs_cosinesimilarities.txt")
  
  sharedmultilayer<-c()
  for(a in 1:length(files)){
    # a<-1
    tt<-fread(files[a],stringsAsFactors = F,sep="\t")
    tt<-tt[which(tt[[5]]<=0.05)]
    tt<-tt[which(tt[[3]]>0)]
    ind1<-c() ; ind2<-c()
    for(b in 1:length(commondis)){
      # b<-1
      ind1<-c(ind1,which(tt[[1]]==commondis[b])) ; ind2<-c(ind2,which(tt[[2]]==commondis[b]))
    }
    comunes<-intersect(ind1,ind2)
    if(length(comunes)>0){
      subnet<-cbind(tt[comunes,1:3],gsub("_.+","",gsub("Similarities/Similarities_with_fdr/ICD10_","",files[a])))
      colnames(subnet)<-c("ICD10_1","ICD10_2","weight","layer")
      sharedmultilayer<-rbind(sharedmultilayer,subnet)
    }
  }
  write.table(sharedmultilayer,"Plots/multilayer_shared_nodes.txt",quote=F,sep="\t",row.names=F)
  
  sharedmultilayer2<-cbind(sharedmultilayer,sharedmultilayer$layer)
  colnames(sharedmultilayer2)<-c("SourceNode","TargetNode","weight","SourceLayer","TargetLayer")
  write.table(sharedmultilayer2,"Plots/multilayer_shared_nodes2.txt",quote=F,sep="\t",row.names=F)
  
  
  # epi<-fread("Epidemiology/ICD10_epidemiology_denmark.txt",stringsAsFactors = F)
  # epi<-epi[which(epi$RR>1.01)]
  # ind1<-c() ; ind2<-c()
  # for(b in 1:length(commondis)){
  #   # b<-1
  #   ind1<-c(ind1,which(epi[[1]]==commondis[b])) ; ind2<-c(ind2,which(epi[[2]]==commondis[b]))
  # }
  # comunes<-intersect(ind1,ind2)
  # epi<-epi[comunes]
  # for(a in 1:length(sharedmultilayer$ICD10_1)){
  #   a<-1
  #   
  # }
  
  
}

if(args[1]=="Dervic_comorbidities"){
  gexf_file <- "Austria_Epidemiology/4.Graphs-gexffiles/Graph_Male_ICD_Year_2003.gexf"  # replace with your actual file path
  g <- read.gexf(gexf_file)
}

## Department seminar example ##
# mm<-fread("WorkFiles/ICD10_Omics/ICD10_microbiome.txt",stringsAsFactors = F,sep="\t")
# setkey(mm,"icd10")
# F20
# G20
# 
# tt<-rbind(cbind("F20",unique(intersect(mm["F20"]$microbiome,mm["G20"]$microbiome))),
#           cbind("G20",unique(intersect(mm["F20"]$microbiome,mm["G20"]$microbiome))))
# 
# write.table(tt,"Network_microbioma_F20_G20.txt",quote=F,sep="\t",row.names=F)
# 
# 
# f20u<-mm["F20"]$microbiome[which(mm["F20"]$fsd>0)]
# f20d<-mm["F20"]$microbiome[which(mm["F20"]$fsd<0)]
# 
# g20u<-mm["G20"]$microbiome[which(mm["G20"]$fsd>0)]
# g20d<-mm["G20"]$microbiome[which(mm["G20"]$fsd<0)]
# 
# 
# intersect(f20u,g20u)
# intersect(f20d,g20d)
# 
# 
# unique(intersect(mm["F20"]$microbiome,mm["G20"]$microbiome))






