#### With this code we are going to generate the multiomic layers (disease-disease similarities) for the BSC-IGC collaboration ####
## Load the needed libraries ##
library("data.table")
# library(VennDiagram)
# library(RColorBrewer)
# library("igraph")
# library("UpSetR")
# library("dendextend")
library(lsa)
library("geometry")
# library("ggplot2")
# library("gplots")
# library(GEOquery)
# library(metapro)

args = commandArgs(trailingOnly=TRUE)

# print(args[1])
# print(args[2])
## Creathe the needed files ##

if("Summaries"%in%list.files()==F){dir.create("Summaries")}
if("Interlayer"%in%list.files("Similarities/")==F){dir.create("Similarities/Interlayer")}
if("Intermediate"%in%list.files("Data//")==F){dir.create("Data/Intermediate")}
if("OriginalInterlayerAssociations"%in%list.files("Similarities/")==FALSE){dir.create("Similarities/OriginalInterlayerAssociations")}
if("RandomInterlayerAssociations"%in%list.files("Similarities/")==FALSE){dir.create("Similarities/RandomInterlayerAssociations")}
#### Build all the interlayers ####
## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
if(args[1]=="buildinterlayers"){
  #### Connect miRNA to gene ####
  if(args[2]=="miRNAgene"){
    if(args[3]=="preparetables"){
      ## Read files ##
      disease_mirna<-fread("WorkFiles/ICD10_Omics/ICD10_miRNA.txt",stringsAsFactors = F,sep="\t") ; disease_mirna$mirna<-tolower(disease_mirna$mirna)
      disease_gene<-fread("WorkFiles/ICD10_Omics/ICD10_gene.txt",stringsAsFactors = F,sep="\t") ; disease_gene$gene<-toupper(disease_gene$gene)
      mirna_gene_associations<-fread("Data/Interlayers/miRTarBase_MTI.csv",stringsAsFactors = F,sep=";") ; mirna_gene_associations$miRNA<-tolower(mirna_gene_associations$miRNA) ; mirna_gene_associations$`Target Gene`<-toupper(mirna_gene_associations$`Target Gene`)
      ## Select only interactions between human miRNAs and genes ##
      mirna_gene_associations<-mirna_gene_associations[which(mirna_gene_associations$`Species (miRNA)`=="Homo sapiens")] ; mirna_gene_associations<-mirna_gene_associations[which(mirna_gene_associations$`Species (Target Gene)`=="Homo sapiens")]
      mirna_gene_associations<-mirna_gene_associations[,c(2,4)]
      colnames(mirna_gene_associations)<-c("mirna","gene")
      beggining<-c(length(unique(disease_mirna$icd10)),length(unique(disease_mirna$mirna)),length(disease_mirna$icd10),
                   length(unique(disease_gene$icd10)),length(unique(disease_gene$gene)),length(disease_gene$icd10),
                   length(mirna_gene_associations$mirna))
      ## Get the miRNAs common to both files
      mirnas_comunes<-intersect(disease_mirna$mirna,mirna_gene_associations$mirna)
      ## Only 434 out of 1,465 miRNAs in the disease network have an interaction with at least one gene ## 29.34%
      length(intersect(disease_mirna$mirna,mirna_gene_associations$mirna))/length(unique(disease_mirna$mirna))
      setkey(mirna_gene_associations,"mirna")
      mirna_gene_associations<-mirna_gene_associations[mirnas_comunes]
      ## Get the genes common to both files
      genes_comunes<-intersect(disease_gene$gene,mirna_gene_associations$gene)
      ## Only 6096 out of 7,360 genes in the disease network have an interaction with at least one miRNA ## 83.13%
      length(intersect(disease_gene$gene,mirna_gene_associations$gene))/length(unique(disease_gene$gene))
      setkey(mirna_gene_associations,"gene")
      mirna_gene_associations<-mirna_gene_associations[genes_comunes]
      write.table(mirna_gene_associations,"WorkFiles/ICD10_Omics/mirna_gene_associations.txt",quote=F,sep="\t",row.names=F)
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
      write.table(disease_mirna,"Data/Intermediate/mirna_gene_ICD10miRNA.txt",quote=F,sep="\t",row.names=F)
      write.table(disease_gene,"Data/Intermediate/mirna_gene_ICD10gene.txt",quote=F,sep="\t",row.names=F)
      mirgene<-matrix(ncol = length(unique(mirna_gene_associations$mirna)),nrow = length(unique(mirna_gene_associations$gene)),0)
      colnames(mirgene)<-unique(mirna_gene_associations$mirna) ; rownames(mirgene)<-unique(mirna_gene_associations$gene)
      for(a in 1:length(mirna_gene_associations$mirna)){mirgene[mirna_gene_associations$gene[a],mirna_gene_associations$mirna[a]]<-1}
      write.table(mirgene,"Data/Intermediate/mirna_gene_variable_associations.txt",quote=F,sep="\t")
    }
    if(args[3]=="original"){
      ## load miRNA-gene associations and transform them into a matrix ##
      mirna_gene_associations<-fread("WorkFiles/ICD10_Omics/mirna_gene_associations.txt",stringsAsFactors = F,sep="\t")
      mirgene<-matrix(ncol = length(unique(mirna_gene_associations$mirna)),nrow = length(unique(mirna_gene_associations$gene)),0)
      colnames(mirgene)<-unique(mirna_gene_associations$mirna) ; rownames(mirgene)<-unique(mirna_gene_associations$gene)
      for(a in 1:length(mirna_gene_associations$mirna)){mirgene[mirna_gene_associations$gene[a],mirna_gene_associations$mirna[a]]<-1}
      ## load disease-mirna and disease-gene associations ##
      disease_mirna<-fread("Data/Intermediate/mirna_gene_ICD10miRNA.txt",stringsAsFactors = F,sep="\t")
      disease_gene<-fread("Data/Intermediate/mirna_gene_ICD10gene.txt",stringsAsFactors = F,sep="\t")
      setkey(disease_mirna,"icd10") ; setkey(disease_gene,"icd10")
      genediseases<-unique(disease_gene$icd10) ; mirnadiseases<-unique(disease_mirna$icd10)
      ## Create a table where we are going to annotate all the results ##
      ncompas<-length(genediseases)*length(mirnadiseases)
      thetable<-matrix(ncol = 3,nrow = ncompas, NA)
      ## Calculate the score ##
      ini<-Sys.time()
      contador<-c()
      for(a in 1:length(mirnadiseases)){
        # a<-1
        mirna<-unique(disease_mirna[mirnadiseases[a]]$mirna)
        for(b in 1:length(genediseases)){
          # b<-1
          contador<-sum(contador,1)
          genes<-unique(disease_gene[genediseases[b]]$gene)
          thetable[contador,]<-c(mirnadiseases[a],genediseases[b],sum(mirgene[genes,mirna]))
        }
        print(paste(round((a/length(mirnadiseases))*100,2),"%",sep=""))
      }
      fin<-Sys.time()
      fin-ini
      colnames(thetable)<-c("icd10_miRNA","icd10_gene","shared_interactions")
      write.table(thetable,"Similarities/OriginalInterlayerAssociations/ICD10_mirna_gene_associations.txt",quote=F,sep="\t",row.names=F)
      ## Create the directory for the random associations ##
      if("ICD10_mirna_gene_associations"%in%list.files("Similarities/RandomInterlayerAssociations/")==FALSE){dir.create("Similarities/RandomInterlayerAssociations/ICD10_mirna_gene_associations")}
    }
    if(args[3]=="random"){
      # input<-"mirna_gene_variable_associations.txt"
      # layer1<-"mirna_gene_ICD10miRNA.txt"
      # layer2<-"mirna_gene_ICD10gene.txt"
      # saveoutput<-"ICD10_mirna_gene_associations"
      # semilla<-1
      # randomnum<-100
      ## Define the function ##
      random_distances<-function(input,layer1,layer2,saveoutput,semilla,randomnum){
        layer1name<-gsub(".txt","",gsub(".+ICD10","",layer1)) ; layer2name<-gsub(".txt","",gsub(".+ICD10","",layer2))
        set.seed(semilla)
        l1_l2_associations<-fread(paste("WorkFiles/ICD10_Omics/",layer1name,"_",layer2name,"_associations.txt",sep=""),stringsAsFactors = F,sep="\t")
        varmat2<-matrix(ncol = length(unique(l1_l2_associations[[1]])),nrow = length(unique(l1_l2_associations[[2]])),0)
        colnames(varmat2)<-unique(l1_l2_associations[[1]]) ; rownames(varmat2)<-unique(l1_l2_associations[[2]])
        for(a in 1:length(l1_l2_associations[[1]])){varmat2[l1_l2_associations[[2]][a],l1_l2_associations[[1]][a]]<-1}
        for(a in 1:randomnum){
          ini<-Sys.time()
          ## Sample variable names ##
          rvarmat<-varmat2
          colnames(rvarmat)<-sample(colnames(rvarmat)) ; rownames(rvarmat)<-sample(rownames(rvarmat))
          
          ## load disease-mirna and disease-gene associations ##
          disease_layer1<-fread(paste("Data/Intermediate/",layer1,sep=""),stringsAsFactors = F,sep="\t")
          disease_layer2<-fread(paste("Data/Intermediate/",layer2,sep=""),stringsAsFactors = F,sep="\t")
          setkey(disease_layer1,"icd10") ; setkey(disease_layer2,"icd10")
          layer2diseases<-unique(disease_layer2$icd10) ; layer1diseases<-unique(disease_layer1$icd10)
          
          ## Create a table where we are going to annotate all the results ##
          ncompas<-length(layer1diseases)*length(layer2diseases)
          rthetable<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
          ## Calculate the score ##
          ini<-Sys.time()
          contador<-c()
          for(d in 1:length(layer1diseases)){
            # d<-1
            layer1v<-unique(disease_layer1[layer1diseases[d]][[2]])
            for(b in 1:length(layer2diseases)){
              # b<-1
              contador<-sum(contador,1)
              layer2v<-unique(disease_layer2[layer2diseases[b]][[2]])
              rthetable[contador,]<-c(layer1diseases[d],layer2diseases[b],sum(rvarmat[layer2v,layer1v]))
            }
            print(paste(round((d/length(layer1diseases))*100,2),"%",sep=""))
          }
          colnames(rthetable)<-c(paste("icd10_",layer1name,sep=""),paste("icd10_",layer2name,sep=""),"shared_interactions")
          write.table(rthetable,paste("Similarities/RandomInterlayerAssociations/",saveoutput,"/Seed_",semilla,"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names=F)
          fin<-Sys.time()
          print(paste("Seed ",semilla," round ",a," time: ",fin-ini,sep=""))
        }
        final<-Sys.time()
        print(final-inicio)
      }
      random_distances(args[4],args[5],args[6],args[7],args[8],args[9])
    }
  }
  #### Connect genes to drug ####
  if(args[2]=="druggene"){
    if(args[3]=="preparetables"){
      ## Read files ##
      disease_drug<-fread("WorkFiles/ICD10_Omics/ICD10_drug.txt",stringsAsFactors = F,sep="\t") ; disease_drug$drug<-tolower(disease_drug$drug)
      disease_gene<-fread("WorkFiles/ICD10_Omics/ICD10_gene.txt",stringsAsFactors = F,sep="\t") ; disease_gene$gene<-toupper(disease_gene$gene)
      drug_gene_associations<-fread("Data/Interlayers/CTD_chem_gene_ixns.csv",stringsAsFactors = F) ; drug_gene_associations<-drug_gene_associations[which(drug_gene_associations$Organism=="Homo sapiens")]
      drug_gene_associations<-drug_gene_associations[,c(1,4,10,11,9)]
      drug_gene_associations$ChemicalName<-tolower(drug_gene_associations$ChemicalName) ; drug_gene_associations$GeneSymbol<-toupper(drug_gene_associations$GeneSymbol)
      colnames(drug_gene_associations)<-c("drug","gene","interaction","pmid","description")
      beggining<-c(length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
                   length(unique(disease_gene$icd10)),length(unique(disease_gene$gene)),length(disease_gene$icd10),
                   length(drug_gene_associations$drug))
      
      ## Get the drugs common to both files
      drugs_comunes<-intersect(disease_drug$drug,drug_gene_associations$drug)
      ## Only 3,196 out of 5,351 drugs in the disease network have an interaction with at least one gene ## 57,443%
      length(intersect(disease_drug$drug,drug_gene_associations$drug))/length(unique(disease_drug$drug))
      setkey(drug_gene_associations,"drug")
      drug_gene_associations<-drug_gene_associations[drugs_comunes]
      ## Get the genes common to both files
      genes_comunes<-intersect(disease_gene$gene,drug_gene_associations$gene)
      ## Only 7,275 out of 7,336 genes in the disease network have an interaction with at least one drug ## 99,16%
      length(intersect(disease_gene$gene,drug_gene_associations$gene))/length(unique(disease_gene$gene))
      setkey(drug_gene_associations,"gene")
      drug_gene_associations<-drug_gene_associations[genes_comunes]
      write.table(drug_gene_associations[,1:2],"WorkFiles/ICD10_Omics/drug_gene_associations.txt",quote=F,sep="\t",row.names=F)
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
      ## Save the files that will be used to calculate the interlayer connections
      write.table(disease_drug,"Data/Intermediate/drug_gene_ICD10drug.txt",quote=F,sep="\t",row.names=F)
      write.table(disease_gene,"Data/Intermediate/drug_gene_ICD10gene.txt",quote=F,sep="\t",row.names=F)
      drugene<-matrix(ncol = length(unique(drug_gene_associations$drug)),nrow = length(unique(drug_gene_associations$gene)),0)
      colnames(drugene)<-unique(drug_gene_associations$drug) ; rownames(drugene)<-unique(drug_gene_associations$gene)
      for(a in 1:length(drug_gene_associations$drug)){drugene[drug_gene_associations$gene[a],drug_gene_associations$drug[a]]<-1}
      write.table(drugene,"Data/Intermediate/drug_gene_variable_associations.txt",quote=F,sep="\t")
    }
    if(args[3]=="original"){
      ## load miRNA-gene associations and transform them into a matrix ##
      drug_gene_associations<-fread("WorkFiles/ICD10_Omics/drug_gene_associations.txt",stringsAsFactors = F,sep="\t")
      drugene<-matrix(ncol = length(unique(drug_gene_associations$drug)),nrow = length(unique(drug_gene_associations$gene)),0)
      colnames(drugene)<-unique(drug_gene_associations$drug) ; rownames(drugene)<-unique(drug_gene_associations$gene)
      for(a in 1:length(drug_gene_associations$drug)){drugene[drug_gene_associations$gene[a],drug_gene_associations$drug[a]]<-1}
      ## load disease-drug and disease-gene associations ##
      disease_drug<-fread("Data/Intermediate/drug_gene_ICD10drug.txt",stringsAsFactors = F,sep="\t")
      disease_gene<-fread("Data/Intermediate/drug_gene_ICD10gene.txt",stringsAsFactors = F,sep="\t")
      setkey(disease_drug,"icd10") ; setkey(disease_gene,"icd10")
      genediseases<-unique(disease_gene$icd10) ; drugdiseases<-unique(disease_drug$icd10)
      ## Create a table where we are going to annotate all the results ##
      ncompas<-length(genediseases)*length(drugdiseases)
      thetable<-matrix(ncol = 3,nrow = ncompas, NA)
      ## Calculate the score ##
      ini<-Sys.time()
      contador<-c()
      for(a in 1:length(drugdiseases)){
        # a<-1
        drug<-unique(disease_drug[drugdiseases[a]]$drug)
        for(b in 1:length(genediseases)){
          # b<-1
          contador<-sum(contador,1)
          genes<-unique(disease_gene[genediseases[b]]$gene)
          thetable[contador,]<-c(drugdiseases[a],genediseases[b],sum(drugene[genes,drug]))
        }
        print(paste(round((a/length(drugdiseases))*100,2),"%",sep=""))
      }
      fin<-Sys.time()
      fin-ini
      colnames(thetable)<-c("icd10_drug","icd10_gene","shared_interactions")
      write.table(thetable,"Similarities/OriginalInterlayerAssociations/ICD10_drug_gene_associations.txt",quote=F,sep="\t",row.names=F)
      ## Create the directory for the random associations ##
      if("ICD10_drug_gene_associations"%in%list.files("Similarities/RandomInterlayerAssociations/")==FALSE){dir.create("Similarities/RandomInterlayerAssociations/ICD10_drug_gene_associations")}
    }
    if(args[3]=="random"){
      # input<-"drug_gene_variable_associations.txt"
      # layer1<-"drug_gene_ICD10drug.txt"
      # layer2<-"drug_gene_ICD10gene.txt"
      # saveoutput<-"ICD10_drug_gene_associations"
      # semilla<-1
      # randomnum<-100
      ## Define the function ##
      random_distances<-function(input,layer1,layer2,saveoutput,semilla,randomnum){
        layer1name<-gsub(".txt","",gsub(".+ICD10","",layer1)) ; layer2name<-gsub(".txt","",gsub(".+ICD10","",layer2))
        set.seed(semilla)
        l1_l2_associations<-fread(paste("WorkFiles/ICD10_Omics/",layer1name,"_",layer2name,"_associations.txt",sep=""),stringsAsFactors = F,sep="\t")
        varmat2<-matrix(ncol = length(unique(l1_l2_associations[[1]])),nrow = length(unique(l1_l2_associations[[2]])),0)
        colnames(varmat2)<-unique(l1_l2_associations[[1]]) ; rownames(varmat2)<-unique(l1_l2_associations[[2]])
        for(a in 1:length(l1_l2_associations[[1]])){varmat2[l1_l2_associations[[2]][a],l1_l2_associations[[1]][a]]<-1}
        for(a in 1:randomnum){
          ini<-Sys.time()
          ## Sample variable names ##
          rvarmat<-varmat2
          colnames(rvarmat)<-sample(colnames(rvarmat)) ; rownames(rvarmat)<-sample(rownames(rvarmat))
          
          ## load disease-drug and disease-gene associations ##
          disease_layer1<-fread(paste("Data/Intermediate/",layer1,sep=""),stringsAsFactors = F,sep="\t")
          disease_layer2<-fread(paste("Data/Intermediate/",layer2,sep=""),stringsAsFactors = F,sep="\t")
          setkey(disease_layer1,"icd10") ; setkey(disease_layer2,"icd10")
          layer2diseases<-unique(disease_layer2$icd10) ; layer1diseases<-unique(disease_layer1$icd10)
          
          ## Create a table where we are going to annotate all the results ##
          ncompas<-length(layer1diseases)*length(layer2diseases)
          rthetable<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
          ## Calculate the score ##
          ini<-Sys.time()
          contador<-c()
          for(d in 1:length(layer1diseases)){
            # d<-1
            layer1v<-unique(disease_layer1[layer1diseases[d]][[2]])
            for(b in 1:length(layer2diseases)){
              # b<-1
              contador<-sum(contador,1)
              layer2v<-unique(disease_layer2[layer2diseases[b]][[2]])
              rthetable[contador,]<-c(layer1diseases[d],layer2diseases[b],sum(rvarmat[layer2v,layer1v]))
            }
            print(paste(round((d/length(layer1diseases))*100,2),"%",sep=""))
          }
          colnames(rthetable)<-c(paste("icd10_",layer1name,sep=""),paste("icd10_",layer2name,sep=""),"shared_interactions")
          write.table(rthetable,paste("Similarities/RandomInterlayerAssociations/",saveoutput,"/Seed_",semilla,"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names=F)
          fin<-Sys.time()
          print(paste("Seed ",semilla," round ",a," time: ",fin-ini,sep=""))
        }
        final<-Sys.time()
        print(final-inicio)
      }
      random_distances(args[4],args[5],args[6],args[7],args[8],args[9])
    }
  }
  #### Connect mirna to drug ####
  if(args[2]=="miRNAdrug"){
    if(args[3]=="preparetables"){
      ## Read files ##
      disease_mirna<-fread("WorkFiles/ICD10_Omics/ICD10_miRNA.txt",stringsAsFactors = F,sep="\t") ; disease_mirna$mirna<-toupper(disease_mirna$mirna)
      disease_drug<-fread("WorkFiles/ICD10_Omics/ICD10_drug.txt",stringsAsFactors = F,sep="\t") ; disease_drug$drug<-tolower(disease_drug$drug)
      ## CTD ##
      ## @@@ ##
      drug_mirna_associations<-fread("Data/Interlayers/CTD_chem_gene_ixns.csv",stringsAsFactors = F,sep=",") ; drug_mirna_associations<-drug_mirna_associations[which(drug_mirna_associations$Organism=="Homo sapiens")]
      drug_mirna_associations<-drug_mirna_associations[,c(1,4,10,11,9)]
      colnames(drug_mirna_associations)<-c("drug","mirna","interaction","pmid","description")
      ## Select the drug-mirna associations ##
      drug_mirna_ctd_associations<-drug_mirna_associations[grep("mir",drug_mirna_associations$mirna,ignore.case = T)]
      drug_mirna_ctd_associations$mirna<-paste("HSA-MIR-",gsub("MIR","",drug_mirna_ctd_associations$mirna),sep="")
      drug_mirna_ctd_associations$drug<-tolower(drug_mirna_ctd_associations$drug)
      ## SM2miR3 ##
      ## @ @ @ @ ##
      drug_mirna_associations<-fread("Data/Interlayers/SM2miR3.csv",stringsAsFactors = F,sep=";")
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
      drug_mirna_associations$mirna<-toupper(drug_mirna_associations$mirna)
      if(length(which(duplicated(drug_mirna_associations)))>0){drug_mirna_associations<-drug_mirna_associations[-which(duplicated(drug_mirna_associations)),]}
      beggining<-c(length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
                   length(unique(disease_mirna$icd10)),length(unique(disease_mirna$mirna)),length(disease_mirna$icd10),
                   length(drug_mirna_associations$drug))
      
      ## Get the drugs common to both files
      drugs_comunes<-intersect(disease_drug$drug,drug_mirna_associations$drug)
      ## Only 258 out of 5,564 drugs in the disease network have an interaction with at least one gene ## 4,63%
      length(intersect(disease_drug$drug,drug_mirna_associations$drug))/length(unique(disease_drug$drug))
      setkey(drug_mirna_associations,"drug")
      drug_mirna_associations<-drug_mirna_associations[drugs_comunes]
      ## Get the genes common to both files
      mirnas_comunes<-intersect(disease_mirna$mirna,drug_mirna_associations$mirna)
      ## Only 765 out of 1,479 mirnas in the disease network have an interaction with at least one drug ## 51,72%
      length(intersect(disease_mirna$mirna,drug_mirna_associations$mirna))/length(unique(disease_mirna$mirna))
      setkey(drug_mirna_associations,"mirna")
      drug_mirna_associations<-drug_mirna_associations[mirnas_comunes]
      write.table(drug_mirna_associations,"WorkFiles/ICD10_Omics/drug_miRNA_associations.txt",quote=F,sep="\t",row.names=F)
      ## Get the diseases associated to mirnas and drugs that have at least one association to each other (mirna - drug association)
      commonmirnas<-unique(drug_mirna_associations$mirna) ; commondrugs<-unique(drug_mirna_associations$drug)
      setkey(disease_drug,"drug") ; setkey(disease_mirna,"mirna")
      disease_mirna<-disease_mirna[commonmirnas] ; disease_drug<-disease_drug[commondrugs]
      ## At the end we have interactions between 430 drugs and 3,818 mirnas ##
      ending<-c(length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
                length(unique(disease_mirna$icd10)),length(unique(disease_mirna$mirna)),length(disease_mirna$icd10),
                length(drug_mirna_associations$drug))
      association_summary<-cbind(beggining,ending)
      rownames(association_summary)<-c("ICD10(drug)","drugs","ICD10-drug associations","ICD10(miRNA)","miRNAs","ICD10-miRNA associations","drug-miRNA associations")
      association_summary<-cbind(association_summary,(association_summary[,2]/association_summary[,1])*100)
      colnames(association_summary)<-c("original","post_filtering","percentage_kept")
      write.table(association_summary,"Summaries/Summary_drug_miRNA_associations.txt",quote=F,sep="\t")
      ## Save the files that will be used to calculate the interlayer connections
      write.table(disease_drug,"Data/Intermediate/drug_mirna_ICD10drug.txt",quote=F,sep="\t",row.names=F)
      write.table(disease_mirna,"Data/Intermediate/drug_mirna_ICD10miRNA.txt",quote=F,sep="\t",row.names=F)
      drugmirna<-matrix(ncol = length(unique(drug_mirna_associations$drug)),nrow = length(unique(drug_mirna_associations$mirna)),0)
      colnames(drugmirna)<-unique(drug_mirna_associations$drug) ; rownames(drugmirna)<-unique(drug_mirna_associations$mirna)
      for(a in 1:length(drug_mirna_associations$drug)){drugmirna[drug_mirna_associations$mirna[a],drug_mirna_associations$drug[a]]<-1}
      write.table(drugmirna,"Data/Intermediate/drug_mirna_variable_associations.txt",quote=F,sep="\t")
    }
    if(args[3]=="original"){
      ## load miRNA-gene associations and transform them into a matrix ##
      drug_mirna_associations<-fread("WorkFiles/ICD10_Omics/drug_miRNA_associations.txt",stringsAsFactors = F,sep="\t")
      drugmirna<-matrix(ncol = length(unique(drug_mirna_associations$drug)),nrow = length(unique(drug_mirna_associations$mirna)),0)
      colnames(drugmirna)<-unique(drug_mirna_associations$drug) ; rownames(drugmirna)<-unique(drug_mirna_associations$mirna)
      for(a in 1:length(drug_mirna_associations$drug)){drugmirna[drug_mirna_associations$mirna[a],drug_mirna_associations$drug[a]]<-1}
      ## load disease-drug and disease-gene associations ##
      disease_drug<-fread("Data/Intermediate/drug_mirna_ICD10drug.txt",stringsAsFactors = F,sep="\t")
      disease_mirna<-fread("Data/Intermediate/drug_mirna_ICD10miRNA.txt",stringsAsFactors = F,sep="\t")
      setkey(disease_drug,"icd10") ; setkey(disease_mirna,"icd10")
      mirnadiseases<-unique(disease_mirna$icd10) ; drugdiseases<-unique(disease_drug$icd10)
      ## Create a table where we are going to annotate all the results ##
      ncompas<-length(mirnadiseases)*length(drugdiseases)
      thetable<-matrix(ncol = 3,nrow = ncompas, NA)
      ## Calculate the score ##
      ini<-Sys.time()
      contador<-c()
      for(a in 1:length(drugdiseases)){
        # a<-1
        drug<-unique(disease_drug[drugdiseases[a]]$drug)
        for(b in 1:length(mirnadiseases)){
          # b<-1
          contador<-sum(contador,1)
          mirnas<-unique(disease_mirna[mirnadiseases[b]]$mirna)
          thetable[contador,]<-c(drugdiseases[a],mirnadiseases[b],sum(drugmirna[mirnas,drug]))
        }
        print(paste(round((a/length(drugdiseases))*100,2),"%",sep=""))
      }
      fin<-Sys.time()
      fin-ini
      colnames(thetable)<-c("icd10_drug","icd10_mirna","shared_interactions")
      write.table(thetable,"Similarities/OriginalInterlayerAssociations/ICD10_drug_mirna_associations.txt",quote=F,sep="\t",row.names=F)
      ## Create the directory for the random associations ##
      if("ICD10_drug_mirna_associations"%in%list.files("Similarities/RandomInterlayerAssociations/")==FALSE){dir.create("Similarities/RandomInterlayerAssociations/ICD10_drug_mirna_associations")}
    }
    if(args[3]=="random"){
      # input<-"drug_mirna_variable_associations.txt"
      # layer1<-"drug_mirna_ICD10drug.txt"
      # layer2<-"drug_mirna_ICD10miRNA.txt"
      # saveoutput<-"ICD10_drug_mirna_associations"
      # semilla<-1
      # randomnum<-100
      ## Define the function ##
      random_distances<-function(input,layer1,layer2,saveoutput,semilla,randomnum){
        layer1name<-gsub(".txt","",gsub(".+ICD10","",layer1)) ; layer2name<-gsub(".txt","",gsub(".+ICD10","",layer2))
        set.seed(semilla)
        l1_l2_associations<-fread(paste("WorkFiles/ICD10_Omics/",layer1name,"_",layer2name,"_associations.txt",sep=""),stringsAsFactors = F,sep="\t")
        varmat2<-matrix(ncol = length(unique(l1_l2_associations[[1]])),nrow = length(unique(l1_l2_associations[[2]])),0)
        colnames(varmat2)<-unique(l1_l2_associations[[1]]) ; rownames(varmat2)<-unique(l1_l2_associations[[2]])
        for(a in 1:length(l1_l2_associations[[1]])){varmat2[l1_l2_associations[[2]][a],l1_l2_associations[[1]][a]]<-1}
        for(a in 1:randomnum){
          ini<-Sys.time()
          ## Sample variable names ##
          rvarmat<-varmat2
          colnames(rvarmat)<-sample(colnames(rvarmat)) ; rownames(rvarmat)<-sample(rownames(rvarmat))
          
          ## load disease-drug and disease-gene associations ##
          disease_layer1<-fread(paste("Data/Intermediate/",layer1,sep=""),stringsAsFactors = F,sep="\t")
          disease_layer2<-fread(paste("Data/Intermediate/",layer2,sep=""),stringsAsFactors = F,sep="\t")
          setkey(disease_layer1,"icd10") ; setkey(disease_layer2,"icd10")
          layer2diseases<-unique(disease_layer2$icd10) ; layer1diseases<-unique(disease_layer1$icd10)
          
          ## Create a table where we are going to annotate all the results ##
          ncompas<-length(layer1diseases)*length(layer2diseases)
          rthetable<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
          ## Calculate the score ##
          ini<-Sys.time()
          contador<-c()
          for(d in 1:length(layer1diseases)){
            # d<-1
            layer1v<-unique(disease_layer1[layer1diseases[d]][[2]])
            for(b in 1:length(layer2diseases)){
              # b<-1
              contador<-sum(contador,1)
              layer2v<-unique(disease_layer2[layer2diseases[b]][[2]])
              rthetable[contador,]<-c(layer1diseases[d],layer2diseases[b],sum(rvarmat[layer1v,layer2v]))
            }
            print(paste(round((d/length(layer1diseases))*100,2),"%",sep=""))
          }
          colnames(rthetable)<-c(paste("icd10_",layer1name,sep=""),paste("icd10_",layer2name,sep=""),"shared_interactions")
          write.table(rthetable,paste("Similarities/RandomInterlayerAssociations/",saveoutput,"/Seed_",semilla,"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names=F)
          fin<-Sys.time()
          print(paste("Seed ",semilla," round ",a," time: ",fin-ini,sep=""))
        }
        final<-Sys.time()
        print(final-inicio)
      }
      random_distances(args[4],args[5],args[6],args[7],args[8],args[9])
    }
  }
  #### Connect genes to symptoms ####
  if(args[2]=="symptomgene"){
    if(args[3]=="preparetables"){
      if(length(which(list.files("Data/Interlayers/")=="SympGAN_symptom_gene_associations.txt"))==0){
        symptom<-fread("WorkFiles/ICD10_Omics/ICD10_symptom.txt",stringsAsFactors = F,sep="\t")
        symptoms<-unique(symptom$symptom)
        ## load symptom - gene associations ##
        sympgenec<-fread("Data/Interlayers/SympGAN/curated_symptom_gene_associations.tsv",stringsAsFactors = F)
        sympgenel<-fread("Data/Interlayers/SympGAN/literature_symptom_gene_associations.tsv",stringsAsFactors = F)
        ## get the CUI for each symptom ##
        mrconso<-fread("MRCONSO.RRF",stringsAsFactors = F,sep="|")
        mrconso$V15<-tolower(mrconso$V15)
        inds<-c()
        for(a in 1:length(symptoms)){
          cuales<-which(mrconso$V15==tolower(symptoms[a]))
          if(length(cuales)>0){inds<-c(inds,paste(unique(mrconso$V1[cuales]),collapse = ";"))}
          if(length(cuales)==0){inds<-c(inds,N)}
          print(paste(round((a/length(symptoms))*100,2),"%",sep=""))
        }
        cursympgen<-c()
        bibsympgen<-c()
        withcsymp<-c() ; withpsymp<-c()
        for(a in 1:length(inds)){
          # a<-1
          cui<-strsplit(inds[a],";")[[1]]
          indicesc<-c() ; indicesl<-c()
          for(b in 1:length(cui)){
            indicesc<-c(indicesc,which(sympgenec$Symptom_CUI==cui[b]))
            indicesl<-c(indicesl,which(sympgenel$Symptom_CUI==cui[b]))
          }
          if(length(indicesc)>0){cursympgen<-rbind(cursympgen,cbind(symptoms[a],sympgenec[unique(indicesc)]));withcsymp<-c(withcsymp,symptoms[a])}
          if(length(indicesl)>0){bibsympgen<-rbind(bibsympgen,cbind(symptoms[a],sympgenel[unique(indicesl)]));withpsymp<-c(withpsymp,symptoms[a])}
          print(paste(round((a/length(inds))*100,2),"%",sep=""))
        }
        if(length(which(duplicated(cursympgen[,c(1,4)])))>0){cursympgen<-cursympgen[-which(duplicated(cursympgen[,c(1,4)])),]}
        colnames(cursympgen)[1]<-"symptom_umls"
        write.table(cursympgen,"Data/Interlayers/SympGAN_symptom_gene_associations.txt",quote=F,sep="\t",row.names = F)
      }
      ## Read files ##
      disease_symptom<-fread("WorkFiles/ICD10_Omics/ICD10_symptom.txt",stringsAsFactors = F,sep="\t") ; disease_symptom$symptom<-tolower(disease_symptom$symptom)
      disease_gene<-fread("WorkFiles/ICD10_Omics/ICD10_gene.txt",stringsAsFactors = F,sep="\t") ; disease_gene$gene<-toupper(disease_gene$gene)
      symptom_gene_associations<-fread("Data/Interlayers/SympGAN_symptom_gene_associations.txt",stringsAsFactors = F,sep="\t") ; symptom_gene_associations$symptom_umls<-tolower(symptom_gene_associations$symptom_umls) ; symptom_gene_associations$Gene_Symbol<-toupper(symptom_gene_associations$Gene_Symbol)
      ## Select the symptom-gene associations ##
      symptom_gene_associations<-symptom_gene_associations[,c(1,4)]
      colnames(symptom_gene_associations)<-c("symptom","gene")
      symptom_gene_associations$symptom<-tolower(symptom_gene_associations$symptom)
      beggining<-c(length(unique(disease_symptom$icd10)),length(unique(disease_symptom$symptom)),length(disease_symptom$icd10),
                   length(unique(disease_gene$icd10)),length(unique(disease_gene$gene)),length(disease_gene$icd10),
                   length(symptom_gene_associations$symptom))
      ## Get the symptoms common to both files
      symptoms_comunes<-intersect(disease_symptom$symptom,symptom_gene_associations$symptom)
      ## Only 159 out of 303 symptoms in the disease network have an interaction with at least one gene ## 52.47%
      length(intersect(disease_symptom$symptom,symptom_gene_associations$symptom))/length(unique(disease_symptom$symptom))
      setkey(symptom_gene_associations,"symptom")
      symptom_gene_associations<-symptom_gene_associations[symptoms_comunes]
      ## Get the genes common to both files
      genes_comunes<-intersect(disease_gene$gene,symptom_gene_associations$gene)
      ## Only 1,954 out of 7,336 genes in the disease network have an interaction with at least one symptom ## 26.63%
      length(intersect(disease_gene$gene,symptom_gene_associations$gene))/length(unique(disease_gene$gene))
      setkey(symptom_gene_associations,"gene")
      symptom_gene_associations<-symptom_gene_associations[genes_comunes]
      write.table(symptom_gene_associations,"WorkFiles/ICD10_Omics/symptom_gene_associations.txt",quote=F,sep="\t",row.names=F)
      ## Get the diseases associated to genes and symptoms that have at least one association to each other (gene - symptom association)
      commongenes<-unique(symptom_gene_associations$gene) ; commonsymptoms<-unique(symptom_gene_associations$symptom)
      setkey(disease_symptom,"symptom") ; setkey(disease_gene,"gene")
      disease_gene<-disease_gene[commongenes] ; disease_symptom<-disease_symptom[commonsymptoms]
      ## At the end we have interactions between 430 symptoms and 3,818 genes ##
      ending<-c(length(unique(disease_symptom$icd10)),length(unique(disease_symptom$symptom)),length(disease_symptom$icd10),
                length(unique(disease_gene$icd10)),length(unique(disease_gene$gene)),length(disease_gene$icd10),
                length(symptom_gene_associations$symptom))
      association_summary<-cbind(beggining,ending)
      rownames(association_summary)<-c("ICD10(symptom)","symptoms","ICD10-symptom associations","ICD10(gene)","Genes","ICD10-gene associations","symptom-gene associations")
      association_summary<-cbind(association_summary,(association_summary[,2]/association_summary[,1])*100)
      colnames(association_summary)<-c("original","post_filtering","percentage_kept")
      write.table(association_summary,"Summaries/Summary_symptom_gene_associations.txt",quote=F,sep="\t")
      write.table(disease_symptom,"Data/Intermediate/symptom_gene_ICD10symptom.txt",quote=F,sep="\t",row.names=F)
      write.table(disease_gene,"Data/Intermediate/symptom_gene_ICD10gene.txt",quote=F,sep="\t",row.names=F)
      sympgene<-matrix(ncol = length(unique(symptom_gene_associations$symptom)),nrow = length(unique(symptom_gene_associations$gene)),0)
      colnames(sympgene)<-unique(symptom_gene_associations$symptom) ; rownames(sympgene)<-unique(symptom_gene_associations$gene)
      for(a in 1:length(symptom_gene_associations$symptom)){sympgene[symptom_gene_associations$gene[a],symptom_gene_associations$symptom[a]]<-1}
      write.table(sympgene,"Data/Intermediate/symptom_gene_variable_associations.txt",quote=F,sep="\t")
    }
    if(args[3]=="original"){
      ## load miRNA-gene associations and transform them into a matrix ##
      symptom_gene_associations<-fread("WorkFiles/ICD10_Omics/symptom_gene_associations.txt",stringsAsFactors = F,sep="\t")
      sympgene<-matrix(ncol = length(unique(symptom_gene_associations$symptom)),nrow = length(unique(symptom_gene_associations$gene)),0)
      colnames(sympgene)<-unique(symptom_gene_associations$symptom) ; rownames(sympgene)<-unique(symptom_gene_associations$gene)
      for(a in 1:length(symptom_gene_associations$symptom)){sympgene[symptom_gene_associations$gene[a],symptom_gene_associations$symptom[a]]<-1}
      ## load disease-symptom and disease-gene associations ##
      disease_symptom<-fread("Data/Intermediate/symptom_gene_ICD10symptom.txt",stringsAsFactors = F,sep="\t")
      disease_gene<-fread("Data/Intermediate/symptom_gene_ICD10gene.txt",stringsAsFactors = F,sep="\t")
      setkey(disease_symptom,"icd10") ; setkey(disease_gene,"icd10")
      genediseases<-unique(disease_gene$icd10) ; symptomdiseases<-unique(disease_symptom$icd10)
      ## Create a table where we are going to annotate all the results ##
      ncompas<-length(genediseases)*length(symptomdiseases)
      thetable<-matrix(ncol = 3,nrow = ncompas, NA)
      ## Calculate the score ##
      ini<-Sys.time()
      contador<-c()
      for(a in 1:length(symptomdiseases)){
        # a<-1
        symptom<-unique(disease_symptom[symptomdiseases[a]]$symptom)
        for(b in 1:length(genediseases)){
          # b<-1
          contador<-sum(contador,1)
          genes<-unique(disease_gene[genediseases[b]]$gene)
          thetable[contador,]<-c(symptomdiseases[a],genediseases[b],sum(sympgene[genes,symptom]))
        }
        print(paste(round((a/length(symptomdiseases))*100,2),"%",sep=""))
      }
      fin<-Sys.time()
      fin-ini
      colnames(thetable)<-c("icd10_symptom","icd10_gene","shared_interactions")
      write.table(thetable,"Similarities/OriginalInterlayerAssociations/ICD10_symptom_gene_associations.txt",quote=F,sep="\t",row.names=F)
      ## Create the directory for the random associations ##
      if("ICD10_symptom_gene_associations"%in%list.files("Similarities/RandomInterlayerAssociations/")==FALSE){dir.create("Similarities/RandomInterlayerAssociations/ICD10_symptom_gene_associations")}
    }
    if(args[3]=="random"){
      # input<-"symptom_gene_variable_associations.txt"
      # layer1<-"symptom_gene_ICD10symptom.txt"
      # layer2<-"symptom_gene_ICD10gene.txt"
      # saveoutput<-"ICD10_symptom_gene_associations"
      # semilla<-1
      # randomnum<-100
      ## Define the function ##
      random_distances<-function(input,layer1,layer2,saveoutput,semilla,randomnum){
        layer1name<-gsub(".txt","",gsub(".+ICD10","",layer1)) ; layer2name<-gsub(".txt","",gsub(".+ICD10","",layer2))
        set.seed(semilla)
        l1_l2_associations<-fread(paste("WorkFiles/ICD10_Omics/",layer1name,"_",layer2name,"_associations.txt",sep=""),stringsAsFactors = F,sep="\t")
        varmat2<-matrix(ncol = length(unique(l1_l2_associations[[1]])),nrow = length(unique(l1_l2_associations[[2]])),0)
        colnames(varmat2)<-unique(l1_l2_associations[[1]]) ; rownames(varmat2)<-unique(l1_l2_associations[[2]])
        for(a in 1:length(l1_l2_associations[[1]])){varmat2[l1_l2_associations[[2]][a],l1_l2_associations[[1]][a]]<-1}
        for(a in 1:randomnum){
          ini<-Sys.time()
          ## Sample variable names ##
          rvarmat<-varmat2
          colnames(rvarmat)<-sample(colnames(rvarmat)) ; rownames(rvarmat)<-sample(rownames(rvarmat))
          
          ## load disease-symptom and disease-gene associations ##
          disease_layer1<-fread(paste("Data/Intermediate/",layer1,sep=""),stringsAsFactors = F,sep="\t")
          disease_layer2<-fread(paste("Data/Intermediate/",layer2,sep=""),stringsAsFactors = F,sep="\t")
          setkey(disease_layer1,"icd10") ; setkey(disease_layer2,"icd10")
          layer2diseases<-unique(disease_layer2$icd10) ; layer1diseases<-unique(disease_layer1$icd10)
          
          ## Create a table where we are going to annotate all the results ##
          ncompas<-length(layer1diseases)*length(layer2diseases)
          rthetable<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
          ## Calculate the score ##
          ini<-Sys.time()
          contador<-c()
          for(d in 1:length(layer1diseases)){
            # d<-1
            layer1v<-unique(disease_layer1[layer1diseases[d]][[2]])
            for(b in 1:length(layer2diseases)){
              # b<-1
              contador<-sum(contador,1)
              layer2v<-unique(disease_layer2[layer2diseases[b]][[2]])
              rthetable[contador,]<-c(layer1diseases[d],layer2diseases[b],sum(rvarmat[layer2v,layer1v]))
            }
            print(paste(round((d/length(layer1diseases))*100,2),"%",sep=""))
          }
          colnames(rthetable)<-c(paste("icd10_",layer1name,sep=""),paste("icd10_",layer2name,sep=""),"shared_interactions")
          write.table(rthetable,paste("Similarities/RandomInterlayerAssociations/",saveoutput,"/Seed_",semilla,"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names=F)
          fin<-Sys.time()
          print(paste("Seed ",semilla," round ",a," time: ",fin-ini,sep=""))
        }
        final<-Sys.time()
        print(final-inicio)
      }
      random_distances(args[4],args[5],args[6],args[7],args[8],args[9])
    }
  }
  #### Connect drugs to symptoms - side effects ####
  if(args[2]=="symptomdrug_se"){
    if(args[3]=="preparetables"){
      ## Read files ##
      disease_symptom<-fread("WorkFiles/ICD10_Omics/ICD10_symptom.txt",stringsAsFactors = F,sep="\t") ; disease_symptom$symptom<-tolower(disease_symptom$symptom)
      disease_drug<-fread("WorkFiles/ICD10_Omics/ICD10_drug.txt",stringsAsFactors = F,sep="\t") ; disease_drug$drug<-tolower(disease_drug$drug)
      symptom_drug_associations<-fread("Data/Interlayers/SIDER/meddra_all_se.tsv",stringsAsFactors = F,sep="\t",header=F)
      ## transform drug codes into drug names ##
      drugname<-fread("Data/Interlayers/SIDER/drug_names.tsv",stringsAsFactors = F,sep="\t",header=F)
      drugnames<-drugname$V2 ; names(drugnames)<-drugname$V1
      symptom_drug_associations$V1<-as.character(drugnames[symptom_drug_associations$V1])
      symptom_drug_associations<-symptom_drug_associations[,c(6,1)]
      colnames(symptom_drug_associations)<-c("symptom","drug")
      symptom_drug_associations$symptom<-tolower(symptom_drug_associations$symptom)
      symptom_drug_associations$drug<-tolower(symptom_drug_associations$drug)
      beggining<-c(length(unique(disease_symptom$icd10)),length(unique(disease_symptom$symptom)),length(disease_symptom$icd10),
                   length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
                   length(symptom_drug_associations$symptom))
      ## Get the symptoms common to both files
      symptoms_comunes<-intersect(disease_symptom$symptom,symptom_drug_associations$symptom)
      ## Only 133 out of 303 symptoms in the disease network have an interaction with at least one drug ## 43.89%
      length(intersect(disease_symptom$symptom,symptom_drug_associations$symptom))/length(unique(disease_symptom$symptom))
      setkey(symptom_drug_associations,"symptom")
      symptom_drug_associations<-symptom_drug_associations[symptoms_comunes]
      ## Get the drugs common to both files
      drugs_comunes<-intersect(disease_drug$drug,symptom_drug_associations$drug)
      ## Only 829 out of 5,564 drugs in the disease network have an interaction with at least one symptom ## 14.89%
      length(intersect(disease_drug$drug,symptom_drug_associations$drug))/length(unique(disease_drug$drug))
      setkey(symptom_drug_associations,"drug")
      symptom_drug_associations<-symptom_drug_associations[drugs_comunes]
      write.table(symptom_drug_associations,"WorkFiles/ICD10_Omics/symptom_drug_associations_se.txt",quote=F,sep="\t",row.names=F)
      ## Get the diseases associated to drugs and symptoms that have at least one association to each other (drug - symptom association)
      commondrugs<-unique(symptom_drug_associations$drug) ; commonsymptoms<-unique(symptom_drug_associations$symptom)
      setkey(disease_symptom,"symptom") ; setkey(disease_drug,"drug")
      disease_drug<-disease_drug[commondrugs] ; disease_symptom<-disease_symptom[commonsymptoms]
      ## At the end we have interactions between 430 symptoms and 3,818 drugs ##
      ending<-c(length(unique(disease_symptom$icd10)),length(unique(disease_symptom$symptom)),length(disease_symptom$icd10),
                length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
                length(symptom_drug_associations$symptom))
      association_summary<-cbind(beggining,ending)
      rownames(association_summary)<-c("ICD10(symptom)","symptoms","ICD10-symptom associations","ICD10(drug)","drugs","ICD10-drug associations","symptom-drug associations")
      association_summary<-cbind(association_summary,(association_summary[,2]/association_summary[,1])*100)
      colnames(association_summary)<-c("original","post_filtering","percentage_kept")
      write.table(association_summary,"Summaries/Summary_symptom_drug_associations_se.txt",quote=F,sep="\t")
      write.table(disease_symptom,"Data/Intermediate/symptom_drug_ICD10symptom_se.txt",quote=F,sep="\t",row.names=F)
      write.table(disease_drug,"Data/Intermediate/symptom_drug_ICD10drug_se.txt",quote=F,sep="\t",row.names=F)
      sympdrug<-matrix(ncol = length(unique(symptom_drug_associations$symptom)),nrow = length(unique(symptom_drug_associations$drug)),0)
      colnames(sympdrug)<-unique(symptom_drug_associations$symptom) ; rownames(sympdrug)<-unique(symptom_drug_associations$drug)
      for(a in 1:length(symptom_drug_associations$symptom)){sympdrug[symptom_drug_associations$drug[a],symptom_drug_associations$symptom[a]]<-1}
      write.table(sympdrug,"Data/Intermediate/symptom_drug_variable_associations_se.txt",quote=F,sep="\t")
    }
    if(args[3]=="original"){
      ## load miRNA-drug associations and transform them into a matrix ##
      symptom_drug_associations<-fread("WorkFiles/ICD10_Omics/symptom_drug_associations_se.txt",stringsAsFactors = F,sep="\t")
      sympdrug<-matrix(ncol = length(unique(symptom_drug_associations$symptom)),nrow = length(unique(symptom_drug_associations$drug)),0)
      colnames(sympdrug)<-unique(symptom_drug_associations$symptom) ; rownames(sympdrug)<-unique(symptom_drug_associations$drug)
      for(a in 1:length(symptom_drug_associations$symptom)){sympdrug[symptom_drug_associations$drug[a],symptom_drug_associations$symptom[a]]<-1}
      ## load disease-symptom and disease-drug associations ##
      disease_symptom<-fread("Data/Intermediate/symptom_drug_ICD10symptom_se.txt",stringsAsFactors = F,sep="\t")
      disease_drug<-fread("Data/Intermediate/symptom_drug_ICD10drug_se.txt",stringsAsFactors = F,sep="\t")
      setkey(disease_symptom,"icd10") ; setkey(disease_drug,"icd10")
      drugdiseases<-unique(disease_drug$icd10) ; symptomdiseases<-unique(disease_symptom$icd10)
      ## Create a table where we are going to annotate all the results ##
      ncompas<-length(drugdiseases)*length(symptomdiseases)
      thetable<-matrix(ncol = 3,nrow = ncompas, NA)
      ## Calculate the score ##
      ini<-Sys.time()
      contador<-c()
      for(a in 1:length(symptomdiseases)){
        # a<-1
        symptom<-unique(disease_symptom[symptomdiseases[a]]$symptom)
        for(b in 1:length(drugdiseases)){
          # b<-1
          contador<-sum(contador,1)
          drugs<-unique(disease_drug[drugdiseases[b]]$drug)
          thetable[contador,]<-c(symptomdiseases[a],drugdiseases[b],sum(sympdrug[drugs,symptom]))
        }
        print(paste(round((a/length(symptomdiseases))*100,2),"%",sep=""))
      }
      fin<-Sys.time()
      fin-ini
      colnames(thetable)<-c("icd10_symptom","icd10_drug","shared_interactions")
      write.table(thetable,"Similarities/OriginalInterlayerAssociations/ICD10_symptom_drug_associations_se.txt",quote=F,sep="\t",row.names=F)
      ## Create the directory for the random associations ##
      if("ICD10_symptom_drug_associations_se"%in%list.files("Similarities/RandomInterlayerAssociations/")==FALSE){dir.create("Similarities/RandomInterlayerAssociations/ICD10_symptom_drug_associations_se")}
    }
    if(args[3]=="random"){
      # input<-"symptom_drug_variable_associations.txt"
      # layer1<-"symptom_drug_ICD10symptom_se.txt"
      # layer2<-"symptom_drug_ICD10drug_se.txt"
      # saveoutput<-"ICD10_symptom_drug_associations_se"
      # semilla<-1
      # randomnum<-100
      ## Define the function ##
      random_distances<-function(input,layer1,layer2,saveoutput,semilla,randomnum){
        layer1name<-gsub("_se.txt","",gsub(".+ICD10","",layer1)) ; layer2name<-gsub("_se.txt","",gsub(".+ICD10","",layer2))
        set.seed(semilla)
        l1_l2_associations<-fread(paste("WorkFiles/ICD10_Omics/",layer1name,"_",layer2name,"_associations_se.txt",sep=""),stringsAsFactors = F,sep="\t")
        varmat2<-matrix(ncol = length(unique(l1_l2_associations[[1]])),nrow = length(unique(l1_l2_associations[[2]])),0)
        colnames(varmat2)<-unique(l1_l2_associations[[1]]) ; rownames(varmat2)<-unique(l1_l2_associations[[2]])
        for(a in 1:length(l1_l2_associations[[1]])){varmat2[l1_l2_associations[[2]][a],l1_l2_associations[[1]][a]]<-1}
        for(a in 1:randomnum){
          ini<-Sys.time()
          ## Sample variable names ##
          rvarmat<-varmat2
          colnames(rvarmat)<-sample(colnames(rvarmat)) ; rownames(rvarmat)<-sample(rownames(rvarmat))
          
          ## load disease-symptom and disease-drug associations ##
          disease_layer1<-fread(paste("Data/Intermediate/",layer1,sep=""),stringsAsFactors = F,sep="\t")
          disease_layer2<-fread(paste("Data/Intermediate/",layer2,sep=""),stringsAsFactors = F,sep="\t")
          setkey(disease_layer1,"icd10") ; setkey(disease_layer2,"icd10")
          layer2diseases<-unique(disease_layer2$icd10) ; layer1diseases<-unique(disease_layer1$icd10)
          
          ## Create a table where we are going to annotate all the results ##
          ncompas<-length(layer1diseases)*length(layer2diseases)
          rthetable<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
          ## Calculate the score ##
          ini<-Sys.time()
          contador<-c()
          for(d in 1:length(layer1diseases)){
            # d<-1
            layer1v<-unique(disease_layer1[layer1diseases[d]][[2]])
            for(b in 1:length(layer2diseases)){
              # b<-1
              contador<-sum(contador,1)
              layer2v<-unique(disease_layer2[layer2diseases[b]][[2]])
              rthetable[contador,]<-c(layer1diseases[d],layer2diseases[b],sum(rvarmat[layer2v,layer1v]))
            }
            print(paste(round((d/length(layer1diseases))*100,2),"%",sep=""))
          }
          colnames(rthetable)<-c(paste("icd10_",layer1name,sep=""),paste("icd10_",layer2name,sep=""),"shared_interactions")
          write.table(rthetable,paste("Similarities/RandomInterlayerAssociations/",saveoutput,"/Seed_",semilla,"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names=F)
          fin<-Sys.time()
          print(paste("Seed ",semilla," round ",a," time: ",fin-ini,sep=""))
        }
        final<-Sys.time()
        print(final-inicio)
      }
      random_distances(args[4],args[5],args[6],args[7],args[8],args[9])
    }
  }
  #### Connect drugs to symptoms - indications ####
  if(args[2]=="symptomdrug_ind"){
    if(args[3]=="preparetables"){
      ## Read files ##
      disease_symptom<-fread("WorkFiles/ICD10_Omics/ICD10_symptom.txt",stringsAsFactors = F,sep="\t") ; disease_symptom$symptom<-tolower(disease_symptom$symptom)
      disease_drug<-fread("WorkFiles/ICD10_Omics/ICD10_drug.txt",stringsAsFactors = F,sep="\t") ; disease_drug$drug<-tolower(disease_drug$drug)
      symptom_drug_associations<-fread("Data/Interlayers/SIDER/meddra_all_indications.tsv",stringsAsFactors = F,sep="\t",header=F)
      ## transform drug codes into drug names ##
      drugname<-fread("Data/Interlayers/SIDER/drug_names.tsv",stringsAsFactors = F,sep="\t",header=F)
      drugnames<-drugname$V2 ; names(drugnames)<-drugname$V1
      symptom_drug_associations$V1<-as.character(drugnames[symptom_drug_associations$V1])
      symptom_drug_associations<-symptom_drug_associations[,c(7,1)]
      colnames(symptom_drug_associations)<-c("symptom","drug")
      symptom_drug_associations$symptom<-tolower(symptom_drug_associations$symptom)
      symptom_drug_associations$drug<-tolower(symptom_drug_associations$drug)
      beggining<-c(length(unique(disease_symptom$icd10)),length(unique(disease_symptom$symptom)),length(disease_symptom$icd10),
                   length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
                   length(symptom_drug_associations$symptom))
      ## Get the symptoms common to both files
      symptoms_comunes<-intersect(disease_symptom$symptom,symptom_drug_associations$symptom)
      ## Only 88 out of 303 symptoms in the disease network have an interaction with at least one drug ## 29.04%
      length(intersect(disease_symptom$symptom,symptom_drug_associations$symptom))/length(unique(disease_symptom$symptom))
      setkey(symptom_drug_associations,"symptom")
      symptom_drug_associations<-symptom_drug_associations[symptoms_comunes]
      ## Get the drugs common to both files
      drugs_comunes<-intersect(disease_drug$drug,symptom_drug_associations$drug)
      ## Only 829 out of 5,564 drugs in the disease network have an interaction with at least one symptom ## 14.89%
      length(intersect(disease_drug$drug,symptom_drug_associations$drug))/length(unique(disease_drug$drug))
      setkey(symptom_drug_associations,"drug")
      symptom_drug_associations<-symptom_drug_associations[drugs_comunes]
      write.table(symptom_drug_associations,"WorkFiles/ICD10_Omics/symptom_drug_associations_ind.txt",quote=F,sep="\t",row.names=F)
      ## Get the diseases associated to drugs and symptoms that have at least one association to each other (drug - symptom association)
      commondrugs<-unique(symptom_drug_associations$drug) ; commonsymptoms<-unique(symptom_drug_associations$symptom)
      setkey(disease_symptom,"symptom") ; setkey(disease_drug,"drug")
      disease_drug<-disease_drug[commondrugs] ; disease_symptom<-disease_symptom[commonsymptoms]
      ## At the end we have interactions between 430 symptoms and 3,818 drugs ##
      ending<-c(length(unique(disease_symptom$icd10)),length(unique(disease_symptom$symptom)),length(disease_symptom$icd10),
                length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
                length(symptom_drug_associations$symptom))
      association_summary<-cbind(beggining,ending)
      rownames(association_summary)<-c("ICD10(symptom)","symptoms","ICD10-symptom associations","ICD10(drug)","drugs","ICD10-drug associations","symptom-drug associations")
      association_summary<-cbind(association_summary,(association_summary[,2]/association_summary[,1])*100)
      colnames(association_summary)<-c("original","post_filtering","percentage_kept")
      write.table(association_summary,"Summaries/Summary_symptom_drug_associations_ind.txt",quote=F,sep="\t")
      write.table(disease_symptom,"Data/Intermediate/symptom_drug_ICD10symptom_ind.txt",quote=F,sep="\t",row.names=F)
      write.table(disease_drug,"Data/Intermediate/symptom_drug_ICD10drug_ind.txt",quote=F,sep="\t",row.names=F)
      sympdrug<-matrix(ncol = length(unique(symptom_drug_associations$symptom)),nrow = length(unique(symptom_drug_associations$drug)),0)
      colnames(sympdrug)<-unique(symptom_drug_associations$symptom) ; rownames(sympdrug)<-unique(symptom_drug_associations$drug)
      for(a in 1:length(symptom_drug_associations$symptom)){sympdrug[symptom_drug_associations$drug[a],symptom_drug_associations$symptom[a]]<-1}
      write.table(sympdrug,"Data/Intermediate/symptom_drug_variable_associations_ind.txt",quote=F,sep="\t")
    }
    if(args[3]=="original"){
      ## load miRNA-drug associations and transform them into a matrix ##
      symptom_drug_associations<-fread("WorkFiles/ICD10_Omics/symptom_drug_associations_ind.txt",stringsAsFactors = F,sep="\t")
      sympdrug<-matrix(ncol = length(unique(symptom_drug_associations$symptom)),nrow = length(unique(symptom_drug_associations$drug)),0)
      colnames(sympdrug)<-unique(symptom_drug_associations$symptom) ; rownames(sympdrug)<-unique(symptom_drug_associations$drug)
      for(a in 1:length(symptom_drug_associations$symptom)){sympdrug[symptom_drug_associations$drug[a],symptom_drug_associations$symptom[a]]<-1}
      ## load disease-symptom and disease-drug associations ##
      disease_symptom<-fread("Data/Intermediate/symptom_drug_ICD10symptom_ind.txt",stringsAsFactors = F,sep="\t")
      disease_drug<-fread("Data/Intermediate/symptom_drug_ICD10drug_ind.txt",stringsAsFactors = F,sep="\t")
      setkey(disease_symptom,"icd10") ; setkey(disease_drug,"icd10")
      drugdiseases<-unique(disease_drug$icd10) ; symptomdiseases<-unique(disease_symptom$icd10)
      ## Create a table where we are going to annotate all the results ##
      ncompas<-length(drugdiseases)*length(symptomdiseases)
      thetable<-matrix(ncol = 3,nrow = ncompas, NA)
      ## Calculate the score ##
      ini<-Sys.time()
      contador<-c()
      for(a in 1:length(symptomdiseases)){
        # a<-1
        symptom<-unique(disease_symptom[symptomdiseases[a]]$symptom)
        for(b in 1:length(drugdiseases)){
          # b<-1
          contador<-sum(contador,1)
          drugs<-unique(disease_drug[drugdiseases[b]]$drug)
          thetable[contador,]<-c(symptomdiseases[a],drugdiseases[b],sum(sympdrug[drugs,symptom]))
        }
        print(paste(round((a/length(symptomdiseases))*100,2),"%",sep=""))
      }
      fin<-Sys.time()
      fin-ini
      colnames(thetable)<-c("icd10_symptom","icd10_drug","shared_interactions")
      write.table(thetable,"Similarities/OriginalInterlayerAssociations/ICD10_symptom_drug_associations_ind.txt",quote=F,sep="\t",row.names=F)
      ## Create the directory for the random associations ##
      if("ICD10_symptom_drug_associations_ind"%in%list.files("Similarities/RandomInterlayerAssociations/")==FALSE){dir.create("Similarities/RandomInterlayerAssociations/ICD10_symptom_drug_associations_ind")}
    }
    if(args[3]=="random"){
      # input<-"symptom_drug_variable_associations.txt"
      # layer1<-"symptom_drug_ICD10symptom_ind.txt"
      # layer2<-"symptom_drug_ICD10drug_ind.txt"
      # saveoutput<-"ICD10_symptom_drug_associations_ind"
      # semilla<-1
      # randomnum<-100
      ## Define the function ##
      random_distances<-function(input,layer1,layer2,saveoutput,semilla,randomnum){
        layer1name<-gsub("_ind.txt","",gsub(".+ICD10","",layer1)) ; layer2name<-gsub("_ind.txt","",gsub(".+ICD10","",layer2))
        set.seed(semilla)
        l1_l2_associations<-fread(paste("WorkFiles/ICD10_Omics/",layer1name,"_",layer2name,"_associations_ind.txt",sep=""),stringsAsFactors = F,sep="\t")
        varmat2<-matrix(ncol = length(unique(l1_l2_associations[[1]])),nrow = length(unique(l1_l2_associations[[2]])),0)
        colnames(varmat2)<-unique(l1_l2_associations[[1]]) ; rownames(varmat2)<-unique(l1_l2_associations[[2]])
        for(a in 1:length(l1_l2_associations[[1]])){varmat2[l1_l2_associations[[2]][a],l1_l2_associations[[1]][a]]<-1}
        for(a in 1:randomnum){
          ini<-Sys.time()
          ## Sample variable names ##
          rvarmat<-varmat2
          colnames(rvarmat)<-sample(colnames(rvarmat)) ; rownames(rvarmat)<-sample(rownames(rvarmat))
          
          ## load disease-symptom and disease-drug associations ##
          disease_layer1<-fread(paste("Data/Intermediate/",layer1,sep=""),stringsAsFactors = F,sep="\t")
          disease_layer2<-fread(paste("Data/Intermediate/",layer2,sep=""),stringsAsFactors = F,sep="\t")
          setkey(disease_layer1,"icd10") ; setkey(disease_layer2,"icd10")
          layer2diseases<-unique(disease_layer2$icd10) ; layer1diseases<-unique(disease_layer1$icd10)
          
          ## Create a table where we are going to annotate all the results ##
          ncompas<-length(layer1diseases)*length(layer2diseases)
          rthetable<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
          ## Calculate the score ##
          ini<-Sys.time()
          contador<-c()
          for(d in 1:length(layer1diseases)){
            # d<-1
            layer1v<-unique(disease_layer1[layer1diseases[d]][[2]])
            for(b in 1:length(layer2diseases)){
              # b<-1
              contador<-sum(contador,1)
              layer2v<-unique(disease_layer2[layer2diseases[b]][[2]])
              rthetable[contador,]<-c(layer1diseases[d],layer2diseases[b],sum(rvarmat[layer2v,layer1v]))
            }
            print(paste(round((d/length(layer1diseases))*100,2),"%",sep=""))
          }
          colnames(rthetable)<-c(paste("icd10_",layer1name,sep=""),paste("icd10_",layer2name,sep=""),"shared_interactions")
          write.table(rthetable,paste("Similarities/RandomInterlayerAssociations/",saveoutput,"/Seed_",semilla,"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names=F)
          fin<-Sys.time()
          print(paste("Seed ",semilla," round ",a," time: ",fin-ini,sep=""))
        }
        final<-Sys.time()
        print(final-inicio)
      }
      random_distances(args[4],args[5],args[6],args[7],args[8],args[9])
    }
  }
  #### Connect genes to metabolites ####
  if(args[2]=="metabolitegene"){
    if(length(which(list.files("Data/Interlayers/")=="gene_metabolite_associations.txt"))==0){
      tt<-fread("Data/Metabolome_HMDB_20220712/hmdb_proteins.xml", sep="*",header=F)
      geneids<-c(grep("gene_name",tt$V1),length(tt$V1))
      genenames<-tt$V1[geneids]
      gene_metabolite<-c()
      for(a in 1:(length(geneids)-1)){
        # a<-1
        subnet<-tt$V1[geneids[a]:(geneids[a+1]-1)]
        gene<-subnet[1]
        metabolites<-gsub("<.+","",gsub("<accession>","",subnet[grep("<accession>",subnet)]))
        metabolites<-metabolites[grep("HMDB",metabolites)]
        gene_metabolite<-rbind(gene_metabolite,cbind(gene,metabolites))
        print(paste(round((a/(length(geneids)-1))*100,2),"%",sep=""))
      }
      if(length(grep("HMDBP",gene_metabolite[,2]))>0){gene_metabolite<-gene_metabolite[-grep("HMDBP",gene_metabolite[,2]),]}
      genename<-gsub("^<gene_name>","",gene_metabolite[,1])
      genename<-gsub("</gene_name>","",genename)
      gene_metabolite[,1]<-genename
      if(length(which(gene_metabolite[,1]==""))>0){gene_metabolite<-gene_metabolite[-which(gene_metabolite[,1]==""),]}
      colnames(gene_metabolite)<-c("gene","metabolome")
      if(length(which(gene_metabolite[,1]=="<gene_name/>"))>0){gene_metabolite<-gene_metabolite[-which(gene_metabolite[,1]=="<gene_name/>"),]}
      write.table(gene_metabolite,"Data/Interlayers/gene_metabolite_associations.txt",quote=F,sep="\t",row.names=F)
      ## transform the HMDB IDs into names ##
      metabtoid<-fread("Data/Metabolome_HMDB_20220712/disease_metabolite_table.txt",stringsAsFactors = F,sep="\t")
      if(length(which(duplicated(metabtoid[,2:3])))>0){metabtoid<-metabtoid[-which(duplicated(metabtoid[,2:3])),]}
      fmetabtoid<-metabtoid$Metabolite ; names(fmetabtoid)<-metabtoid$HMDB_ID
      gene_metabolite<-fread("Data/Interlayers/gene_metabolite_associations.txt",stringsAsFactors = F)
      gene_metabolite$metabolome<-as.character(fmetabtoid[gene_metabolite$metabolome])
      if(length(which(is.na(gene_metabolite$metabolome)))>0){gene_metabolite<-gene_metabolite[-which(is.na(gene_metabolite$metabolome)),]}
      write.table(gene_metabolite,"Data/Interlayers/gene_metabolite_associations.txt",quote=F,sep="\t",row.names=F)
    }
    if(args[3]=="preparetables"){
      ## Read files ##
      disease_metabolome<-fread("WorkFiles/ICD10_Omics/ICD10_metabolome.txt",stringsAsFactors = F,sep="\t") ; disease_metabolome$metabolome<-tolower(disease_metabolome$metabolome)
      disease_gene<-fread("WorkFiles/ICD10_Omics/ICD10_gene.txt",stringsAsFactors = F,sep="\t") ; disease_gene$gene<-toupper(disease_gene$gene)
      metabolome_gene_associations<-fread("Data/Interlayers/gene_metabolite_associations.txt",stringsAsFactors = F)
      metabolome_gene_associations$gene<-toupper(metabolome_gene_associations$gene)
      metabolome_gene_associations$metabolome<-tolower(metabolome_gene_associations$metabolome)
      beggining<-c(length(unique(disease_gene$icd10)),length(unique(disease_gene$gene)),length(disease_gene$icd10),
                   length(unique(disease_metabolome$icd10)),length(unique(disease_metabolome$metabolome)),length(disease_metabolome$icd10),
                   length(metabolome_gene_associations$gene))
      ## Get the metabolomes common to both files
      metabolomes_comunes<-intersect(disease_metabolome$metabolome,metabolome_gene_associations$metabolome)
      ## Only 1,123 out of 21,338 metabolomes in the disease network have an interaction with at least one gene ## 5.26%
      length(intersect(disease_metabolome$metabolome,metabolome_gene_associations$metabolome))/length(unique(disease_metabolome$metabolome))
      setkey(metabolome_gene_associations,"metabolome")
      metabolome_gene_associations<-metabolome_gene_associations[metabolomes_comunes]
      ## Get the genes common to both files
      genes_comunes<-intersect(disease_gene$gene,metabolome_gene_associations$gene)
      ## Only 2,106 out of 7,336 genes in the disease network have an interaction with at least one metabolome ## 28.71%
      length(intersect(disease_gene$gene,metabolome_gene_associations$gene))/length(unique(disease_gene$gene))
      setkey(metabolome_gene_associations,"gene")
      metabolome_gene_associations<-metabolome_gene_associations[genes_comunes]
      write.table(metabolome_gene_associations,"WorkFiles/ICD10_Omics/metabolome_gene_associations.txt",quote=F,sep="\t",row.names=F)
      ## Get the diseases associated to genes and metabolomes that have at least one association to each other (gene - metabolome association)
      commongenes<-unique(metabolome_gene_associations$gene) ; commonmetabolomes<-unique(metabolome_gene_associations$metabolome)
      setkey(disease_metabolome,"metabolome") ; setkey(disease_gene,"gene")
      disease_gene<-disease_gene[commongenes] ; disease_metabolome<-disease_metabolome[commonmetabolomes]
      ## At the end we have interactions between 430 metabolomes and 3,818 genes ##
      ending<-c(length(unique(disease_metabolome$icd10)),length(unique(disease_metabolome$metabolome)),length(disease_metabolome$icd10),
                length(unique(disease_gene$icd10)),length(unique(disease_gene$gene)),length(disease_gene$icd10),
                length(metabolome_gene_associations$metabolome))
      association_summary<-cbind(beggining,ending)
      rownames(association_summary)<-c("ICD10(metabolome)","metabolomes","ICD10-metabolome associations","ICD10(gene)","Genes","ICD10-gene associations","metabolome-gene associations")
      association_summary<-cbind(association_summary,(association_summary[,2]/association_summary[,1])*100)
      colnames(association_summary)<-c("original","post_filtering","percentage_kept")
      write.table(association_summary,"Summaries/Summary_metabolome_gene_associations.txt",quote=F,sep="\t")
      ## Save the files that will be used to calculate the interlayer connections
      write.table(disease_metabolome,"Data/Intermediate/metabolome_gene_ICD10metabolome.txt",quote=F,sep="\t",row.names=F)
      write.table(disease_gene,"Data/Intermediate/metabolome_gene_ICD10gene.txt",quote=F,sep="\t",row.names=F)
      metabolomeene<-matrix(ncol = length(unique(metabolome_gene_associations$metabolome)),nrow = length(unique(metabolome_gene_associations$gene)),0)
      colnames(metabolomeene)<-unique(metabolome_gene_associations$metabolome) ; rownames(metabolomeene)<-unique(metabolome_gene_associations$gene)
      for(a in 1:length(metabolome_gene_associations$metabolome)){metabolomeene[metabolome_gene_associations$gene[a],metabolome_gene_associations$metabolome[a]]<-1}
      write.table(metabolomeene,"Data/Intermediate/metabolome_gene_variable_associations.txt",quote=F,sep="\t")
    }
    if(args[3]=="original"){
      ## load metabolome-gene associations and transform them into a matrix ##
      metabolome_gene_associations<-fread("WorkFiles/ICD10_Omics/metabolome_gene_associations.txt",stringsAsFactors = F,sep="\t")
      mirgene<-matrix(ncol = length(unique(metabolome_gene_associations$metabolome)),nrow = length(unique(metabolome_gene_associations$gene)),0)
      colnames(mirgene)<-unique(metabolome_gene_associations$metabolome) ; rownames(mirgene)<-unique(metabolome_gene_associations$gene)
      for(a in 1:length(metabolome_gene_associations$metabolome)){mirgene[metabolome_gene_associations$gene[a],metabolome_gene_associations$metabolome[a]]<-1}
      ## load disease-metabolome and disease-gene associations ##
      disease_metabolome<-fread("Data/Intermediate/metabolome_gene_ICD10metabolome.txt",stringsAsFactors = F,sep="\t")
      disease_gene<-fread("Data/Intermediate/metabolome_gene_ICD10gene.txt",stringsAsFactors = F,sep="\t")
      setkey(disease_metabolome,"icd10") ; setkey(disease_gene,"icd10")
      genediseases<-unique(disease_gene$icd10) ; metabolomediseases<-unique(disease_metabolome$icd10)
      ## Create a table where we are going to annotate all the results ##
      ncompas<-length(genediseases)*length(metabolomediseases)
      thetable<-matrix(ncol = 3,nrow = ncompas, NA)
      ## Calculate the score ##
      ini<-Sys.time()
      contador<-c()
      for(a in 1:length(metabolomediseases)){
        # a<-1
        metabolome<-unique(disease_metabolome[metabolomediseases[a]]$metabolome)
        for(b in 1:length(genediseases)){
          # b<-1
          contador<-sum(contador,1)
          genes<-unique(disease_gene[genediseases[b]]$gene)
          thetable[contador,]<-c(metabolomediseases[a],genediseases[b],sum(mirgene[genes,metabolome]))
        }
        print(paste(round((a/length(metabolomediseases))*100,2),"%",sep=""))
      }
      fin<-Sys.time()
      fin-ini
      colnames(thetable)<-c("icd10_metabolome","icd10_gene","shared_interactions")
      write.table(thetable,"Similarities/OriginalInterlayerAssociations/ICD10_metabolome_gene_associations.txt",quote=F,sep="\t",row.names=F)
      ## Create the directory for the random associations ##
      if("ICD10_metabolome_gene_associations"%in%list.files("Similarities/RandomInterlayerAssociations/")==FALSE){dir.create("Similarities/RandomInterlayerAssociations/ICD10_metabolome_gene_associations")}
    }
    if(args[3]=="random"){
      # input<-"metabolome_gene_variable_associations.txt"
      # layer1<-"metabolome_gene_ICD10metabolome.txt"
      # layer2<-"metabolome_gene_ICD10gene.txt"
      # saveoutput<-"ICD10_metabolome_gene_associations"
      # semilla<-1
      # randomnum<-100
      ## Define the function ##
      random_distances<-function(input,layer1,layer2,saveoutput,semilla,randomnum){
        layer1name<-gsub(".txt","",gsub(".+ICD10","",layer1)) ; layer2name<-gsub(".txt","",gsub(".+ICD10","",layer2))
        set.seed(semilla)
        l1_l2_associations<-fread(paste("WorkFiles/ICD10_Omics/",layer1name,"_",layer2name,"_associations.txt",sep=""),stringsAsFactors = F,sep="\t")
        varmat2<-matrix(ncol = length(unique(l1_l2_associations[[1]])),nrow = length(unique(l1_l2_associations[[2]])),0)
        colnames(varmat2)<-unique(l1_l2_associations[[1]]) ; rownames(varmat2)<-unique(l1_l2_associations[[2]])
        for(a in 1:length(l1_l2_associations[[1]])){varmat2[l1_l2_associations[[2]][a],l1_l2_associations[[1]][a]]<-1}
        for(a in 1:randomnum){
          ini<-Sys.time()
          ## Sample variable names ##
          rvarmat<-varmat2
          colnames(rvarmat)<-sample(colnames(rvarmat)) ; rownames(rvarmat)<-sample(rownames(rvarmat))
          
          ## load disease-metabolome and disease-gene associations ##
          disease_layer1<-fread(paste("Data/Intermediate/",layer1,sep=""),stringsAsFactors = F,sep="\t")
          disease_layer2<-fread(paste("Data/Intermediate/",layer2,sep=""),stringsAsFactors = F,sep="\t")
          setkey(disease_layer1,"icd10") ; setkey(disease_layer2,"icd10")
          layer2diseases<-unique(disease_layer2$icd10) ; layer1diseases<-unique(disease_layer1$icd10)
          
          ## Create a table where we are going to annotate all the results ##
          ncompas<-length(layer1diseases)*length(layer2diseases)
          rthetable<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
          ## Calculate the score ##
          ini<-Sys.time()
          contador<-c()
          for(d in 1:length(layer1diseases)){
            # d<-1
            layer1v<-unique(disease_layer1[layer1diseases[d]][[2]])
            for(b in 1:length(layer2diseases)){
              # b<-1
              contador<-sum(contador,1)
              layer2v<-unique(disease_layer2[layer2diseases[b]][[2]])
              rthetable[contador,]<-c(layer1diseases[d],layer2diseases[b],sum(rvarmat[layer1v,layer2v]))
            }
            print(paste(round((d/length(layer1diseases))*100,2),"%",sep=""))
          }
          colnames(rthetable)<-c(paste("icd10_",layer1name,sep=""),paste("icd10_",layer2name,sep=""),"shared_interactions")
          write.table(rthetable,paste("Similarities/RandomInterlayerAssociations/",saveoutput,"/Seed_",semilla,"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names=F)
          fin<-Sys.time()
          print(paste("Seed ",semilla," round ",a," time: ",fin-ini,sep=""))
        }
        final<-Sys.time()
        print(final-inicio)
      }
      random_distances(args[4],args[5],args[6],args[7],args[8],args[9])
    }
  }
  #### Connect microbes to drugs ####
  if(args[2]=="microbedrug"){
    if(args[3]=="preparetables"){
      if(args[4]=="Microbes_metabolize_substances"){
        ## Read files ##
        disease_microbiome<-fread("WorkFiles/ICD10_Omics/ICD10_microbiome.txt",stringsAsFactors = F,sep="\t") ; disease_microbiome$microbiome<-tolower(disease_microbiome$microbiome)
        disease_drug<-fread("WorkFiles/ICD10_Omics/ICD10_drug.txt",stringsAsFactors = F,sep="\t") ; disease_drug$drug<-tolower(disease_drug$drug)
        drug_microbiome_associations<-fread("Data/Interlayers/MASI_20240725/MASI_v1.0_download_microbeSubstanceInteractionRecords_ver20200928.txt",stringsAsFactors = F)
        drug_microbiome_associations<-drug_microbiome_associations[which(drug_microbiome_associations$Interaction_Category=="Microbes metabolize substances")]
        drug_microbiome_associations<-drug_microbiome_associations[,c(5,7,8,9,12,22,24)]
        colnames(drug_microbiome_associations)[1:2]<-c("microbiome","drug")
        drug_microbiome_associations$drug<-tolower(drug_microbiome_associations$drug)
        drug_microbiome_associations$microbiome<-tolower(drug_microbiome_associations$microbiome)
        beggining<-c(length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
                     length(unique(disease_microbiome$icd10)),length(unique(disease_microbiome$microbiome)),length(disease_microbiome$icd10),
                     length(drug_microbiome_associations$drug))
        ## Get the drugs common to both files
        drugs_comunes<-intersect(disease_drug$drug,drug_microbiome_associations$drug)
        ## Only 609 out of 5,564 drugs in the disease network have an interaction with at least one gene ## 10,94%
        length(intersect(disease_drug$drug,drug_microbiome_associations$drug))/length(unique(disease_drug$drug))
        setkey(drug_microbiome_associations,"drug")
        drug_microbiome_associations<-drug_microbiome_associations[drugs_comunes]
        ## Get the genes common to both files
        microbiomes_comunes<-intersect(disease_microbiome$microbiome,drug_microbiome_associations$microbiome)
        ## Only 182 out of 833 microbiomes in the disease network have an interaction with at least one drug ## 21.84%
        length(intersect(disease_microbiome$microbiome,drug_microbiome_associations$microbiome))/length(unique(disease_microbiome$microbiome))
        setkey(drug_microbiome_associations,"microbiome")
        drug_microbiome_associations<-drug_microbiome_associations[microbiomes_comunes]
        write.table(drug_microbiome_associations,"WorkFiles/ICD10_Omics/drug_microbiome_associations_microbaffectdrug.txt",quote=F,sep="\t",row.names=F)
        ## Get the diseases associated to microbiomes and drugs that have at least one association to each other (microbiome - drug association)
        commonmicrobiomes<-unique(drug_microbiome_associations$microbiome) ; commondrugs<-unique(drug_microbiome_associations$drug)
        setkey(disease_drug,"drug") ; setkey(disease_microbiome,"microbiome")
        disease_microbiome<-disease_microbiome[commonmicrobiomes] ; disease_drug<-disease_drug[commondrugs]
        ## At the end we have interactions between 430 drugs and 3,818 microbiomes ##
        ending<-c(length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
                  length(unique(disease_microbiome$icd10)),length(unique(disease_microbiome$microbiome)),length(disease_microbiome$icd10),
                  length(drug_microbiome_associations$drug))
        association_summary<-cbind(beggining,ending)
        rownames(association_summary)<-c("ICD10(drug)","drugs","ICD10-drug associations","ICD10(microbiome)","microbiomes","ICD10-microbiome associations","drug-microbiome associations")
        association_summary<-cbind(association_summary,(association_summary[,2]/association_summary[,1])*100)
        colnames(association_summary)<-c("original","post_filtering","percentage_kept")
        write.table(association_summary,"Summaries/Summary_drug_microbiome_associations_microbaffectdrug.txt",quote=F,sep="\t")
        ## Save the files that will be used to calculate the interlayer connections
        write.table(disease_drug,"Data/Intermediate/drug_microbiome_ICD10drug_microbaffectdrug.txt",quote=F,sep="\t",row.names=F)
        write.table(disease_microbiome,"Data/Intermediate/drug_microbiome_ICD10microbiome_microbaffectdrug.txt",quote=F,sep="\t",row.names=F)
        drugmicrobiome<-matrix(ncol = length(unique(drug_microbiome_associations$drug)),nrow = length(unique(drug_microbiome_associations$microbiome)),0)
        colnames(drugmicrobiome)<-unique(drug_microbiome_associations$drug) ; rownames(drugmicrobiome)<-unique(drug_microbiome_associations$microbiome)
        for(a in 1:length(drug_microbiome_associations$drug)){drugmicrobiome[drug_microbiome_associations$microbiome[a],drug_microbiome_associations$drug[a]]<-1}
        write.table(drugmicrobiome,"Data/Intermediate/drug_microbiome_variable_associations_microbaffectdrug.txt",quote=F,sep="\t")
      }
      if(args[4]=="Substances_alter_microbe_abundance"){
        ## Read files ##
        disease_microbiome<-fread("WorkFiles/ICD10_Omics/ICD10_microbiome.txt",stringsAsFactors = F,sep="\t") ; disease_microbiome$microbiome<-tolower(disease_microbiome$microbiome)
        disease_drug<-fread("WorkFiles/ICD10_Omics/ICD10_drug.txt",stringsAsFactors = F,sep="\t") ; disease_drug$drug<-tolower(disease_drug$drug)
        drug_microbiome_associations<-fread("Data/Interlayers/MASI_20240725/MASI_v1.0_download_microbeSubstanceInteractionRecords_ver20200928.txt",stringsAsFactors = F)
        drug_microbiome_associations<-drug_microbiome_associations[which(drug_microbiome_associations$Interaction_Category=="Substances alter microbe abundance")]
        drug_microbiome_associations<-drug_microbiome_associations[,c(5,7,8,9,12,22,24)]
        colnames(drug_microbiome_associations)[1:2]<-c("microbiome","drug")
        drug_microbiome_associations$drug<-tolower(drug_microbiome_associations$drug)
        drug_microbiome_associations$microbiome<-tolower(drug_microbiome_associations$microbiome)
        beggining<-c(length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
                     length(unique(disease_microbiome$icd10)),length(unique(disease_microbiome$microbiome)),length(disease_microbiome$icd10),
                     length(drug_microbiome_associations$drug))
        ## Get the drugs common to both files
        drugs_comunes<-intersect(disease_drug$drug,drug_microbiome_associations$drug)
        ## Only 609 out of 5,564 drugs in the disease network have an interaction with at least one gene ## 10,94%
        length(intersect(disease_drug$drug,drug_microbiome_associations$drug))/length(unique(disease_drug$drug))
        setkey(drug_microbiome_associations,"drug")
        drug_microbiome_associations<-drug_microbiome_associations[drugs_comunes]
        ## Get the genes common to both files
        microbiomes_comunes<-intersect(disease_microbiome$microbiome,drug_microbiome_associations$microbiome)
        ## Only 182 out of 833 microbiomes in the disease network have an interaction with at least one drug ## 21.84%
        length(intersect(disease_microbiome$microbiome,drug_microbiome_associations$microbiome))/length(unique(disease_microbiome$microbiome))
        setkey(drug_microbiome_associations,"microbiome")
        drug_microbiome_associations<-drug_microbiome_associations[microbiomes_comunes]
        write.table(drug_microbiome_associations,"WorkFiles/ICD10_Omics/drug_microbiome_associations_drugaffectmicrob.txt",quote=F,sep="\t",row.names=F)
        ## Get the diseases associated to microbiomes and drugs that have at least one association to each other (microbiome - drug association)
        commonmicrobiomes<-unique(drug_microbiome_associations$microbiome) ; commondrugs<-unique(drug_microbiome_associations$drug)
        setkey(disease_drug,"drug") ; setkey(disease_microbiome,"microbiome")
        disease_microbiome<-disease_microbiome[commonmicrobiomes] ; disease_drug<-disease_drug[commondrugs]
        ## At the end we have interactions between 430 drugs and 3,818 microbiomes ##
        ending<-c(length(unique(disease_drug$icd10)),length(unique(disease_drug$drug)),length(disease_drug$icd10),
                  length(unique(disease_microbiome$icd10)),length(unique(disease_microbiome$microbiome)),length(disease_microbiome$icd10),
                  length(drug_microbiome_associations$drug))
        association_summary<-cbind(beggining,ending)
        rownames(association_summary)<-c("ICD10(drug)","drugs","ICD10-drug associations","ICD10(microbiome)","microbiomes","ICD10-microbiome associations","drug-microbiome associations")
        association_summary<-cbind(association_summary,(association_summary[,2]/association_summary[,1])*100)
        colnames(association_summary)<-c("original","post_filtering","percentage_kept")
        write.table(association_summary,"Summaries/Summary_drug_microbiome_associations_drugaffectmicrob.txt",quote=F,sep="\t")
        ## Save the files that will be used to calculate the interlayer connections
        write.table(disease_drug,"Data/Intermediate/drug_microbiome_ICD10drug_drugaffectmicrob.txt",quote=F,sep="\t",row.names=F)
        write.table(disease_microbiome,"Data/Intermediate/drug_microbiome_ICD10microbiome_drugaffectmicrob.txt",quote=F,sep="\t",row.names=F)
        drugmicrobiome<-matrix(ncol = length(unique(drug_microbiome_associations$drug)),nrow = length(unique(drug_microbiome_associations$microbiome)),0)
        colnames(drugmicrobiome)<-unique(drug_microbiome_associations$drug) ; rownames(drugmicrobiome)<-unique(drug_microbiome_associations$microbiome)
        for(a in 1:length(drug_microbiome_associations$drug)){drugmicrobiome[drug_microbiome_associations$microbiome[a],drug_microbiome_associations$drug[a]]<-1}
        write.table(drugmicrobiome,"Data/Intermediate/drug_microbiome_variable_associations_drugaffectmicrob.txt",quote=F,sep="\t")
      }
    }
    if(args[3]=="original"){
      if(args[4]=="Microbes_metabolize_substances"){
        ## load microbiome-gene associations and transform them into a matrix ##
        drug_microbiome_associations<-fread("WorkFiles/ICD10_Omics/drug_microbiome_associations_microbaffectdrug.txt",stringsAsFactors = F,sep="\t")
        drugmicrobiome<-matrix(ncol = length(unique(drug_microbiome_associations$drug)),nrow = length(unique(drug_microbiome_associations$microbiome)),0)
        colnames(drugmicrobiome)<-unique(drug_microbiome_associations$drug) ; rownames(drugmicrobiome)<-unique(drug_microbiome_associations$microbiome)
        for(a in 1:length(drug_microbiome_associations$drug)){drugmicrobiome[drug_microbiome_associations$microbiome[a],drug_microbiome_associations$drug[a]]<-1}
        ## load disease-drug and disease-gene associations ##
        disease_drug<-fread("Data/Intermediate/drug_microbiome_ICD10drug_microbaffectdrug.txt",stringsAsFactors = F,sep="\t")
        disease_microbiome<-fread("Data/Intermediate/drug_microbiome_ICD10microbiome_microbaffectdrug.txt",stringsAsFactors = F,sep="\t")
        setkey(disease_drug,"icd10") ; setkey(disease_microbiome,"icd10")
        microbiomediseases<-unique(disease_microbiome$icd10) ; drugdiseases<-unique(disease_drug$icd10)
        ## Create a table where we are going to annotate all the results ##
        ncompas<-length(microbiomediseases)*length(drugdiseases)
        thetable<-matrix(ncol = 3,nrow = ncompas, NA)
        ## Calculate the score ##
        ini<-Sys.time()
        contador<-c()
        for(a in 1:length(drugdiseases)){
          # a<-1
          drug<-unique(disease_drug[drugdiseases[a]]$drug)
          for(b in 1:length(microbiomediseases)){
            # b<-1
            contador<-sum(contador,1)
            microbiomes<-unique(disease_microbiome[microbiomediseases[b]]$microbiome)
            thetable[contador,]<-c(drugdiseases[a],microbiomediseases[b],sum(drugmicrobiome[microbiomes,drug]))
          }
          print(paste(round((a/length(drugdiseases))*100,2),"%",sep=""))
        }
        fin<-Sys.time()
        fin-ini
        colnames(thetable)<-c("icd10_drug","icd10_microbiome","shared_interactions")
        write.table(thetable,"Similarities/OriginalInterlayerAssociations/ICD10_drug_microbiome_associations_microbaffectdrug.txt",quote=F,sep="\t",row.names=F)
        ## Create the directory for the random associations ##
        if("ICD10_drug_microbiome_associations_microbaffectdrug"%in%list.files("Similarities/RandomInterlayerAssociations/")==FALSE){dir.create("Similarities/RandomInterlayerAssociations/ICD10_drug_microbiome_associations_microbaffectdrug")}
      }
      if(args[4]=="Substances_alter_microbe_abundance"){
        ## load microbiome-gene associations and transform them into a matrix ##
        drug_microbiome_associations<-fread("WorkFiles/ICD10_Omics/drug_microbiome_associations_drugaffectmicrob.txt",stringsAsFactors = F,sep="\t")
        drugmicrobiome<-matrix(ncol = length(unique(drug_microbiome_associations$drug)),nrow = length(unique(drug_microbiome_associations$microbiome)),0)
        colnames(drugmicrobiome)<-unique(drug_microbiome_associations$drug) ; rownames(drugmicrobiome)<-unique(drug_microbiome_associations$microbiome)
        for(a in 1:length(drug_microbiome_associations$drug)){drugmicrobiome[drug_microbiome_associations$microbiome[a],drug_microbiome_associations$drug[a]]<-1}
        ## load disease-drug and disease-gene associations ##
        disease_drug<-fread("Data/Intermediate/drug_microbiome_ICD10drug_drugaffectmicrob.txt",stringsAsFactors = F,sep="\t")
        disease_microbiome<-fread("Data/Intermediate/drug_microbiome_ICD10microbiome_drugaffectmicrob.txt",stringsAsFactors = F,sep="\t")
        setkey(disease_drug,"icd10") ; setkey(disease_microbiome,"icd10")
        microbiomediseases<-unique(disease_microbiome$icd10) ; drugdiseases<-unique(disease_drug$icd10)
        ## Create a table where we are going to annotate all the results ##
        ncompas<-length(microbiomediseases)*length(drugdiseases)
        thetable<-matrix(ncol = 3,nrow = ncompas, NA)
        ## Calculate the score ##
        ini<-Sys.time()
        contador<-c()
        for(a in 1:length(drugdiseases)){
          # a<-1
          drug<-unique(disease_drug[drugdiseases[a]]$drug)
          for(b in 1:length(microbiomediseases)){
            # b<-1
            contador<-sum(contador,1)
            microbiomes<-unique(disease_microbiome[microbiomediseases[b]]$microbiome)
            thetable[contador,]<-c(drugdiseases[a],microbiomediseases[b],sum(drugmicrobiome[microbiomes,drug]))
          }
          print(paste(round((a/length(drugdiseases))*100,2),"%",sep=""))
        }
        fin<-Sys.time()
        fin-ini
        colnames(thetable)<-c("icd10_drug","icd10_microbiome","shared_interactions")
        write.table(thetable,"Similarities/OriginalInterlayerAssociations/ICD10_drug_microbiome_associations_drugaffectmicrob.txt",quote=F,sep="\t",row.names=F)
        ## Create the directory for the random associations ##
        if("ICD10_drug_microbiome_associations_drugaffectmicrob"%in%list.files("Similarities/RandomInterlayerAssociations/")==FALSE){dir.create("Similarities/RandomInterlayerAssociations/ICD10_drug_microbiome_associations_drugaffectmicrob")}
      }
    }
    if(args[3]=="random"){
      # input<-"drug_microbiome_variable_associations_microbaffectdrug.txt"
      # layer1<-"drug_microbiome_ICD10drug_microbaffectdrug.txt"
      # layer2<-"drug_microbiome_ICD10microbiome_microbaffectdrug.txt"
      # saveoutput<-"ICD10_drug_microbiome_associations_microbaffectdrug"
      # assoc<-"microbaffectdrug"
      # semilla<-1
      # randomnum<-100
      ## Define the function ##
      random_distances<-function(input,layer1,layer2,saveoutput,semilla,randomnum,assoc){
        layer1name<-gsub(paste("_",assoc,sep=""),"",gsub(".txt","",gsub(".+ICD10","",layer1))) ; layer2name<-gsub(paste("_",assoc,sep=""),"",gsub(".txt","",gsub(".+ICD10","",layer2)))
        set.seed(semilla)
        l1_l2_associations<-fread(paste("WorkFiles/ICD10_Omics/",layer1name,"_",layer2name,"_associations_",assoc,".txt",sep=""),stringsAsFactors = F,sep="\t")
        varmat2<-matrix(ncol = length(unique(l1_l2_associations[[1]])),nrow = length(unique(l1_l2_associations[[2]])),0)
        colnames(varmat2)<-unique(l1_l2_associations[[1]]) ; rownames(varmat2)<-unique(l1_l2_associations[[2]])
        for(a in 1:length(l1_l2_associations[[1]])){varmat2[l1_l2_associations[[2]][a],l1_l2_associations[[1]][a]]<-1}
        for(a in 1:randomnum){
          ini<-Sys.time()
          ## Sample variable names ##
          rvarmat<-varmat2
          colnames(rvarmat)<-sample(colnames(rvarmat)) ; rownames(rvarmat)<-sample(rownames(rvarmat))
          
          ## load disease-drug and disease-gene associations ##
          disease_layer1<-fread(paste("Data/Intermediate/",layer1,sep=""),stringsAsFactors = F,sep="\t")
          disease_layer2<-fread(paste("Data/Intermediate/",layer2,sep=""),stringsAsFactors = F,sep="\t")
          setkey(disease_layer1,"icd10") ; setkey(disease_layer2,"icd10")
          layer2diseases<-unique(disease_layer2$icd10) ; layer1diseases<-unique(disease_layer1$icd10)
          
          ## Create a table where we are going to annotate all the results ##
          ncompas<-length(layer1diseases)*length(layer2diseases)
          rthetable<-matrix(ncol = 3,nrow = ncompas, NA) ; fuzzyjaccard<-matrix(ncol = 3,nrow = ncompas, NA)
          ## Calculate the score ##
          ini<-Sys.time()
          contador<-c()
          for(d in 1:length(layer1diseases)){
            # d<-1
            layer1v<-unique(disease_layer1[layer1diseases[d]][[2]])
            for(b in 1:length(layer2diseases)){
              # b<-1
              contador<-sum(contador,1)
              layer2v<-unique(disease_layer2[layer2diseases[b]][[2]])
              rthetable[contador,]<-c(layer1diseases[d],layer2diseases[b],sum(rvarmat[layer1v,layer2v]))
            }
            print(paste(round((d/length(layer1diseases))*100,2),"%",sep=""))
          }
          colnames(rthetable)<-c(paste("icd10_",layer1name,sep=""),paste("icd10_",layer2name,sep=""),"shared_interactions")
          write.table(rthetable,paste("Similarities/RandomInterlayerAssociations/",saveoutput,"/Seed_",semilla,"_round_",a,".txt",sep=""),quote=F,sep="\t",row.names=F)
          fin<-Sys.time()
          print(paste("Seed ",semilla," round ",a," time: ",fin-ini,sep=""))
        }
        final<-Sys.time()
        print(final-inicio)
      }
      random_distances(args[4],args[5],args[6],args[7],args[8],args[9],args[10])
    }
  }
}


#### calculate the p-value and remove all the individual files ####
if(args[1]=="calculate_pval"){
  if("Similarities_Interlayer_with_pval"%in%list.files("Similarities")==FALSE){dir.create("Similarities/Similarities_Interlayer_with_pval")}
  # thefile<-"ICD10_drug_gene_associations"
  getpval<-function(thefile){
    originalnetwork<-fread(paste("Similarities/OriginalInterlayerAssociations/",thefile,".txt",sep=""),stringsAsFactors=F,sep="\t")
    ## continue only if we haven't created the pvalue file yet ## 
    if(paste(thefile,".txt",sep="")%in%list.files("Similarities/Similarities_Interlayer_with_pval")==F){
      randomfiles<-list.files(paste("Similarities/RandomInterlayerAssociations/",thefile,sep=""))
      ## continue only if we have 10.000 random associations ##
      if(length(randomfiles)==10000){
        randommatrix<-matrix(ncol=10000,nrow=length(originalnetwork[[1]]),NA)
        for(b in 1:length(randomfiles)){
          # b<-1
          randomnetwork<-fread(paste("Similarities/RandomInterlayerAssociations/",thefile,"/",randomfiles[b],sep=""),stringsAsFactors = F,sep="\t")
          randommatrix[,b]<-randomnetwork[[3]]
          print(paste(round((b/length(randomfiles))*100,2),"%",sep=""))
        }
        resultmatrix<-matrix(ncol=4,nrow=length(originalnetwork[[1]]),NA)
        for(a in 1:length(originalnetwork[[1]])){
          # a<-1
          thepvalup<-length(which(randommatrix[a,]>=originalnetwork[[3]][a]))/length(randommatrix[a,])
          thepvaldown<-length(which(randommatrix[a,]<=originalnetwork[[3]][a]))/length(randommatrix[a,])
          meandistance<-mean(randommatrix[a,])
          sddistance<-sd(randommatrix[a,])
          resultmatrix[a,]<-c(thepvalup,thepvaldown,meandistance,sddistance)
          print(paste(round((a/length(originalnetwork[[1]]))*100,2),"%",sep=""))
        }
        resultmatrix<-cbind(originalnetwork,resultmatrix)
        colnames(resultmatrix)<-c("disease1","disease2","interactions","pval_up","pval_down","mean","sd")
        write.table(resultmatrix,paste("Similarities/Similarities_Interlayer_with_pval/",thefile,".txt",sep=""),quote=F,sep="\t",row.names=F)
      }
      if(length(randomfiles)!=10000){print(paste(thefile,"has only",length(randomfiles),"random files"))}
    }
    if(paste(thefile,".txt",sep="")%in%list.files("Similarities/Similarities_Interlayer_with_pval")==T){print("The file has been already created")}
    print("Finished!")
  }
  ## example: getpval("ICD10_drug_gene_associations")
  ## get the p-values ##
  ## @@ @@ @@@@ @@ @@ ##
  getpval(args[2])
}

#### calculate FDR (Bonferroni correction) and sort associations ####
if(args[1]=="calculate_fdr"){
  # thefiles<-gsub(".txt","",list.files("Similarities/Similarities_Interlayer_with_pval/"))
  # thefile<-thefiles[1]
  getfdr<-function(thefile){
    originalnetwork<-fread(paste("Similarities/Similarities_Interlayer_with_pval/",thefile,".txt",sep=""),stringsAsFactors = F,sep="\t")
    fdr_up<-originalnetwork$pval_up*length(originalnetwork$pval_up)
    if(length(which(fdr_up>1))>0){fdr_up[which(fdr_up>1)]<-1}
    fdr_down<-originalnetwork$pval_down*length(originalnetwork$pval_down)
    if(length(which(fdr_down>1))>0){fdr_down[which(fdr_down>1)]<-1}
    originalnetwork<-cbind(originalnetwork,fdr_up,fdr_down)
    colnames(originalnetwork)[6:7]<-c("random_mean","random_sd")
    originalnetwork<-originalnetwork[,c(1:4,8,5,9,6,7)]
    ## put the correct name to the interactions ##
    ori<-fread(paste("Similarities/OriginalInterlayerAssociations/",thefile,".txt",sep=""),stringsAsFactors = F,sep="\t")
    colnames(originalnetwork)[1:2]<-colnames(ori)[1:2]
    write.table(originalnetwork,paste("Similarities/Similarities_Interlayer_with_fdr/",thefile,".txt",sep=""),quote=F,sep="\t",row.names=F)
    print(paste(thefile," finished!"))
  }
  ## get the FDR ##
  ## @@ @@ @@ @@ ##
  files<-list.files("Similarities/Similarities_Interlayer_with_pval/")
  dones<-list.files("Similarities/Similarities_Interlayer_with_fdr/")
  ## run it only for those networks that have not been generated yet ##
  files<-gsub(".txt","",setdiff(files,dones))
  for(a in 1:length(files)){getfdr(files[a]) ; print(paste(round((a/length(files))*100,2),"%",sep=""))}
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



#### Build an heterogeneous network putting together all the information used ####
if(args[1]=="Heterogeneous_network"){
  ## Disease - variable association ##
  a1<-fread("WorkFiles/ICD10_Omics/ICD10_drugchembl.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(a1)<-c("n1","n2")
  a2<-fread("WorkFiles/ICD10_Omics/ICD10_mksymptom.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(a2)<-c("n1","n2")
  a3<-fread("WorkFiles/ICD10_Omics/ICD10_gene.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(a3)<-c("n1","n2")
  a4<-fread("WorkFiles/ICD10_Omics/ICD10_microbiome.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(a4)<-c("n1","n2")
  a5<-fread("WorkFiles/ICD10_Omics/ICD10_metabolome.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(a5)<-c("n1","n2")
  a6<-fread("WorkFiles/ICD10_Omics/ICD10_methylome.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(a6)<-c("n1","n2")
  a7<-fread("WorkFiles/ICD10_Omics/ICD10_miRNA.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(a7)<-c("n1","n2")
  ## inter layer ##
  e1<-fread("WorkFiles/ICD10_Omics/drug_microbiome_associations_microbaffectdrug.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(e1)<-c("n1","n2")
  e2<-fread("WorkFiles/ICD10_Omics/drug_microbiome_associations_drugaffectmicrob.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(e2)<-c("n1","n2")
  e3<-fread("WorkFiles/ICD10_Omics/drug_miRNA_associations.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(e3)<-c("n1","n2")
  e4<-fread("WorkFiles/ICD10_Omics/drug_gene_associations.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(e4)<-c("n1","n2")
  e5<-fread("WorkFiles/ICD10_Omics/symptom_drug_associations_se.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(e5)<-c("n1","n2")
  e6<-fread("WorkFiles/ICD10_Omics/symptom_drug_associations_ind.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(e6)<-c("n1","n2")
  e7<-fread("WorkFiles/ICD10_Omics/symptom_gene_associations.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(e7)<-c("n1","n2")
  e8<-fread("WorkFiles/ICD10_Omics/mirna_gene_associations.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(e8)<-c("n1","n2")
  e9<-fread("WorkFiles/ICD10_Omics/metabolome_gene_associations.txt",stringsAsFactors = F,sep="\t")[,1:2] ; colnames(e9)<-c("n1","n2")
  
  
  
  todas<-rbind(a1,a2,a3,a4,a5,a6,a7,e1,e2,e3,e4,e5,e6,e7,e8,e9)
  length(unique(toupper(c(todas$n1,todas$n2))))
}

























