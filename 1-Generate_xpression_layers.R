                                             ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
                                             ## Generate the disease-disease expression networks based on RNAseq and Microarray data ##
                                             ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@  @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##


## Load al the required packages
source("MetaDE/R/individual_analysis03282012.r")
source("MetaDE/R/meta_analysis03282012.r")
source("MetaDE/R/utility03282012.r")
library("affy");library(frma);library("limma");library("genefilter");require(graphics);library("hgu133plus2frmavecs");library("massiR");library("hgu133plus2.db")
library("cluster");library("NbClust");library(e1071);library("DESeq2")
## Create required directories
if ("ExpressionSet_Together"%in%list.files("Data/") == FALSE){dir.create("Data/ExpressionSet_Together")}
if ("ExpressionGenes_Together"%in%list.files("Data/") == FALSE){dir.create("Data/ExpressionGenes_Together")}
if ("Data/ExpressionSetfiles"%in%list.files("Data/") == FALSE){dir.create("Data/ExpressionSetfiles")}
if("RNAseq_genes_not_expressed_per_sample"%in%list.files("Data/")==FALSE){dir.create("Data/RNAseq_genes_not_expressed_per_sample/")}
if("RNAseq_normalized"%in%list.files("Data/")==FALSE){dir.create("Data/RNAseq_normalized/")}
if("RNAseq_not_expressed_genes_removed"%in%list.files("Data/")==FALSE){dir.create("Data/RNAseq_not_expressed_genes_removed/")}
set.seed(1)

                                                                                  ## @ @ @@ @ @ ##
                                                                                  ## Microarray ##
                                                                                  ## @ @ @@ @ @ ##

## Normalize data
load("Data/Remove_Symbols_u133plus2.Rdata") # remo
load("Data/Gene_Symbols_u133plus2.Rdata") # simbolos
tabla<-read.csv2("Data/Microarray_and_RNAseq_names_and_mesh_terms_v2.txt",stringsAsFactors = F,sep="\t")
tabla<-tabla[-which(tabla$MESH_term==""),]
rawdatadir<-"./Microarray_raw_data/"
rawdatafiles<-list.files(rawdatadir)
for (a in 1:length(rawdatafiles)){
  inicio<-Sys.time()
  ## print information regarding the dataset that is going to be analyzed
  print(paste(rawdatafiles[a],"Control:",length(list.files(paste("./",rawdatadir,"/",rawdatafiles[a],"/Control/",sep=""))),"Case:",
              length(list.files(paste("./",rawdatadir,"/",rawdatafiles[a],"/Case/",sep=""))),sep="      "))
  control<-list.files(paste("./",rawdatadir,"/",rawdatafiles[a],"/Control",sep=""))
  path_con<-paste(rawdatadir,rawdatafiles[a],"/Control/",control,sep="")
  case<-list.files(paste("./",rawdatadir,"/",rawdatafiles[a],"/Case",sep=""))
  path_ca<-paste(rawdatadir,rawdatafiles[a],"/Case/",case,sep="")
  path_read<-c(path_con, path_ca)
  controls<-cbind(control,rep("control",length(control)))
  cases<-cbind(case,rep("case",length(case)))
  targets<-rbind(controls, cases)
  affyBatch<-ReadAffy(filenames=path_read)
  expSet<-frma(affyBatch)
  save(expSet,file=paste("Data/ExpressionSetfiles/",rawdatafiles[a],".Rdata",sep=""))
  comoda<-t(as.data.frame(expSet))
  ## save the normalized expression matrix (with probes)
  save(comoda,file=paste("Data/ExpressionSet_Together/",rawdatafiles[a],".Rdata",sep=""))
  comoda1<-comoda[-length(comoda[,1]),]
  ## remove probes not associated to gene symbols
  comoda2<-comoda1[-remo,]
  coluna<-colnames(comoda2)
  genes<-unique(simbolos)
  expression<-c()
  ## convert the probes' expression matrix into a gene expression matrix, median values are calculated when seveeral probes refer to the same gene symbol
  for (c in 1:length(genes)){
    com<-which(simbolos==as.character(genes[c]))
    if (length(com)>1){iq<-apply(comoda[com,],2,median)}
    if (length(com)==1){iq<-comoda[com,]}
    expression<-rbind(expression,iq)
  }
  colnames(expression)<-coluna
  rownames(expression)<-genes
  columna<-c(control,case)
  ## change samples names
  newpatientnames<-c(paste("Control_",1:length(control),sep=""),paste("Patient_",1:length(case),sep=""))
  colnames(expression)<-newpatientnames
  controls<-cbind(paste("Control_",1:length(control),sep=""),rep("control",length(control)))
  cases<-cbind(paste("Patient_",1:length(case),sep=""),rep("case",length(case)))
  targets<-rbind(controls, cases)
  ## save gene expression matrix with the new names
  save(expression,file=paste("Data/ExpressionGenes_Together/",rawdatafiles[a],".Rdata",sep=""))
}
print("Gene expression has been normalized!!")
## Write the normalize matrix as a txt
if ("Microarray_normalized"%in%list.files("Data/") == FALSE){dir.create("Data/Microarray_normalized")}
files<-list.files("Data/Microarray_ExpressionGenes_Together/")
for(a in 1:length(files)){
  load(paste("Data/Microarray_ExpressionGenes_Together/",files[a],sep=""))
  queson<-rep(1,length(colnames(expression))) ; queson[grep("Control",colnames(expression))]<-0
  expression2<-rbind(queson,expression)
  write.table(expression2,paste("Data/Microarray_normalized/",gsub("Rdata","txt",files[a]),sep=""),quote=F,sep="\t")
}
## Conduct differential expression analyses
analyzethem<-intersect(gsub("_.+","",list.files("Microarray_raw_data/")),tabla$NameInTheStudy)
nombres<-c();for(a in 1:length(analyzethem)){nombres<-c(nombres,which(tabla$NameInTheStudy==analyzethem[a]))}
correctnames<-tabla[nombres,]
normalized<-list.files("Data/Microarray_normalized/")
if("Microarray_metanalyzed"%in%list.files("Data/")==FALSE){dir.create("Data/Microarray_metanalyzed/")}
for(a in 1:length(correctnames[,1])){
  if(length(grep(paste("^",correctnames[a,1],"_",sep=""),normalized))==1){
    cojo<-grep(paste("^",correctnames[a,1],"_",sep=""),normalized)
    mydata.raw<-MetaDE.Read(gsub(".txt","",paste("Data/Microarray_normalized/",normalized[cojo],sep="")),via="txt",skip=1,matched=TRUE,log=TRUE)
    indanal<-ind.analysis(mydata.raw,ind.method = "regt",nperm = 100,tail="abs") 
    resul<-cbind(indanal$stat,indanal$p) ; colnames(resul)<-c("tvalue","pvalue")
    write.table(resul,paste("Data/Microarray_metanalyzed/",correctnames[a,2],".txt",sep=""),quote=F,sep="\t")
  }
  if(length(grep(paste("^",correctnames[a,1],"_",sep=""),normalized))>1){
    cojo<-grep(paste("^",correctnames[a,1],"_",sep=""),normalized)
    mydata.raw<-MetaDE.Read(gsub(".txt","",paste("Data/Microarray_normalized/",normalized[cojo],sep="")),via="txt",skip=rep(1,length(cojo)),matched=TRUE,log=TRUE)
    resul<-ind.cal.ES(mydata.raw,paired=rep(FALSE,length(cojo)),nperm=1000)
    resul2<-MetaDE.ES(resul,meta.method="REM")
    resul3<-cbind(resul2$zval,resul2$FDR) ; colnames(resul3)<-c("zvalue","FDR")
    write.table(resul3,paste("Data/Microarray_metanalyzed/",correctnames[a,2],".txt",sep=""),quote=F,sep="\t")
  }
  print(paste(a," out of ",length(correctnames[,1])," finished!"))
}
## Generate a list with all the diseases and their corresponding significantly differentially up- and down-regulated genes (FDR 0.05)
ficheros<-list.files("Data/Microarray_metanalyzed/")
sdegs<-list()
for(a in 1:length(ficheros)){
  tabla<-read.csv2(paste("Data/Microarray_metanalyzed/",ficheros[a],sep=""),stringsAsFactors = F,sep="\t")
  intersect(which(as.numeric(tabla[,1])>0),which(as.numeric(tabla[,2])<=0.05))
  up<-rownames(tabla)[intersect(which(as.numeric(tabla[,1])>0),which(as.numeric(tabla[,2])<=0.05))]
  down<-rownames(tabla)[intersect(which(as.numeric(tabla[,1])<0),which(as.numeric(tabla[,2])<=0.05))]
  sdegs[[gsub(".txt","",ficheros[a])]]$up<-up
  sdegs[[gsub(".txt","",ficheros[a])]]$down<-down
}
## Calculate overlaps between the significantly up- and down-regulated genes using fisher exact tests
tt<-read.csv2(paste("Data/Microarray_metanalyzed/",ficheros[a],sep=""),stringsAsFactors = F,sep="\t")
background<-length(tt[,1])
fisin<-Sys.time()
todos<-c()
toditos<-c()
for (a in 1:(length(names(sdegs))-1)){
  positi1<-sdegs[[a]]$up ; negati1<-sdegs[[a]]$down ; pega<-c() ; pegados<-c()
  for (b in (a+1):length(names(sdegs))){
    positi2<-sdegs[[b]]$up ; negati2<-sdegs[[b]]$down
    interspd <- length(intersect(positi1,positi2))
    kkpd <- matrix(c(interspd,length(positi1)-interspd,length(positi2)-interspd,background+interspd-length(positi1)-length(positi2)),nrow=2,ncol=2)
    fispd<-fisher.test(kkpd,alternative="greater") $p.value
    intersnd <- length(intersect(negati1,negati2))
    kknd <- matrix(c(intersnd,length(negati1)-intersnd,length(negati2)-intersnd,background+intersnd-length(negati1)-length(negati2)),nrow=2,ncol=2)
    fisnd<-fisher.test(kknd,alternative="greater") $p.value
    interspi <- length(intersect(positi1,negati2))
    kkpi <- matrix(c(interspi,length(positi1)-interspi,length(negati2)-interspi,background+interspi-length(positi1)-length(negati2)),nrow=2,ncol=2)
    fispi<-fisher.test(kkpi,alternative="greater") $p.value
    intersni <- length(intersect(negati1,positi2))
    kkni <- matrix(c(intersni,length(negati1)-intersni,length(positi2)-intersni,background+intersni-length(negati1)-length(positi2)),nrow=2,ncol=2)
    fisni<-fisher.test(kkni,alternative="greater") $p.value
    junt<-c(names(sdegs)[a],names(sdegs)[b],fispi,fispd,fisni,fisnd)
    jun<-t(junt)
    write.table(jun,"Data/Microarray_fisher_overlaps.txt",quote=F,sep="\t",row.names=F,col.names=F,append = T)
  }
  print(a/(length(names(sdegs))-1))
}
print("Finished!!") ; fisfin<-Sys.time() ; print(fisfin-fisin)
## Add a positive link when there is a significant overlap between the genes significantly deregulated in the same direction in both diseases (up-up and down-down)
## Add a negative link when there is a significant overlap between the genes significantly up-regulated in one disease and down-regulated in the other one (and vice versa)
fis<-read.csv("Data/Microarray_fisher_overlaps.txt",stringsAsFactors = F,sep="\t",header=F)
threshold<-0.05
negs<-list("Pos_Neg"=which((fis[,3]*length(fis[,1]))<=threshold),"Neg_Pos"=which((fis[,5]*length(fis[,1]))<=threshold),"Pos_Pos"=which((fis[,4]*length(fis[,1]))>threshold),"Neg_Neg"=which((fis[,6]*length(fis[,1]))>threshold))
poss<-list("Pos_Neg"=which((fis[,3]*length(fis[,1]))>threshold),"Neg_Pos"=which((fis[,5]*length(fis[,1]))>threshold),"Pos_Pos"=which((fis[,4]*length(fis[,1]))<=threshold),"Neg_Neg"=which((fis[,6]*length(fis[,1]))<=threshold))
interactions<-rep(0,length(fis[,1])) ; interactions[Reduce(intersect,negs)]<--1 ; interactions[Reduce(intersect,poss)]<-1
fisher_table<-cbind(fis[,1:2],interactions) ; colnames(fisher_table)<-c("Disease_1","Disease_2","Interaction")
write.table(fisher_table,"Data/Microarray_complete_network.txt",sep="\t",quote=F,row.names=F)
network<-fisher_table[which(fisher_table[,3]!=0),]
## Save the generated network
write.table(network,"Data/Microarray_1s_and_-1s.txt",sep="\t",quote=F,row.names=F)


                                                                                     ## @@  @@ ##
                                                                                     ## RNAseq ##
                                                                                     ## @@  @@ ##

## Normalize data
rnaseqs<-list.files("RNAseq_raw_data/")
for(a in 1:length(rnaseqs)){
  countdata<-read.csv2(paste("RNAseq_raw_data/",rnaseqs[a],sep=""),stringsAsFactors = F,sep="\t")
  ## Remove genes not expressed in at least 25% of the samples ##
  ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
  ceros<-list();long<-c();todos<-c()
  for(t in 1:length(countdata[1,])){ceros[[t]]<-rownames(countdata)[(which(countdata[,t]==0))];long<-c(long,length(which(countdata[,t]==0)));todos<-c(todos,rownames(countdata)[(which(countdata[,t]==0))])}
  tabla<-as.data.frame(table(todos))[order(as.data.frame(table(todos))[,2],decreasing = T),]
  pdf(file=paste("Data/RNAseq_genes_not_expressed_per_sample/",gsub(".txt","",rnaseqs[a]),".pdf",sep=""))
  hist(tabla[,2],breaks=length(countdata[1,]),main=gsub(".txt","",rnaseqs[a]),xlab="Sample numbers",cex.axis=0.8)
  dev.off()
  ## Remove those genes not expressed in 20% or more of the samples 
  quita<-as.character(tabla[which(tabla[,2]>=round(length(countdata[1,])*0.2)),1])
  quitind<-c();for(b in 1:length(quita)){quitind<-c(quitind,which(rownames(countdata)==quita[b]))}
  countdata2<-countdata[-quitind,]
  write.table(countdata2,paste("Data/RNAseq_not_expressed_genes_removed/",rnaseqs[a],sep=""),quote=F,sep="\t")
  ## Lets see which is the difference removing and keeping the non-expressed genes from the counts ##
  condition <- factor(c(rep("exp", length(grep("Patient",colnames(countdata)))), rep("ctl", length(grep("Control",colnames(countdata))))))
  coldata <- data.frame(row.names=colnames(countdata), condition)
  queson<-rep("0",length(colnames(countdata))) ; queson[grep("Patient",colnames(countdata))]<-1
  ## Removing the genes with zero expression
  dds <- DESeqDataSetFromMatrix(countData=countdata2, colData=coldata, design=~condition)
  normalized2<-estimateSizeFactors(dds,type="ratio")
  withoutzeros<-counts(normalized2,normalized=TRUE)
  mwithoutzeros<-rbind(queson,withoutzeros)
  write.table(mwithoutzeros,paste("Data/RNAseq_normalized/",rnaseqs[a],sep=""),quote=F,sep="\t")
  print(paste(gsub(".txt","",rnaseqs[a]),"  -->  ",length(quitind)," transcripts removed",sep=""))
}
tabla<-read.csv2("Data/Microarray_and_RNAseq_names_and_mesh_terms_v2.txt",stringsAsFactors = F,sep="\t")
tabla<-tabla[-which(tabla$MESH_term==""),]
## Conduct differential expression analyses
analyzethem<-intersect(gsub("_.+","",list.files("RNAseq_raw_data/")),tabla$NameInTheStudy)
nombres<-c();for(a in 1:length(analyzethem)){nombres<-c(nombres,which(tabla$NameInTheStudy==analyzethem[a]))}
correctnames<-tabla[nombres,]
normalized<-list.files("Data/RNAseq_normalized/")
unnormalized<-list.files("Data/RNAseq_not_expressed_genes_removed/")
if("RNAseq_metanalyzed"%in%list.files("Data/")==FALSE){dir.create("Data/RNAseq_metanalyzed/")}
for(a in 1:length(correctnames[,1])){
  ## If there is only one study analyzing the corresponding disease, we conduct differential expression analysis using DESeq2
  if(length(grep(paste("^",correctnames[a,1],"_",sep=""),normalized))==1){
    cojo<-grep(paste("^",correctnames[a,1],"_",sep=""),normalized)
    mydata.raw<-MetaDE.Read(gsub(".txt","",paste("Data/RNAseq_normalized/",normalized[cojo],sep="")),via="txt",skip=1,matched=TRUE,log=TRUE)
    indanal<-ind.analysis(mydata.raw,ind.method = "regt",nperm = 100,tail="abs") 
    resul<-cbind(indanal$stat,indanal$p) ; colnames(resul)<-c("tvalue","pvalue")
    write.table(resul,paste("Data/RNAseq_metanalyzed/",correctnames[a,2],".txt",sep=""),quote=F,sep="\t")
  }
  ## If there are more than one study analyzing the corresponding disease, we conduct random effects model meta-analyses
  if(length(grep(paste("^",correctnames[a,1],"_",sep=""),normalized))>1){
    cojo<-grep(paste("^",correctnames[a,1],"_",sep=""),normalized)
    mydata.raw<-MetaDE.Read(gsub(".txt","",paste("Data/RNAseq_normalized/",normalized[cojo],sep="")),via="txt",skip=rep(1,length(cojo)),matched=TRUE,log=TRUE)
    ## In this part we add the merge step for selecting only those genes analyzed by all the studies
    mydata.raw<-MetaDE.merge(mydata.raw)
    resul<-ind.cal.ES(mydata.raw,paired=rep(FALSE,length(cojo)),nperm=1000)
    resul2<-MetaDE.ES(resul,meta.method="REM")
    resul3<-cbind(resul2$zval,resul2$FDR) ; colnames(resul3)<-c("zvalue","FDR")
    write.table(resul3,paste("Data/RNAseq_metanalyzed/",correctnames[a,2],".txt",sep=""),quote=F,sep="\t")
  }
  print(paste(a," out of ",length(correctnames[,1])," finished!"))
}
## Generate a list with all the diseases and their corresponding significantly differentially up- and down-regulated genes (FDR 0.05)
ficheros<-list.files("Data/RNAseq_metanalyzed/")
sdegs<-list()
for(a in 1:length(ficheros)){
  tabla<-read.csv2(paste("Data/RNAseq_metanalyzed/",ficheros[a],sep=""),stringsAsFactors = F,sep="\t")
  intersect(which(as.numeric(tabla[,1])>0),which(as.numeric(tabla[,2])<=0.05))
  up<-rownames(tabla)[intersect(which(as.numeric(tabla[,1])>0),which(as.numeric(tabla[,2])<=0.05))]
  down<-rownames(tabla)[intersect(which(as.numeric(tabla[,1])<0),which(as.numeric(tabla[,2])<=0.05))]
  sdegs[[gsub(".txt","",ficheros[a])]]$up<-up
  sdegs[[gsub(".txt","",ficheros[a])]]$down<-down
}
## Calculate overlaps between the significantly up- and down-regulated genes using fisher exact tests
fisin<-Sys.time()
todos<-c()
toditos<-c()
for (a in 1:(length(names(sdegs))-1)){
  tt1<-read.csv2(paste("Data/RNAseq_metanalyzed/",ficheros[a],sep=""),stringsAsFactors = F,sep="\t") ; gen1<-rownames(tt1)
  for (b in (a+1):length(names(sdegs))){
    tt2<-read.csv2(paste("Data/RNAseq_metanalyzed/",ficheros[b],sep=""),stringsAsFactors = F,sep="\t") ; gen2<-rownames(tt2)
    genes<-intersect(gen1,gen2) ; background<-length(genes)
    positi1<-intersect(sdegs[[a]]$up,genes) ; negati1<-intersect(sdegs[[a]]$down,genes)
    positi2<-intersect(sdegs[[b]]$up,genes) ; negati2<-intersect(sdegs[[b]]$down,genes)
    interspd <- length(intersect(positi1,positi2))
    kkpd <- matrix(c(interspd,length(positi1)-interspd,length(positi2)-interspd,background+interspd-length(positi1)-length(positi2)),nrow=2,ncol=2)
    fispd<-fisher.test(kkpd,alternative="greater") $p.value
    intersnd <- length(intersect(negati1,negati2))
    kknd <- matrix(c(intersnd,length(negati1)-intersnd,length(negati2)-intersnd,background+intersnd-length(negati1)-length(negati2)),nrow=2,ncol=2)
    fisnd<-fisher.test(kknd,alternative="greater") $p.value
    interspi <- length(intersect(positi1,negati2))
    kkpi <- matrix(c(interspi,length(positi1)-interspi,length(negati2)-interspi,background+interspi-length(positi1)-length(negati2)),nrow=2,ncol=2)
    fispi<-fisher.test(kkpi,alternative="greater") $p.value
    intersni <- length(intersect(negati1,positi2))
    kkni <- matrix(c(intersni,length(negati1)-intersni,length(positi2)-intersni,background+intersni-length(negati1)-length(positi2)),nrow=2,ncol=2)
    fisni<-fisher.test(kkni,alternative="greater") $p.value
    junt<-c(names(sdegs)[a],names(sdegs)[b],fispi,fispd,fisni,fisnd)
    jun<-t(junt)
    write.table(jun,"Data/RNAseq_fisher_overlaps.txt",quote=F,sep="\t",row.names=F,col.names=F,append = T)
  }
  print(a/(length(names(sdegs))-1))
}
print("Finished!!") ; fisfin<-Sys.time() ; print(fisfin-fisin)
## Add a positive link when there is a significant overlap between the genes significantly deregulated in the same direction in both diseases (up-up and down-down)
## Add a negative link when there is a significant overlap between the genes significantly up-regulated in one disease and down-regulated in the other one (and vice versa)
fis<-read.csv("Data/RNAseq_fisher_overlaps.txt",stringsAsFactors = F,sep="\t",header=F)
threshold<-0.05
negs<-list("Pos_Neg"=which((fis[,3]*length(fis[,1]))<=threshold),"Neg_Pos"=which((fis[,5]*length(fis[,1]))<=threshold),"Pos_Pos"=which((fis[,4]*length(fis[,1]))>threshold),"Neg_Neg"=which((fis[,6]*length(fis[,1]))>threshold))
poss<-list("Pos_Neg"=which((fis[,3]*length(fis[,1]))>threshold),"Neg_Pos"=which((fis[,5]*length(fis[,1]))>threshold),"Pos_Pos"=which((fis[,4]*length(fis[,1]))<=threshold),"Neg_Neg"=which((fis[,6]*length(fis[,1]))<=threshold))
interactions<-rep(0,length(fis[,1])) ; interactions[Reduce(intersect,negs)]<--1 ; interactions[Reduce(intersect,poss)]<-1
fisher_table<-cbind(fis[,1:2],interactions) ; colnames(fisher_table)<-c("Disease_1","Disease_2","Interaction")
write.table(fisher_table,"Data/RNAseq_complete_network.txt",sep="\t",quote=F,row.names=F)
network<-fisher_table[which(fisher_table[,3]!=0),]
## Save the generated network
write.table(network,"Data/RNAseq_1s_and_-1s.txt",sep="\t",quote=F,row.names=F)