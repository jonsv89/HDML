## Code for conducting the multiplex analysis to study comorbidity between COVID-19 and neurodegenerative diseases ##
## This code has been developed by Jon Sanchez-Valle (BSC) ##

## Load the needed libraries ##
suppressPackageStartupMessages({
  library("data.table")
  library("igraph")
  library("UpSetR")
  library("corrplot")
  library("CmmD")
  library("RWRtoolkit") ## information on how to run the analyses: https://github.com/dkainer/RWRtoolkit?tab=readme-ov-file
  library("RCy3")
  library("lsa")
  library("dplyr")
  library("ggplot2")
  library("dendextend")
  library("circlize")
  library(e1071)
  library("sigclust2")
  library("pvclust")
  library("parallelDist")
  library("circlize")
  # library(readr)
  # library(stringr)
  # library(reshape)
  # library(reshape2)
  # library(tidyverse)
  # library("RandomWalkRestartMH")
  # library(multinet)
  # library(parallel)
  library("forestplot")
  library(RColorBrewer)
  library("geometry")
  library(ggpattern)
  require(gridExtra)
  library(VennDiagram)
  library("dendextend")
  library(ggcorrplot)
})

## Provide and argument
args = commandArgs(trailingOnly=TRUE)

#### Conduct the clustering analysis using CmmD ####
if(args[1]=="CmmD"){
  if("CmmD_out"%in%list.files()==FALSE){dir.create("CmmD_out")}
  if(args[2]=="get_clusters_CTD"){
    if("NoNAsforMolti"%in%list.files("Networks/")==FALSE){dir.create("Networks/NoNAsforMolti")}
    redes<-paste("Networks/Intralayer/",list.files("Networks/Intralayer/"),sep="")
    redes<-redes[c(9,10,2:7)]
    for(a in 1:length(redes)){
      tt<-fread(redes[a],stringsAsFactors = F,sep="\t",header=F)
      if(length(which(is.na(tt[[2]])))>0){tt<-tt[-which(is.na(tt[[2]])),]}
      tt<-tt[,1:2]
      write.table(tt,gsub("Networks/Intralayer","Networks/NoNAsforMolti",redes[a]),quote=F,sep="\t",row.names=F,col.names = F)
    }
    # First, we generate a vector with the paths where our network files are saved: 
    netsctd<-paste("Networks/NoNAsforMolti/",list.files("Networks/NoNAsforMolti//"),sep="")
    # In order to run MolTi and analyze the outputs given by the program, we call the function CmmD():
    # If we want all the nodes of the network to be returned in the final output:
    outcmmd<-CmmD(nodelist=NULL, input_layers= netsctd, resolution_start= 0, resolution_end= 100, interval= 0.5, destfile_community_analysis= "CmmD_out/",distmethod= "hamming",threads = 4)
    ## get the matrix ##
    gencommat<-outcmmd$gene_community_matrix
    write.table(gencommat,"CmmD_out/Multiplex_CTD.txt",quote=F,sep="\t")
  }
  if(args[2]=="bootstrap_CTD"){
    set.seed(1)
    tt<-fread("CmmD_out/Multiplex_CTD.txt",stringsAsFactors = F,sep="\t")
    diseases<-tt$V1
    matriz<-as.matrix(tt[,-1]) ; rownames(matriz)<-diseases
    class(matriz)<-"numeric"
    matriz<-matriz[,order(as.numeric(gsub(".csv","",colnames(matriz))),decreasing=F)]
    
    ## Identify the limit of the resolution parameter ##
    ## plot the change in the number of clusters as we change the resolution parameter ##
    numbers<-c()
    meansize<-c()
    for(a in 1:length(matriz[1,])){
      clustersize<-table(as.character(matriz[,a]))
      numbers<-c(numbers,length(clustersize))
      meansize<-c(meansize,mean(as.numeric(clustersize)))
      print(paste(round((a/length(matriz[1,]))*100,3),"%",sep=""))
    }
    resoluciones<-as.numeric(gsub(".csv","",colnames(matriz)))
    
    ## With all the resolution parameters 
    derivative <- diff(meansize) / diff(resoluciones)
    tolerance <- 0.05
    plateau_index <- which(abs(derivative) < tolerance)[1]
    resolucion <- resoluciones[plateau_index + 1]  # +1 because diff() shortens the vector by 1
    
    ## Plot meansize ##
    ## @@ @@ @ @@ @@ ##
    plot(numbers, meansize, type = "b", col = "blue", pch = 19,
         ylab = "Average Community Size", xlab = "number of communities")
    # Add derivative on a secondary axis
    par(new = TRUE)
    plot(numbers[-1], derivative, type = "b", col = "red", pch = 17, axes = FALSE, xlab = "", ylab = "")
    axis(side = 4, col = "red", col.axis = "red")
    mtext("Derivative", side = 4, line = 3, col = "red")
    abline(v=numbers[which(resoluciones==resolucion)])
    
    matriz<-matriz[,1:which(colnames(matriz)==paste(resolucion,".csv",sep=""))]
    
    ## Define the distance matrix ##
    hammingfun <- function(x) { parDist(x,method = "hamming") }
    pvvs<-pvclust(matriz,method.dist = hammingfun,method.hclust = "single")
    pvvc<-pvclust(matriz,method.dist = hammingfun,method.hclust = "complete")
    pvva<-pvclust(matriz,method.dist = hammingfun,method.hclust = "average")
    
    lista<-list("single"==pvvs,"complete"=pvvc,"average"=pvva)
    save(lista,file="pvclust_results_nonas.RData")
    
    load("pvclust_results.RData")
    plot(lista$complete)
    pvrect(lista$complete, alpha=0.95)
    
    plot(lista$average)
    pvrect(lista$average, alpha=0.95)
    
  }
}

if(args[1]=="create_directories"){
  if("MultiplexAnalysis"%in%list.files()==FALSE){dir.create("MultiplexAnalysis")}
  if("Networks"%in%list.files("MultiplexAnalysis/")==FALSE){dir.create("MultiplexAnalysis/Networks")}
  if("Results/"%in%list.files("MultiplexAnalysis/")==FALSE){dir.create("MultiplexAnalysis/Results")}
  if("Monolayers"%in%list.files("MultiplexAnalysis")==FALSE){dir.create("MultiplexAnalysis/Monolayers")}
  if("MultiplexICD10s_ctd"%in%list.files("MultiplexAnalysis")==FALSE){dir.create("MultiplexAnalysis/MultiplexICD10s_ctd")}
  if("MultiplexICD10s_chembl"%in%list.files("MultiplexAnalysis")==FALSE){dir.create("MultiplexAnalysis/MultiplexICD10s_chembl")}
  if("MultiplexICD10s_chemblall"%in%list.files("MultiplexAnalysis")==FALSE){dir.create("MultiplexAnalysis/MultiplexICD10s_chemblall")}
  if("Objects"%in%list.files("MultiplexAnalysis/")==FALSE){dir.create("MultiplexAnalysis/Objects")}
}

if(args[1]=="prepare_layers"){
  ## remove NAs from the networks ##
  redes<-paste("Networks/Intralayer/",list.files("Networks/Intralayer/"),sep="")
  redes<-redes[-9]
  for(a in 1:length(redes)){
    tt<-fread(redes[a],stringsAsFactors = F,sep="\t",header=F)
    if(length(which(is.na(tt[[2]])))>0){tt<-tt[-which(is.na(tt[[2]])),]}
    tt<-tt[,1:2] ; colnames(tt)<-c("node1","node2")
    write.table(tt,gsub("Networks/Intralayer","MultiplexAnalysis/Networks",redes[a]),quote=F,sep="\t",row.names=F)
  }
  ## Save the file with the list of layers composing the multilayer ##
  redes<-paste("MultiplexAnalysis/Networks/",list.files("MultiplexAnalysis/Networks/"),sep="")
  ## save one file per monolayer ##
  for(a in 1:length(redes)){
    subtab<-cbind(redes[a],gsub("_.+","",gsub(".txt","",gsub(".+works/","",redes[a]))))
    write.table(subtab,paste("MultiplexAnalysis/Monolayers/",gsub(".txt",".tsv",gsub(".+works/","",redes[a])),sep=""),quote=F,sep="\t",col.names = F,row.names = F)
  }
  ## single-layers ##
  ## @@ @@ @ @@ @@ ##
  ini0<-Sys.time()
  files<-list.files("MultiplexAnalysis/Monolayers/")
  for(a in 1:length(files)){
    RWR_make_multiplex(flist = paste("MultiplexAnalysis/Monolayers/",files[a],sep=""),
                       delta=0.5,
                       output = gsub(".tsv",".Rdata",paste("MultiplexAnalysis/Objects/",files[a],sep="")))
  }
  fin0<-Sys.time()
}
#### drug layer - CTD ####
if(args[1]=="multiplex_analysis_ctd"){
  if(args[2]=="create_multiplex"){
    if(args[3]=="write_layer_lists"){
      ## Save the file with the list of layers composing the multilayer ##
      redes<-paste("MultiplexAnalysis/Networks/",list.files("MultiplexAnalysis/Networks/"),sep="")[c(9,10,3:8)]
      redestab<-cbind(redes,gsub("_.+","",gsub(".txt","",gsub(".+works/","",redes))))
      write.table(redestab,"MultiplexAnalysis/layer_list_ctd.tsv",quote=F,sep="\t",col.names = F,row.names = F)
      ## Create a list for each ICD10 to be able to run the RWR ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ ##
      ## get the entire list of nodes ##
      icds<-c() ; for(a in 1:length(redestab[,1])){aa<-fread(redestab[a,1],stringsAsFactors = F,sep="\t") ; icds<-unique(c(icds,aa[[1]],aa[[2]]))}
      if(length(which(is.na(icds)))>0){icds<-icds[-which(is.na(icds))]}
      icds<-cbind(paste("run_",icds,sep=""),icds)
      colnames(icds)<-c("V1","V2")
      ## Save one file per icd10 ##
      for(a in 1:length(icds[,1])){icdtab<-icds[a,] ; write.table(t(icdtab),paste("MultiplexAnalysis/MultiplexICD10s_ctd/",as.character(icds[a,2]),".tsv",sep=""),quote=F,sep="\t",row.names=F,col.names = F)}
    }
    if(args[3]=="create_the_object"){
      ## multiplex - CTD ## - done in 17,8 secs
      ## @@ @@ @ @ @@ @@ ##
      ini1<-Sys.time()
      RWR_make_multiplex(flist = "MultiplexAnalysis/layer_list_ctd.tsv",
                         delta=0.5,
                         output = "MultiplexAnalysis/Objects/multiplex_ctd.Rdata")
      fin1<-Sys.time()
    }
  }
  if(args[2]=="run_RWR_icd_by_icd"){
    if("DifProf_ctd"%in%list.files("MultiplexAnalysis/")==FALSE){dir.create("MultiplexAnalysis/DifProf_ctd")}
    ini<-Sys.time()
    ## load the multiplex
    string.interactions.fp <- "./MultiplexAnalysis/Objects/multiplex_ctd.Rdata"
    ## indicate the output file
    outdir.path <- 'MultiplexAnalysis/DifProf_ctd/'
    semillas<-list.files("MultiplexAnalysis/MultiplexICD10s_ctd/")
    for(a in 1:length(semillas)){
      ## load the geneset ##
      geneset.path <- paste('MultiplexAnalysis/MultiplexICD10s_ctd/',semillas[a],sep="")
      ## run the RWR ##
      RWRtoolkit::RWR_LOE(
        data= string.interactions.fp,
        seed_geneset= geneset.path,
        tau=rep(1,8),
        outdir= outdir.path)
      print(paste(round((a/length(semillas))*100,3),"%",sep=))
    }
    fin<-Sys.time()
    print(fin-ini)
  }
  if(args[2]=="individual_diffusion_profiles_to_single_matrix"){
    files<-list.files("MultiplexAnalysis/DifProf_ctd/")
    icds<-gsub("_default.+","",gsub("RWR-LOE_run_","",files))
    matriz<-matrix(nrow = length(files),ncol=length(files),NA)
    colnames(matriz)<-icds ; rownames(matriz)<-icds
    numeros<-rep(0,length(icds))
    for(a in 1:length(files)){
      # a<-1
      tabla<-fread(paste("MultiplexAnalysis/DifProf_ctd/",files[a],sep=""),stringsAsFactors = F,sep = "\t")
      setkey(tabla,"NodeNames")
      tabla<-tabla[intersect(icds,tabla$NodeNames)]
      matriz[tabla$NodeNames,a]<-as.numeric(tabla$Score)
      numeros[a]<-length(intersect(icds,tabla$NodeNames))
      print(paste(round((a/length(files))*100,3),"%",sep=""))
    }
    save(matriz,file="MultiplexAnalysis/Multiplex_CTD_diffusion_profile_matrix.Rdata")
    ## seed nodes are in columns, so, for getting the profile of a single disease we have to calculate the mean of the disease-associated icds in the columns
  }
  if(args[2]=="cluster_diseases_1"){
    load(file="MultiplexAnalysis/Multiplex_CTD_diffusion_profile_matrix.Rdata")
    clusters<-pvclust(matriz,method.dist = "euclidean",method.hclust = "ward.D2")
    save(clusters,file="pvclust_diffussion_profiles_CTD.RData")
  }
  if(args[2]=="cluster_diseases_2"){
    load("pvclust_diffussion_profiles_CTD.RData")
    pdf(file="MultiplexAnalysis/Results/Circle_plot_with_bars_CTD.pdf")
    # 1. Prepare dendrogram and significant clusters
    dend <- as.dendrogram(clusters)
    labels <- clusters$hclust$labels
    
    # Identify significant clusters (AU p-value >= 0.95)
    signif_clusters <- pvpick(clusters, alpha = 0.95, pv = "au", type = "geq")$clusters
    
    n_signif <- length(signif_clusters)
    cluster_palette <- rainbow(n_signif)
    
    # Map each label to its cluster color or grey if not significant
    label_colors <- rep("grey60", length(labels)) # default: grey
    bar_colors   <- rep("white",   length(labels)) # default: white
    cluster_ids  <- rep(NA, length(labels)) # to store cluster number for each label
    for (i in seq_along(signif_clusters)) {
      idx <- labels %in% signif_clusters[[i]]
      label_colors[idx] <- cluster_palette[i]
      bar_colors[idx]   <- cluster_palette[i]
      cluster_ids[idx]  <- i
    }
    
    # 2. Set all branches to black, labels as above
    dend <- set(dend, "branches_col", "black")
    labels_colors(dend) <- label_colors[order.dendrogram(dend)]
    
    # 3. Order bar colors and cluster ids to match dendrogram
    ordered_labels <- labels(dend)
    bar_colors <- bar_colors[match(ordered_labels, labels)]
    cluster_ids <- cluster_ids[match(ordered_labels, labels)]
    
    n <- length(ordered_labels)
    circos.clear()
    # circos.par(cell.padding = c(0, 0, 0, 0))
    circos.par(canvas.xlim = c(-1.1, 1.1), canvas.ylim = c(-1.1, 1.1), cell.padding = c(0, 0, 0, 0))
    circos.initialize(sectors = "a", xlim = c(0, n))
    
    # OUTERMOST: cluster numbers
    circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                           panel.fun = function(x, y) {
                             # Example for manual adjustment:
                             for (i in seq_len(n_signif)) {
                               inds <- which(cluster_ids == i)
                               if (length(inds) > 0) {
                                 mean_pos <- mean(inds - 0.5)
                                 ypos<-1.1
                                 # Manually adjust for specific clusters if needed
                                 if(i==15 || i==28 || i==158 || i==121 || i==132 || i==163 || i==112 || i==41 || i==146) ypos<-ypos+1.2
                                 if(i==122 || i==149 || i==128 || i==126 || i==114 || i==123 || i==129 || i==143) ypos<-ypos+1.2
                                 if(i==109 || i==57 || i==77 || i==111 || i==100 || i==86 || i==97 || i==58) ypos<-ypos+1.2
                                 if(i==25 || i==8 || i==75 || i==34 || i==74 || i==50 || i==39 || i==17) ypos<-ypos+1.2
                                 if(i==12 || i==29 || i==110 || i==137 || i==99 || i==118 || i==56 || i==93) ypos<-ypos+1.2
                                 if(i==104 || i==107 || i==85 || i==98 || i==156 || i==92 || i==153 || i==68 || i==125) ypos<-ypos+1.2
                                 circos.text(mean_pos, ypos, labels = as.character(i), 
                                             facing = "clockwise", niceFacing = TRUE,
                                             col = cluster_palette[i], cex = 0.6, adj = c(0.5, 0.5), font = 2)
                               }
                             }
                           }
                           
    )
    
    # OUTER: colored bar/ring
    circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                           panel.fun = function(x, y) {
                             for (i in seq_len(n)) {
                               circos.rect(i-1, 0, i, 1, cex=0.3, col = bar_colors[i], border = NA)
                             }
                           }
    )
    
    # MIDDLE: labels
    circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.01,
                           panel.fun = function(x, y) {
                             for (i in seq_len(n)) {
                               circos.text(i-0.5, 0, ordered_labels[i], adj = c(0, 0.5),
                                           facing = "clockwise", niceFacing = TRUE,
                                           col = label_colors[match(ordered_labels[i], labels)], cex = 0.12)
                             }
                           }
    )
    
    # INNER: dendrogram (all branches black)
    max_height <- attr(dend, "height")
    circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA, track.height = 0.4,
                           panel.fun = function(x, y) {
                             circos.dendrogram(dend, max_height = max_height)
                           }
    )
    dev.off()
  }
  if(args[2]=="cluster_diseases_larger_than_6"){
    load("pvclust_diffussion_profiles_CTD.RData")
    pdf(file="MultiplexAnalysis/Results/Circle_plot_with_bars_CTD_atleast6.pdf")
      # 1. Prepare dendrogram and significant clusters
      dend <- as.dendrogram(clusters)
      labels <- clusters$hclust$labels
      
      # Identify significant clusters (AU p-value >= 0.95)
      signif_clusters <- pvpick(clusters, alpha = 0.95, pv = "au", type = "geq", max.only = FALSE)$clusters
      # Filter to keep only clusters with ≥6 elements
      signif_clusters <- signif_clusters[sapply(signif_clusters, length) >= 6]
      
      n_signif <- length(signif_clusters)
      cluster_palette <- rainbow(n_signif)
      
      # Map each label to its cluster color or grey if not significant
      label_colors <- rep("grey60", length(labels)) # default: grey
      bar_colors   <- rep("white",   length(labels)) # default: white
      cluster_ids  <- rep(NA, length(labels)) # to store cluster number for each label
      for (i in seq_along(signif_clusters)) {
        idx <- labels %in% signif_clusters[[i]]
        label_colors[idx] <- cluster_palette[i]
        bar_colors[idx]   <- cluster_palette[i]
        cluster_ids[idx]  <- i
      }
      
      # 2. Set all branches to black, labels as above
      dend <- set(dend, "branches_col", "black")
      labels_colors(dend) <- label_colors[order.dendrogram(dend)]
      
      # 3. Order bar colors and cluster ids to match dendrogram
      ordered_labels <- labels(dend)
      bar_colors <- bar_colors[match(ordered_labels, labels)]
      cluster_ids <- cluster_ids[match(ordered_labels, labels)]
      
      n <- length(ordered_labels)
      circos.clear()
      # circos.par(cell.padding = c(0, 0, 0, 0))
      circos.par(canvas.xlim = c(-1.1, 1.1), canvas.ylim = c(-1.1, 1.1), cell.padding = c(0, 0, 0, 0))
      circos.initialize(sectors = "a", xlim = c(0, n))
      
      # OUTERMOST: cluster numbers
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                             panel.fun = function(x, y) {
                               # Example for manual adjustment:
                               for (i in seq_len(n_signif)) {
                                 inds <- which(cluster_ids == i)
                                 if (length(inds) > 0) {
                                   mean_pos <- mean(inds - 0.5)
                                   ypos<-1.1
                                   circos.text(mean_pos, ypos, labels = as.character(i), 
                                               facing = "clockwise", niceFacing = TRUE,
                                               col = cluster_palette[i], cex = 0.6, adj = c(0.5, 0.5), font = 2)
                                 }
                               }
                             }
                             
      )
      
      # OUTER: colored bar/ring
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                             panel.fun = function(x, y) {
                               for (i in seq_len(n)) {
                                 circos.rect(i-1, 0, i, 1, cex=0.3, col = bar_colors[i], border = NA)
                               }
                             }
      )
      
      # MIDDLE: labels
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.01,
                             panel.fun = function(x, y) {
                               for (i in seq_len(n)) {
                                 circos.text(i-0.5, 0, ordered_labels[i], adj = c(0, 0.5),
                                             facing = "clockwise", niceFacing = TRUE,
                                             col = label_colors[match(ordered_labels[i], labels)], cex = 0.12)
                               }
                             }
      )
      
      # INNER: dendrogram (all branches black)
      max_height <- attr(dend, "height")
      circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA, track.height = 0.4,
                             panel.fun = function(x, y) {
                               circos.dendrogram(dend, max_height = max_height)
                             }
      )
    dev.off()
  }
  if(args[2]=="describing_clusters"){
    load("pvclust_diffussion_profiles_CTD.RData")
    sclusters <- pvpick(clusters, alpha = 0.95, pv = "au", type = "geq")$clusters
    #### 2 - Plot the size of the clusters and the heterogeneity in terms of diseases composing them ####
    ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
    tipos<-c()
    sizes<-c()
    for(a in 1:168){
      tipos<-c(tipos,length(table(substr(sclusters$clusters[[a]],1,1))))
      sizes<-c(sizes,length(sclusters$clusters[[a]]))
      if(length(which(sclusters$clusters[[a]]=="G30"))>0){print(paste("G30",a))}
      if(length(which(sclusters$clusters[[a]]=="U07"))>0){print(paste("U07",a))}
    }
    dev.off()
    hist(sizes,xlim=c(0,20),breaks = 20)
    diseases<-sclusters$clusters[[135]]
    
    #### 3 - Plot some subnetworks representing the connections among the studies diseases in the different layers ####
    ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
    ## load the colors for each node ##
    dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F,sep="\t")
    colores<-dictionary$color ; names(colores)<-dictionary$icd10
    ## load the networks and induce the subgraphs ##
    redes<-paste("MultiplexAnalysis/Networks/",list.files("MultiplexAnalysis/Networks/"),sep="")[c(7,8,1:6)]
    symp<-fread(redes[1],stringsAsFactors = F,sep="\t") ; gsymp<-graph_from_data_frame(symp,directed=F)
    trans<-fread(redes[2],stringsAsFactors = F,sep="\t") ; gtrans<-graph_from_data_frame(trans,directed=F)
    drug<-fread(redes[3],stringsAsFactors = F,sep="\t") ; gdrug<-graph_from_data_frame(drug,directed=F)
    gene<-fread(redes[4],stringsAsFactors = F,sep="\t") ; ggene<-graph_from_data_frame(gene,directed=F)
    metab<-fread(redes[5],stringsAsFactors = F,sep="\t") ; gmetab<-graph_from_data_frame(metab,directed=F)
    methy<-fread(redes[6],stringsAsFactors = F,sep="\t") ; gmethy<-graph_from_data_frame(methy,directed=F)
    micro<-fread(redes[7],stringsAsFactors = F,sep="\t") ; gmicro<-graph_from_data_frame(micro,directed=F)
    mirna<-fread(redes[8],stringsAsFactors = F,sep="\t") ; gmirna<-graph_from_data_frame(mirna,directed=F)
    
    ## induce subgraph ##
    # nodes_of_interest<-sclusters$clusters[[135]]
    nodes_of_interest<-sclusters$clusters[[148]]
    # Specify nodes of interest by name
    if(length(which(is.na(nodes_of_interest)))>0){nodes_of_interest<-nodes_of_interest[-which(is.na(nodes_of_interest))]}
    
    selnodessymp<-V(gsymp)[name %in% nodes_of_interest]
    selnodestrans<-V(gtrans)[name %in% nodes_of_interest]
    selnodesdrug<-V(gdrug)[name %in% nodes_of_interest]
    selnodesgene<-V(ggene)[name %in% nodes_of_interest]
    selnodesmetab<-V(gmetab)[name %in% nodes_of_interest]
    selnodesmethy<-V(gmethy)[name %in% nodes_of_interest]
    selnodesmicro<-V(gmicro)[name %in% nodes_of_interest]
    selnodesmirna<-V(gmirna)[name %in% nodes_of_interest]
    
    sgsymp <- induced_subgraph(gsymp, selnodessymp)
    sgtrans <- induced_subgraph(gtrans, selnodestrans)
    sgdrug <- induced_subgraph(gdrug, selnodesdrug)
    sggene <- induced_subgraph(ggene, selnodesgene)
    sgmetab <- induced_subgraph(gmetab, selnodesmetab)
    sgmethy <- induced_subgraph(gmethy, selnodesmethy)
    sgmicro <- induced_subgraph(gmicro, selnodesmicro)
    sgmirna <- induced_subgraph(gmirna, selnodesmirna)
    
    ## Add colors to the nodes ##
    V(sgsymp)$color<-as.character(colores[as_ids(V(sgsymp))])
    V(sgtrans)$color<-as.character(colores[as_ids(V(sgtrans))])
    V(sgdrug)$color<-as.character(colores[as_ids(V(sgdrug))])
    V(sggene)$color<-as.character(colores[as_ids(V(sggene))])
    V(sgmetab)$color<-as.character(colores[as_ids(V(sgmetab))])
    V(sgmethy)$color<-as.character(colores[as_ids(V(sgmethy))])
    V(sgmicro)$color<-as.character(colores[as_ids(V(sgmicro))])
    V(sgmirna)$color<-as.character(colores[as_ids(V(sgmirna))])
    
    ## Plot the networks
    plot(sgsymp, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgtrans, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgdrug, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sggene, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmetab, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmethy, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmicro, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmirna, vertex.label.degree = pi/2, vertex.label.dist = 2)
    
    #### Calculate the number of edges within the cluster compared to the ones outside the clusters ####
  }
  if(args[2]=="get_shortest_paths"){
    # Can be run from anywhere so long as RWRtoolkit is installed. 
    string.interactions.fp <- "MultiplexAnalysis/Objects/multiplex_ctd.Rdata"
    source.geneset.path <- "MultiplexAnalysis/MultiplexICD10s_ctd/U07.tsv"
    target.geneset.path <- "MultiplexAnalysis/MultiplexICD10s_ctd/G30.tsv"
    outdir.path <- './MultiplexAnalysis/Results/'
    ## Get shortes paths
    shortpath<-RWRtoolkit::RWR_ShortestPaths(
      data = string.interactions.fp,
      source_geneset = source.geneset.path,
      target_geneset = target.geneset.path,
      outdir = outdir.path
    )
    #### do it layer by layer ####
    ## Microbiome ##
    microbiome<-fread("MultiplexAnalysis/Networks/microbiome.txt",stringsAsFactors = F,sep="\t")
    gmicrobiome<-graph_from_data_frame(microbiome,directed = F)
    microbiome1<-all_simple_paths(gmicrobiome, "U07", "G30", cutoff=1)
    microbiome2<-all_simple_paths(gmicrobiome, "U07", "G30", cutoff=2)
    microbiome3<-all_simple_paths(gmicrobiome, "U07", "G30", cutoff=3)
    ## Symptoms ##
    symptom<-fread("MultiplexAnalysis/Networks/symptom_2025.txt",stringsAsFactors = F,sep="\t")
    gsymptom<-graph_from_data_frame(symptom,directed = F)
    symptom1<-all_simple_paths(gsymptom, "U07", "G30", cutoff=1)
    symptom2<-all_simple_paths(gsymptom, "U07", "G30", cutoff=2)
    symptom3<-all_simple_paths(gsymptom, "U07", "G30", cutoff=3)
    ## Metabolome ##
    metabolome<-fread("MultiplexAnalysis/Networks/metabolome.txt",stringsAsFactors = F,sep="\t")
    gmetabolome<-graph_from_data_frame(metabolome,directed = F)
    metabolome1<-all_simple_paths(gmetabolome, "U07", "G30", cutoff=1)
    metabolome2<-all_simple_paths(gmetabolome, "U07", "G30", cutoff=2)
    metabolome3<-all_simple_paths(gmetabolome, "U07", "G30", cutoff=3)
    ## Drugs ##
    drug<-fread("MultiplexAnalysis/Networks/drugs_ctd.txt",stringsAsFactors = F,sep="\t")
    gdrug<-graph_from_data_frame(drug,directed = F)
    drug1<-all_simple_paths(gdrug, "U07", "G30", cutoff=1)
    drug2<-all_simple_paths(gdrug, "U07", "G30", cutoff=2)
    drug3<-all_simple_paths(gdrug, "U07", "G30", cutoff=3)
    ## Genes ##
    gene<-fread("MultiplexAnalysis/Networks/gene.txt",stringsAsFactors = F,sep="\t")
    ggene<-graph_from_data_frame(gene,directed = F)
    gene1<-all_simple_paths(ggene, "U07", "G30", cutoff=1)
    gene2<-all_simple_paths(ggene, "U07", "G30", cutoff=2)
    gene3<-all_simple_paths(ggene, "U07", "G30", cutoff=3)
    ## miRNA ##
    mirna<-fread("MultiplexAnalysis/Networks/miRNA.txt",stringsAsFactors = F,sep="\t")
    gmirna<-graph_from_data_frame(mirna,directed = F)
    mirna1<-all_simple_paths(gmirna, "U07", "G30", cutoff=1)
    mirna2<-all_simple_paths(gmirna, "U07", "G30", cutoff=2)
    mirna3<-all_simple_paths(gmirna, "U07", "G30", cutoff=3)
    
    ## create the table ##
    u07g30net<-c()
    for(a in 1:length(microbiome2)){
      u07g30net<-rbind(u07g30net,c(as_ids(microbiome2[[a]][1:2]),"microbiome"),c(as_ids(microbiome2[[a]][2:3]),"microbiome"))
    }
    for(a in 1:length(symptom2)){
      u07g30net<-rbind(u07g30net,c(as_ids(symptom2[[a]][1:2]),"symptom"),c(as_ids(symptom2[[a]][2:3]),"symptom"))
    }
    for(a in 1:length(drug2)){
      u07g30net<-rbind(u07g30net,c(as_ids(drug2[[a]][1:2]),"drug"),c(as_ids(drug2[[a]][2:3]),"drug"))
    }
    u07g30net<-rbind(u07g30net,c(as_ids(gene1[[1]]),"gene"))
    u07g30net<-rbind(u07g30net,c(as_ids(mirna1[[1]]),"miRNA"))
    u07g30net<-rbind(u07g30net,c(as_ids(metabolome1[[1]]),"metabolome"))
    colnames(u07g30net)<-c("node1","node2","layer")
    write.table(u07g30net,"MultiplexAnalysis/Results/shortest_path_net_ctd_U07_G30.txt",quote=F,sep="\t",row.names=F)
    
    dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F,sep="\t")
    colores<-dictionary$color ; names(colores)<-dictionary$icd10
    u07g30net<-fread("MultiplexAnalysis/Results/shortest_path_net_ctd_U07_G30.txt",stringsAsFactors = F,sep="\t")
    diseases<-unique(c(u07g30net$node1,u07g30net$node2))
    nodes<-cbind(diseases,as.character(colores[diseases]))
    colnames(nodes)<-c("nodes","colores")
    write.table(nodes,"MultiplexAnalysis/Results/nodes_U07_G30_ctd.txt",quote=F,sep="\t",row.names=F)
    
    # coloresoriginales<-c("#FDBF6F","dodgerblue2","#E31A1C","yellow4","palegreen2","khaki2","gold1","#FF7F00","#CAB2D6",
    #                      "#6A3D9A","deeppink1","orchid1","brown","skyblue2","maroon","steelblue4","darkturquoise","green4",
    #                      "blue1","#FB9A99","darkorange4")
    # newcolors<-c("#FDBF6F","#4184E7","#E31A1C","#8C8B2A","#A6EC9A","#EDE792","#F8D949","#FF7F00","#CAB2D6",
    #              "#6A3D9A","#EB3891","#EF8AF4","#983430","#8DBEEA","#A33A60","#426388","#5DCBCF","#3D8925",
    #              "#0F01F5","#FB9A99","#824817")
    # names(newcolors)<-coloresoriginales
    # nodes[,2]<-as.character(newcolors[nodes[,2]])
    # dictionary$color<-as.character(newcolors[dictionary$color])
    # write.table(dictionary,"DiseaseName_ICD10_dictionary.txt",quote=F,sep="\t",row.names=F)
  }
  if(args[2]=="calculate_similarities_between_diseases"){
    load("MultiplexAnalysis/Multiplex_CTD_diffusion_profile_matrix.Rdata") # matriz
    ## seed nodes are in columns, so, for getting the profile of a single disease we have to calculate the mean of the disease-associated icds in the columns
    ##substitute NA by 1 in the diagonal 
    for(a in 1:length(matriz[,1])){matriz[a,a]<-1}
    tmatriz<-t(matriz)
    cosinemat<-cosine(tmatriz)
    u07ranks<-sort(cosinemat[,"U07"],decreasing = T)[-1]
    ## position of the prediction ##
    which(names(u07ranks)=="G30")
    valores<-c()
    for(a in 1:(length(cosinemat[,1])-1)){
      valores<-c(valores,cosinemat[a,(a+1):length(cosinemat[,1])])
    }
    plot(density(valores),xlim=c(0,0.0002))
    abline(v=u07ranks[17])
    pval<-length(which(valores>=as.numeric(u07ranks[17])))/length(valores)
    pvalchanges<-c()
    for(a in 1:length(u07ranks)){
      pvalchanges<-c(pvalchanges,length(which(valores>=as.numeric(u07ranks[a])))/length(valores))
    }
    hist(valores,breaks = 100)
    
    hist(log(valores))
    abline(v=log(u07ranks[17]))
  }
}

#### drug layer - ChEMBL ####
if(args[1]=="multiplex_analysis_chembl"){
  if(args[2]=="create_multiplex"){
    if(args[3]=="write_layer_lists"){
      ## Save the file with the list of layers composing the multilayer ##
      redes<-paste("MultiplexAnalysis/Networks/",list.files("MultiplexAnalysis/Networks/"),sep="")[c(9,10,1,4:8)]
      redestab<-cbind(redes,gsub("_.+","",gsub(".txt","",gsub(".+works/","",redes))))
      write.table(redestab,"MultiplexAnalysis/layer_list_chembl.tsv",quote=F,sep="\t",col.names = F,row.names = F)
      ## Create a list for each ICD10 to be able to run the RWR ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ ##
      ## get the entire list of nodes ##
      icds<-c() ; for(a in 1:length(redestab[,1])){aa<-fread(redestab[a,1],stringsAsFactors = F,sep="\t") ; icds<-unique(c(icds,aa[[1]],aa[[2]]))}
      if(length(which(is.na(icds)))>0){icds<-icds[-which(is.na(icds))]}
      icds<-cbind(paste("run_",icds,sep=""),icds)
      colnames(icds)<-c("V1","V2")
      ## Save one file per icd10 ##
      for(a in 1:length(icds[,1])){icdtab<-icds[a,] ; write.table(t(icdtab),paste("MultiplexAnalysis/MultiplexICD10s_chembl/",as.character(icds[a,2]),".tsv",sep=""),quote=F,sep="\t",row.names=F,col.names = F)}
    }
    if(args[3]=="create_the_object"){
      ## multiplex - CTD ## - done in 17,8 secs
      ## @@ @@ @ @ @@ @@ ##
      ini1<-Sys.time()
      RWR_make_multiplex(flist = "MultiplexAnalysis/layer_list_chembl.tsv",
                         delta=0.5,
                         output = "MultiplexAnalysis/Objects/multiplex_chembl.Rdata")
      fin1<-Sys.time()
    }
  }
  if(args[2]=="run_RWR_icd_by_icd"){
    if("DifProf_chembl"%in%list.files("MultiplexAnalysis/")==FALSE){dir.create("MultiplexAnalysis/DifProf_chembl")}
    ini<-Sys.time()
    ## load the multiplex
    string.interactions.fp <- "./MultiplexAnalysis/Objects/multiplex_chembl.Rdata"
    ## indicate the output file
    outdir.path <- 'MultiplexAnalysis/DifProf_chembl/'
    semillas<-list.files("MultiplexAnalysis/MultiplexICD10s_chembl/")
    for(a in 1:length(semillas)){
      ## load the geneset ##
      geneset.path <- paste('MultiplexAnalysis/MultiplexICD10s_chembl/',semillas[a],sep="")
      ## run the RWR ##
      RWRtoolkit::RWR_LOE(
        data= string.interactions.fp,
        seed_geneset= geneset.path,
        tau=rep(1,8),
        outdir= outdir.path)
      print(paste(round((a/length(semillas))*100,3),"%",sep=))
    }
    fin<-Sys.time()
    print(fin-ini)
  }
  if(args[2]=="individual_diffusion_profiles_to_single_matrix"){
    files<-list.files("MultiplexAnalysis/DifProf_chembl/")
    icds<-gsub("_default.+","",gsub("RWR-LOE_run_","",files))
    matriz<-matrix(nrow = length(files),ncol=length(files),NA)
    colnames(matriz)<-icds ; rownames(matriz)<-icds
    numeros<-rep(0,length(icds))
    for(a in 1:length(files)){
      # a<-1
      tabla<-fread(paste("MultiplexAnalysis/DifProf_chembl/",files[a],sep=""),stringsAsFactors = F,sep = "\t")
      setkey(tabla,"NodeNames")
      tabla<-tabla[intersect(icds,tabla$NodeNames)]
      matriz[tabla$NodeNames,a]<-as.numeric(tabla$Score)
      numeros[a]<-length(intersect(icds,tabla$NodeNames))
      print(paste(round((a/length(files))*100,3),"%",sep=""))
    }
    save(matriz,file="MultiplexAnalysis/Multiplex_ChEMBL_diffusion_profile_matrix.Rdata")
    ## seed nodes are in columns, so, for getting the profile of a single disease we have to calculate the mean of the disease-associated icds in the columns
  }
  if(args[2]=="cluster_diseases_1"){
    load(file="MultiplexAnalysis/Multiplex_ChEMBL_diffusion_profile_matrix.Rdata")
    clusters<-pvclust(matriz,method.dist = "euclidean",method.hclust = "ward.D2")
    save(clusters,file="pvclust_diffussion_profiles_ChEMBL.RData")
  }
  if(args[2]=="cluster_diseases_2"){
    load("pvclust_diffussion_profiles_ChEMBL.RData")
    pdf(file="MultiplexAnalysis/Results/Circle_plot_with_bars_ChEMBL.pdf")
      # 1. Prepare dendrogram and significant clusters
      dend <- as.dendrogram(clusters)
      labels <- clusters$hclust$labels
      
      # Identify significant clusters (AU p-value >= 0.95)
      signif_clusters <- pvpick(clusters, alpha = 0.95, pv = "au", type = "geq")$clusters
      
      n_signif <- length(signif_clusters)
      cluster_palette <- rainbow(n_signif)
      
      # Map each label to its cluster color or grey if not significant
      label_colors <- rep("grey60", length(labels)) # default: grey
      bar_colors   <- rep("white",   length(labels)) # default: white
      cluster_ids  <- rep(NA, length(labels)) # to store cluster number for each label
      for (i in seq_along(signif_clusters)) {
        idx <- labels %in% signif_clusters[[i]]
        label_colors[idx] <- cluster_palette[i]
        bar_colors[idx]   <- cluster_palette[i]
        cluster_ids[idx]  <- i
      }
      
      # 2. Set all branches to black, labels as above
      dend <- set(dend, "branches_col", "black")
      labels_colors(dend) <- label_colors[order.dendrogram(dend)]
      
      # 3. Order bar colors and cluster ids to match dendrogram
      ordered_labels <- labels(dend)
      bar_colors <- bar_colors[match(ordered_labels, labels)]
      cluster_ids <- cluster_ids[match(ordered_labels, labels)]
      
      n <- length(ordered_labels)
      circos.clear()
      # circos.par(cell.padding = c(0, 0, 0, 0))
      circos.par(canvas.xlim = c(-1.1, 1.1), canvas.ylim = c(-1.1, 1.1), cell.padding = c(0, 0, 0, 0))
      circos.initialize(sectors = "a", xlim = c(0, n))
      
      # OUTERMOST: cluster numbers
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                             panel.fun = function(x, y) {
                               # Example for manual adjustment:
                               for (i in seq_len(n_signif)) {
                                 inds <- which(cluster_ids == i)
                                 if (length(inds) > 0) {
                                   mean_pos <- mean(inds - 0.5)
                                   ypos<-1.1
                                   # Manually adjust for specific clusters if needed
                                   if(i==15 || i==28 || i==158 || i==121 || i==132 || i==163 || i==112 || i==41 || i==146) ypos<-ypos+1.2
                                   if(i==122 || i==149 || i==128 || i==126 || i==114 || i==123 || i==129 || i==143) ypos<-ypos+1.2
                                   if(i==109 || i==57 || i==77 || i==111 || i==100 || i==86 || i==97 || i==58) ypos<-ypos+1.2
                                   if(i==25 || i==8 || i==75 || i==34 || i==74 || i==50 || i==39 || i==17) ypos<-ypos+1.2
                                   if(i==12 || i==29 || i==110 || i==137 || i==99 || i==118 || i==56 || i==93) ypos<-ypos+1.2
                                   if(i==104 || i==107 || i==85 || i==98 || i==156 || i==92 || i==153 || i==68 || i==125) ypos<-ypos+1.2
                                   circos.text(mean_pos, ypos, labels = as.character(i), 
                                               facing = "clockwise", niceFacing = TRUE,
                                               col = cluster_palette[i], cex = 0.6, adj = c(0.5, 0.5), font = 2)
                                 }
                               }
                             }
                             
      )
      
      # OUTER: colored bar/ring
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                             panel.fun = function(x, y) {
                               for (i in seq_len(n)) {
                                 circos.rect(i-1, 0, i, 1, cex=0.3, col = bar_colors[i], border = NA)
                               }
                             }
      )
      
      # MIDDLE: labels
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.01,
                             panel.fun = function(x, y) {
                               for (i in seq_len(n)) {
                                 circos.text(i-0.5, 0, ordered_labels[i], adj = c(0, 0.5),
                                             facing = "clockwise", niceFacing = TRUE,
                                             col = label_colors[match(ordered_labels[i], labels)], cex = 0.12)
                               }
                             }
      )
      
      # INNER: dendrogram (all branches black)
      max_height <- attr(dend, "height")
      circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA, track.height = 0.4,
                             panel.fun = function(x, y) {
                               circos.dendrogram(dend, max_height = max_height)
                             }
      )
    dev.off()
  }
  if(args[2]=="cluster_diseases_larger_than_6"){
    load("pvclust_diffussion_profiles_ChEMBL.RData")
    pdf(file="MultiplexAnalysis/Results/Circle_plot_with_bars_ChEMBL_atleast6.pdf")
      # 1. Prepare dendrogram and significant clusters
      dend <- as.dendrogram(clusters)
      labels <- clusters$hclust$labels
      
      # Identify significant clusters (AU p-value >= 0.95)
      signif_clusters <- pvpick(clusters, alpha = 0.95, pv = "au", type = "geq", max.only = FALSE)$clusters
      # Filter to keep only clusters with ≥6 elements
      signif_clusters <- signif_clusters[sapply(signif_clusters, length) >= 6]
      
      n_signif <- length(signif_clusters)
      cluster_palette <- rainbow(n_signif)
      
      # Map each label to its cluster color or grey if not significant
      label_colors <- rep("grey60", length(labels)) # default: grey
      bar_colors   <- rep("white",   length(labels)) # default: white
      cluster_ids  <- rep(NA, length(labels)) # to store cluster number for each label
      for (i in seq_along(signif_clusters)) {
        idx <- labels %in% signif_clusters[[i]]
        label_colors[idx] <- cluster_palette[i]
        bar_colors[idx]   <- cluster_palette[i]
        cluster_ids[idx]  <- i
      }
      
      # 2. Set all branches to black, labels as above
      dend <- set(dend, "branches_col", "black")
      labels_colors(dend) <- label_colors[order.dendrogram(dend)]
      
      # 3. Order bar colors and cluster ids to match dendrogram
      ordered_labels <- labels(dend)
      bar_colors <- bar_colors[match(ordered_labels, labels)]
      cluster_ids <- cluster_ids[match(ordered_labels, labels)]
      
      n <- length(ordered_labels)
      circos.clear()
      # circos.par(cell.padding = c(0, 0, 0, 0))
      circos.par(canvas.xlim = c(-1.1, 1.1), canvas.ylim = c(-1.1, 1.1), cell.padding = c(0, 0, 0, 0))
      circos.initialize(sectors = "a", xlim = c(0, n))
      
      # OUTERMOST: cluster numbers
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                             panel.fun = function(x, y) {
                               # Example for manual adjustment:
                               for (i in seq_len(n_signif)) {
                                 inds <- which(cluster_ids == i)
                                 if (length(inds) > 0) {
                                   mean_pos <- mean(inds - 0.5)
                                   ypos<-1.1
                                   circos.text(mean_pos, ypos, labels = as.character(i), 
                                               facing = "clockwise", niceFacing = TRUE,
                                               col = cluster_palette[i], cex = 0.6, adj = c(0.5, 0.5), font = 2)
                                 }
                               }
                             }
                             
      )
      
      # OUTER: colored bar/ring
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                             panel.fun = function(x, y) {
                               for (i in seq_len(n)) {
                                 circos.rect(i-1, 0, i, 1, cex=0.3, col = bar_colors[i], border = NA)
                               }
                             }
      )
      
      # MIDDLE: labels
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.01,
                             panel.fun = function(x, y) {
                               for (i in seq_len(n)) {
                                 circos.text(i-0.5, 0, ordered_labels[i], adj = c(0, 0.5),
                                             facing = "clockwise", niceFacing = TRUE,
                                             col = label_colors[match(ordered_labels[i], labels)], cex = 0.12)
                               }
                             }
      )
      
      # INNER: dendrogram (all branches black)
      max_height <- attr(dend, "height")
      circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA, track.height = 0.4,
                             panel.fun = function(x, y) {
                               circos.dendrogram(dend, max_height = max_height)
                             }
      )
    dev.off()
  }
  if(args[2]=="describing_clusters"){
    load("pvclust_diffussion_profiles_ChEMBL.RData")
    sclusters <- pvpick(clusters, alpha = 0.95, pv = "au", type = "geq")$clusters
    #### 2 - Plot the size of the clusters and the heterogeneity in terms of diseases composing them ####
    ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
    tipos<-c()
    sizes<-c()
    for(a in 1:168){
      tipos<-c(tipos,length(table(substr(sclusters$clusters[[a]],1,1))))
      sizes<-c(sizes,length(sclusters$clusters[[a]]))
      if(length(which(sclusters$clusters[[a]]=="G30"))>0){print(paste("G30",a))}
      if(length(which(sclusters$clusters[[a]]=="U07"))>0){print(paste("U07",a))}
    }
    dev.off()
    hist(sizes,xlim=c(0,20),breaks = 20)
    diseases<-sclusters$clusters[[135]]
    
    #### 3 - Plot some subnetworks representing the connections among the studies diseases in the different layers ####
    ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
    ## load the colors for each node ##
    dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F,sep="\t")
    colores<-dictionary$color ; names(colores)<-dictionary$icd10
    ## load the networks and induce the subgraphs ##
    redes<-paste("MultiplexAnalysis/Networks/",list.files("MultiplexAnalysis/Networks/"),sep="")[c(7,8,1:6)]
    symp<-fread(redes[1],stringsAsFactors = F,sep="\t") ; gsymp<-graph_from_data_frame(symp,directed=F)
    trans<-fread(redes[2],stringsAsFactors = F,sep="\t") ; gtrans<-graph_from_data_frame(trans,directed=F)
    drug<-fread(redes[3],stringsAsFactors = F,sep="\t") ; gdrug<-graph_from_data_frame(drug,directed=F)
    gene<-fread(redes[4],stringsAsFactors = F,sep="\t") ; ggene<-graph_from_data_frame(gene,directed=F)
    metab<-fread(redes[5],stringsAsFactors = F,sep="\t") ; gmetab<-graph_from_data_frame(metab,directed=F)
    methy<-fread(redes[6],stringsAsFactors = F,sep="\t") ; gmethy<-graph_from_data_frame(methy,directed=F)
    micro<-fread(redes[7],stringsAsFactors = F,sep="\t") ; gmicro<-graph_from_data_frame(micro,directed=F)
    mirna<-fread(redes[8],stringsAsFactors = F,sep="\t") ; gmirna<-graph_from_data_frame(mirna,directed=F)
    
    ## induce subgraph ##
    # nodes_of_interest<-sclusters$clusters[[135]]
    nodes_of_interest<-sclusters$clusters[[148]]
    # Specify nodes of interest by name
    if(length(which(is.na(nodes_of_interest)))>0){nodes_of_interest<-nodes_of_interest[-which(is.na(nodes_of_interest))]}
    
    selnodessymp<-V(gsymp)[name %in% nodes_of_interest]
    selnodestrans<-V(gtrans)[name %in% nodes_of_interest]
    selnodesdrug<-V(gdrug)[name %in% nodes_of_interest]
    selnodesgene<-V(ggene)[name %in% nodes_of_interest]
    selnodesmetab<-V(gmetab)[name %in% nodes_of_interest]
    selnodesmethy<-V(gmethy)[name %in% nodes_of_interest]
    selnodesmicro<-V(gmicro)[name %in% nodes_of_interest]
    selnodesmirna<-V(gmirna)[name %in% nodes_of_interest]
    
    sgsymp <- induced_subgraph(gsymp, selnodessymp)
    sgtrans <- induced_subgraph(gtrans, selnodestrans)
    sgdrug <- induced_subgraph(gdrug, selnodesdrug)
    sggene <- induced_subgraph(ggene, selnodesgene)
    sgmetab <- induced_subgraph(gmetab, selnodesmetab)
    sgmethy <- induced_subgraph(gmethy, selnodesmethy)
    sgmicro <- induced_subgraph(gmicro, selnodesmicro)
    sgmirna <- induced_subgraph(gmirna, selnodesmirna)
    
    ## Add colors to the nodes ##
    V(sgsymp)$color<-as.character(colores[as_ids(V(sgsymp))])
    V(sgtrans)$color<-as.character(colores[as_ids(V(sgtrans))])
    V(sgdrug)$color<-as.character(colores[as_ids(V(sgdrug))])
    V(sggene)$color<-as.character(colores[as_ids(V(sggene))])
    V(sgmetab)$color<-as.character(colores[as_ids(V(sgmetab))])
    V(sgmethy)$color<-as.character(colores[as_ids(V(sgmethy))])
    V(sgmicro)$color<-as.character(colores[as_ids(V(sgmicro))])
    V(sgmirna)$color<-as.character(colores[as_ids(V(sgmirna))])
    
    ## Plot the networks
    plot(sgsymp, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgtrans, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgdrug, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sggene, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmetab, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmethy, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmicro, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmirna, vertex.label.degree = pi/2, vertex.label.dist = 2)
    
    #### Calculate the number of edges within the cluster compared to the ones outside the clusters ####
  }
  if(args[2]=="get_shortest_paths"){
    # Can be run from anywhere so long as RWRtoolkit is installed. 
    string.interactions.fp <- "MultiplexAnalysis/Objects/multiplex_chembl.Rdata"
    source.geneset.path <- "MultiplexAnalysis/MultiplexICD10s_chembl/U07.tsv"
    target.geneset.path <- "MultiplexAnalysis/MultiplexICD10s_chembl/G30.tsv"
    outdir.path <- './MultiplexAnalysis/Results/'
    ## Get shortes paths
    shortpath<-RWRtoolkit::RWR_ShortestPaths(
      data = string.interactions.fp,
      source_geneset = source.geneset.path,
      target_geneset = target.geneset.path,
      outdir = outdir.path
    )
    #### do it layer by layer ####
    ## Microbiome ##
    microbiome<-fread("MultiplexAnalysis/Networks/microbiome.txt",stringsAsFactors = F,sep="\t")
    gmicrobiome<-graph_from_data_frame(microbiome,directed = F)
    microbiome1<-all_simple_paths(gmicrobiome, "U07", "G30", cutoff=1)
    microbiome2<-all_simple_paths(gmicrobiome, "U07", "G30", cutoff=2)
    microbiome3<-all_simple_paths(gmicrobiome, "U07", "G30", cutoff=3)
    ## Symptoms ##
    symptom<-fread("MultiplexAnalysis/Networks/symptom_2025.txt",stringsAsFactors = F,sep="\t")
    gsymptom<-graph_from_data_frame(symptom,directed = F)
    symptom1<-all_simple_paths(gsymptom, "U07", "G30", cutoff=1)
    symptom2<-all_simple_paths(gsymptom, "U07", "G30", cutoff=2)
    symptom3<-all_simple_paths(gsymptom, "U07", "G30", cutoff=3)
    ## Metabolome ##
    metabolome<-fread("MultiplexAnalysis/Networks/metabolome.txt",stringsAsFactors = F,sep="\t")
    gmetabolome<-graph_from_data_frame(metabolome,directed = F)
    metabolome1<-all_simple_paths(gmetabolome, "U07", "G30", cutoff=1)
    metabolome2<-all_simple_paths(gmetabolome, "U07", "G30", cutoff=2)
    metabolome3<-all_simple_paths(gmetabolome, "U07", "G30", cutoff=3)
    ## Drugs ##
    drug<-fread("MultiplexAnalysis/Networks/drugs_chembl.txt",stringsAsFactors = F,sep="\t")
    gdrug<-graph_from_data_frame(drug,directed = F)
    drug1<-all_simple_paths(gdrug, "U07", "G30", cutoff=1)
    drug2<-all_simple_paths(gdrug, "U07", "G30", cutoff=2)
    drug3<-all_simple_paths(gdrug, "U07", "G30", cutoff=3)
    ## Genes ##
    gene<-fread("MultiplexAnalysis/Networks/gene.txt",stringsAsFactors = F,sep="\t")
    ggene<-graph_from_data_frame(gene,directed = F)
    gene1<-all_simple_paths(ggene, "U07", "G30", cutoff=1)
    gene2<-all_simple_paths(ggene, "U07", "G30", cutoff=2)
    gene3<-all_simple_paths(ggene, "U07", "G30", cutoff=3)
    ## miRNA ##
    mirna<-fread("MultiplexAnalysis/Networks/miRNA.txt",stringsAsFactors = F,sep="\t")
    gmirna<-graph_from_data_frame(mirna,directed = F)
    mirna1<-all_simple_paths(gmirna, "U07", "G30", cutoff=1)
    mirna2<-all_simple_paths(gmirna, "U07", "G30", cutoff=2)
    mirna3<-all_simple_paths(gmirna, "U07", "G30", cutoff=3)
    
    ## create the table ##
    u07g30net<-c()
    for(a in 1:length(microbiome2)){
      u07g30net<-rbind(u07g30net,c(as_ids(microbiome2[[a]][1:2]),"microbiome"),c(as_ids(microbiome2[[a]][2:3]),"microbiome"))
    }
    for(a in 1:length(symptom2)){
      u07g30net<-rbind(u07g30net,c(as_ids(symptom2[[a]][1:2]),"symptom"),c(as_ids(symptom2[[a]][2:3]),"symptom"))
    }
    for(a in 1:length(drug2)){
      u07g30net<-rbind(u07g30net,c(as_ids(drug2[[a]][1:2]),"drug"),c(as_ids(drug2[[a]][2:3]),"drug"))
    }
    u07g30net<-rbind(u07g30net,c(as_ids(gene1[[1]]),"gene"))
    u07g30net<-rbind(u07g30net,c(as_ids(mirna1[[1]]),"miRNA"))
    u07g30net<-rbind(u07g30net,c(as_ids(metabolome1[[1]]),"metabolome"))
    colnames(u07g30net)<-c("node1","node2","layer")
    write.table(u07g30net,"MultiplexAnalysis/Results/shortest_path_net_chembl_U07_G30.txt",quote=F,sep="\t",row.names=F)
    
    dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F,sep="\t")
    colores<-dictionary$color ; names(colores)<-dictionary$icd10
    u07g30net<-fread("MultiplexAnalysis/Results/shortest_path_net_chembl_U07_G30.txt",stringsAsFactors = F,sep="\t")
    diseases<-unique(c(u07g30net$node1,u07g30net$node2))
    nodes<-cbind(diseases,as.character(colores[diseases]))
    colnames(nodes)<-c("nodes","colores")
    write.table(nodes,"MultiplexAnalysis/Results/nodes_U07_G30_chembl.txt",quote=F,sep="\t",row.names=F)
  }
  if(args[2]=="calculate_similarities_between_diseases"){
    load("MultiplexAnalysis/Multiplex_ChEMBL_diffusion_profile_matrix.Rdata") # matriz
    ## seed nodes are in columns, so, for getting the profile of a single disease we have to calculate the mean of the disease-associated icds in the columns
    ##substitute NA by 1 in the diagonal 
    for(a in 1:length(matriz[,1])){matriz[a,a]<-1}
    tmatriz<-t(matriz)
    cosinemat<-cosine(tmatriz)
    u07ranks<-sort(cosinemat[,"U07"],decreasing = T)[-1]
    ## position of the prediction ##
    which(names(u07ranks)=="G30")
    valores<-c()
    for(a in 1:(length(cosinemat[,1])-1)){
      valores<-c(valores,cosinemat[a,(a+1):length(cosinemat[,1])])
    }
    plot(density(valores),xlim=c(0,0.0002))
    abline(v=u07ranks[17])
    pval<-length(which(valores>=as.numeric(u07ranks[17])))/length(valores)
    pvalchanges<-c()
    for(a in 1:length(u07ranks)){
      pvalchanges<-c(pvalchanges,length(which(valores>=as.numeric(u07ranks[a])))/length(valores))
    }
    hist(valores,breaks = 100)
    
    hist(log(valores))
    abline(v=log(u07ranks[17]))
  }
}

#### drug layer - ChEMBL-all ####
if(args[1]=="multiplex_analysis_chemblall"){
  if(args[2]=="create_multiplex"){
    if(args[3]=="write_layer_lists"){
      ## Save the file with the list of layers composing the multilayer ##
      redes<-paste("MultiplexAnalysis/Networks/",list.files("MultiplexAnalysis/Networks/"),sep="")[c(9,10,2,4:8)]
      redestab<-cbind(redes,gsub("_.+","",gsub(".txt","",gsub(".+works/","",redes))))
      write.table(redestab,"MultiplexAnalysis/layer_list_chemblall.tsv",quote=F,sep="\t",col.names = F,row.names = F)
      ## Create a list for each ICD10 to be able to run the RWR ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ @@ @@ @@ @@ @@ @@ @@ @@ ##
      ## get the entire list of nodes ##
      icds<-c() ; for(a in 1:length(redestab[,1])){aa<-fread(redestab[a,1],stringsAsFactors = F,sep="\t") ; icds<-unique(c(icds,aa[[1]],aa[[2]]))}
      if(length(which(is.na(icds)))>0){icds<-icds[-which(is.na(icds))]}
      icds<-cbind(paste("run_",icds,sep=""),icds)
      colnames(icds)<-c("V1","V2")
      ## Save one file per icd10 ##
      for(a in 1:length(icds[,1])){icdtab<-icds[a,] ; write.table(t(icdtab),paste("MultiplexAnalysis/MultiplexICD10s_chemblall/",as.character(icds[a,2]),".tsv",sep=""),quote=F,sep="\t",row.names=F,col.names = F)}
    }
    if(args[3]=="create_the_object"){
      ## multiplex - CTD ## - done in 17,8 secs
      ## @@ @@ @ @ @@ @@ ##
      ini1<-Sys.time()
      RWR_make_multiplex(flist = "MultiplexAnalysis/layer_list_chemblall.tsv",
                         delta=0.5,
                         output = "MultiplexAnalysis/Objects/multiplex_chemblall.Rdata")
      fin1<-Sys.time()
    }
  }
  if(args[2]=="run_RWR_icd_by_icd"){
    if("DifProf_chemblall"%in%list.files("MultiplexAnalysis/")==FALSE){dir.create("MultiplexAnalysis/DifProf_chemblall")}
    ini<-Sys.time()
    ## load the multiplex
    string.interactions.fp <- "./MultiplexAnalysis/Objects/multiplex_chemblall.Rdata"
    ## indicate the output file
    outdir.path <- 'MultiplexAnalysis/DifProf_chemblall/'
    semillas<-list.files("MultiplexAnalysis/MultiplexICD10s_chemblall/")
    for(a in 1:length(semillas)){
      ## load the geneset ##
      geneset.path <- paste('MultiplexAnalysis/MultiplexICD10s_chemblall/',semillas[a],sep="")
      ## run the RWR ##
      RWRtoolkit::RWR_LOE(
        data= string.interactions.fp,
        seed_geneset= geneset.path,
        tau=rep(1,8),
        outdir= outdir.path)
      print(paste(round((a/length(semillas))*100,3),"%",sep=))
    }
    fin<-Sys.time()
    print(fin-ini)
  }
  if(args[2]=="individual_diffusion_profiles_to_single_matrix"){
    files<-list.files("MultiplexAnalysis/DifProf_chemblall/")
    icds<-gsub("_default.+","",gsub("RWR-LOE_run_","",files))
    matriz<-matrix(nrow = length(files),ncol=length(files),NA)
    colnames(matriz)<-icds ; rownames(matriz)<-icds
    numeros<-rep(0,length(icds))
    for(a in 1:length(files)){
      # a<-1
      tabla<-fread(paste("MultiplexAnalysis/DifProf_chemblall/",files[a],sep=""),stringsAsFactors = F,sep = "\t")
      setkey(tabla,"NodeNames")
      tabla<-tabla[intersect(icds,tabla$NodeNames)]
      matriz[tabla$NodeNames,a]<-as.numeric(tabla$Score)
      numeros[a]<-length(intersect(icds,tabla$NodeNames))
      print(paste(round((a/length(files))*100,3),"%",sep=""))
    }
    save(matriz,file="MultiplexAnalysis/Multiplex_ChEMBLall_diffusion_profile_matrix.Rdata")
    ## seed nodes are in columns, so, for getting the profile of a single disease we have to calculate the mean of the disease-associated icds in the columns
  }
  if(args[2]=="cluster_diseases_1"){
    load(file="MultiplexAnalysis/Multiplex_ChEMBLall_diffusion_profile_matrix.Rdata")
    clusters<-pvclust(matriz,method.dist = "euclidean",method.hclust = "ward.D2")
    save(clusters,file="pvclust_diffussion_profiles_ChEMBLall.RData")
  }
  if(args[2]=="cluster_diseases_2"){
    load("pvclust_diffussion_profiles_ChEMBLall.RData")
    pdf(file="MultiplexAnalysis/Results/Circle_plot_with_bars_ChEMBLall.pdf")
      # 1. Prepare dendrogram and significant clusters
      dend <- as.dendrogram(clusters)
      labels <- clusters$hclust$labels
      
      # Identify significant clusters (AU p-value >= 0.95)
      signif_clusters <- pvpick(clusters, alpha = 0.95, pv = "au", type = "geq")$clusters
      
      n_signif <- length(signif_clusters)
      cluster_palette <- rainbow(n_signif)
      
      # Map each label to its cluster color or grey if not significant
      label_colors <- rep("grey60", length(labels)) # default: grey
      bar_colors   <- rep("white",   length(labels)) # default: white
      cluster_ids  <- rep(NA, length(labels)) # to store cluster number for each label
      for (i in seq_along(signif_clusters)) {
        idx <- labels %in% signif_clusters[[i]]
        label_colors[idx] <- cluster_palette[i]
        bar_colors[idx]   <- cluster_palette[i]
        cluster_ids[idx]  <- i
      }
      
      # 2. Set all branches to black, labels as above
      dend <- set(dend, "branches_col", "black")
      labels_colors(dend) <- label_colors[order.dendrogram(dend)]
      
      # 3. Order bar colors and cluster ids to match dendrogram
      ordered_labels <- labels(dend)
      bar_colors <- bar_colors[match(ordered_labels, labels)]
      cluster_ids <- cluster_ids[match(ordered_labels, labels)]
      
      n <- length(ordered_labels)
      circos.clear()
      # circos.par(cell.padding = c(0, 0, 0, 0))
      circos.par(canvas.xlim = c(-1.1, 1.1), canvas.ylim = c(-1.1, 1.1), cell.padding = c(0, 0, 0, 0))
      circos.initialize(sectors = "a", xlim = c(0, n))
      
      # OUTERMOST: cluster numbers
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                             panel.fun = function(x, y) {
                               # Example for manual adjustment:
                               for (i in seq_len(n_signif)) {
                                 inds <- which(cluster_ids == i)
                                 if (length(inds) > 0) {
                                   mean_pos <- mean(inds - 0.5)
                                   ypos<-1.1
                                   # Manually adjust for specific clusters if needed
                                   if(i==26 || i==87 || i==45 || i==145 || i==119 || i==92 || i==138 || i==89 || i==148) ypos<-ypos+1.2
                                   if(i==143 || i==125 || i==133 || i==98 || i==131 || i==62 || i==90 || i==71) ypos<-ypos+1.2
                                   if(i==132 || i==118 || i==114 || i==120 || i==116 || i==110 || i==99 || i==56) ypos<-ypos+1.2
                                   if(i==115 || i==61 || i==25 || i==109 || i==91 || i==82 || i==49 || i==73) ypos<-ypos+1.2
                                   if(i==24 || i==5 || i==53 || i==9 || i==19 || i==105 || i==104 || i==86) ypos<-ypos+1.2
                                   if(i==81 || i==74 || i==69 || i==43) ypos<-ypos+1.2
                                   circos.text(mean_pos, ypos, labels = as.character(i), 
                                               facing = "clockwise", niceFacing = TRUE,
                                               col = cluster_palette[i], cex = 0.6, adj = c(0.5, 0.5), font = 2)
                                 }
                               }
                             }
                             
      )
      
      # OUTER: colored bar/ring
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                             panel.fun = function(x, y) {
                               for (i in seq_len(n)) {
                                 circos.rect(i-1, 0, i, 1, cex=0.3, col = bar_colors[i], border = NA)
                               }
                             }
      )
      
      # MIDDLE: labels
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.01,
                             panel.fun = function(x, y) {
                               for (i in seq_len(n)) {
                                 circos.text(i-0.5, 0, ordered_labels[i], adj = c(0, 0.5),
                                             facing = "clockwise", niceFacing = TRUE,
                                             col = label_colors[match(ordered_labels[i], labels)], cex = 0.12)
                               }
                             }
      )
      
      # INNER: dendrogram (all branches black)
      max_height <- attr(dend, "height")
      circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA, track.height = 0.4,
                             panel.fun = function(x, y) {
                               circos.dendrogram(dend, max_height = max_height)
                             }
      )
    dev.off()
  }
  if(args[2]=="cluster_diseases_larger_than_6"){
    load("pvclust_diffussion_profiles_ChEMBLall.RData")
    pdf(file="MultiplexAnalysis/Results/Circle_plot_with_bars_ChEMBLall_atleast6.pdf")
      # 1. Prepare dendrogram and significant clusters
      dend <- as.dendrogram(clusters)
      labels <- clusters$hclust$labels
      
      # Identify significant clusters (AU p-value >= 0.95)
      signif_clusters <- pvpick(clusters, alpha = 0.95, pv = "au", type = "geq", max.only = FALSE)$clusters
      # Filter to keep only clusters with ≥6 elements
      signif_clusters <- signif_clusters[sapply(signif_clusters, length) >= 6]
      
      n_signif <- length(signif_clusters)
      cluster_palette <- rainbow(n_signif)
      
      # Map each label to its cluster color or grey if not significant
      label_colors <- rep("grey60", length(labels)) # default: grey
      bar_colors   <- rep("white",   length(labels)) # default: white
      cluster_ids  <- rep(NA, length(labels)) # to store cluster number for each label
      for (i in seq_along(signif_clusters)) {
        idx <- labels %in% signif_clusters[[i]]
        label_colors[idx] <- cluster_palette[i]
        bar_colors[idx]   <- cluster_palette[i]
        cluster_ids[idx]  <- i
      }
      
      # 2. Set all branches to black, labels as above
      dend <- set(dend, "branches_col", "black")
      labels_colors(dend) <- label_colors[order.dendrogram(dend)]
      
      # 3. Order bar colors and cluster ids to match dendrogram
      ordered_labels <- labels(dend)
      bar_colors <- bar_colors[match(ordered_labels, labels)]
      cluster_ids <- cluster_ids[match(ordered_labels, labels)]
      
      n <- length(ordered_labels)
      circos.clear()
      # circos.par(cell.padding = c(0, 0, 0, 0))
      circos.par(canvas.xlim = c(-1.1, 1.1), canvas.ylim = c(-1.1, 1.1), cell.padding = c(0, 0, 0, 0))
      circos.initialize(sectors = "a", xlim = c(0, n))
      
      # OUTERMOST: cluster numbers
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                             panel.fun = function(x, y) {
                               # Example for manual adjustment:
                               for (i in seq_len(n_signif)) {
                                 inds <- which(cluster_ids == i)
                                 if (length(inds) > 0) {
                                   mean_pos <- mean(inds - 0.5)
                                   ypos<-1.1
                                   circos.text(mean_pos, ypos, labels = as.character(i), 
                                               facing = "clockwise", niceFacing = TRUE,
                                               col = cluster_palette[i], cex = 0.6, adj = c(0.5, 0.5), font = 2)
                                 }
                               }
                             }
                             
      )
      
      # OUTER: colored bar/ring
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.07,
                             panel.fun = function(x, y) {
                               for (i in seq_len(n)) {
                                 circos.rect(i-1, 0, i, 1, cex=0.3, col = bar_colors[i], border = NA)
                               }
                             }
      )
      
      # MIDDLE: labels
      circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.01,
                             panel.fun = function(x, y) {
                               for (i in seq_len(n)) {
                                 circos.text(i-0.5, 0, ordered_labels[i], adj = c(0, 0.5),
                                             facing = "clockwise", niceFacing = TRUE,
                                             col = label_colors[match(ordered_labels[i], labels)], cex = 0.12)
                               }
                             }
      )
      
      # INNER: dendrogram (all branches black)
      max_height <- attr(dend, "height")
      circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA, track.height = 0.4,
                             panel.fun = function(x, y) {
                               circos.dendrogram(dend, max_height = max_height)
                             }
      )
    dev.off()
  }
  if(args[2]=="describing_clusters"){
    load("pvclust_diffussion_profiles_ChEMBLall.RData")
    sclusters <- pvpick(clusters, alpha = 0.95, pv = "au", type = "geq")$clusters
    #### 2 - Plot the size of the clusters and the heterogeneity in terms of diseases composing them ####
    ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
    tipos<-c()
    sizes<-c()
    for(a in 1:168){
      tipos<-c(tipos,length(table(substr(sclusters$clusters[[a]],1,1))))
      sizes<-c(sizes,length(sclusters$clusters[[a]]))
      if(length(which(sclusters$clusters[[a]]=="G30"))>0){print(paste("G30",a))}
      if(length(which(sclusters$clusters[[a]]=="U07"))>0){print(paste("U07",a))}
    }
    dev.off()
    hist(sizes,xlim=c(0,20),breaks = 20)
    diseases<-sclusters$clusters[[135]]
    
    #### 3 - Plot some subnetworks representing the connections among the studies diseases in the different layers ####
    ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
    ## load the colors for each node ##
    dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F,sep="\t")
    colores<-dictionary$color ; names(colores)<-dictionary$icd10
    ## load the networks and induce the subgraphs ##
    redes<-paste("MultiplexAnalysis/Networks/",list.files("MultiplexAnalysis/Networks/"),sep="")[c(7,8,1:6)]
    symp<-fread(redes[1],stringsAsFactors = F,sep="\t") ; gsymp<-graph_from_data_frame(symp,directed=F)
    trans<-fread(redes[2],stringsAsFactors = F,sep="\t") ; gtrans<-graph_from_data_frame(trans,directed=F)
    drug<-fread(redes[3],stringsAsFactors = F,sep="\t") ; gdrug<-graph_from_data_frame(drug,directed=F)
    gene<-fread(redes[4],stringsAsFactors = F,sep="\t") ; ggene<-graph_from_data_frame(gene,directed=F)
    metab<-fread(redes[5],stringsAsFactors = F,sep="\t") ; gmetab<-graph_from_data_frame(metab,directed=F)
    methy<-fread(redes[6],stringsAsFactors = F,sep="\t") ; gmethy<-graph_from_data_frame(methy,directed=F)
    micro<-fread(redes[7],stringsAsFactors = F,sep="\t") ; gmicro<-graph_from_data_frame(micro,directed=F)
    mirna<-fread(redes[8],stringsAsFactors = F,sep="\t") ; gmirna<-graph_from_data_frame(mirna,directed=F)
    
    ## induce subgraph ##
    # nodes_of_interest<-sclusters$clusters[[135]]
    nodes_of_interest<-sclusters$clusters[[148]]
    # Specify nodes of interest by name
    if(length(which(is.na(nodes_of_interest)))>0){nodes_of_interest<-nodes_of_interest[-which(is.na(nodes_of_interest))]}
    
    selnodessymp<-V(gsymp)[name %in% nodes_of_interest]
    selnodestrans<-V(gtrans)[name %in% nodes_of_interest]
    selnodesdrug<-V(gdrug)[name %in% nodes_of_interest]
    selnodesgene<-V(ggene)[name %in% nodes_of_interest]
    selnodesmetab<-V(gmetab)[name %in% nodes_of_interest]
    selnodesmethy<-V(gmethy)[name %in% nodes_of_interest]
    selnodesmicro<-V(gmicro)[name %in% nodes_of_interest]
    selnodesmirna<-V(gmirna)[name %in% nodes_of_interest]
    
    sgsymp <- induced_subgraph(gsymp, selnodessymp)
    sgtrans <- induced_subgraph(gtrans, selnodestrans)
    sgdrug <- induced_subgraph(gdrug, selnodesdrug)
    sggene <- induced_subgraph(ggene, selnodesgene)
    sgmetab <- induced_subgraph(gmetab, selnodesmetab)
    sgmethy <- induced_subgraph(gmethy, selnodesmethy)
    sgmicro <- induced_subgraph(gmicro, selnodesmicro)
    sgmirna <- induced_subgraph(gmirna, selnodesmirna)
    
    ## Add colors to the nodes ##
    V(sgsymp)$color<-as.character(colores[as_ids(V(sgsymp))])
    V(sgtrans)$color<-as.character(colores[as_ids(V(sgtrans))])
    V(sgdrug)$color<-as.character(colores[as_ids(V(sgdrug))])
    V(sggene)$color<-as.character(colores[as_ids(V(sggene))])
    V(sgmetab)$color<-as.character(colores[as_ids(V(sgmetab))])
    V(sgmethy)$color<-as.character(colores[as_ids(V(sgmethy))])
    V(sgmicro)$color<-as.character(colores[as_ids(V(sgmicro))])
    V(sgmirna)$color<-as.character(colores[as_ids(V(sgmirna))])
    
    ## Plot the networks
    plot(sgsymp, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgtrans, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgdrug, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sggene, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmetab, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmethy, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmicro, vertex.label.degree = pi/2, vertex.label.dist = 2)
    plot(sgmirna, vertex.label.degree = pi/2, vertex.label.dist = 2)
    
    #### Calculate the number of edges within the cluster compared to the ones outside the clusters ####
  }
  if(args[2]=="get_shortest_paths"){
    # Can be run from anywhere so long as RWRtoolkit is installed. 
    string.interactions.fp <- "MultiplexAnalysis/Objects/multiplex_chemblall.Rdata"
    source.geneset.path <- "MultiplexAnalysis/MultiplexICD10s_chemblall/U07.tsv"
    target.geneset.path <- "MultiplexAnalysis/MultiplexICD10s_chemblall/G30.tsv"
    outdir.path <- './MultiplexAnalysis/Results/'
    ## Get shortes paths
    shortpath<-RWRtoolkit::RWR_ShortestPaths(
      data = string.interactions.fp,
      source_geneset = source.geneset.path,
      target_geneset = target.geneset.path,
      outdir = outdir.path
    )
    #### do it layer by layer ####
    ## Microbiome ##
    microbiome<-fread("MultiplexAnalysis/Networks/microbiome.txt",stringsAsFactors = F,sep="\t")
    gmicrobiome<-graph_from_data_frame(microbiome,directed = F)
    microbiome1<-all_simple_paths(gmicrobiome, "U07", "G30", cutoff=1)
    microbiome2<-all_simple_paths(gmicrobiome, "U07", "G30", cutoff=2)
    microbiome3<-all_simple_paths(gmicrobiome, "U07", "G30", cutoff=3)
    ## Symptoms ##
    symptom<-fread("MultiplexAnalysis/Networks/symptom_2025.txt",stringsAsFactors = F,sep="\t")
    gsymptom<-graph_from_data_frame(symptom,directed = F)
    symptom1<-all_simple_paths(gsymptom, "U07", "G30", cutoff=1)
    symptom2<-all_simple_paths(gsymptom, "U07", "G30", cutoff=2)
    symptom3<-all_simple_paths(gsymptom, "U07", "G30", cutoff=3)
    ## Metabolome ##
    metabolome<-fread("MultiplexAnalysis/Networks/metabolome.txt",stringsAsFactors = F,sep="\t")
    gmetabolome<-graph_from_data_frame(metabolome,directed = F)
    metabolome1<-all_simple_paths(gmetabolome, "U07", "G30", cutoff=1)
    metabolome2<-all_simple_paths(gmetabolome, "U07", "G30", cutoff=2)
    metabolome3<-all_simple_paths(gmetabolome, "U07", "G30", cutoff=3)
    ## Drugs ##
    drug<-fread("MultiplexAnalysis/Networks/drugs_chembl.txt",stringsAsFactors = F,sep="\t")
    gdrug<-graph_from_data_frame(drug,directed = F)
    drug1<-all_simple_paths(gdrug, "U07", "G30", cutoff=1)
    drug2<-all_simple_paths(gdrug, "U07", "G30", cutoff=2)
    drug3<-all_simple_paths(gdrug, "U07", "G30", cutoff=3)
    ## Genes ##
    gene<-fread("MultiplexAnalysis/Networks/gene.txt",stringsAsFactors = F,sep="\t")
    ggene<-graph_from_data_frame(gene,directed = F)
    gene1<-all_simple_paths(ggene, "U07", "G30", cutoff=1)
    gene2<-all_simple_paths(ggene, "U07", "G30", cutoff=2)
    gene3<-all_simple_paths(ggene, "U07", "G30", cutoff=3)
    ## miRNA ##
    mirna<-fread("MultiplexAnalysis/Networks/miRNA.txt",stringsAsFactors = F,sep="\t")
    gmirna<-graph_from_data_frame(mirna,directed = F)
    mirna1<-all_simple_paths(gmirna, "U07", "G30", cutoff=1)
    mirna2<-all_simple_paths(gmirna, "U07", "G30", cutoff=2)
    mirna3<-all_simple_paths(gmirna, "U07", "G30", cutoff=3)
    
    ## create the table ##
    u07g30net<-c()
    for(a in 1:length(microbiome2)){
      u07g30net<-rbind(u07g30net,c(as_ids(microbiome2[[a]][1:2]),"microbiome"),c(as_ids(microbiome2[[a]][2:3]),"microbiome"))
    }
    for(a in 1:length(symptom2)){
      u07g30net<-rbind(u07g30net,c(as_ids(symptom2[[a]][1:2]),"symptom"),c(as_ids(symptom2[[a]][2:3]),"symptom"))
    }
    for(a in 1:length(drug2)){
      u07g30net<-rbind(u07g30net,c(as_ids(drug2[[a]][1:2]),"drug"),c(as_ids(drug2[[a]][2:3]),"drug"))
    }
    u07g30net<-rbind(u07g30net,c(as_ids(gene1[[1]]),"gene"))
    u07g30net<-rbind(u07g30net,c(as_ids(mirna1[[1]]),"miRNA"))
    u07g30net<-rbind(u07g30net,c(as_ids(metabolome1[[1]]),"metabolome"))
    colnames(u07g30net)<-c("node1","node2","layer")
    write.table(u07g30net,"MultiplexAnalysis/Results/shortest_path_net_chemblall_U07_G30.txt",quote=F,sep="\t",row.names=F)
    
    dictionary<-fread("DiseaseName_ICD10_dictionary.txt",stringsAsFactors = F,sep="\t")
    colores<-dictionary$color ; names(colores)<-dictionary$icd10
    u07g30net<-fread("MultiplexAnalysis/Results/shortest_path_net_chemblall_U07_G30.txt",stringsAsFactors = F,sep="\t")
    diseases<-unique(c(u07g30net$node1,u07g30net$node2))
    nodes<-cbind(diseases,as.character(colores[diseases]))
    colnames(nodes)<-c("nodes","colores")
    write.table(nodes,"MultiplexAnalysis/Results/nodes_U07_G30_chemblall.txt",quote=F,sep="\t",row.names=F)
  }
  if(args[2]=="calculate_similarities_between_diseases"){
    load("MultiplexAnalysis/Multiplex_ChEMBLall_diffusion_profile_matrix.Rdata") # matriz
    ## seed nodes are in columns, so, for getting the profile of a single disease we have to calculate the mean of the disease-associated icds in the columns
    ##substitute NA by 1 in the diagonal 
    for(a in 1:length(matriz[,1])){matriz[a,a]<-1}
    tmatriz<-t(matriz)
    cosinemat<-cosine(tmatriz)
    u07ranks<-sort(cosinemat[,"U07"],decreasing = T)[-1]
    ## position of the prediction ##
    which(names(u07ranks)=="G30")
    valores<-c()
    for(a in 1:(length(cosinemat[,1])-1)){
      valores<-c(valores,cosinemat[a,(a+1):length(cosinemat[,1])])
    }
    plot(density(valores),xlim=c(0,0.0002))
    abline(v=u07ranks[17])
    pval<-length(which(valores>=as.numeric(u07ranks[17])))/length(valores)
    pvalchanges<-c()
    for(a in 1:length(u07ranks)){
      pvalchanges<-c(pvalchanges,length(which(valores>=as.numeric(u07ranks[a])))/length(valores))
    }
    hist(valores,breaks = 100)
    
    hist(log(valores))
    abline(v=log(u07ranks[17]))
  }
}





