#### Xiaodong Feng 202202 ####
# Contains the scripts for using MetGUI.

setwd("D:/github/MetGUI") # the dir
roxygen2::roxygenise()
devtools::test()
devtools::install_github("xiaodfeng/DynamicTol")
.libPaths()
.libPaths("d:/PackageR/4.1/") # change the default package directory
# .libPaths(.libPaths()[1])
# Sys.setenv(JAVA_HOME='c:/Program Files/Java/jre1.8.0_191')

#### Packages ####
## Install packages from official sources
InstallRequiredPackages <- function(CRAN.packages=c(), bioconductor.packages=c()) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  for (p in bioconductor.packages) {
    if (!require(p, character.only=T)) {
      BiocManager::install(p) 
    }
    library(p,character.only=T)
  }
  
  for (p in CRAN.packages) {	
    if (!require(p, character.only=T)) { 	
      install.packages(p, dependencies = TRUE) 	
    }	
    library(p, character.only=T)
  }
}
CRAN.packages <- c('gWidgets2RGtk2','gWidgets2')
bioconductor.packages <- c(
  'blandr',
  'CAMERA','ComplexHeatmap',
  'data.table','dplyr','dendextend',
  'EnhancedVolcano',
  'fmcsR',
  'ggplot2','gtools',
  'imager',
  'MSnbase','msPurity','MetaboAnalystR','magrittr','MsBackendMgf','MetaboAnnotation',
  'parallel','pls','purrr','plyr',
  'RGtk2','RColorBrewer','readxl','RSQLite',
  'Spectra','stringr',
  'usethis',
  'VennDiagram',
  'writexl',
  'xcms','xlsx'
  )

InstallRequiredPackages(CRAN.packages, bioconductor.packages)

## Install packages from other sources
# install.packages("https://cran.microsoft.com/snapshot/2021-12-15/bin/windows/contrib/4.1/RGtk2_2.20.36.2.zip", repos=NULL) # https://github.com/BroVic/RQDAassist/issues/4
# devtools::install_github("computational-metabolomics/msPurity")
# remove.packages('msPurity')
# detach("package:msPurity")
# BiocManager::install("moments")
packageVersion("msPurity")
# install.packages("ComplexHeatmap")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.13")
# BiocManager::valid()
# BiocManager::available("xcms")
# BiocManager::install("msPurity",update = F)
# library(circlize)
# library(ChemmineR)
# library(ChemmineOB)
# library(cowplot)
# library(ggrepel)
# library(imager)
# library(multtest)
# library(mzR)
# library(metfRag)
# library(msdata)
# library(pander)
# library(parallel)
# library(pheatmap)
# library(Rcpp)
# library(rJava)
# detach("package:MetaboAnalystR")
# remove.packages(c("xcms"))
# metanr_packages <- function(){
#   
#   metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn","httr","qs")
#   
#   list_installed <- installed.packages()
#   
#   new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
#   
#   if(length(new_pkgs)!=0){
#     
#     if (!requireNamespace("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#     BiocManager::install(new_pkgs)
#     print(c(new_pkgs, " packages added..."))
#   }
#   
#   if((length(new_pkgs)<1)){
#     print("No new packages added...")
#   }
# }
# metanr_packages()
# install("d:/PackageR/MetaboAnalystR_3.0.3.1")
# install("d:/github/DynamicXCMS/DynamicXCMS3.16",update = F)
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
# BiocManager::install("fmcsR",update=FALSE)
# install("d:/software/MetaboAnalystR_3.0.3.1")

#### Functions ####
source('d:/github/dynamic/R/Functions_Linux.R')
F_readMSData<-function(Dir,Group){
  cdfs <- dir(Dir, full.names = TRUE,recursive = TRUE)
  print(cdfs)
  sample_name<-sub(basename(cdfs), pattern = ".mzML",replacement = "", fixed = TRUE)
  pd <- data.frame(sample_name = sample_name,
                   sample_group = Group,
                   stringsAsFactors = FALSE)
  print(pd)
  Raw <- MSnbase::readMSData(files = cdfs, pdata = new("NAnnotatedDataFrame", pd),mode = "onDisk")
  return(Raw)
}
F_readMSDataMix12<-function(Dir){
  myFiles <- dir(Dir, full.names = TRUE,recursive = TRUE) %>% mixedsort(.)
  sample_name<-sub(basename(myFiles), pattern = ".mzML",replacement = "", fixed = TRUE)
  pd <- data.frame(sample_name = sample_name,
                   sample_group = sample_name,
                   stringsAsFactors = FALSE)
  Raw <- readMSData(files = myFiles, pdata = new("NAnnotatedDataFrame", pd),mode = "onDisk")
  return(Raw)
}



F_findChromPeaks<-function(Raw,cwp,ParallelNum=3){
  ## Parallel
  # BPPARAM  = SnowParam(detectCores()-ParallelNum, progressbar = TRUE)
  ## Check .mzML to set proper noise threshold, prefilter, and peakwidth
  xdata <- findChromPeaks(Raw, param = cwp) #,BPPARAM = BPPARAM
  return(xdata)
}
F_adjustRtime<-function(xdata,param){
  # BPPARAM  = SnowParam(detectCores()-ParallelNum, progressbar = TRUE)
  # xdata <- dropAdjustedRtime(xdata) 
  xdata <- adjustRtime(xdata, param = param) 
  return(xdata)
}
F_adjustRtimeEffect<-function(Raw,xdata,Name){
  ## Base peak chromatograms before RT alignment
  bpis <- chromatogram(Raw, aggregationFun = "max")
  png(file=paste("RT alignment before bpis",Name,".png",seq=""),width=700,height=480) #res=300,
  plot(bpis, col = group_colors[bpis$sample_group])#
  legend("topright",legend=unique(bpis$sample_group),col=unique(group_colors[bpis$sample_group]), lty=1, cex=0.4,box.lty=0) #
  title(sub=paste("RT alignment before bpis",Name,seq=""))
  dev.off()
  # Base peak chromatograms after RT alignment
  bpis_adj <- chromatogram(xdata, aggregationFun = "max") 
  return(bpis_adj)
}
F_GroupFill<-function(xdata,pdp){
  xdata <- groupChromPeaks(xdata, param = pdp)
  xdata <- fillChromPeaks(xdata)
  return(xdata)
}




F_boxplotmSet<-function(mSet,group_colors){
  NewNameNoB<-mSet[["dataSet"]][["prenorm.smpl.nms"]]
  NewGroupNoB<-mSet[["dataSet"]][["cls"]]
  names(group_colors)<-unique(NewGroupNoB)
  # Boxplot after MedianAll adjustment
  Intensity<-do.call(rbind,mSet[["dataSet"]][["row.norm"]])
  MedianAll<-median(do.call(rbind,mSet[["dataSet"]][["proc"]]),na.rm = T)
  png(file=paste("Boxplot MedianAll adjustment.png",seq=""),width=700,height=480) #res=300,
  boxplot(log2(Intensity*MedianAll),ylab = expression(intensity), main = "Median adjustment",cex.axis=0.5, varwidth = TRUE, 
          col = group_colors[NewGroupNoB],names=NewNameNoB)
  legend("topright",legend=unique(NewGroupNoB),
         col=unique(group_colors[NewGroupNoB]), lty=1, cex=0.5,box.lty=0) 
  dev.off()
  # Boxplot after scaling
  Intensity<-do.call(rbind,mSet[["dataSet"]][["norm"]])
  png(file=paste("Boxplot scaling adjustment.png",seq=""),width=700,height=480) #res=300,
  boxplot(Intensity,ylab = expression(intensity), main = "Scaling adjustment",cex.axis=0.5, varwidth = TRUE, 
          col = group_colors[NewGroupNoB],names=NewNameNoB)
  legend("topright",legend=unique(NewGroupNoB),
         col=unique(group_colors[NewGroupNoB]), lty=1, cex=0.5,box.lty=0) 
  dev.off()
  
}

mSet<-F_PlotHeatMap2(mSet, paste0(Name,TopFeature,"_heatmap2_"), "png", 300, width=NA, 
                     smplDist="euclidean",clstDist = "average",colors ="bwm",
                     viewOpt ="overview", hiRes=T, sortInx=1, useSigFeature =TopFeature, drawBorder =FALSE)

F_PlotHeatMap2 <- function (mSetObj = NA, imgName, format = "png", dpi = 72, 
                            width = NA, smplDist = "pearson", clstDist = "average", 
                            colors = "bwm", viewOpt = "overview", hiRes = FALSE, 
                            sortInx = 1, useSigFeature, drawBorder) 
{
  # mSetObj=mSet; imgName=paste0(Name,TopFeature,"_heatmap2_");format='png';dpi=300;
  # width=NA;smplDist="euclidean";clstDist = "average";colors ="bwm";
  # viewOpt ="overview";hiRes=T; sortInx=1; useSigFeature =TopFeature; drawBorder =FALSE
  # mSetObj <- .get.mSet(mSetObj)
  facA <- mSetObj$dataSet$facA
  facB <- mSetObj$dataSet$facB
  if (sortInx == 1) {
    ordInx <- order(facA, facB)
  }
  else {
    ordInx <- order(facB, facA)
  }
  new.facA <- facA[ordInx]
  new.facB <- facB[ordInx]
  data <- mSetObj$dataSet$norm[ordInx, ]
  # if (useSigFeature) {
  #   hits <- colnames(data) %in% rownames(mSetObj$analSet$aov2$sig.mat)
  #   data <- mSetObj$dataSet$norm[ordInx, hits]
  # }
  Sigmat <- data.frame(mSetObj$analSet$aov2$sig.mat)
  # names(Sigmat)
  # rownames(Sigmat)
  setorder(Sigmat,'Interaction.raw.p.')
  # useSigFeature <- 100
  SigmatSel <- Sigmat[1:useSigFeature,]
  hits <- colnames(data) %in% rownames(SigmatSel)
  data <- mSetObj$dataSet$norm[ordInx, hits]
  hc.dat <- as.matrix(data)
  colnames(hc.dat) <- substr(colnames(data), 1, 18)
  if (colors == "gbr") {
    colors <- colorRampPalette(c("green", "black", 
                                 "red"), space = "rgb")(256)
  }
  else if (colors == "heat") {
    colors <- heat.colors(256)
  }
  else if (colors == "topo") {
    colors <- topo.colors(256)
  }
  else if (colors == "gray") {
    colors <- colorRampPalette(c("grey90", "grey10"), 
                               space = "rgb")(256)
  }
  else {
    colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, 
                                                            "RdBu"))(256))
  }
  if (drawBorder) {
    border.col <- "grey60"
  }
  else {
    border.col <- NA
  }
  imgName = paste(imgName, "dpi", dpi, ".", format, 
                  sep = "")
  if (viewOpt == "overview") {
    if (is.na(width)) {
      w <- 9
    }
    else if (width == 0) {
      w <- 7.2
    }
    else {
      w <- 7.2
    }
    mSetObj$imgSet$htmaptwo <- imgName
    h <- w
  }
  else {
    if (is.na(width)) {
      minW <- 650
      myW <- nrow(hc.dat) * 11 + 150
      if (myW < minW) {
        myW <- minW
      }
      w <- round(myW/72, 2)
    }
    else if (width == 0) {
      w <- 7.2
    }
    else {
      w <- 7.2
    }
    if (ncol(hc.dat) > 100) {
      myH <- ncol(hc.dat) * 12 + 120
    }
    else if (ncol(hc.dat) > 50) {
      myH <- ncol(hc.dat) * 12 + 60
    }
    else {
      myH <- ncol(hc.dat) * 12 + 20
    }
    mSetObj$imgSet$htmaptwo <- imgName
    h <- round(myH/72, 2)
  }
  if (format == "pdf") {
    pdf(file = imgName, width = w, height = h, bg = "white", 
        onefile = FALSE)
  }
  else {
    Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, 
                 width = w, height = h, type = format, bg = "white")
  }
  annotation <- data.frame(new.facB, new.facA)
  colnames(annotation) <- c(mSetObj$dataSet$facB.lbl, mSetObj$dataSet$facA.lbl)
  rownames(annotation) <- rownames(hc.dat)
  pheatmap::pheatmap(t(hc.dat), annotation = annotation, fontsize = 8, 
                     fontsize_row = 8, clustering_distance_rows = smplDist, 
                     clustering_distance_cols = smplDist, clustering_method = clstDist, 
                     border_color = border.col, cluster_rows = TRUE, cluster_cols = FALSE, 
                     scale = "row", color = colors)
  dev.off()
  mSetObj$analSet$htmap2 <- list(dist.par = smplDist, clust.par = clstDist)
  return(mSetObj)
}

F_EICLoopPureAA<-function(DT,massI,addI,mzI,rtI,Raw){
  group_colors <- c("Orange","Black")#brewer.pal(2, "Set3")
  names(group_colors) <- unique(Raw$sample_group)
  F_EICSingle<-function(x){ 
    ExactMass<-as.numeric(x[[massI]])
    Adduct<-x[[addI]]
    mz<-as.numeric(x[[mzI]])
    rt<-as.numeric(x[[rtI]])
    mzValue<-c(mz-0.02,mz+0.02)
    rtValue<-c(rt-30,rt+30)
    EIC<-chromatogram(Raw, aggregationFun = "sum",mz = mzValue, rt = rtValue)
    name<-paste("mass",ExactMass,Adduct,"mz",mz,"rt",rt,"Pure AA.png",seq="")
    png(file=name,res=200,width = 1000,height = 1000)
    plot(EIC, col = group_colors[Raw$sample_group],peakType = "none")
    legend("topright",legend=unique(Raw$sample_group),col=unique(group_colors[Raw$sample_group]), lty=1, cex=0.8,box.lty=0) 
    title(sub=paste("mass",ExactMass,Adduct,"mz",mz,"rt",rt,seq=""))
    dev.off()
  }
  apply(DT,1,F_EICSingle)
}
F_ReplaceMZRT<-function(l,lmz,lrt,r,rmz,rrt,tolmz = 0.001,tolrt = 5){ 
  l<-as.data.frame(l)
  r<-as.data.frame(r)
  for (i in 1:nrow(r)){
    l[, lmz][(abs(r[, rmz][i] - l[,lmz])<tolmz)&(abs(r[, rrt][i] - l[,lrt])<tolrt)] <- r[,rmz][i]
    l[, lrt][(abs(r[, rmz][i] - l[,lmz])<tolmz)&(abs(r[, rrt][i] - l[,lrt])<tolrt)] <- r[,rrt][i]
  }
  return(data.table(l))
}
F_ReplaceMZRTRange<-function(l,lmz,lmzmin,lmzmax,lrt,lrtmin,lrtmax,r,rmz,rrt,tolmz = 0.001,tolrt = 5){ 
  l<-as.data.frame(l)
  r<-as.data.frame(r)
  for (i in 1:nrow(r)){
    l[, lmz][(r[, rmz][i]>(l[,lmzmin]-tolmz ))&(r[, rmz][i]<(l[,lmzmax]+tolmz ))&(r[, rrt][i]>l[,lrtmin]-tolrt)&(r[, rrt][i]<l[,lrtmax]+tolrt)] <- r[,rmz][i]
    l[, lrt][(r[, rmz][i]>(l[,lmzmin]-tolmz ))&(r[, rmz][i]<(l[,lmzmax]+tolmz ))&(r[, rrt][i]>l[,lrtmin]-tolrt)&(r[, rrt][i]<l[,lrtmax]+tolrt)] <- r[,rrt][i]
  }
  return(data.table(l))
}
F_AddStatistics<-function(MCRep,Feature,DTNoB){
  ## According to the mz and rt in MC, find the statistics in Feature
  MCRepFT<-left_join(MCRep,Feature,by=c("mz","rt")) %>% data.table(.)
  ## According to the Feature index, find the quantification in DTNoB
  DTNoB<-DTNoB[-1,]
  setnames(DTNoB,"Sample","rownames")
  MCRepFT<-left_join(MCRepFT,DTNoB,by="rownames") %>% data.table(.)
  return(MCRepFT)
}

F_Boxplot<-function(Intensity,Name){
  boxplot(Intensity,ylab = expression(intensity), main = Name,cex.axis=0.5, varwidth = TRUE, 
          col = group_colors[xdata1_NoB$sample_group],
          names=NewNameNoB)
  legend("topright",legend=unique(xdata1_NoB$sample_group),
         col=unique(group_colors[xdata1_NoB$sample_group]), lty=1, cex=0.5,box.lty=0) 
}

F_BoxplotNeg<-function(Intensity,Name){
  boxplot(Intensity,ylab = expression(intensity), main = Name,cex.axis=0.5, varwidth = TRUE, 
          col = group_colors[NewGroupNoB],
          names=NewNameNoB)
  legend("topright",legend=unique(NewGroupNoB),
         col=unique(group_colors[NewGroupNoB]), lty=1, cex=0.5,box.lty=0) 
}
#group_colors[as.matrix(NewGroupNoB)]

F_Plot24HeatMap<-function(cormat){
  sample_number<-c(1:18)
  colnames(cormat) <- rownames(cormat) <- sample_number
  ## Define which phenodata columns should be highlighted in the plot
  ann <- data.frame(group = Raw_data_MS1_NoB$sample_group)
  rownames(ann) <- sample_number
  ## Perform the Hierarchical  cluster analysis
  group_colors <- c("coral1","coral2","coral3","cyan1","cyan2","cyan3")#brewer.pal(2, "Set3")
  names(group_colors) <- unique(Raw_data_MS1_NoB$sample_group)
  bk <- seq(0,1,by=0.02) # Set fixed indicater
  pheatmap(cormat, annotation = ann,clustering_method="average",annotation_color = list(group = group_colors),
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(bk)),
           #color=c(colorRampPalette(colors=c("blue","white"))((length(bk))/2),colorRampPalette(color=c("white","red"))((length(bk))/2)),
           legend_breaks=seq(0,1,0.1),breaks=bk) #Set fixed indicator
}
F_EICAA<-function(mass,rt,raw_data1){ # plot EIC according to the selected mz for AA
  mass<-as.numeric(mass)
  rt<-as.numeric(rt)
  rtValue<-c(rt-30,rt+30)
  MH<-mass+1.0073
  EIC<-chromatogram(raw_data1, aggregationFun = "sum",mz = c(MH-0.01,MH+0.01),rt=rtValue)#
  jpeg(file=paste("Mass",mass,"RT",rt,"M+H",MH,".jpeg",seq=""))
  plot(EIC, col = group_colors[raw_data1$sample_group])#
  legend("topright",legend=unique(raw_data1$sample_group),col=group_colors[unique(raw_data1$sample_group)], lty=1, cex=0.8,box.lty=0) #
  title(sub=paste("Mass",mass,"RT",rt,"M+H",MH,seq=""))
  dev.off()
  
  MNa<-mass+22.9892
  EIC<-chromatogram(raw_data1, aggregationFun = "sum",mz = c(MNa-0.01,MNa+0.01),rt=rtValue)#
  jpeg(file=paste("Mass",mass,"RT",rt,"M+Na",MNa,".jpeg",seq=""))
  plot(EIC, col = group_colors[raw_data1$sample_group])#
  legend("topright",legend=unique(raw_data1$sample_group),col=group_colors[unique(raw_data1$sample_group)], lty=1, cex=0.8,box.lty=0) #
  title(sub=paste("Mass",mass,"RT",rt,"M+Na",MNa,seq=""))
  dev.off()
  MK<-mass+38.9632
  EIC<-chromatogram(raw_data1, aggregationFun = "sum",mz = c(MK-0.01,MK+0.01),rt=rtValue)#
  jpeg(file=paste("Mass",mass,"RT",rt,"M+K",MK,".jpeg",seq=""))
  plot(EIC, col = group_colors[raw_data1$sample_group])#
  legend("topright",legend=unique(raw_data1$sample_group),col=group_colors[unique(raw_data1$sample_group)], lty=1, cex=0.8,box.lty=0) #
  title(sub=paste("Mass",mass,"RT",rt,"M+K",MK,seq=""))
  dev.off()
  MNH4<-mass+18.0338
  EIC<-chromatogram(raw_data1, aggregationFun = "sum",mz = c(MNH4-0.01,MNH4+0.01),rt=rtValue)#
  jpeg(file=paste("Mass",mass,"RT",rt,"M+NH4",MNH4,".jpeg",seq=""))
  plot(EIC, col = group_colors[raw_data1$sample_group])#
  legend("topright",legend=unique(raw_data1$sample_group),col=group_colors[unique(raw_data1$sample_group)], lty=1, cex=0.8,box.lty=0) #
  title(sub=paste("Mass",mass,"RT",rt,"M+NH4",MNH4,seq=""))
  dev.off()
}

F_VolcanoPlot<-function(DT,Name,FC=0.6,P=0.05,NSAA='pink',SigAA='red2',NSMet='lightblue',SigMet='royalblue'){ 
  # Customize colour
  keyvals.colour <- ifelse(
    (abs(DT$log2FoldChange) > FC)&(abs(DT$pvalue) < P)&(is.na(DT$ShortName)), SigMet, # 
    ifelse((abs(DT$log2FoldChange) > FC)&(abs(DT$pvalue) < P)&(!is.na(DT$ShortName)), SigAA,
           ifelse((!is.na(DT$ShortName)), NSAA,
                  NSMet))) 
  keyvals.colour[is.na(keyvals.colour)] <-NSMet
  names(keyvals.colour)[keyvals.colour == NSAA] <- 'AAs NS'
  names(keyvals.colour)[keyvals.colour == SigAA] <- 'AAs Significant'
  names(keyvals.colour)[keyvals.colour == NSMet] <- 'Metabolites NS'
  names(keyvals.colour)[keyvals.colour == SigMet] <- 'Metabolites Significant'
  # Customize shape
  res<-boxplot.stats(log2(DT$Mean))
  LowerHinge<-res$stats[2]
  UpperHinge<-res$stats[4]
  Log2<-log2(DT$Mean)
  keyvals.shape <- ifelse(
    Log2 <= LowerHinge, 46,
    ifelse((Log2 > LowerHinge)&(Log2 < UpperHinge), 20,
           19))
  keyvals.shape[is.na(keyvals.shape)] <- 46
  names(keyvals.shape)[keyvals.shape == 46] <- paste(" 0 < log2(Int) <=",round(LowerHinge,digits = 0)) 
  names(keyvals.shape)[keyvals.shape == 20] <- paste(round(LowerHinge,digits = 0),"< log2(Int) <=",round(UpperHinge,digits = 0)) 
  names(keyvals.shape)[keyvals.shape == 19] <- paste(round(UpperHinge,digits = 0),"< log2(Int) < Inf") 
  EnhancedVolcano(DT,x = "log2FoldChange",y = "pvalue",
                  pCutoff = 0.05,FCcutoff = 0.6,
                  xlim = c(min(DT[,log2FoldChange], na.rm=TRUE), max(DT[,log2FoldChange], na.rm=TRUE)),
                  ylim = c(0, max(-log10(DT[,pvalue]), na.rm=TRUE) + 1),
                  pointSize = 4,
                  # lab = DT$ShortName, 
                  # selectLab = DT[,ShortName][which(names(keyvals.colour) %in% c('AAs Significant'))],
                  #labSize = 4, drawConnectors = TRUE,
                  lab = DT$name, 
                  selectLab = DT[,name][which(names(keyvals.colour) %in% c('Metabolites Significant'))],
                  labSize = 2,drawConnectors = FALSE,
                  shapeCustom = keyvals.shape,
                  colCustom = keyvals.colour,colAlpha = 0.7,
                  title = Name,subtitle = "",
                  caption = "Normal / HYP, log2FC 0.6, pvalue 0.05",
                  legendPosition = "right", legendLabSize = 7,
                  gridlines.major=FALSE,gridlines.minor=FALSE,
                  widthConnectors = 0.5)
  ggsave(file=paste("VolcanoPlot",Name,".jpeg"),dpi=300)
}
# F_VolcanoPlot(Pos,"Positive",FC=0.6,P=0.05,NSAA='pink',SigAA='red2',NSMet = 'white',SigMet = 'white')
# F_VolcanoPlot(Neg,"Negative",FC=0.6,P=0.05,NSAA='pink',SigAA='red2')

F_Bland<-function(x,y,Label,Name){
  Stat<-blandr.statistics(x,y)
  DT<-data.table("Means"=Stat$means,"Differences"=Stat$differences,"Label"=Label)
  ggplot(DT,aes(x=Means,y=Differences,label=Label))+
    geom_point(size=3,alpha=0.7,color="black")+
    geom_text_repel() +
    geom_hline(yintercept=0, color="black")+
    geom_hline(yintercept=Stat$upperLOA, linetype="dashed",color="green")+
    geom_hline(yintercept=Stat$bias, linetype="dashed",color="blue")+
    geom_hline(yintercept=Stat$lowerLOA, linetype="dashed",color="red")+
    labs(title=name)+
    theme_classic() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
          axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),plot.title =element_text(size=17,face="bold"),
          legend.text=element_text(size=12),legend.justification=c(1,0),legend.position=c(1,0),legend.title=element_blank() )
  # xlim(8,24)+ylim(-0.5,0.5)
  ggsave(file=paste("Bland",Name,".jpeg"),dpi=300)
  # return(g_plot)
}
F_Adduct<-function(DT){
  # Negative
  M_H <- data.table(DT,Adduct="M-H") %>% .[,Adductmz:=(monoiso-1.007276)]
  M_HH2O <- data.table(DT,Adduct="M-H2O-H") %>% .[,Adductmz:=(monoiso-19.01839)]
  MNa2H <- data.table(DT,Adduct="M+Na-2H") %>% .[,Adductmz:=(monoiso+20.974666)]
  MCl <- data.table(DT,Adduct="M+Cl") %>% .[,Adductmz:=(monoiso+34.969402)]
  MFAH <- data.table(DT,Adduct="M+FA-H") %>% .[,Adductmz:=(monoiso+44.998201)]
  AdductNeg<-rbind(M_H,M_HH2O,MNa2H,MCl,MFAH)
  AdductNeg[,mode:='Negative']
  # Positive
  MH <- data.table(DT,Adduct="M+H") %>% .[,Adductmz:=(monoiso+1.007276)]
  MNa <- data.table(DT,Adduct="M+Na") %>% .[,Adductmz:=(monoiso+22.989218)]
  MK <- data.table(DT,Adduct="M+K") %>% .[,Adductmz:=(monoiso+38.963158)]
  MNH4 <- data.table(DT,Adduct="M+NH4") %>% .[,Adductmz:=(monoiso+18.033823)]
  M2NaH <- data.table(DT,Adduct="M+2Na-H") %>% .[,Adductmz:=(monoiso+44.971160)]
  AdductPos<-rbind(MH,MNa,MK,MNH4,M2NaH)
  AdductPos[,mode:='Positive']
  # Combine
  Adduct<-rbind(AdductNeg,AdductPos)
  Adduct[,mz:=Adductmz]
  setorder(Adduct,monoiso,Adductmz)
  return(Adduct)
}

F_Align <- function(top, bottom,topindex=33, bottomindex=1,ppm = 5){
  ## mz and intensity extraction
  top <- data.frame(top) 
  bottom <- data.frame(bottom)
  ## find in common
  for (i in 1:nrow(bottom)) 
    top[, topindex][abs(bottom[, bottomindex][i] - top[,topindex]) < ppm*bottom[,bottomindex][i]/1000000] <- bottom[,bottomindex][i]
  ## merge according to mz
  alignment <- data.table(left_join(top, bottom, by = 'mz'))
  # colnames(alignment) <- c("mz","intensity.top", "intensity.bottom")
  # setorder(alignment,-into)
  alignment<-alignment[!is.na(rt)] 
  alignment<-alignment[mode==PeakMode] 
  print(unique(alignment$PRIMARY_NAME))
  print(alignment)
  return(alignment)
}

F_AlignMix12<-function(top,bottom,topindex=33, bottomindex=1,ppm = 7){
  Mix1<-F_Align(top[Mix==1],bottom[sample==1],topindex,bottomindex,ppm)
  Mix2<-F_Align(top[Mix==2],bottom[sample==2],topindex,bottomindex,ppm)
  Mix3<-F_Align(top[Mix==3],bottom[sample==3],topindex,bottomindex,ppm)
  Mix4<-F_Align(top[Mix==4],bottom[sample==4],topindex,bottomindex,ppm)
  Mix5<-F_Align(top[Mix==5],bottom[sample==5],topindex,bottomindex,ppm)
  Mix6<-F_Align(top[Mix==6],bottom[sample==6],topindex,bottomindex,ppm)
  
  Mix7<-F_Align(top[Mix==7],bottom[sample==7],topindex,bottomindex,ppm)
  Mix8<-F_Align(top[Mix==8],bottom[sample==8],topindex,bottomindex,ppm)
  Mix9<-F_Align(top[Mix==9],bottom[sample==9],topindex,bottomindex,ppm)
  Mix10<-F_Align(top[Mix==10],bottom[sample==10],topindex,bottomindex,ppm)
  Mix11<-F_Align(top[Mix==11],bottom[sample==11],topindex,bottomindex,ppm)
  Mix12<-F_Align(top[Mix==12],bottom[sample==12],topindex,bottomindex,ppm)
  
  MatchedPeaks<-rbind(Mix1,Mix2,Mix3,Mix4,Mix5,Mix6,Mix7,Mix8,Mix9,Mix10,Mix11,Mix12)
  MatchedPeaks[,Deltamz:=abs(Adductmz-mz)]
  MatchedPeaks[,Deltappm:=(Deltamz/mz)*1000000]
  return(MatchedPeaks)
}
F_AlignAA<-function(top,bottom,topindex=7, bottomindex=1,ppm = 7){
  MatchedPeaks<-F_Align(top,bottom,topindex,bottomindex,ppm)
  MatchedPeaks[,Deltamz:=abs(Adductmz-mz)]
  MatchedPeaks[,Deltappm:=(Deltamz/mz)*1000000]
  return(MatchedPeaks)
}
F_EICLoop<-function(DT,massI,compoundI,addI,mzI,rtI,Raw){
  group_colors <- brewer.pal(12, "Set3")
  names(group_colors) <- unique(Raw$sample_group)
  F_EICSingle<-function(x){ 
    ExactMass<-as.numeric(x[[massI]])
    Adduct<-x[[addI]]
    Compound<-x[[compoundI]]
    mz<-as.numeric(x[[mzI]])
    rt<-as.numeric(x[[rtI]])
    mzValue<-c(mz-0.02,mz+0.02)
    rtValue<-c(rt-30,rt+30)
    EIC<-chromatogram(Raw, aggregationFun = "sum",mz = mzValue, rt = rtValue)
    name<-paste("mass",ExactMass,Adduct,"mz",mz,"rt",rt,Compound,"mixture.png",seq="")
    png(file=name,res=200,width = 1000,height = 1000)
    plot(EIC, col = group_colors[Raw$sample_group],peakType = "none")
    legend("topright",legend=unique(Raw$sample_group),col=unique(group_colors[Raw$sample_group]), lty=1, cex=0.8,box.lty=0) 
    title(sub=paste("mass",ExactMass,Adduct,"mz",mz,"rt",rt,Compound,seq=""))
    dev.off()
  }
  apply(DT,1,F_EICSingle)
}






#### Install GUI ####
# library(devtools)
# install("d:/github/GUI/MetGUI",upgrade = "never")
# library(MetGUI)
# MetGUI::MetGUI()
# options("guiToolkit"="RGtk2")
# options('gWidgets2RGtk2')
# require('gWidgets2RGtk2')
# library(statTarget)
# statTarget::statTargetGUI()






