#### Packages and Functions ####
# This code is to integrate the codes from both Positive and Negative mode from HILIC
# Use Msconvert to transfer dataset from .raw into .mzML, spread the .mzML dataset into different groups
library(data.table)
library(dplyr)
library(ggplot2)
library(MSnbase)
#### Fig. 3 MS2 distribution Fig. S1 based on PNNL dataset ####
## Prepare PNNL dataset
Path <- 'd:/OneDrive/github/GUI/'
Output <- 'd:/OneDrive/github/GUI/output'
setwd(Output)
PNNL_Pos <- MSnbase::readMgfData(paste0(Path,"/input/Mgf/PNNL-LIPIDS-POSITIVE.mgf"))
PNNL_Neg <- MSnbase::readMgfData(paste0(Path,"/input/Mgf/PNNL-LIPIDS-NEGATIVE.mgf"))
nrow(fData(PNNL_Pos)) #30582 ## Data Information
nrow(fData(PNNL_Neg)) #16142
# write.csv(fData(PNNL_Pos), "fData(PNNL_Pos).csv")
# write.csv(fData(PNNL_Neg), "fData(PNNL_Neg).csv")
Fdata_Pos <- data.table(read.csv("fData(PNNL_Pos).csv"))
Fdata_Neg <- data.table(read.csv("fData(PNNL_Neg).csv"))
Fdata_Pos[(SOURCE_INSTRUMENT == "LC-ESI-CID; Lumos"), ] #538 # LC-ESI-CID; Lumos LC-ESI-HCD; Lumos LC-ESI-CID; Velos LC-ESI-HCD; Velos
Fdata_Pos[(SOURCE_INSTRUMENT == "LC-ESI-HCD; Lumos"), ] #537
Fdata_Pos[(SOURCE_INSTRUMENT == "LC-ESI-CID; Velos"), ] #14754
Fdata_Pos[(SOURCE_INSTRUMENT == "LC-ESI-HCD; Velos"), ] #14753
CIDPos <- Fdata_Pos[(SOURCE_INSTRUMENT == "LC-ESI-CID; Velos"), ]
HCDPos <- Fdata_Pos[(SOURCE_INSTRUMENT == "LC-ESI-HCD; Velos"), ]
CIDNeg <- Fdata_Neg[(SOURCE_INSTRUMENT == "LC-ESI-CID; Velos"), ] #7784
HCDNeg <- Fdata_Neg[(SOURCE_INSTRUMENT == "LC-ESI-HCD; Velos"), ] #7783
F_ExtractGNPSMS2 <- function(x) { # Extract MS2 according to index
  ## Extract all MS2 spectra that were associated with the candidate using the Index
  # print(x)
  Spectra <- GNPSSpectra[[x]]
  PeakList_matrix <- matrix(c(Spectra@mz, Spectra@intensity), ncol = 2, byrow = FALSE)
  DT <- data.table(PeakList_matrix)
  DT[, Index := x]
  # print(PeakList_matrix)
  return(DT)
}
F_SpectraToMS2 <- function(Fdata, Name, PNNL) {
  # OptWeiCIDPos <- F_MSMSDistribution(Fdata=CIDPos, Name='CID Pos', PNNL=PNNL_Pos, binwidth = 10)
  ## Extract Spectra
  GNPSSpectra <<- spectra(PNNL)
  Name <<- Name
  # Cutoff <<- Cutoff
  results <- data.frame()
  results <- apply(Fdata[, 1], 1, F_ExtractGNPSMS2)
  results <- plyr::ldply (results, data.frame) #Transfer results from list to data.frame
  DT <- data.table(results) %>% setnames(., c('V1', 'V2'), c('mz', 'intensity'))
  ## Add the compounds names for unique purpose
  IndexName <- Fdata[, c('X', 'NAME')] %>% setnames(., 'X', 'Index')
  DT <- left_join(DT, IndexName, by = 'Index') %>% data.table(.)
  return(DT)
}
RawMS2CIDPos <- F_SpectraToMS2(Fdata=CIDPos, Name='CID Pos', PNNL=PNNL_Pos)
RawMS2HCDPos <- F_SpectraToMS2(HCDPos, 'HCD Pos', PNNL_Pos)
RawMS2CIDNeg <- F_SpectraToMS2(CIDNeg, 'CID Neg', PNNL_Neg)
RawMS2HCDNeg <- F_SpectraToMS2(HCDNeg, 'HCD Neg', PNNL_Neg)

F_MzIntDis <- function(Name, DT,statsValue) {
  DT[, mzRound := round(mz)]
  setorder(DT, -intensity) # sort intensity from big to small
  MS2 <- unique(DT, by = c('NAME', 'mzRound')) # for the same compound, only keep one scan at the same mzRound
  ## Boxplot visualization and use extreme of the lower whisker as cutoff
  Cutoff <<- boxplot.stats(MS2$intensity)$stats[statsValue] %>% round(.,digits = 1)
  MS2 <- MS2[intensity >= Cutoff]
  ## Plot MSMS distribution
  svg(file = paste(Name,'Cutoff', Cutoff,"Number MS2", nrow(MS2), 'MZ Distribution.svg'))
  print(ggplot(MS2, aes(mz)) + geom_histogram(colour = "black", binwidth = 10) +
          labs(title = Name) + theme_classic()+
          theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),plot.title =element_text(size=17,face="bold"))
  )#+xlim(0,xlimrange)
  dev.off()
  ## Plot Intensity distribution
  svg(file = paste(Name,'Cutoff', Cutoff,"Number MS2", nrow(MS2), 'Intensity Distribution.svg'))
  print(ggplot(MS2, aes(log2(intensity))) +
          geom_histogram(colour = "black", binwidth = 1) +
          labs(title = Name) + theme_classic() + xlim(0,30)+
          theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),plot.title =element_text(size=17,face="bold"))
  )
  dev.off()
}
F_MzIntDis(Name = "CID Pos",DT=RawMS2CIDPos,statsValue= 2)
F_MzIntDis(Name = "HCD Pos",DT=RawMS2HCDPos,statsValue= 2)
F_MzIntDis(Name = "CID Neg",DT=RawMS2CIDNeg,statsValue= 2)
F_MzIntDis(Name = "HCD Neg",DT=RawMS2HCDNeg,statsValue= 2)

#### Fig. S2 mz distribution and intensity distribution scaling factor based on the Hema dataset ####
## Prepare the Hematology dataset
BiocParallel::register(SerialParam()) #disable paralell
dda_2020 <- xcms::featureSpectra(xdata_MixAA_Pos,return.type = "Spectra")
ex_spectra <- combineSpectra(dda_2020,  ppm = 40, peaks = "union",
                             intensityFun = median, mzFun = median,f = dda_2020$feature_id)
ex_spectra <- setBackend(ex_spectra, MsBackendDataFrame())
dda_2020_meta <- spectraData(ex_spectra,spectraVariables(ex_spectra)) %>% as.data.table(.)
dda_2020_meta <- unique(dda_2020_meta,by='feature_id')
F_ExtractMS2 <- function(x) {
  # print('Extract MS2 according to NeutralPrecursorMass')
  ## Extract all MS2 spectra that were associated with the candidate using the ID of Feature (rownames)
  # x[1] <- 'FT0020'
  ex_spectrum<- ex_spectra[ex_spectra$feature_id == x[1]]
  # View(ex_spectrum@backend@spectraData@listData)
  # print('From the only one spectra, extract the mz and intensity for MS2')
  Length <- length(ex_spectrum@backend@spectraData@listData[["mz"]]@listData)
  PeakList_matrix <- matrix(c(0, 0), ncol = 2, byrow = FALSE)
  if (Length > 0) {
    for (id in c(1:Length)) {
      PeakList <- matrix(c(ex_spectrum@backend@spectraData@listData[["mz"]]@listData[[id]]
                           , ex_spectrum@backend@spectraData@listData[["intensity"]]@listData[[id]]),
                         ncol = 2, byrow = FALSE)
      PeakList_matrix <- rbind(PeakList_matrix, PeakList)
    }
  }
  DT <- data.table(PeakList_matrix)
  DT[, Index := x[1]]
  # print(PeakList_matrix)
  return(DT)
}
results <- apply(dda_2020_meta[,'feature_id'], 1, F_ExtractMS2)
results <- ldply (results, data.frame) #Transfer results from list to data.frame
RawMS2Hema <- data.table(results) %>% setnames(., c('V1', 'V2'), c('mz', 'intensity'))
RawMS2Hema[,NAME:=Index]
RawMS2Hema <- RawMS2Hema[intensity>0]
wmzVector <<- c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.3, 2, 3, 4, 5,10)
# Original wmz = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1.3,2,3,4,5,10)
# mz used before 0, 0.9,1, 1.3, 2, 3
wintVector <<- c(.01, .05, .1, .2,.33, .3, .4, .5, .53, .6, .7, .8, .9, 1, 1.2, 2, 3, 5,10)
# Original wint = c(.01,.05,.1,.2,.3,.4,.5,.53,.6,.7,.8,.9,1,2,3,5,10)
# intensity used before 0.33, 0.5, 0.53, 0.6, 1, 1.2

F_ScalingFactor <- function(Name, DT,statsValue,RfoldValue) {
  DT[, mzRound := round(mz)]
  setorder(DT, -intensity) # sort intensity from big to small
  MS2 <- unique(DT, by = c('NAME', 'mzRound')) # for the same compound, only keep one scan at the same mzRound
  ## Boxplot visualization and use extreme of the lower whisker as cutoff
  Cutoff <<- boxplot.stats(MS2$intensity)$stats[statsValue] %>% round(.,digits = 1)
  MS2 <- MS2[intensity >= Cutoff]
  ## Plot MSMS distribution
  svg(file = paste(Name,'Cutoff', Cutoff,"Number MS2", nrow(MS2), 'MZ Distribution.svg'))
  print(ggplot(MS2, aes(mz)) + geom_histogram(colour = "black", binwidth = 10) +
          labs(title = Name) + theme_classic()+
          theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),plot.title =element_text(size=17,face="bold"))
  )#+xlim(0,xlimrange)
  dev.off()
  ## Plot Intensity distribution
  svg(file = paste(Name,'Cutoff', Cutoff,"Number MS2", nrow(MS2), 'Intensity Distribution.svg'))
  print(ggplot(MS2, aes(log2(intensity))) +
          geom_histogram(colour = "black", binwidth = 1) +
          labs(title = Name) + theme_classic() + xlim(0,30)+
          theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),plot.title =element_text(size=17,face="bold"))
  )
  dev.off()
  # Creat a matrix for scaling
  setorder(MS2, NAME) # NAME from small to large
  MS2NAME <- unique(MS2, by = 'NAME')  # unique
  alignment = data.table("mzRound" = 1:round(max(MS2$mz))) # creat rows temple
  for (id in seq(1, nrow(MS2NAME), 1)) {
    # for each unique index, select the related mz and intensity
    # id <-1
    x <- as.character(MS2NAME[id, NAME])
    # print(paste("Row",id,"NAME",x))
    MS2Selected <- MS2[NAME == x]
    setorder(MS2Selected, -intensity) # sort from big to small
    Selected <- MS2Selected[, c('mzRound', 'intensity')]
    Selected<-unique(Selected,by="mzRound") # for multiple entries with same mz, only select the one with max intensity
    setnames(Selected, "intensity", x) # change the column into index
    alignment <- merge(alignment, Selected, by = "mzRound", all.x = TRUE)
  }
  alignment <- alignment[, -"mzRound"] # remove the column mz
  alignment[is.na(alignment)] <- 0 #replace missing value with 0
  write.csv(alignment, paste(Name, "Cutoff", Cutoff, "Number MS2", nrow(MS2), "Matrix.csv"))
  Name <<- Name
  # OptimizedWeight<-opt.weight(rdata = as.data.frame(alignment),R = 10,fold = 10,plot = T,verbose = T)
  OptimizedWeight <- opt.weight(rdata = as.data.frame(alignment), R = RfoldValue, fold = RfoldValue, plot = T, verbose = T)
  return(OptimizedWeight)
}
OptHema <- F_ScalingFactor(Name = "K562", DT = RawMS2Hema,statsValue= 1, RfoldValue = 1)

## for testing purpose
# range(MS2CID$Norm)
# max(MS2HCD$intensity)
# MatrixCID<-read.csv('MatrixCID.csv')
# MS2<-MS2[mz<1000]
load("d:/github/Identification/Optimal_weight_factors/subWebNIST.rda")
write.csv(webnist,'webnist.csv')
# webnist.test <- read.csv('d:/github/Identification/webnist.test.csv')
# wt.ran <- opt.weight(rdata = webnist.test,plot = T,R = 10,fold = 10,verbose = T)
# MS2Neg<-F_MSMSDistribution(PNNL=PNNL_Neg,Fdata=Fdata_Neg[1:100],Cutoff=0.1,Name='Negative Density',binwidth=9,xlimrange=1000)
# MatrixNeg<-F_MatrixForScaling(MS2=MS2Neg)
# wt.ran.Neg<-F_opt.weight(MatrixNeg)

#### Fig. 4 Demo analysis based on the hematology ####
#### ImportData using readMSData command from packagee MSnbase
## Load the dataset for the pure AA
# Raw_PureAA_Neg<-F_readMSData(Dir = "d:/VirtualMachineDisk/Hematology/HilicNegative/DDA/mzMLCentroid/AAStandard",
#                              Group= c(rep("AAstandard", 3),rep("BLANK", 1)))
# Raw_PureAA_Pos<-F_readMSData(Dir = "D:/VirtualMachineDisk/Hematology/Centroid/AAStandard",
#                              Group= c(rep("AAstandard", 3),rep("BLANK", 1)))
# Raw_PureAA_Pos22<-F_readMSData(Dir = "D:/VirtualMachineDisk/Hematology/Centroid/AAStandard22",
#                              Group= c(rep("AAstandard", 3),rep("BLANK", 1)))
## Load the dataset for the AA mixture
# Raw_MixAA_Neg<-F_readMSData(Dir = "d:/VirtualMachineDisk/Hematology/HilicNegative/DDA/mzMLCentroid/NormalVSHYP",
#                             Group= c(rep("BLANK", 3),rep("HYP01", 3),rep("HYP02", 3),rep("HYP03", 3),rep("Normal01", 3),rep("Normal02", 3),rep("Normal03", 3)))
Raw_MixAA_Pos<-F_readMSData(Dir = "D:/VirtualMachineDisk/Hematology/Centroid/NormalVSHYP",
              Group= c(rep("HYP01", 3),rep("HYP02", 3),rep("HYP03", 3),rep("Normal01", 3),rep("Normal02", 3),rep("Normal03", 3)))

# Raw_MixAA_Pos22<-F_readMSData(Dir = "D:/VirtualMachineDisk/Hematology/Centroid/NormalVSHYP22",
#                             Group= c(rep("BLANK", 3),rep("HYP", 3),rep("Norm", 3)))
#### Perform peak detection ####
## Calculate the proper constant A
MRP <- 30000 # Velos, 30000 for MS1 and 7500 for MS2 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3673017/
RefMZ <- 400 #  Velos is mz reference 400 while Q-Exactive is 200.
A <- 1 / (MRP * (RefMZ^0.5))
print(A)
B <- A / 2.35482
print(B)
# cwpNeg <- CentWaveParam(ppm=20 ,peakwidth=c(8,30),snthresh = 6,noise=200,prefilter=c(2, 1000),firstBaselineCheck = FALSE,integrate=2)
# cwpNeg<- CentWaveParam(A = 7.077682e-07, ppm = 1, Instrument = 2, integrate=2,
#                         peakwidth=c(8,30),snthresh = 6,noise=200,prefilter=c(2, 1000),firstBaselineCheck = FALSE)
# xdata_MixAA_Neg<-F_findChromPeaks(Raw=Raw_MixAA_Neg,cwp =cwpNeg)
# xdata_PureAA_Neg<-F_findChromPeaks(Raw=Raw_PureAA_Neg,cwp =cwpNeg)

# cwpPos <- CentWaveParam(ppm=20 ,peakwidth=c(2,30),snthresh = 6,noise=200,prefilter=c(2, 1000),firstBaselineCheck = FALSE,integrate=2)
cwpPos<- CentWaveParam(A = 7.077682e-07, ppm = 1, Instrument = 2,
                       peakwidth=c(2,10),snthresh = 5,noise=100,prefilter=c(3, 500),firstBaselineCheck = FALSE)
xdata_MixAA_Pos <- F_findChromPeaks(Raw = Raw_MixAA_Pos, cwp = cwpPos)
xdata_PureAA_Pos <- F_findChromPeaks(Raw = Raw_PureAA_Pos, cwp = cwpPos)

# cwpPos22 <- CentWaveParam(A = 7.077682e-07, ppm = 1, Instrument = 2,
#   peakwidth = c(2, 10), snthresh = 5, noise = 100, prefilter = c(3, 500),firstBaselineCheck = FALSE)
#
# xdata_MixAA_Pos22 <- F_findChromPeaks(Raw = Raw_MixAA_Pos22, cwp = cwpPos22)
# xdata_PureAA_Pos22 <- F_findChromPeaks(Raw = Raw_PureAA_Pos22, cwp = cwpPos22)

#### Retention time alignment ####
register(SerialParam()) #disable paralell
# paramNeg<-ObiwarpParam(subset = 13:21, subsetAdjust = "average") # use normal condition as reference
# xdata_MixAA_Neg <- F_adjustRtime(xdata_MixAA_Neg,param = paramNeg)
# xdata_PureAA_Neg <- F_adjustRtime(xdata_PureAA_Neg,param = ObiwarpParam(subset = 1:3, subsetAdjust = "average"))

paramPos<-ObiwarpParam(subsetAdjust = "average") # use normal condition as reference
xdata_MixAA_Pos <- F_adjustRtime(xdata_MixAA_Pos,param = paramPos)
xdata_PureAA_Pos <- F_adjustRtime(xdata_PureAA_Pos,param = ObiwarpParam(subset = 1:3, subsetAdjust = "average"))

# xdata_MixAA_Pos22$sample_name
# xdata_MixAA_Pos22 <- F_adjustRtime(xdata_MixAA_Pos22,param = ObiwarpParam())
# xdata_PureAA_Pos22$sample_name
# xdata_PureAA_Pos22 <- F_adjustRtime(xdata_PureAA_Pos22,param =ObiwarpParam(subset = 1:3, subsetAdjust = "average"))

## See the retention time alignment effect
F_adjustRtimeEffect<-function(Raw,xdata,Name){
  ## Base peak chromatograms before RT alignment
  bpis <- chromatogram(Raw, aggregationFun = "max")
  ## Base peak chromatograms after RT alignment
  bpis_adj <- chromatogram(xdata, aggregationFun = "max", include = "none")
  ## Plot before and after
  # png(file=paste("RT alignment effect",Name,".png",seq=""),width=700,height=480) #res=300,
  svg(file=paste("RT alignment effect",Name,".svg",seq="")) #res=300,,width=700,height=480
  # par(mfrow = c(2, 1), mar = c(5, 4.2, 1, 2))
  # plot(bpis, col = group_colors[bpis$sample_group])#
  # legend("topright",legend=unique(bpis$sample_group),col=unique(group_colors[bpis$sample_group]), lty=1, cex=0.4,box.lty=0) #
  # title(sub=paste("RT alignment before bpis",Name,seq=""))
  plot(bpis_adj, col = group_colors[bpis_adj$sample_group],
       cex.axis=1.5,cex.lab = 2,cex.main = 2,cex.sub = 1)
  legend("topright",legend=unique(bpis$sample_group),
         col=unique(group_colors[bpis$sample_group]), lty=1, cex=1,box.lty=0) #
  # title(sub=paste("RT alignment after bpis",Name,seq=""))
  dev.off()
}
group_colors <- c("coral1","coral2","coral3","cyan1","cyan2","cyan3")
names(group_colors) <- unique(xdata_MixAA_Pos$sample_group)
F_adjustRtimeEffect(Raw = Raw_MixAA_Pos, xdata = xdata_MixAA_Pos, Name = "K562_2020_Pos")

# group_colors <- c("black","coral1","cyan1")
# names(group_colors) <- unique(xdata_MixAA_Pos22$sample_group)
# F_adjustRtimeEffect(Raw = Raw_MixAA_Pos22, xdata = xdata_MixAA_Pos22, Name = "K562_2022_Pos")
#### Grouping, fill missing ####
# xdata_MixAA_Neg<-F_GroupFill(xdata_MixAA_Neg,
#                              pdp =PeakDensityParam(sampleGroups = xdata_MixAA_Neg$sample_group,bw=5,binSize=0.04))
xdata_MixAA_Pos<-F_GroupFill(xdata_MixAA_Pos,
                             pdp =PeakDensityParam(sampleGroups = xdata_MixAA_Pos$sample_group, bw=10,binSize=0.02, minFraction = 0.6, minSamples = 2))
# xdata_MixAA_Pos22<-F_GroupFill(xdata_MixAA_Pos22,
#                              pdp =PeakDensityParam(sampleGroups = xdata_MixAA_Pos22$sample_group, bw=10,binSize=0.02, minFraction = 0.6, minSamples = 2))
#### Boxplot ####
F_boxplotxdata<-function(xdata,Cell,group_colors,Keywords){ ## Boxplot after peak detection per-sample peak intensities (in log2 scale)
  # names(group_colors) <- unique(xdata$sample_group)
  ints <- split(log2(chromPeaks(xdata)[, "into"]),
                f = chromPeaks(xdata)[, "sample"])
  # ints[ints>0]
  # png(file=paste("Boxplot xdata peak intensity.png",seq=""),res=100,width=700,height=480) #
  svg(file=paste(Cell,"Boxplot xdata peak intensity.svg",seq="")) #
  boxplot(ints,ylab = "log2 intensity", main = Cell, varwidth = TRUE,
          col = group_colors[xdata$sample_group],names=xdata$sample_name,xaxt = "n",
          cex.axis=1.5,cex.lab = 2,cex.main = 2,cex.sub = 1)
  # ## Draw the labels on x-axis
  # text(x = 1:length(ints),
  #      ## Move labels to just below bottom of chart.
  #      y = par("usr")[3] - 0.45,
  #      ## Use names from the data list.
  #      labels = sub(xdata$sample_name, pattern = Keywords, replacement = "", perl = TRUE),
  #      ## Change the clipping region.
  #      xpd = NA,
  #      ## Rotate the labels by 35 degrees.
  #      srt = 35,
  #      ## Adjust the labels to almost 100% right-justified.
  #      adj = 0.965,
  #      ## Increase label size.
  #      cex = 1.2)
  legend("topright",legend=unique(xdata$sample_group),col=unique(group_colors[xdata$sample_group]), lty=1, cex=1,box.lty=0) #
  dev.off()
}

# F_boxplotxdata(xdata_NoB=filterFile(xdata_MixAA_Neg, file = c(4:21)),
#                Name=NewName[5:22],
#                group_colors = c("coral1","coral2","coral3","cyan1","cyan2","cyan3"))
# xdata_MixAA_Pos_NoB <- filterFile(xdata_MixAA_Pos, file = c(4:21))
group_colors <- c("coral1","coral2","coral3","cyan1","cyan2","cyan3")
names(group_colors) <- unique(xdata_MixAA_Pos$sample_group)
F_boxplotxdata(xdata=xdata_MixAA_Pos, Cell = 'Peak intensity', Keywords = "ACN.*",group_colors = group_colors)

# xdata_MixAA_Pos22_NoB <- filterFile(xdata_MixAA_Pos22, file = c(4:9))
# group_colors <- c("black","coral1","cyan1")
# names(group_colors) <- unique(xdata_MixAA_Pos22$sample_group)
# F_boxplotxdata(xdata=xdata_MixAA_Pos22_NoB, Cell = '2022', Keywords = "ACN.*",group_colors = group_colors)
# F_boxplotmSet(mSetNeg,group_colors = c("coral1","cyan1"))
# F_boxplotmSet(mSetPos,group_colors = c("coral1","cyan1"))

#### Add mz, rt and isotope information ####
F_AddMzRtIso <- function(xdata, Cell, Mode) {
  ## Export quantification using exportMetaboAnalyst function
  exportMetaboAnalyst(xdata,
                      file = paste0(Cell, "exportMetab.csv", seq = ""), label = xdata$sample_group,
                      value = "into", method = "maxint", msLevel = 1,filled=TRUE
  )
  exportMetab <- data.table(read.csv(paste0(Cell, "exportMetab.csv", seq = "")))
  print(names(exportMetab))
  ## Add mz, rt and isotope information
  if (length(unique(msLevel(xdata))) != 1) {
    xdata <- filterMsLevel(xdata, msLevel = 1)
  } # Use MS1 as input of Camera
  xset <- as(xdata, "xcmsSet")
  xsa <- xsAnnotate(xset, polarity = Mode) # create an CAMERA object
  xsag <- groupFWHM(xsa, intval = "into") # group based on RT
  xsaC <- groupCorr(xsag, intval = "into") # Verify grouping
  xsaI <- findIsotopes(xsaC, intval = "into") # Annotate isotopes
  peaklist <- data.table(getPeaklist(xsaI, intval = "into")) # use the peak area to represent the intensity
  IsoFilter <- str_replace_all(peaklist$isotopes, pattern = ".*M\\].*", replacement = "")
  DT_peaklist <- cbind(peaklist[, c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "isotopes", "pcgroup")],
                       "IsoFilter" = IsoFilter, "mode" = Mode,
                       "Sample" = xdata@msFeatureData[["featureDefinitions"]]@rownames)
  ## Filter based on the isotope annotation
  exportMetab <- left_join(exportMetab, DT_peaklist, by = "Sample") %>% data.table(.)
  Exp <- exportMetab[IsoFilter == "" | is.na(IsoFilter)]
  ## Change names of the exported table for mSet
  # head(Exp)
  # Exp[Sample=='Label',]
  Sample<-sub(names(Exp), pattern = ".mzML",replacement = "",perl = TRUE)
  setnames(Exp,names(Exp),Sample) # change the name of the data table
  print(names(Exp))
  ## Add two row as Phenotype and time
  Phenotype<-sub("Label", "Label", Exp[Sample=='Label',]) %>% sub("[0-9]+", "", .) %>%
    sub("Normal", "Norm", .)# update the Phenotype
  Phe <-rbind(as.list(Phenotype),Exp[Sample!='Label',])
  ## Export all the results except blank and NA
  inxLabel <- Phe[Sample=='Label',] %in% c('Label')
  inxNorm <- Phe[Sample=='Label',] %in% c('Norm')
  inxHYP <- Phe[Sample=='Label',] %in% c('HYP')
  PheExport <- cbind(Phe[,..inxLabel],Phe[,..inxHYP],Phe[,..inxNorm])
  write.csv(PheExport,paste0(Cell,"Phe.csv"),row.names = FALSE)
  return(Exp)
}
setwd(Path)
Exp_2020 <- F_AddMzRtIso(xdata = xdata_MixAA_Pos, Cell = "2020", Mode = "positive")
# Exp_2022 <- F_AddMzRtIso(xdata = xdata_MixAA_Pos22, Cell = "2022", Mode = "positive")
#### InitDataObjects ####
.get.mSet <- function(mSetObj=NA){
    return(mSetObj)
}
.set.mSet <- function(mSetObj=NA){
  return(mSetObj);
}
Plot3D <- function(x, y = NULL, z = NULL, color = par("col"), pch = NULL,
                   main = NULL, sub = NULL, xlim = NULL, ylim = NULL, zlim = NULL,
                   xlab = NULL, ylab = NULL, zlab = NULL, scale.y = 1, angle = 40,
                   axis = TRUE, tick.marks = TRUE, label.tick.marks = TRUE,
                   x.ticklabs = NULL, y.ticklabs = NULL, z.ticklabs = NULL,
                   y.margin.add = 0, grid = TRUE, box = TRUE, lab = par("lab"),
                   lab.z = mean(lab[1:2]), type = "p", highlight.3d = FALSE,
                   mar = c(5, 3, 4, 3) + 0.1, col.axis = par("col.axis"),
                   col.grid = "grey", col.lab = par("col.lab"), cex.symbols = par("cex"),
                   cex.axis = 0.8 * par("cex.axis"), cex.lab = par("cex.lab"),
                   font.axis = par("font.axis"), font.lab = par("font.lab"),
                   lty.axis = par("lty"), lty.grid = 2, lty.hide = 1,
                   lty.hplot = par("lty"), log = "", ...)
# log not yet implemented
{
  ## Uwe Ligges <ligges@statistik.tu-dortmund.de>,
  ## http://www.statistik.tu-dortmund.de/~ligges
  ##
  ## For MANY ideas and improvements thanks to Martin Maechler!!!
  ## Parts of the help files are stolen from the standard plotting functions in R.

  mem.par <- par(mar = mar)
  x.scal <- y.scal <- z.scal <- 1
  xlabel <- if (!missing(x)) deparse(substitute(x))
  ylabel <- if (!missing(y)) deparse(substitute(y))
  zlabel <- if (!missing(z)) deparse(substitute(z))

  ## color as part of `x' (data.frame or list):
  if(!is.null(d <- dim(x)) && (length(d) == 2) && (d[2] >= 4))
    color <- x[,4]
  else if(is.list(x) && !is.null(x$color))
    color <- x$color

  ## convert 'anything' -> vector
  xyz <- xyz.coords(x=x, y=y, z=z, xlab=xlabel, ylab=ylabel, zlab=zlabel,
                    log=log)
  if(is.null(xlab)) { xlab <- xyz$xlab; if(is.null(xlab)) xlab <- "" }
  if(is.null(ylab)) { ylab <- xyz$ylab; if(is.null(ylab)) ylab <- "" }
  if(is.null(zlab)) { zlab <- xyz$zlab; if(is.null(zlab)) zlab <- "" }

  if(length(color) == 1)
    color <- rep(color, length(xyz$x))
  else if(length(color) != length(xyz$x))
    stop("length(color) ", "must be equal length(x) or 1")

  angle <- (angle %% 360) / 90
  yz.f <- scale.y * abs(if(angle < 1) angle else if(angle > 3) angle - 4 else 2 - angle)
  yx.f <- scale.y * (if(angle < 2) 1 - angle else angle - 3)
  if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
    temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
    temp <- xlab;  xlab <- ylab;   ylab <- temp
    temp <- xlim;  xlim <- ylim;   ylim <- temp
  }
  angle.1 <- (1 < angle && angle < 2) || angle > 3
  angle.2 <- 1 <= angle && angle <= 3
  dat <- cbind(as.data.frame(xyz[c("x","y","z")]), col = color)

  n <- nrow(dat);
  y.range <- range(dat$y[is.finite(dat$y)])

  ### 3D-highlighting / colors / sort by y
  if(type == "p" || type == "h") {
    y.ord <- rev(order(dat$y))
    dat <- dat[y.ord, ]
    if(length(pch) > 1)
      if(length(pch) != length(y.ord))
        stop("length(pch) ", "must be equal length(x) or 1")
    else pch <- pch[y.ord]
    daty <- dat$y
    daty[!is.finite(daty)] <- mean(daty[is.finite(daty)])
    if(highlight.3d && !(all(diff(daty) == 0)))
      dat$col <- rgb(seq(0, 1, length = n) * (y.range[2] - daty) / diff(y.range), g=0, b=0)
  }

  ### optim. axis scaling
  p.lab <- par("lab")
  ## Y
  y.range <- range(dat$y[is.finite(dat$y)], ylim)
  y.prty <- pretty(y.range, n = lab[2],
                   min.n = max(1, min(.5 * lab[2], p.lab[2])))
  y.scal <- round(diff(y.prty[1:2]), digits = 12)
  y.add <- min(y.prty)
  dat$y <- (dat$y - y.add) / y.scal
  y.max <- (max(y.prty) - y.add) / y.scal

  x.range <- range(dat$x[is.finite(dat$x)], xlim)
  x.prty <- pretty(x.range, n = lab[1],
                   min.n = max(1, min(.5 * lab[1], p.lab[1])))
  x.scal <- round(diff(x.prty[1:2]), digits = 12)
  dat$x <- dat$x / x.scal
  x.range <- range(x.prty) / x.scal
  x.max <- ceiling(x.range[2])
  x.min <-   floor(x.range[1])
  if(!is.null(xlim)) {
    x.max <- max(x.max, ceiling(xlim[2] / x.scal))
    x.min <- min(x.min,   floor(xlim[1] / x.scal))
  }
  x.range <- range(x.min, x.max)
  ## Z
  z.range <- range(dat$z[is.finite(dat$z)], zlim)
  z.prty <- pretty(z.range, n = lab.z,
                   min.n = max(1, min(.5 * lab.z, p.lab[2])))
  z.scal <- round(diff(z.prty[1:2]), digits = 12)
  dat$z <- dat$z / z.scal
  z.range <- range(z.prty) / z.scal
  z.max <- ceiling(z.range[2])
  z.min <-   floor(z.range[1])
  if(!is.null(zlim)) {
    z.max <- max(z.max, ceiling(zlim[2] / z.scal))
    z.min <- min(z.min,   floor(zlim[1] / z.scal))
  }
  z.range <- range(z.min, z.max)

  ### init graphics
  plot.new()
  if(angle.2) {x1 <- x.min + yx.f * y.max; x2 <- x.max}
  else        {x1 <- x.min; x2 <- x.max + yx.f * y.max}
  plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
  temp <- strwidth(format(rev(y.prty))[1], cex = cex.axis/par("cex"))
  if(angle.2) x1 <- x1 - temp - y.margin.add
  else        x2 <- x2 + temp + y.margin.add
  plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
  if(angle > 2) par("usr" = par("usr")[c(2, 1, 3:4)])
  usr <- par("usr") # we have to remind it for use in closures
  title(main, sub, ...)

  ### draw axis, tick marks, labels, grid, ...
  xx <- if(angle.2) c(x.min, x.max) else c(x.max, x.min)
  if(grid) {
    ## grids
    ###################
    # XY wall
    i <- x.min:x.max;
    segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + z.min,
             col = col.grid, lty = lty.grid);

    i <- 0:y.max;
    segments(x.min + (i * yx.f), i * yz.f + z.min,
             x.max + (i * yx.f), i * yz.f + z.min,
             col = col.grid, lty = lty.grid);

    ######################
    # XZ wall
    # verticle lines
    temp <- yx.f * y.max;
    temp1 <- yz.f * y.max;
    i <- (x.min + temp):(x.max + temp);
    segments(i, z.min + temp1, i, z.max + temp1,
             col = col.grid, lty = lty.grid);

    # horizontal lines
    i <- (z.min + temp1):(z.max + temp1);
    segments(x.min + temp, i, x.max + temp, i,
             col = col.grid, lty = lty.grid)


    ##################
    # YZ wall
    # horizontal lines
    i <- xx[2]:x.min;
    mm <- z.min:z.max;
    segments(i, mm, i + temp, mm + temp1,
             col = col.grid, lty = lty.grid);
    # verticle lines
    i <- 0:y.max;
    segments(x.min + (i * yx.f), i * yz.f + z.min,
             xx[2] + (i * yx.f), i * yz.f + z.max,
             col = col.grid, lty = lty.grid)


    # make the axis into solid line
    segments(x.min, z.min, x.min + (yx.f * y.max), yz.f * y.max + z.min,
             col = col.grid, lty = lty.hide);
    segments(x.max, z.min, x.max + (yx.f * y.max), yz.f * y.max + z.min,
             col = col.axis, lty = lty.hide);
    segments(x.min + (y.max * yx.f), y.max * yz.f + z.min,
             x.max + (y.max* yx.f), y.max * yz.f + z.min,
             col = col.grid, lty = lty.hide);
    segments(x.min + temp, z.min + temp1, x.min + temp, z.max + temp1,
             col = col.grid, lty = lty.hide);
    segments(x.max + temp, z.min + temp1, x.max + temp, z.max + temp1,
             col = col.axis, lty = lty.hide);
    segments(x.min + temp, z.max + temp1, x.max + temp, z.max + temp1,
             col = col.axis, lty = lty.hide);
    segments(xx[2], z.max, xx[2] + temp, z.max + temp1,
             col = col.axis, lty = lty.hide);
  }
  if(axis) {
    if(tick.marks) { ## tick marks
      xtl <- (z.max - z.min) * (tcl <- -par("tcl")) / 50
      ztl <- (x.max - x.min) * tcl / 50
      mysegs <- function(x0,y0, x1,y1)
        segments(x0,y0, x1,y1, col=col.axis, lty=lty.axis)
      ## Y
      i.y <- 0:y.max
      mysegs(yx.f * i.y - ztl + xx[1], yz.f * i.y + z.min,
             yx.f * i.y + ztl + xx[1], yz.f * i.y + z.min)
      ## X
      i.x <- x.min:x.max
      mysegs(i.x, -xtl + z.min, i.x, xtl + z.min)
      ## Z
      i.z <- z.min:z.max
      mysegs(-ztl + xx[2], i.z, ztl + xx[2], i.z)

      if(label.tick.marks) { ## label tick marks
        las <- par("las")
        mytext <- function(labels, side, at, ...)
          mtext(text = labels, side = side, at = at, line = -.5,
                col=col.lab, cex=cex.axis, font=font.lab, ...)
        ## X
        if(is.null(x.ticklabs))
          x.ticklabs <- format(i.x * x.scal)
        mytext(x.ticklabs, side = 1, at = i.x)
        ## Z
        if(is.null(z.ticklabs))
          z.ticklabs <- format(i.z * z.scal)
        mytext(z.ticklabs, side = if(angle.1) 4 else 2, at = i.z,
               adj = if(0 < las && las < 3) 1 else NA)
        ## Y
        temp <- if(angle > 2) rev(i.y) else i.y ## turn y-labels around
        if(is.null(y.ticklabs))
          y.ticklabs <- format(y.prty)
        else if (angle > 2)
          y.ticklabs <- rev(y.ticklabs)
        text(i.y * yx.f + xx[1],
             i.y * yz.f + z.min, y.ticklabs,
             pos=if(angle.1) 2 else 4, offset=1,
             col=col.lab, cex=cex.axis/par("cex"), font=font.lab)
      }
    }

    ## axis and labels

    mytext2 <- function(lab, side, line, at)
      mtext(lab, side = side, line = line, at = at, col = col.lab,
            cex = cex.lab, font = font.axis, las = 0)
    ## X
    lines(c(x.min, x.max), c(z.min, z.min), col = col.axis, lty = lty.axis)
    mytext2(xlab, 1, line = 1.5, at = mean(x.range))
    ## Y
    lines(xx[1] + c(0, y.max * yx.f), c(z.min, y.max * yz.f + z.min),
          col = col.axis, lty = lty.axis)
    mytext2(ylab, if(angle.1) 2 else 4, line= 0.5, at = z.min + y.max * yz.f)

    ## Z
    lines(xx[c(2,2)], c(z.min, z.max), col = col.axis, lty = lty.axis)
    mytext2(zlab, if(angle.1) 4 else 2, line= 1.5, at = mean(z.range))

  }

  ### plot points
  x <- dat$x + (dat$y * yx.f)
  z <- dat$z + (dat$y * yz.f)
  col <- as.character(dat$col)
  if(type == "h") {
    z2 <- dat$y * yz.f + z.min
    segments(x, z, x, z2, col = col, cex = cex.symbols, lty = lty.hplot, ...)
    points(x, z, type = "p", col = col, pch = pch, cex = cex.symbols, ...)
  }
  else points(x, z, type = type, col = col, pch = pch, cex = cex.symbols, ...)

  ### box-lines in front of points (overlay)
  if(axis && box) {
    lines(c(x.min, x.max), c(z.max, z.max),
          col = col.axis, lty = lty.axis)
    lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max,
          col = col.axis, lty = lty.axis)
    lines(xx[c(1,1)], c(z.min, z.max), col = col.axis, lty = lty.axis)
  }


  # par(mem.par) # we MUST NOT set the margins back
  ### Return Function Object
  ob <- ls() ## remove all unused objects from the result's enviroment:
  rm(list = ob[!ob %in% c("angle", "mar", "usr", "x.scal", "y.scal", "z.scal", "yx.f",
                          "yz.f", "y.add", "z.min", "z.max", "x.min", "x.max", "y.max",
                          "x.prty", "y.prty", "z.prty")])
  rm(ob)
  invisible(list(
    xyz.convert = function(x, y=NULL, z=NULL) {
      xyz <- xyz.coords(x, y, z)
      if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
        temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
      }
      y <- (xyz$y - y.add) / y.scal
      return(list(x = xyz$x / x.scal + yx.f * y,
                  y = xyz$z / z.scal + yz.f * y))
    },
    points3d = function(x, y = NULL, z = NULL, type = "p", ...) {
      xyz <- xyz.coords(x, y, z)
      if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
        temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
      }
      y2 <- (xyz$y - y.add) / y.scal
      x <- xyz$x / x.scal + yx.f * y2
      y <- xyz$z / z.scal + yz.f * y2
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      if(type == "h") {
        y2 <- z.min + yz.f * y2
        segments(x, y, x, y2, ...)
        points(x, y, type = "p", ...)
      }
      else points(x, y, type = type, ...)
    },
    plane3d = function(Intercept, x.coef = NULL, y.coef = NULL,
                       lty = "dashed", lty.box = NULL, ...){
      if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
      if(is.null(lty.box)) lty.box <- lty
      if(is.null(x.coef) && length(Intercept) == 3){
        x.coef <- Intercept[if(angle > 2) 3 else 2]
        y.coef <- Intercept[if(angle > 2) 2 else 3]
        Intercept <- Intercept[1]
      }
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      x <- x.min:x.max
      ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
      x.coef <- x.coef * x.scal
      z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
      z2 <- (Intercept + x * x.coef +
               (y.max * y.scal + y.add) * y.coef) / z.scal
      segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
      y <- 0:y.max
      ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
      y.coef <- (y * y.scal + y.add) * y.coef
      z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
      z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
      segments(x.min + y * yx.f, z1 + y * yz.f,
               x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
    },

    wall3d = function(Intercept, x.coef = NULL, y.coef = NULL,
                      lty = "dashed", lty.box = NULL, ...){
      if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
      if(is.null(lty.box)) lty.box <- lty
      if(is.null(x.coef) && length(Intercept) == 3){
        x.coef <- Intercept[if(angle > 2) 3 else 2]
        y.coef <- Intercept[if(angle > 2) 2 else 3]
        Intercept <- Intercept[1]
      }
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      x <- x.min:x.max
      ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
      x.coef <- x.coef * x.scal
      z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
      z2 <- (Intercept + x * x.coef +
               (y.max * y.scal + y.add) * y.coef) / z.scal
      segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
      y <- 0:y.max
      ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
      y.coef <- (y * y.scal + y.add) * y.coef
      z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
      z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
      segments(x.min + y * yx.f, z1 + y * yz.f,
               x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
    },
    box3d = function(...){
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      lines(c(x.min, x.max), c(z.max, z.max), ...)
      lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max, ...)
      lines(c(0, y.max * yx.f) + x.min, c(0, y.max * yz.f) + z.max, ...)
      lines(c(x.max, x.max), c(z.min, z.max), ...)
      lines(c(x.min, x.min), c(z.min, z.max), ...)
      lines(c(x.min, x.max), c(z.min, z.min), ...)
    }
  ))
}
PlotPLS3DScoreImg<-function(mSetObj=NA, imgName, format="svg", dpi=72, width=NA, inx1, inx2, inx3, angl){

  mSetObj <- .get.mSet(mSetObj);

  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  if(is.na(width)){
    w <- 9;
  }else if(width == 0){
    w <- 7.2;
  }else{
    w <- width;
  }
  h <- w;

  mSetObj$imgSet$pls.score3d <- imgName;

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  par(mar=c(5,5,3,3));

  xlabel <- paste("Component", inx1, "(", round(100*mSetObj$analSet$plsr$Xvar[inx1]/mSetObj$analSet$plsr$Xtotvar,1), "%)");
  ylabel <- paste("Component", inx2, "(", round(100*mSetObj$analSet$plsr$Xvar[inx2]/mSetObj$analSet$plsr$Xtotvar,1), "%)");
  zlabel <- paste("Component", inx3, "(", round(100*mSetObj$analSet$plsr$Xvar[inx3]/mSetObj$analSet$plsr$Xtotvar,1), "%)");

  # cols <- GetColorSchema(mSetObj);
  cols <<- c(rep('coral1',9),rep('cyan1',9));
  legend.nm <- unique(as.character(mSetObj$dataSet$cls));
  uniq.cols <- unique(cols);
  pchs <- as.numeric(mSetObj$dataSet$cls)+1;
  uniq.pchs <- unique(pchs);
  Plot3D(mSetObj$analSet$plsr$score[,inx1], mSetObj$analSet$plsr$score[,inx2], mSetObj$analSet$plsr$score[,inx3], xlab= xlabel, ylab=ylabel,
         zlab=zlabel, angle =angl, color=cols, pch=pchs, box=F,cex.lab = 2,cex.axis = 1.5);
  legend("topleft", legend = legend.nm, pch=uniq.pchs, col=uniq.cols);
  dev.off();


    # # 3D View using plotly
    # if(length(uniq.pchs) > 3){
    #   col <- RColorBrewer::brewer.pal(length(uniq.pchs), "Set3")
    # }else{
    #   col <- c("#1972A4", "#FF7070")
    # }
    # p <- plot_ly(x = mSetObj$analSet$plsr$score[, inx1], y = mSetObj$analSet$plsr$score[, inx2], z = mSetObj$analSet$plsr$score[, inx3],
    #              color = mSetObj$dataSet$cls, colors = col)
    # p <- add_markers(p, sizes = 5)
    # p <- layout(p, scene = list(xaxis = list(title = xlabel),
    #                             yaxis = list(title = ylabel),
    #                             zaxis = list(title = zlabel)))
    #
    # mSetObj$imgSet$plsda.3d <- p;
    # print("The Interactive 3D PLS-DA plot has been created, please find it in mSet$imgSet$plsda.3d.")
    #
    #
  return(.set.mSet(mSetObj));
}
F_InitDataObjects<-function(FileName, Exp){
  rm(mSet)
  mSet <- InitDataObjects("pktable", "stat", FALSE) # is a peak table ("pktable") and that statistical analysis will be performed ("stat"),data is not paired
  ## read in the processed data (created above)
  mSet<-Read.TextData(mSet, FileName, "colts", "disc")
  mSet$msgSet$read.msg ## To view messages from the data import and processing
  ncol(mSet[["dataSet"]][["orig"]])
  ## Filter out missing value
  mSet <- SanityCheckData(mSet) #check function evaluates the accuracy of sample and class labels, data structure, deals with non-numeric values, removes columns that are constant across all samples (variance = 0), and by default replaces missing values with half of the original minimal positive value in your dataset.
  mSet <- RemoveMissingPercent(mSet, percent=0.2) #Remove features containing a user-defined % cut-off of missing values
  # ncol(mSet[["dataSet"]][["preproc"]])
  ## Fill in missing value
  mSet <- ImputeVar(mSet, method="knn") #replace missing value
  print(mSet$msgSet$replace.msg)
  ## Prepare data for normalization
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, rowNorm="MedianNorm", transNorm = "NULL", scaleNorm = "AutoNorm") #  ,"LogNorm","AutoNorm"
  Name <- sub(FileName, pattern = ".csv",replacement = "",perl = TRUE)
  mSet <- PlotNormSummary(mSet, paste0(Name,"_norm_"), "png", 300, width=NA) #Feature View of before and after normalization
  mSet <- PlotSampleNormSummary(mSet, paste0(Name,"_snorm_"), "png", 300, width=NA) #Sample View of before and after normalization
  ## PCA
  mSet<-PCA.Anal(mSet)
  # mSet<-PlotPCAPairSummary(mSet, paste0(Name,"_pca_pair_"), "png", 300, width=NA, 5)
  # mSet<-PlotPCAScree(mSet, paste0(Name,"_pca_scree_"), "png", 300, width=NA, 5)
  mSet<-PlotPCA2DScore(mSet, paste0(Name,"_pca_score2d_"), "png", 300, width=NA, pcx=1,pcy=2)
  # mSet<-PlotPCALoading(mSet, paste0(Name,"_pca_loading_"), format="png", dpi=300, width=NA, 1,2)
  # mSet<-PlotPCABiplot(mSet, paste0(Name,"_pca_biplot_"), "png", 300, width=NA, 1,2)
  mSet<-PlotPCA3DScoreImg(mSet, paste0(Name,"_pca_score3d_"), "png", 300, width=NA,inx1 = 1,inx2=2,inx3=3,angl=40)
  #mSet$imgSet$pca.3d
  ## Create an Interactive PCA Plot
  # mSet<-iPCA.Anal(mSet, paste0(Name,"_ipca_3d.json"))
  # mSet$imgSet$time$score3d  # View the interactive scores plot
  ## PLS-DA
  mSet<-PLSR.Anal(mSet, reg=TRUE) # Perform PLS-DA using Orthogonal scores algorithm
  # mSet<-PlotPLSPairSummary(mSet, paste0(Name,"_pls_pair_"), "png", 300, width=NA, 5) #5 means the number of principal components
  mSet<-PlotPLS2DScore(mSet, paste0(Name,"_pls_score2d_"), "png", 300, width=NA, inx1=1,inx2=2) # component 1 VS component 2, show labels
  mSet<-PlotPLS3DScoreImg(mSet, paste0(Name,"pls_score3d_"), "png", 300, width=NA, inx1=1,inx2=2,inx3=3, angl=40) # 3 dimension Component 1, 2 and 3
  ## Volcano plot
  # mSet<-Volcano.Anal(mSet, threshp=0.05,fcthresh=2,cmpType=1,equal.var=FALSE,pval.type='raw') # used for plotting
  # mSet<-PlotVolcano(mSet, imgName=paste0(Name,"Volcano"),plotLbl=0,"png", dpi=500, width=NA)
  # mSet<-Volcano.Anal(mSet, threshp=1,fcthresh=1,pval.type="fdr",cmpType=1,paired=F,equal.var=F,nonpar=F) # used for mummichog
  # mSet<-FC.Anal.unpaired(mSet, 2.0, 0)
  ## Calculate the raw pvalue
  mSet<-Volcano.Anal(mSet, threshp=1,fcthresh=1,cmpType=0,pval.type='raw') # left/right
  Volcano_raw <- data.table(read.csv("volcano.csv")) # load in the dataset generated in the last step
  Volcano_raw <- setnames(Volcano_raw,"X","Sample") %>% .[,-'X.log10.p.']
  ## Calculate the adjusted pvalue
  mSet<-Volcano.Anal(mSet, threshp=1,fcthresh=1,cmpType=0,pval.type='fdr')
  Volcano_fdr<-data.table(read.csv("volcano.csv")) # load in the dataset generated in the last step
  Volcano_fdr <- setnames(Volcano_fdr,"X","Sample") %>% .[,-c('FC','log2.FC.','X.log10.p.')] #
  ## Combine the raw and adjusted pvalue
  Volcano <- left_join(Volcano_raw, Volcano_fdr, by='Sample')
  Volcano$raw.pval<-as.numeric(Volcano$raw.pval)
  Volcano$p.ajusted<-as.numeric(Volcano$p.ajusted)
  Volcano$log2.FC.<-as.numeric(Volcano$log2.FC.)
  Volcano <- left_join(Volcano,Exp,by='Sample')
  ## Output all the results for follow-up pathway enrichment
  write.csv(Volcano,paste0(Name,"_Volcano.csv"),row.names = FALSE)
  return(Volcano)
}
setwd('./output')
Volcano <- F_InitDataObjects(FileName = "2020Phe.csv", Exp = Exp_2020 )
Volcano_Check <- F_CheckMH20AA(DT=Volcano,MH20AA)
View(Volcano_Check)
#### MS peaks to Pathways mummichog ####
F_Mummichog <- function(FT, Cell, Exp){
  ## reform the volcano results
  # write.csv(PheTimeUniExp,paste0(Cell,"PheTimeUniExp.csv"),row.names = FALSE)
  MumInputUp <- FT[log2.FC.>0, c("mz", "mode", "raw.pval","log2.FC.", "rt")] %>%
    setnames(.,c("mz","raw.pval","log2.FC.","rt"),c("m.z","p.value","t.score","r.t"))
  MumInputUp <- unique(MumInputUp,by=c("m.z","r.t"))
  write.csv(MumInputUp,paste0(Cell,"MumInputUp.txt"),row.names = FALSE)
  MumInputDown <- FT[log2.FC.<0, c("mz", "mode", "raw.pval","log2.FC.", "rt")] %>%
    setnames(.,c("mz","raw.pval","log2.FC.","rt"),c("m.z","p.value","t.score","r.t"))
  MumInputDown <- unique(MumInputDown,by=c("m.z","r.t"))
  write.csv(MumInputDown,paste0(Cell,"MumInputDown.txt"),row.names = FALSE)
  F_mMum <- function(Input, Label){
    ## Initialize mMum
    rm(mMum)
    mMum <- InitDataObjects("mass_all", "mummichog", FALSE) # Create objects for storing processed data from the MS peaks to pathways module
    SetPeakFormat("mprt") # Set the format of the peak list, contains m.z, p.value, rt, and t.score
    mMum <- UpdateInstrumentParameters(mMum,10, "positive",force_primary_ion="yes",rt_frac = 0.02) # Set parameters for analysis, in this case the mass accuracy is set to 5 ppm, the mode of the MS instrument is negative
    mMum<-Read.PeakListData(mMum, paste0(Cell, Input)) #only support .txt format
    mMum <- SanityCheckMummichogData(mMum) # Sanity check of the uploaded data
    ## Now customize adducts
    #add.vec <- c("M [1+]","M+H [1+]","M+2H [2+]","M+3H [3+]","M+Na [1+]","M+H+Na [2+]","M+K [1+]","M+H2O+H [1+]","M-H2O+H [1+]","M-H4O2+H [1+]","M(C13)+H [1+]","M(C13)+2H [2+]","M(C13)+3H [3+]","M(Cl37)+H [1+]","M-NH3+H [1+]","M-CO+H [1+]","M-CO2+H [1+]","M-HCOOH+H [1+]","M+HCOONa [1+]","M-HCOONa+H [1+]","M+NaCl [1+]","M-C3H4O2+H [1+]","M+HCOOK [1+]","M-HCOOK+H [1+]","M-H [1-]","M-2H [2-]","M-H2O-H [1-]","M-H+O [1-]","M+K-2H [1-]","M+Na-2H [1- ]","M+Cl [1-]","M+Cl37 [1-]","M+HCOO [1-]","M+CH3COO [1-]")
    add.vec <- c("M+H [1+]","M+Na [1+]","M-H2O+H [1+]")
    mMum<-Setup.AdductData(mMum, add.vec);
    mMum<-PerformAdductMapping(mMum, "mixed")
    # Perform the mummichog algorithm, First set the algorithm to be used to mummichog, then. This function may take sometime for processing, and will output the pathway-results and the compound matching tables in your working directory
    mMum<-SetPeakEnrichMethod(mMum, "mum", "v2") #in this case the pathway library is from the human MFN model.
    mMum<-SetMummichogPval(mMum, 0.01)  #set the p-value cutoff for mummichog analysis, also use 0.15
    mMum<-PerformPSEA(mMum, "hsa_mfn", permNum = 100,libVersion = "current") #hsa_mfn library,permutations to 1000 bsu_kegg for bacteria
    mMum<-PlotPeaks2Paths(mMum, paste0(Cell,Label,"peaks_to_paths"), "png", 300, width=NA) #The color and size of each circle corresponds to its p-value and enrichment factor,respectively. The enrichment factor of a pathway is calculated as the ratio between the number of significant pathway hits and the expected number of compound hits within the pathway.
    head(mMum$mummi.resmat) # To view the results of the pathway analysis
    ## Network view of meta-analysis
    mMum<-PlotMSPeaksCpdEcpdNetwork(mMum, "cpd", 0.25, "static", "fr", "YlOrRd", 3.5)
    mMum<-PlotMSPeaksCpdEcpdNetwork(mMum, "ec", 0.25, "static", "fr", "YlOrRd", 3.5)
    file.rename(c("mummichog_cpd_network.png","mummichog_ec_network.png"),
                c(paste0(Cell,Label,"_mummichog_cpd_network.png"),paste0(Cell,Label,"_mummichog_ec_network.png")))
    ## Load the pathway and compound files into the final result excel
    mummichog_pathway_enrichment <- read.csv("mummichog_pathway_enrichment.csv")
    write.csv(mummichog_pathway_enrichment,paste0(Cell,Label,"_mummichog_pathway_enrichment.csv"),row.names = FALSE)
    ## Process then load the mummichog_matched_compound_all.csv file
    matched_compound<-data.table(read.csv("mummichog_matched_compound_all.csv")) %>% setorder(.,Query.Mass,Retention.Time,Mass.Diff)
    setnames(matched_compound,c("Query.Mass","Matched.Compound","Matched.Form","Retention.Time"),c("mz","kegg_id","adduct","rt"))
    compound_db<-data.table(readRDS(paste0(Path,"compound_db.rds")))  %>% unique(.,by="kegg_id") # only keep 1 entry for one kegg id
    matched_compound<-left_join(matched_compound,compound_db,by=("kegg_id")) %>% data.table(.)
    ## Calculate mass according to adduct type
    matched_compound[adduct=="M+H [1+]",mass:=mz-1.007276]
    matched_compound[adduct=="M-H2O+H [1+]",mass:=mz+17.00381]
    matched_compound[adduct=="M+Na [1+]",mass:=mz-22.989218]
    matched_compound[,massRound:=round(mass,digits = 2)]
    setcolorder(matched_compound, c('mz','rt',colnames(matched_compound)[!(colnames(matched_compound) %in% c('mz','rt'))]))
    ## Add quantification value to the matched compound
    # names(matched_compound)
    # names(Exp)
    setcolorder(FT, c('mz','mzmin','mzmax','rt','rtmin','rtmax',colnames(FT)[!(colnames(FT) %in% c('mz','mzmin','mzmax','rt','rtmin','rtmax'))]))
    FTRep<-F_ReplaceMZRTRange(FT,lmz=1,lmzmin=2,lmzmax=3,lrt=4,lrtmin = 5,lrtmax = 6,
                                   matched_compound,rmz=1,rrt=2,tolmz=0,tolrt=0)
    matched_compound_stat<-left_join(matched_compound,FTRep,by=c("mz","rt")) %>% data.table(.)
    # matched_compound_stat<-left_join(matched_compound_stat,Volcano,by=c("Sample")) %>% data.table(.)
    write.csv(matched_compound_stat,paste0(Cell,Label,"_matched_compound_stat.csv"),row.names = FALSE)
  }
  F_mMum(Input = "MumInputUp.txt", Label = "Up")
  F_mMum(Input = "MumInputDown.txt", Label = "Down")
  return(mMum)
}
setwd(Path)
#### Exclude the overlaping compounds between up or down regulated pathway ####
## Find the overlaping compounds in common
MumUp <- data.table(read_excel(paste0(Path,"Reference/MetGUIResults0723.xlsx"),sheet="MumUp"))
MumDown <- data.table(read_excel(paste0(Path,"Reference/MetGUIResults0723.xlsx"),sheet="MumDown"))
MumOver <- rbind(MumUp[kegg_id %in% MumDown$kegg_id],MumDown[kegg_id %in% MumUp$kegg_id])
setorder(MumOver,kegg_id) # 366
# write.csv(MumOver,'MumOver.csv')
## Check the MumOver.csv files by hand with 50 already checked overlaps, Copy the checked MumOver.csv file into OverlapCheck sheet of the MetGUIResults0723.xlsx
OverlapCheck <- data.table(read_excel(paste0(Path,"Reference/MetGUIResults0723.xlsx"),sheet="OverlapCheck"))
names(OverlapCheck) # 106
## Exclude the feature that are not selected during the checking process
FTToExclude <- MumOver[!Sample %in% OverlapCheck$Sample]
FTToExclude[Sample=='FT0406']
FTToExclude[Empirical.Compound=='EC00051'] # used to check whether or not the selection is correct
unique(FTToExclude,by='Sample')
VolcanoChecked <- Volcano[!Sample %in% FTToExclude$Sample]
VolcanoChecked[Sample=='FT0406']
## Pathway enrichment
mMum_2020 <- F_Mummichog(FT=VolcanoChecked,Cell = '2020', Exp = Exp_2020)
## Double check whether or not still overlap
MumUp <- data.table(read.csv('d:/github/MetGUI/output/2020Up_matched_compound_stat.csv'))
MumDown <- data.table(read.csv('d:/github/MetGUI/output/2020Down_matched_compound_stat.csv'))
MumOver <- rbind(ComUp[kegg_id %in% ComDown$kegg_id],ComDown[kegg_id %in% ComUp$kegg_id])
MumOver # if MumOver is empty, then it it fine
MumAll <- rbind(MumUp,MumDown)
MumAll[Empirical.Compound=='EC0001']
## Add the extra column indicating the enriched compounds and the kegg index
EnriUp <- data.table(read.csv('d:/github/MetGUI/output/2020Up_mummichog_pathway_enrichment.csv'))
EnriDown <- data.table(read.csv('d:/github/MetGUI/output/2020Down_mummichog_pathway_enrichment.csv'))
F_AddKeggName<-function(Enri, Mum){
  Strsplit<-strsplit(Enri$EC.Hits,";",fixed=TRUE)
  for (id in 1:length(Strsplit)) {
    # id <- 1
    ECSel <- Strsplit[[id]]
    ListKegg <- list()
    ListName <- list()
    for (j in 1:length(ECSel)) {
      # j <- 1
      ListKegg[[j]] <- toString(unique(Mum[Empirical.Compound==ECSel[[j]]]$kegg_id))
      ListName[[j]] <- toString(unique(Mum[Empirical.Compound==ECSel[[j]]]$name))
    }
    Enri[id,StrKegg:=paste(ListKegg, collapse = "; ")]
    Enri[id,StrName:=paste(ListName, collapse = "; ")]
  }
  return(Enri)
}
EnriUpKegg <- F_AddKeggName(Enri = EnriUp,Mum=MumUp)
EnriDownKegg <- F_AddKeggName(Enri = EnriDown,Mum=MumDown)
EnriUpKegg[Pathway.Number=='P20']$StrName




#### Create the consensus library spectra ####
## Link the library MoNA-export-LC-MS-MS_Spectra
MoNAcon <- DBI::dbConnect(RSQLite::SQLite(), "d:/github/MetGUI/input/MoNA-export-LC-MS-MS_Spectra.sqlite")
library_spectra_meta <- MoNAcon %>% dplyr::tbl("library_spectra_meta") %>% dplyr::collect() %>% as.data.table(.)
metab_compound <- MoNAcon %>% dplyr::tbl("metab_compound") %>% dplyr::collect() %>% as.data.table(.)
library_spectra <- MoNAcon %>% dplyr::tbl("library_spectra") %>% dplyr::collect() %>% as.data.table(.)

## Reduce the library by restrict the instrument
unique(library_spectra_meta$instrument)
MetaL <- library_spectra_meta[is.na(instrument) | instrument %in%
         c('LTQ Orbitrap XL Thermo Scientific','Q Exactive Plus Orbitrap Thermo Scientific',
           'Q Exactive Orbitrap Thermo Scientific','Q-Exactive + Thermo Scientific',
           'Q-Exactive Orbitrap Thermo Scientific','Q-Exactive Orbitrap Thermo Scientific',
           'Q Exactive Thermo Fisher Scientific','Q Exactive Orbitrap (Thermo Scientific)',
           'Orbitrap','Thermo Q Exactive HF','LTQ Orbitrap XL, Thermo Scientfic; HP-1100 HPLC, Agilent',
           'LTQ Orbitrap XL, Thermo Scientfic','Q-Exactive Plus','LTQ Orbitrap Velos Thermo Scientific',
           'Q Exactive Plus Orbitrap Thermo Scientific')]
MetaL[abs(precursor_mz-205.0972)<0.01]$name
MetaL[,inchikey_14:=sub(inchikey_id,pattern = "-.*",replacement = "",perl = TRUE)] # only use top 14 characters
MetaL[,inchikey_14_precursor_type:=paste0(inchikey_14,'_',precursor_type)]
MetaL[,inchikey_14_precursor_type_instrument:=paste0(inchikey_14,'_',precursor_type,'_',instrument)]
MetaL[,inchikey_14_precursor_mz:=paste0(inchikey_14,'_',precursor_mz)]
unique(MetaL$polarity)
MetaL <- MetaL[polarity!='N']
nrow(MetaL) # 38540 raw spectra from Massbank
length(unique(MetaL$inchikey_14_precursor_type)) #10238
length(unique(MetaL$inchikey_14_precursor_type_instrument)) #11074
length(unique(MetaL$inchikey_14_precursor_mz)) #10978
## Start creating conSensus library_spectra
library_spectra_Sensus <- data.table()

for (inchi in unique(MetaL$inchikey_14_precursor_mz)) {
  ## Extract the meta id related to each unique inchikey
  # inchi <- 'MTCFGRXMJLQNBG_[M+H]+_LTQ Orbitrap XL, Thermo Scientfic; HP-1100 HPLC, Agilent'
  # inchi <- 'MTCFGRXMJLQNBG_106.04987'
  print(inchi)
  Selected <- MetaL[inchikey_14_precursor_mz==inchi]
  ## Extract MS2 spectra related to each meta id
  ListMS2 <- list()
  for (j in 1:nrow(Selected)){
    # j <- 1
    SelectedMS2 <- library_spectra[library_spectra_meta_id==Selected[j,]$id]
    ListMS2[[j]]<- SelectedMS2[,c('mz','i')] %>% setnames(.,'i','intensity') %>% as.matrix(.)
  }
  ## Combine peaks
  Combined <- combinePeaks(ListMS2,ppm = 10, peaks ='union') #, peaks ='intersect',minProp = 0.1
  # Combined
  ## Add meta information
  CombinedMeta <- as.data.table(Combined) %>% setnames(.,'intensity','i') %>%
    .[,library_spectra_meta_id:=min(Selected$id)] %>% .[,inchikey_14_precursor_mz:=inchi]
  library_spectra_Sensus<- rbind(library_spectra_Sensus, CombinedMeta)
}
unique(library_spectra_Sensus,by='inchikey_14_precursor_mz')
# library_spectra_meta_Sensus$inchikey_14_precursor_mz

## Creat Sensus.sqlite
con_Sensus <- DBI::dbConnect(RSQLite::SQLite(),'d:/github/MetGUI/input/Sensus0710.sqlite')
l_dbPthValue <- 'd:/github/MetGUI/input/Sensus0710.sqlite'
## library_spectra_meta
library_spectra_meta_Sensus <- MetaL[,-c('resolution')] %>%
  left_join(., metab_compound[,-c('name')], by='inchikey_id')
setorder(library_spectra_meta_Sensus,id) # small to big
library_spectra_meta_Sensus <- unique(library_spectra_meta_Sensus, by= 'inchikey_14_precursor_mz')
DBI::dbWriteTable(con_Sensus, name = "library_spectra_meta", value = library_spectra_meta_Sensus, overwrite = T)
## Add library_spectra_source
library_spectra_source_Sensus <- data.frame(id=1,
                                            name=paste('ConSensus Database',  format(Sys.time(), "%Y-%m-%d-%I%M%S"), sep='-'),
                                            parsing_software=paste('DBI::dbWriteTable'))
DBI::dbWriteTable(con_Sensus, name='library_spectra_source',
                  value=library_spectra_source_Sensus, overwrite = T)
## Add metab_compound
metab_compound_Sensus <- metab_compound[inchikey_id %in% MetaL$inchikey_id, ]
DBI::dbWriteTable(con_Sensus, name='metab_compound',value=metab_compound_Sensus, overwrite = T)

DBI::dbWriteTable(con_Sensus, name='library_spectra', value=library_spectra_Sensus, overwrite=T)
DBI::dbDisconnect(con_Sensus)

#### Create the consensus query spectra ####
## Purity assesments and linking fragmentation to XCMS features
cdfs <- dir("D:/VirtualMachineDisk/Hematology/Centroid/NormalVSHYP", full.names = TRUE,recursive = TRUE)
cdfs
pa <- purityA(cdfs,ilim = 0,isotopes = FALSE)
# View(pa)
onDisk <- xdata_MixAA_Pos
# onDisk <- xdata_MixAA_Pos
# if (length(unique(msLevel(onDisk))) != 1) {
#   onDisk <- filterMsLevel(onDisk, msLevel = 1)
# } # to use MS1
pa_ft <- frag4feature(pa, onDisk,ppm = 580) # link the MS2 with MS1
View(pa_ft)
pa_ft_fil <- filterFragSpectra(pa_ft) # filter the fragmentation spectra , allfrag = TRUE
View(pa_ft_fil)
pa_ft_fil_ave <- averageAllFragSpectra(pa_ft_fil) # treat the inter and intra fragmentation scans the same
# View(pa_ft_fil_ave)

# q_dbPthValue <- createDatabase(pa_ft_fil_ave, onDisk, dbName = "q_dbPth0708.sqlite")
q_dbPthValue <- 'd:/github/MetGUI/output/q_dbPth0708.sqlite'
qon <- DBI::dbConnect(RSQLite::SQLite(), q_dbPthValue)
## Extract from the query database the "c_peak_groups", which will be used for subset.
c_peak_groups <- qon %>% dplyr::tbl("c_peak_groups") %>% dplyr::collect() %>% as.data.table(.)
names(c_peak_groups)
PeakGroups <- c_peak_groups[,c('mz','mzmin','mzmax','rt','rtmin','rtmax','grpid','grp_name')]
nrow(PeakGroups) # 2192 consensus spectra
# View(PeakGroups)
## Extract the query with the MS2 peaks
s_peaks <- qon %>% dplyr::tbl("s_peaks") %>% dplyr::collect() %>% as.data.table(.)
## Extract the query with the related xcms groups
s_peak_meta <- qon %>% dplyr::tbl("s_peak_meta") %>% dplyr::collect() %>% as.data.table(.)
names(s_peak_meta)
nrow(s_peak_meta)
s_peak_meta_rbind <-s_peak_meta
# View(s_peak_meta)
PeakMeta <- s_peak_meta[,c('precursor_mz','retention_time','grpid','pid','spectrum_type','precursorIntensity')]
## change the spectrum_type and grpid for s_peak_meta file
s_peak_meta_withgrpid <- s_peak_meta[!is.na(grpid)]
s_peak_meta_nogrpid <- s_peak_meta[is.na(grpid)]
setorder(PeakMeta,precursorIntensity)

PeakGroupsRep<-F_ReplaceMZRTRange(PeakGroups,lmz=1,lmzmin=2,lmzmax=3,lrt=4,lrtmin = 5,lrtmax = 6,
                                PeakMeta,rmz=1,rrt=2,tolmz=0.07,tolrt=1)
setnames(PeakGroupsRep,c('grpid','mz','rt'),c('Grpid','precursor_mz','retention_time'))
# setnames(PeakGroupsRep,c('mz','rt'),c('precursor_mz','retention_time'))
s_peak_meta_nogrpid_join <- left_join(s_peak_meta_nogrpid,PeakGroupsRep,by=c('precursor_mz','retention_time'))
s_peak_meta_nogrpid_join[mzmin>0, spectrum_type:='all']
s_peak_meta_nogrpid_join[mzmin>0, grpid:=Grpid]
s_peak_meta_nogrpid_join[mzmin>0]
s_peak_meta_rbind <- rbind(s_peak_meta_withgrpid,
                           s_peak_meta_nogrpid_join[,-c('mzmin','mzmax','rtmin','rtmax','Grpid','grp_name')])
nrow(s_peak_meta_rbind) # 51458 raw spectra from 18 experimental dataset
## Check with the 20 AAs
s_peak_meta_rbind$grpid <- as.numeric(s_peak_meta_rbind$grpid)
s_peak_meta_rbind_PeakGroups <- left_join(s_peak_meta_rbind[spectrum_type=='all'],PeakGroups,by='grpid')
s_peak_meta_rbind_PeakGroups_Check <- F_CheckMH20AA(DT=s_peak_meta_rbind_PeakGroups,MH20AA)
View(s_peak_meta_rbind_PeakGroups_Check)
# View(c_peak_groups)
## Update the s_peak_meta file in the database
DBI::dbWriteTable(qon, name='s_peak_meta',value=s_peak_meta_rbind, overwrite = T)
DBI::dbDisconnect(qon)

#### MS2 spectra matching based on the modified msPurity package ####
## Library search without weighted intensity
setwd(paste0(Path,'output/LibrarySearch/0710'))
MRPValue <<- 7500
# MRPValue <<- 17500
RefMZValue <<- 400
# RefMZValue <<- 200
q_pidValue <- s_peak_meta_rbind[!is.na(grpid)]$pid
l_dbPthValue0710 <- 'd:/github/MetGUI/input/Sensus0710.sqlite'
l_dbPthValue0708 <- 'd:/github/MetGUI/input/Sensus0708.sqlite'
l_dbPthValueSensus <- 'd:/github/MetGUI/input/Sensus.sqlite'
for (id in unique(q_pidValue)) {
  # id <- 44036
  spectralMatching(q_dbPth=q_dbPthValue, l_dbPth=l_dbPthValue0710,
                   l_ppmPrec = 20, q_pids= id, mztol = NA)
}
## Library search with weighted intensity
# setwd(paste0(Path,'output/LibrarySearch/mz2int0.4')) # decoy also used the weighted intensity
setwd(paste0(Path,'output/LibrarySearch/mz2int0.4_0827')) # decoy do not use the weighted intensity
for (id in unique(q_pidValue)) {
  # id <- 339
  spectralMatching(q_dbPth=q_dbPthValue, l_dbPth=l_dbPthValue0710,
                   l_ppmPrec = 20, q_pids= id,  mztol = NA,raW = 0.4,mzW = 2)
}
## For testing purpose using the 20 AAs
# q_pidAA <-c(51120,7696,44036,9621,30286,51162,51163,51170,51173,51182,10757,13381,1226,5031,931,2196,51208,955)
# for (id in unique(q_pidAA)) {
#   # id <- 44036
#   spectralMatching(q_dbPth=q_dbPthValue, l_dbPth=l_dbPthValue0710,
#                    l_ppmPrec = 20, q_pids= id, mztol = NA)
# }


## Merge the files together, count candidates

Dir <- paste0(Path, "output/LibrarySearch/0710")
setwd(Dir)
MergedDT <- F_MergeLibrarySearch(Dir= paste0(Dir))

Dir2.4 <- paste0(Path, "output/LibrarySearch/mz2int0.4_0827")
setwd(Dir2.4)
MergedDT2.4 <- F_MergeLibrarySearch(Dir= paste0(Dir2.4))

## Add meta and stat
setwd(Path)
F_AddMetaStat <- function(DT,Meta,Vol){
  ## Add meta data
  setnames(Meta,c('pid','grp_name'),c('query_qpid','Sample'),skip_absent = TRUE)
  Meta<- Meta[,c('query_qpid','Sample')] #,'mz','mzmin','mzmax','rt','rtmin','rtmax'
  DT <- left_join(DT,Meta,by='query_qpid') # add the meta dat
  DT <- left_join(DT,Vol,by='Sample') # Add the statistics
  DT <- DT[!is.na(FC)] # filtration based on the Volcano NA, as these are the features after isotop filtration
  return(DT)
}
SingleWhole <- F_AddMetaStat(DT=MergedDT,Meta=s_peak_meta_rbind_PeakGroups,Vol=Volcano)
SingleWhole2.4 <- F_AddMetaStat(DT=MergedDT2.4,Meta=s_peak_meta_rbind_PeakGroups,Vol=Volcano)
# write.csv(SingleWhole,'SingleWhole.csv')
# write.csv(SingleWhole2.4,'SingleWhole2.4.csv')
## Check the identification using the 20 AAs
MH20AA<-data.table(read_excel(paste0(Path,"Reference/AAMass.xlsx"),sheet="M+H"))
F_CheckMH20AA <- function(DT, MH20AA){
  NameCol <- c('mz','mzmin','mzmax','rt','rtmin','rtmax')
  DT <- unique(DT,by=c('mz','rt'))
  setcolorder(DT, c(NameCol,colnames(DT)[!(colnames(DT) %in% NameCol )]))
  NameCol <- c('mz','rt')
  setcolorder(MH20AA, c(NameCol,colnames(MH20AA)[!(colnames(MH20AA) %in% NameCol )]))
  DTRep<-F_ReplaceMZRTRange(DT,lmz=1,lmzmin=2,lmzmax=3,lrt=4,lrtmin = 5,lrtmax = 6,
                            MH20AA,rmz=1,rrt=2,tolmz=0.02,tolrt=10)
  DTRep_MH20AA <- left_join(DTRep,MH20AA,by=c("mz","rt")) %>% data.table(.) %>%
    setorder(.,mz) %>% .[!is.na(letter3)]
  NameCol <- c('letter3','Monoisotopic','Adduct')
  setcolorder(DTRep_MH20AA , c(NameCol,colnames(DTRep_MH20AA )[!(colnames(DTRep_MH20AA ) %in% NameCol )]))
  return(DTRep_MH20AA)
}
SingleWhole[abs(library_precursor_mz - 106.0499 )<0.02]
SingleWhole2.4[abs(library_precursor_mz - 106.0499 )<0.02]
SingleWhole_Check <- F_CheckMH20AA(DT=SingleWhole,MH20AA)
View(SingleWhole_Check)
SingleWhole[Sample=='FT0054']
# write.csv(SingleWhole_Check, 'SingleWhole_Check.csv',row.names = F)
PeakGroups[grpid==54] # Check the reason for the missing 20 AAs
s_peaks[pid==24620]
s_peak_meta_rbind[grpid==337]
library_spectra_Sensus[library_spectra_meta_id==39160]
library_spectra_Sensus[inchikey_14_precursor_mz=='DHMQDGOQFOQNFH_76.0394']

## Precursor filtration
F_PrecursorFilt <- function(DT, MatchCutoff=1){
  print(paste0('before filtration ',nrow(DT)))
  DT[, MZDiff := abs(library_precursor_mz - query_precursor_mz)]
  DT[, MS1Tol := F_CalPMMT(library_precursor_mz, 30000, 400)]
  # DT <- DT[MZDiff < MS1Tol]
  DT <- DT[Match >= MatchCutoff]
  print(paste0('after match numbers filtration ',nrow(DT)))
  # F_CSV(DT)
  # View(DT)
  return(DT)
}
SingleWhole_Pre <- F_PrecursorFilt(DT=SingleWhole,MatchCutoff = 2)
SingleWhole_Pre2.4 <- F_PrecursorFilt(DT=SingleWhole2.4,MatchCutoff = 2)
nrow(SingleWhole_Pre) # 1772    candidates
nrow(SingleWhole_Pre2.4)

SingleWhole_Pre[Sample=='FT0054']
# unique(SingleWhole$query_qpid)
# s_peak_meta_rbind_PeakGroups[!is.na(grpid)]
# write.csv(SingleWhole_Pre, 'SingleWhole_Pre.csv',row.names = F)
SingleWhole_Pre_Check <- F_CheckMH20AA(DT=SingleWhole_Pre,MH20AA) ## Check whether or not the 20 AAs are identified
View(SingleWhole_Pre_Check)
write.csv(SingleWhole_Pre_Check,'SingleWhole_Pre_Check.csv')
#### FDR and q-value estimation ####
F_FDRCutoff <- function(DT, MatchCutoff=2,TopCutoff=10, NumCutoff=1000,CsvName){
  FDRCutoff <- 0.05
  print(paste('Before match cutoff filtration',MatchCutoff,nrow(DT)))
  DT <- DT[Match>=MatchCutoff]
  print(paste('After match cutoff filtration',MatchCutoff,nrow(DT)))
  setorder(DT,-dpc.xcorr) # decreasing of the score according to dpc
  NumCutoffMin <- min(nrow(DT),NumCutoff)
  DT <- DT[1:NumCutoffMin]
  print(paste('After NumCutoff filtration',NumCutoff,nrow(DT)))
  print('Select Top for each query based on TopCutoff')
  TarDynamic <- setorder(DT[mztol=='NA'], -dpc) %>% .[, head(.SD, TopCutoff), by='query_qpid'] %>%
    .[, Database:='Target'] %>% setnames(.,'dpc','Thres')
  DecDynamic <- setorder(DT[mztol=='NA'], -dpc.decoy) %>% .[, head(.SD, TopCutoff), by='query_qpid'] %>%
    .[, Database:='Decoy'] %>%  setnames(.,'dpc.decoy','Thres')
  print(paste('TopCutoff',TopCutoff,'Before',nrow(DT),'After',nrow(TarDynamic)))
  Dynamic <- rbind(TarDynamic, DecDynamic,fill=TRUE)
  # setorder(Dynamic,-Thres)
  DynamicDis <- Dynamic
  print('Calculate the estimated qValue based on Pep accumulation')
  # Lambda<-getPEPFromScoreLambda(DynamicDis[Database=='Target']$Thres, DynamicDis[Database=='Decoy']$Thres,
  #                               paste0('FDR cutoff prediction DynamicCurve'))

  LambdaDynamic<-F_getPEPFromScoreLambda(DynamicDis[Database=='Target']$Thres, #& Thres>ThresCutoff
                                         DynamicDis[Database=='Decoy']$Thres, paste0('FDR cutoff prediction DynamicCurve'))
  Dynamic[,pep:=sapply(Dynamic$Thres, LambdaDynamic[[1]])]
  Dynamic[, SumPep := cumsum(pep)]
  Dynamic[, Length:=seq_len(length(pep))]
  Dynamic[, qValue_Pep:=rev(cummin( rev(SumPep / Length)))]
  # print('Calculate the estimated qValue based on target decoy')
  Dynamic[, decoy_hit := grepl("Decoy",Dynamic$Database)]
  Dynamic[, NDecoy := cumsum(decoy_hit)]
  Dynamic[, NTarget := Length - NDecoy]
  Dynamic[,qValue_DecoyTarget:=rev(cummin( rev(NDecoy / NTarget)))]
  Colu <- c('Sample','query_precursor_mz','query_rt','Thres','query_qpid','library_lpid',
            'library_entry_name','library_inchikey','library_precursor_type',
            'FC','log2.FC.','raw.pval','p.ajusted','npeaks',
            'dpc.decoy', 'pep','qValue_Pep',
            'mz','mzmin','mzmax','rt','rtmin','rtmax',
            'SumPep','Length')
  setcolorder(Dynamic, c(Colu,colnames(Dynamic)[!(colnames(Dynamic) %in% Colu)]))
  write.csv(Dynamic,paste(CsvName,'DynamicWithQvalue.csv'))
  ## Density plot
  Pep_Thres <- round(min(Dynamic[qValue_Pep <= FDRCutoff]$Thres), digits = 2)
  # FDR_Thres <- round(min(Dynamic[qValue_DecoyTarget <= FDRCutoff]$Thres), digits = 2)
  print(paste("Pep_Thres", Pep_Thres))
  # Score.dpc<-rbind(data.table('Score'=Dynamic[,dpc],'Legend'='Target'),data.table('Score'=Dynamic[,dpc.decoy.Dec],'Legend'='Decoy'))
  g_density <- ggplot(DynamicDis)  +
    ylab('Density') + xlab('Score') + theme_classic() + theme(legend.position="top") + #ylim(0,3000)+
    # geom_density(aes(x = Thres, y = after_stat(count),colour = Database), adjust = 2) +
    geom_density(aes(x = Thres,  colour = Database)) + #..scaled..,
    # geom_vline(xintercept = FDR_Thres, linetype="dotted",color = "orange", size=1.5)+
    # annotate(geom = "label", x = FDR_Thres, y = 1000, label = FDR_Thres,color = "orange")+
    geom_vline(xintercept = Pep_Thres, linetype="dotted",color = "red", size=1.5)+
    annotate(geom = "label", x = Pep_Thres, y = 1, label = Pep_Thres,color = "red")+
    scale_colour_manual( name="Legend",values = c("Target" = "blue","Decoy" = "green"))
  return(g_density)
}
F_FDRCutoffComp <- function(Normal, Weighted,FileName,MatchCutoff,TopCutoff,NumCutoff){
  g_Normal <- F_FDRCutoff(DT = Normal, MatchCutoff,TopCutoff,NumCutoff,CsvName = 'Normal')
  g_Weighted <- F_FDRCutoff(DT = Weighted, MatchCutoff,TopCutoff,NumCutoff,CsvName = 'Weighted')
  grDevices::svg(paste(FileName,'MatchCutoff',MatchCutoff,'TopCutoff',TopCutoff,'NumCutoff',NumCutoff,'.svg')
                 ,width = 14, height = 7) #
  print(plot_grid(g_Normal, g_Weighted,align="h")) #,ncol = 2,nrow = 1
  dev.off()
}
# SingleWhole_PreUni <- unique(SingleWhole_Pre,by=c('query_qpid','library_lpid'))
# SingleWhole_Pre2.4Uni <- unique(SingleWhole_Pre2.4,by=c('query_qpid','library_lpid'))
F_FDRCutoffComp(Normal = SingleWhole_Pre,Weighted = SingleWhole_Pre2.4,
                FileName = 'CompareSingleWhole_Pre',
                MatchCutoff = 2,TopCutoff = 10, NumCutoff = 1000)
## Plot the scatter plot of similarity scores at different cutoff
WeightedQvalue <- data.table(read.csv('Weighted DynamicWithQvalue.csv'))
NormalQvalue <- data.table(read.csv('Normal DynamicWithQvalue.csv'))
for (FDRCutoff in seq(0.01,0.5,0.01)) {
  # FDRCutoff <- 0.1
  # WeightedQvalue <- WeightedQvalue[Database=='Target']
  # NormalQvalue <- NormalQvalue[Database=='Target']
  Weighted_Thres <- round(min(WeightedQvalue[qValue_DecoyTarget <= FDRCutoff]$Thres), digits = 2)
  Normal_Thres <- round(min(NormalQvalue[qValue_DecoyTarget <= FDRCutoff]$Thres), digits = 2)
  print(paste('FDRCutoff',FDRCutoff,'Weighted_Thres',Weighted_Thres,'Normal_Thres',Normal_Thres))
}



SingleWhole_Pre_Top10 <- F_FDRCutoff(DT = SingleWhole_Pre, FileName='SingleWhole_Pre_Top10',
                                     MatchCutoff = 2,TopCutoff = 10, XcorrCutoff = 0.1)
nrow(SingleWhole_Pre_Top10[qValue_Pep<=0.05]) #453 candidates
SingleWhole_Pre_Top10_Check <- F_CheckMH20AA(DT=SingleWhole_Pre_Top10,MH20AA) ## Check whether or not the 20 AAs are identified
View(SingleWhole_Pre_Top10_Check) # only Gly and Cys are not identified

## Select the top candidates for each query
SingleWhole_Pre_Top10_Uni <- setorder(SingleWhole_Pre_Top10, -Thres) %>%
  .[, head(.SD, 1), by='Sample'] %>% .[, head(.SD, 1), by='library_inchikey']
nrow(SingleWhole_Pre_Top10_Uni[qValue_Pep<=0.05]) #88 candidates
SingleWhole_Pre_Top10_Uni_Check <- F_CheckMH20AA(DT=SingleWhole_Pre_Top10_Uni,MH20AA) ## Check whether or not the 20 AAs are identified
View(SingleWhole_Pre_Top10_Uni_Check)


#### Write results ####
setwd(paste0(Path))
MetGUIResults <- list("Volcano" = Volcano, "SingleWhole_Pre" = SingleWhole_Pre,
                      "SingleWhole_Pre_Top10" = SingleWhole_Pre_Top10,
                      "SingleWhole_Pre_Top10_Uni" = SingleWhole_Pre_Top10_Uni,
                      'EnriUp' = EnriUp, 'EnriDown' = EnriDown,
                      'EnriUpKegg' = EnriUpKegg, 'EnriDownKegg' = EnriDownKegg,
                      'MumUp' = MumUp, 'MumDown' = MumDown)
write_xlsx(MetGUIResults, "MetGUIResults.xlsx")
## Check the results of 20 AAs, save the results as Check sheet in MetGUIResults0710.xlsx
## Creat new column called ShortName in the Volcano sheet, label AAs if exist, NA if not
#### MS2 visualization ####
setwd(paste0(Path,'MirrorPlot/Check0710'))
F_SpectrumSingle <- function(x) {
  query_qpid <- x[5]
  query_precursor_mz <- as.numeric(x[2])
  dpc.dynamic <- as.numeric(x[4])
  library_lpid <- x[6]
  library_entry_name <- x[7]
  library_inchikey <- x[8]
  Match <- x[27]
  library_accession <- x[30]
  inchikey_14 <- sub(library_inchikey,pattern = "-.*",replacement = "",perl = TRUE)
  # library_lpid <- 27440
  # inchikey_14 <- gsub(inchikey_14, pattern = "\\[",replacement = "_",perl = TRUE)
  # inchikey_14 <- 'inchikey_14'
  # query_qpid <- 16135 # for Glycine
  # library_lpid <- 43890 # for Glycine
  SpecQ <-s_peaks[pid==query_qpid]
  SpecL <-library_spectra_Sensus[library_spectra_meta_id==library_lpid]
  if (length(SpecQ) < 1 | length(SpecQ) < 1) {
    next
  } else {
    spec.top <- data.table("mz" = SpecQ$mz, "intensity" = SpecQ$i)
    spec.bottom <- data.table("mz" = SpecL$mz, "intensity" = SpecL$i)
    ## Target
    top_tmp <- data.frame(mz = spec.top[, 1], intensity = spec.top[, 2])
    top_tmp$normalized <- (top_tmp$intensity / max(top_tmp$intensity)) * 100
    top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized) # data frame for plotting spectrum

    bottom_tmp <- data.frame(mz = spec.bottom[, 1], intensity = spec.bottom[, 2])
    bottom_tmp$normalized <- (bottom_tmp$intensity / max(bottom_tmp$intensity)) * 100
    bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized) # data frame for plotting spectrum
    AlignDynamic <- F_DynamicMatching(top = top_plot, bottom = bottom_plot,RefMZ = 400,MRP=7500)
    dpD <- dpc.dynamic
    # dpD <- 0.95
    ## Plot
    png(paste('precursor', round(query_precursor_mz, digits = 3),'score', round(dpD, digits = 3),
              'Match',Match,"inchikey_14", inchikey_14,
              # 'library_entry_name ',library_entry_name,
              "query_qpid",  query_qpid,"library_lpid", library_lpid,
              ".png"), width = 1000, height = 1100, res = 150)
    F_PlotMSMS(alignment = AlignDynamic, top_plot = top_plot, bottom_plot = bottom_plot)
    title(paste0(library_entry_name), outer = F,cex.main = 2)
    dev.off()
  }
}

ForMirror <- data.table(read_xlsx(paste0(Path,"MetGUIResults0710.xlsx"), sheet = "Check"))
# ForMirror <- SingleWhole_Pre_Match2_Top %>% unique(.,by=c('query_qpid', 'library_lpid'))
names(ForMirror)
setorder(ForMirror, query_precursor_mz)
View(ForMirror)
for (id in (1:nrow(ForMirror))) {
  # id <- 39
  print(id)
  F_SpectrumSingle(x = as.matrix(ForMirror[id]))
}
## Plot the missing Glycine and cystine
GlyCys <- SingleWhole[Sample %in% c('FT0054','FT0337') ]
setnames(GlyCys,'dpc','Thres')
GlyCys$pep <- 'pep'
GlyCys$qValue_Pep <- 'qValue_Pep'
GlyCys$SumPep <- 'SumPep'
GlyCys$Length <- 'Length'
Colu <- c('Sample','query_precursor_mz','query_rt','Thres','query_qpid','library_lpid',
          'library_entry_name','library_inchikey','library_precursor_type',
          'FC','log2.FC.','raw.pval','p.ajusted','npeaks',
          'dpc.decoy', 'pep','qValue_Pep',
          'mz','mzmin','mzmax','rt','rtmin','rtmax',
          'SumPep','Length')
setcolorder(GlyCys, c(Colu,colnames(GlyCys)[!(colnames(GlyCys) %in% Colu)]))
names(GlyCys)
GlyCys[Sample=='FT0054']
F_PlotMirMissing <- function(query_qpid,library_lpid,Name){
  SpecQ <-s_peaks[pid==query_qpid]
  SpecL <-library_spectra_Sensus[library_spectra_meta_id==library_lpid]
  spec.top <- data.table("mz" = SpecQ$mz, "intensity" = SpecQ$i)
  spec.bottom <- data.table("mz" = SpecL$mz, "intensity" = SpecL$i)
  ## Target
  top_tmp <- data.frame(mz = spec.top[, 1], intensity = spec.top[, 2])
  top_tmp$normalized <- (top_tmp$intensity / max(top_tmp$intensity)) * 100
  top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized) # data frame for plotting spectrum

  bottom_tmp <- data.frame(mz = spec.bottom[, 1], intensity = spec.bottom[, 2])
  bottom_tmp$normalized <- (bottom_tmp$intensity / max(bottom_tmp$intensity)) * 100
  bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized) # data frame for plotting spectrum
  AlignDynamic <- F_DynamicMatching(top = top_plot, bottom = bottom_plot,RefMZ = 400,MRP=7500)
  ## Plot
  png(paste('PlotMirMissing', Name,".png"), width = 1000, height = 1100, res = 150)
  F_PlotMSMS(alignment = AlignDynamic, top_plot = top_plot, bottom_plot = bottom_plot)
  title(paste0('Name = ',Name), outer = F,cex.main = 2)
  dev.off()
}
F_PlotMirMissing(query_qpid=1676,library_lpid=21178,Name='Cysteine')

#### EIC visualization for the 20 AAs ####
## Extract the features that are related to 20 AAs
NameCol <- c('mz','mzmin','mzmax','rt','rtmin','rtmax')
setcolorder(Volcano, c(NameCol,colnames(Volcano)[!(colnames(Volcano) %in% NameCol )]))
MH20AA<-data.table(read_excel(paste0(Path,"Reference/AAMass.xlsx"),sheet="M+H"))
NameCol <- c('mz','rt')
setcolorder(MH20AA, c(NameCol,colnames(MH20AA)[!(colnames(MH20AA) %in% NameCol )]))
VolcanoRep<-F_ReplaceMZRTRange(Volcano,lmz=1,lmzmin=2,lmzmax=3,lrt=4,lrtmin = 5,lrtmax = 6,
                               MH20AA,rmz=1,rrt=2,tolmz=0.02,tolrt=10)
Volcano_MH20AA <- left_join(VolcanoRep,MH20AA,by=c("mz","rt")) %>% data.table(.) %>%
  setorder(.,mz) %>% .[!is.na(letter3)]

NameCol <- c('letter3','Monoisotopic','Adduct')
setcolorder(Volcano_MH20AA, c(NameCol,colnames(Volcano_MH20AA)[!(colnames(Volcano_MH20AA) %in% NameCol )]))
# View(Volcano)
write_xlsx(Volcano_MH20AA, 'Volcano_MH20AA.xlsx') # add into the MetGUIResults.xlsx
## Extract the unique feature according to rt and mz
matched_compound[,mzRound:=round(mz,digits = 2)]
matched_compound[,rtRound:=round(rt,digits = 0)]
MCUni<-unique(matched_compound,by=c("mzRound","rtRound"))
# MCUniNeg<-MCUni[(adduct=="M-H [1-]")|(adduct=="M-H2O-H [1-]")|(adduct=="M-2H [2-]")]
MCUniPos<-MCUni[(adduct=="M+H [1+]")|(adduct=="M-H2O+H [1+]")|(adduct=="M+Na [1+]")]

## EIC of the AA mixture
names(MCUniPos)
F_EICLoopMix<-function(DT,massI,addI,mzI,rtI,Raw){
  # group_colors <- c("Black","coral1","cyan1")
  names(group_colors) <- unique(Raw$sample_group)
  F_EICSingle<-function(x){
    ExactMass<-as.numeric(x[[massI]])
    Adduct<-x[[addI]]
    mz<-as.numeric(x[[mzI]])
    rt<-as.numeric(x[[rtI]])
    mzValue<-c(mz-0.02,mz+0.02)
    rtValue<-c(rt-30,rt+30)
    EIC<-chromatogram(Raw, aggregationFun = "sum",mz = mzValue, rt = rtValue)
    name<-paste("mass",ExactMass,Adduct,"mz",mz,"rt",rt,"mixture.png",seq="")
    png(file=name,res=300,width = 2000,height = 2000)
    plot(EIC, col = group_colors[Raw$sample_group],peakType = "none",
         cex.axis=1.5,cex.lab = 2,cex.main = 2)
    legend("topright",legend=unique(Raw$sample_group),col=unique(group_colors[Raw$sample_group]), lty=1, cex=0.8,box.lty=0)
    # title(sub=paste("mass",ExactMass,Adduct,"mz",mz,"rt",round(rt),seq=""),cex.sub = 2)
    dev.off()
  }
  apply(DT,1,F_EICSingle)
}
115.0633+1.007276
ForMirror[,mass:=mz-1.007276]
names(ForMirror)
F_EICLoopMix(ForMirror[abs(mass-174.11)<0.005],massI = 62,addI = 9,mzI = 18,
             rtI = 21,Raw=xdata_MixAA_Pos)

# F_EICLoopMix(unique(AAVolcano[(massRound==147.05)&(mode=="negative")],by="mass"),massI = 46,addI = 3,mzI = 47,rtI = 48,Raw_MixAA_Neg)
# F_EICLoopMix(unique(AAVolcano[(mode=="positive")],by="mass"),
#              massI = 46,addI = 3,mzI = 47,rtI = 48,xdata_MixAA_Pos)
# MCUniVisualized<-rbind(MCUniNeg,MCUniPos)

## EIC of the pure AAs
names(AAVolcano[(massRound==132.05)&(mode=="negative")])
F_EICLoopPureAA(MCUniPos,massI = 15,addI = 3,mzI = 16,rtI = 17,Raw=xdata_PureAA_Pos22)
# F_EICLoopPureAA(MCUniNeg[(massRound==147.05)&(rtRound==897)],
#              massI = 15,addI = 3,mzI = 16,rtI = 17,xdata_PureAA_Neg)
# F_EICLoopPureAA(AAVolcano[(massRound==147.05)&(mode=="negative")],massI = 46,addI = 3,mzI = 47,rtI = 48,xdata_PureAA_Neg)
F_EICLoopPureAA(AAVolcano[(mode=="positive")],massI = 46,addI = 3,mzI = 47,rtI = 48,xdata_PureAA_Pos)

## EIC of the mixture single
ExactMass<-115.0633
Adduct<-'M+H [1+]'
mz<-116.0706
rt<-518
Raw<-xdata_MixAA_Pos
mzValue<-c(mz-0.02,mz+0.02)
rtValue<-c(rt-30,rt+30)
EIC<-chromatogram(Raw, aggregationFun = "sum",mz = mzValue, rt = rtValue)
name<-paste("mass",ExactMass,Adduct,"mz",mz,"rt",rt,"SingleOutput.png",seq="")
png(file=name,res=200,width = 1000,height = 1000)
plot(EIC, col = group_colors[Raw$sample_group],peakType = "none")
legend("topright",legend=unique(Raw$sample_group),col=unique(group_colors[Raw$sample_group]), lty=1, cex=0.8,box.lty=0)
title(sub=paste("mass",round(ExactMass,digits = 2),Adduct,"mz",round(mz,digits = 2),"rt",rt,seq=""))
dev.off()



#### Volcano plot EnhancedVolcano package####
# https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
setwd(Path)
## Save the matched_compound.xlsx as matched_compound_20AAs.xlsx, add extra sheet named AA Volcano
# PureAA<-data.table(read_xlsx(paste0(Path,"MetGUIResults0710.xlsx"),sheet="Volcano")) %>%
#   .[!is.na(ShortName)] # For pure AAs
# names(PureAA)
# Feature <-data.table(read_xlsx(paste0(Path,"MetGUIResults0710.xlsx"),sheet="Check")) # For features
# Others <- anti_join(Feature,PureAA,by='Sample')
# Merge<-rbind(Others,PureAA,fill = TRUE)
Merge <-data.table(read_xlsx(paste0(Path,"MetGUIResults0710.xlsx"),sheet="Volcano")) # For features
head(Merge)
setnames(Merge,"log2.FC.","log2FoldChange")
setnames(Merge,"raw.pval","pvalue")
# setorder(Merge,mz) # from small to big
# Merge<-unique(Merge,by=c("mz","rt","ShortName"))
# Merge$adduct<-sub(Merge$adduct, pattern = " .*",replacement = "",perl = TRUE) # clean the adduct
names(Merge)
Merge[,13:30]=lapply(Merge[,13:30],type.convert,as.is=TRUE)
Merge[,Mean:=apply(Merge[,13:30],1,mean,na.rm=TRUE)]
# Merge[,ShortName:= letter3]
# DT = Merge

F_VolcanoPlot<-function(DT,Name,FC=0.6,P=0.05,NSAA='pink',SigAA='red2',NSMet='lightblue',SigMet='royalblue'){
  # Customize colour
  keyvals.colour <- ifelse(
    (abs(DT$log2FoldChange) > FC)&(abs(DT$pvalue) < P)&(is.na(DT$ShortName)), SigMet, #
    ifelse((abs(DT$log2FoldChange) > FC)&(abs(DT$pvalue) < P)&(!is.na(DT$ShortName)), SigAA,
           ifelse((!is.na(DT$ShortName)), NSAA,
                  NSMet)))
  keyvals.colour[is.na(keyvals.colour)] <-NSMet
  names(keyvals.colour)[keyvals.colour == NSAA] <- 'NSAA'
  names(keyvals.colour)[keyvals.colour == SigAA] <- 'SigAA'
  names(keyvals.colour)[keyvals.colour == NSMet] <- 'NSMet'
  names(keyvals.colour)[keyvals.colour == SigMet] <- 'SigMet'
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
  names(keyvals.shape)[keyvals.shape == 46] <- paste(" 0 < Int <",round(LowerHinge,digits = 0))
  names(keyvals.shape)[keyvals.shape == 20] <- paste(round(LowerHinge,digits = 0),"< Int <",round(UpperHinge,digits = 0))
  names(keyvals.shape)[keyvals.shape == 19] <- paste(round(UpperHinge,digits = 0),"< Int < Inf")
  EnhancedVolcano(DT,x = "log2FoldChange",y = "pvalue",
                  pCutoff = 0.05,FCcutoff = 0.6,
                  xlim = c(min(DT[,log2FoldChange], na.rm=TRUE), max(DT[,log2FoldChange], na.rm=TRUE)),
                  ylim = c(0, max(-log10(DT[,pvalue]), na.rm=TRUE) + 1),
                  pointSize = 3,
                  # lab = DT$ShortName,
                  # selectLab = DT[,ShortName][which(names(keyvals.colour) %in% c('AAs Significant'))],
                  lab = DT$ShortName,
                  selectLab = DT[,ShortName][which(names(keyvals.colour) %in% c('NSAA','SigAA'))],
                  axisLabSize = 22,
                  labSize = 8,drawConnectors = FALSE, # increase label size 5
                  shapeCustom = keyvals.shape,
                  colCustom = keyvals.colour,colAlpha = 0.5,
                  title = Name,subtitle = "",
                  caption = "HYP / Norm, log2FC 0.6, pvalue 0.05",
                  legendPosition = "", legendLabSize = 6, # legendPosition = "" to get pure volcano without legend
                  gridlines.major=FALSE,gridlines.minor=FALSE,
                  widthConnectors = 0.5)
  ggsave(file=paste("VolcanoPlot",Name,".jpeg"),dpi=300)
}

F_VolcanoPlot(DT=Merge,"Volcano plot",FC=0.6,P=0.05)
#### No use Venn diagram between massbank, mzcloud and 142 priority ####
setwd(paste0(Path,'/Mzcloud'))
## Process massbank results
Massbank <- Iden_2020[precursorMz>0]
setcolorder(Massbank, c('mz','mzmin','mzmax','rt','rtmin','rtmax',colnames(Massbank)[!(colnames(Massbank) %in% c('mz','mzmin','mzmax','rt','rtmin','rtmax'))]))
names(Massbank)
setorder(Massbank, -score)
MassbankUni <- unique(Massbank, by='target_compound_name')
setnames(MassbankUni,'target_formula','Formula')
## Process mzcloud results
mzCloud <- data.table(read_xlsx(paste0(Path,'Mzcloud/mzCloud.xlsx')))
View(mzCloud)
names(mzCloud)
mzCloud$Formula
mzCloud$Formula <- gsub(mzCloud$Formula, pattern = " ",replacement = "",perl = TRUE) # clean the adduct
mzCloud[,mz:= `Molecular Weight`+1.007276]
setorder(mzCloud,-`Best Match`, -`Best Sim. Match`)
mzCloudUni <- unique(mzCloud, by='Formula')
## Process 142 Priority
Priority <- data.table(read_xlsx(paste0(Path,"Reference/List of metabolite priorities 210308.xlsx"), sheet = "FinalClean"))
names(Priority)
head(Priority)
nrow(Priority)


## Venn between massbank, mzcloud and the 142 priority
names(MassbankUni)
names(mzCloudUni)
names(Priority)
ListFormula = list(
  "mzCloud" = mzCloudUni$Formula,
  "Massbank" = MassbankUni$Formula,
  "Priority" = as.character(Priority$Formula)
)
venn.plot <- venn.diagram(
  x = ListFormula,
  disable.logging = TRUE,
  force.unique = FALSE,
  euler.d = TRUE,
  fill = c("cornflowerblue", "green", "yellow"),
  alpha = 0.50,
  filename = 'Venn.png',resolution = 300,
  # cat.pos = c(-20, 0, 20),
  cat.dist = c(0.05, 0.05, 0.02),
  cex = 2.5,
  cat.cex = 2.5,
  reverse = TRUE
);
overlap <- calculate.overlap(ListFormula)
## Find overlap based on molecular formula
# MassbankRep<-F_ReplaceMZRTRange(Massbank,lmz=1,lmzmin=2,lmzmax=3,lrt=4,lrtmin = 5,lrtmax = 6,
#                                 Compounds,rmz=1,rrt=2,tolmz=0,tolrt=0)
# setcolorder(Compounds, c('mz','rt','Molecular Weight',colnames(Compounds)[!(colnames(Compounds) %in% c('mz','rt','Molecular Weight'))]))
# MassbankUni_mzCloudUni_inner <- inner_join(MassbankUni, mzCloudUni, by = c("Formula")) %>% data.table(.)
# MassbankUni_anti <- anti_join(MassbankUni, mzCloudUni, by = c("Formula")) %>% data.table(.)
# mzCloudUni_anti <- anti_join(mzCloudUni,MassbankUni , by = c("Formula")) %>% data.table(.)

#### No use Bland for the rt difference between positive and negative####
F_Bland(Pos$rt,Neg$rt,"AA RT Difference (Seconds) Positive-Negative")
NegAA<-Neg[!is.na(ShortName)]
PosAA<-Pos[!is.na(ShortName)]
NegAA$ShortName
PosAA$ShortName
F_Bland(x=PosAA$rt,y=NegAA$rt,Label=NegAA$ShortName,Name="AA RT Difference (Seconds) Positive-Negative")

#### No use Cluster based on SMILES####
Stat<-data.table(read.xlsx('matched_compound.xlsx','stat'))
names(Stat)
setorder(Stat,p.ajusted) # sort from small to big
Smi<-Stat[!is.na(smiles)]
Smi<-unique(Smi,by=c('rownames','mode'))
Smi<-unique(Smi,by='smiles')
Smi<-unique(Smi,by='name')
smiles<-Smi$smiles
# count(smiles)
# write.csv(DTSmi,'DTSmi.csv')
names(smiles)<-Smi$name
Convert<-convertFormat("SMI","SDF",paste(paste(smiles,names(smiles),sep="\t"),collapse="\n"))
Strsplit<-strsplit(Convert,"\n",fixed=TRUE)
sdfset<-read.SDFset(unlist(Strsplit))
apset <- sdf2ap(sdfset)
apset@ID <- sdfid(sdfset)
dummy <- cmp.cluster(db=apset, cutoff=0, save.distances="distmat.rda", quiet=TRUE)
load("distmat.rda")
# https://stats.stackexchange.com/questions/4062/how-to-plot-a-fan-polar-dendrogram-in-r
hc <- hclust(as.dist(distmat), method="ave")
hc[["labels"]] <- cid(apset) # Assign correct item labels
dend <- as.dendrogram(hc) # create a dendrogram
# cols <- c("#009000", "#FF033E", "#CB410B", "#3B444B", "#007FFF")
dend <- dend %>%  # modify the dendrogram to have some colors in the branches and labels
  set("labels_cex", 0.5) %>%
  set("branches_lwd", 0.5) %>%
  color_branches(k = 13) %>%
  color_labels(k = 13)
svg(file=paste("ClusterSmiles.svg",seq=""),width = 15,height = 15)
circlize_dendrogram(dend, labels_track_height=0.5,dend_track_height = 0.3)
dev.off()
## Add the cluter color into the table
colors <- labels_colors(dend)
Color<-data.table('name'=names(colors),'Color'=as.character(colors))
SmiColor<-left_join(Smi,Color,by=('name')) %>% data.table(.)

#### No use Heatmap ####
names(SmiColor)
Heatmap<-SmiColor[,c("H1.1","H1.2","H2.1","H2.2","H2.3","H3.1","H3.2","H3.3",
                     "N1.1","N1.2","N1.3","N2.1","N2.2","N2.3","N3.1","N3.2","N3.3",
                     'name','p.ajusted','log2.FC.','rownames','Color')]
# setnames(Heatmap,c('name','p.ajusted'),c('Name','pvalue'))
Heatmap[,1:17][]=lapply(Heatmap[,1:17],type.convert,as.is=TRUE) # translate from character into number
Heatmap[,Mean:=apply(Heatmap[,1:17],1,mean)]
# HeatmapSig<-Heatmap[(p.ajusted<0.05)&(abs(log2.FC.)>1),]
# nrow(HeatmapSig[log2.FC.>1]) # increased
# nrow(HeatmapSig[log2.FC.<(-1)]) # decreased
Heatmap<-Heatmap[!is.na(Mean)]
setorder(Heatmap,p.ajusted)
sample_number<-c("N1.1","N1.2","N1.3","N2.1","N2.2","N2.3","N3.1","N3.2","N3.3",
                 "H1.1","H1.2","H2.1","H2.2","H2.3","H3.1","H3.2","H3.3")
F_PlotHeatMap<-function(DT,Cutoff,group,Name,sample_number,Size){
  ## Translate into numbers
  DT<-DT[1:Cutoff]
  df<-DT[,..sample_number]
  df[]=lapply(df,type.convert,as.is=TRUE)
  colnames(df)<- sample_number
  rownames(df)<-DT$name
  ## Define annotation for top
  annotation_col = data.frame(group = factor(group))
  rownames(annotation_col) = sample_number
  group_colors <- c("cyan1","coral1")
  names(group_colors) <- unique(group)
  # ## Define annotation for left
  # annotation_row = data.frame(GO = factor(DT$GO_names))
  # rownames(annotation_row) = DT$Name
  ## Define annotation for right, the text
  right_annotation<-rowAnnotation(width=unit(2, "cm"),
                                  text = anno_text(gt_render(rownames(df),  align_widths = TRUE),
                                                   gp = gpar(fontsize = 8,col=DT$Color)))
  ## Define heatmap colors
  ann_colors = list(group = group_colors) #
  bk <- seq(-3,3,by=1) # Set fixed indicater
  ha1 = HeatmapAnnotation(foo1 = runif(10), bar1 = sample(c("f", "m"), 10, replace = TRUE))
  ## Perform the Hierarchical  cluster analysis
  p<-pheatmap(df,
              na_col = "white",
              treeheight_row = 10,treeheight_col = 10,
              annotation_col = annotation_col,
              # annotation_row = annotation_row,
              annotation_colors = ann_colors,
              clustering_method = "ward.D",
              border_color = "grey",
              main = Name,
              color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(bk)),
              legend_breaks=seq(-3,3,1),
              cellheight = 8, cellwidth = 8,
              # heatmap_legend_param = list(direction = "horizontal"),
              scale = "row",
              right_annotation=right_annotation)
  # p
  svg(file=paste(Name,".svg",seq=""),width = Size,height = Size)
  draw(p,heatmap_legend_side = "left",annotation_legend_side = "left",merge_legend=TRUE)
  dev.off()
}

F_PlotHeatMap(DT=Heatmap,Cutoff=50,group=c(rep("Normal",9),rep("HYP",8)),
              Name="Heatmap of differential compounds",sample_number=sample_number,Size=15)


#### No use MS2 spectra matching based on the Spectra package ####
# https://github.com/jorainer/SpectraTutorials/blob/main/vignettes/Spectra-matching-with-MetaboAnnotation.Rmd
# https://github.com/jorainer/SpectraTutorials/blob/main/vignettes/analyzing-MS-data-from-different-sources-with-Spectra.Rmd
dda_2020 <- xcms::featureSpectra(xdata_MixAA_Pos,return.type = "Spectra")
dda_2022 <- xcms::featureSpectra(xdata_MixAA_Pos22,return.type = "Spectra")

con <- dbConnect(SQLite(), paste0(Path,"Reference/MassbankSql-2021-03.db"))
mbank <- Spectra(con, source = MsBackendMassbankSql())
print(object.size(mbank), units = "MB")

names(Mbank)
head(Mbank)
nrow(Mbank) # 86576 spectra
unique(Mbank, by='inchikey') # 16538 compounds
MbankUni <- unique(Mbank, by='smiles') # 19353 compounds

## Normalize intensities for all MassBank spectra
norm_int <- function(x, ...) {
  #' Define a function to *normalize* the intensities
  maxint <- max(x[, "intensity"], na.rm = TRUE)
  x[, "intensity"] <- 100 * x[, "intensity"] / maxint
  x
}
mbank <- addProcessing(mbank, norm_int)
intensity(mbank[1])
spectraVariables(mbank)
Mbank <- as.data.table(spectraData(mbank))
Mbank[is.na(smiles)] # no entry
low_int <- function(x, ...) {
  x > max(x, na.rm = TRUE) * 0.0001
}

label_fun <- function(x) {
  #' Specifying a function to draw peak labels
  ints <- unlist(intensity(x))
  mzs <- format(unlist(mz(x)), digits = 4)
  mzs[ints < 5] <- ""
  mzs
}


F_CalPMMT <- function(mz, MRP = 17500, RefMZ = 200) {
  A <- 1 / (MRP * (RefMZ^0.5))
  B <- A / 2.35482
  Da <- B * mz^1.5 + mz / 1000000
  Ppm <- 1000000*Da/mz
  print(paste0('Da ', Da))
  print(paste0('Ppm ', Ppm))
}
F_CalPMMT(mz=200,MRP = 7500, RefMZ = 400)

F_compareSpectra <- function(dda,lib,Cell){
  ## Load the volcano results
  DT <- data.table(read.csv(paste0(Path,Cell,'Phe_Volcano.csv')))
  # names(DT)
  sps_all  <- dda[dda$feature_id %in% unique(DT$Sample)] # only use the dda that in the DT
  ## Process dda from the experimental results
  sps_all  <- filterIntensity(sps_all, intensity = low_int)   #' Remove peaks with an intensity below 5% or the spectra's BPC
  sps_all <- sps_all[lengths(sps_all) > 1]   #' Remove spectra with a single peak
  # spectraVariables(sps_all)   #' List all available spectra variables (attributes)
  # head(sps_all$rtime) #' Access the spectras' retention time
  # print(object.size(sps_all), units = "MB")
  sps_all <- addProcessing(sps_all, norm_int) # #' Apply the normalization to the data
  # intensity(sps_all)[[1]] #' Get the intensities after normalization
  ## Build the consensus spectrum of experimental spectra
  ex_spectrum <- Spectra::combineSpectra(sps_all,  ppm = 40, peaks = "intersect", minProp = 0.8,
                                         intensityFun = median, mzFun = median,f = sps_all$feature_id)
  ## Library search in mass bank whole
  prm <- CompareSpectraParam(ppm = 40, requirePrecursor = TRUE,
                             THRESHFUN = function(x) which(x > 0.001))
  mtch <- matchSpectra(ex_spectrum, lib, param = prm)
  mtch <- mtch[whichQuery(mtch)]
  ## Export the related the candidates list
  head(spectraData(mtch,spectraVariables(mtch)))
  res <- spectraData(mtch, c('feature_id','precursorMz','collisionEnergy',
                             'target_pubchem','target_precursorMz','target_spectrum_id',
                             'target_spectrum_name','target_publication',
                             'target_compound_id','target_adduct','target_ionization',
                             'target_fragmentation_mode','target_collision_energy_text',
                             'target_instrument','target_instrument_type',
                             'target_formula','target_exactmass','target_smiles',
                             'target_inchikey','target_cas','target_pubchem',
                             'target_compound_name','score'))
  res <- as.data.table(res) %>% setorder(.,-score)
  ## Exclude the results from different instrument than Orbitrap
  unique(res$feature_id)
  unique(res,by='target_instrument') %>% .[,c('target_instrument','target_instrument_type')]
  unique(res$target_instrument_type)
  res <- res[target_instrument %in% c('Q Exactive Orbitrap (Thermo Scientific)',
                                      'LTQ Orbitrap Velos Thermo Scientific',
                                      'LTQ Orbitrap XL, Thermo Scientfic; HP-1100 HPLC, Agilent',
                                      'Q Exactive Plus Orbitrap Thermo Scientific',
                                      'Q Exactive Orbitrap Thermo Scientific',
                                      'LTQ Orbitrap XL Thermo Scientific',
                                      'Orbitrap Classic, Thermo Scientific',
                                      'LTQ Orbitrap XL, Thermo Scientfic',
                                      'Q-Exactive Orbitrap Thermo Scientific'
  )]
  # write.csv(res,'res.csv')
  res$target_pubchem <- sub(res$target_pubchem, pattern = "CID:",replacement = "",perl = TRUE)
  res$target_adduct<-sub(res$target_adduct, pattern = "\\[",replacement = "",perl = TRUE) %>%
    sub(., pattern = "\\].*",replacement = "",perl = TRUE)
  setnames(res,c('feature_id','target_adduct'),c('Sample','adduct'))
  unique(res$adduct)
  res <- res[adduct=='M+H']
  F <- unique(res,by=c('Sample'))
  Iden <- left_join(DT,F,by=c('Sample'))
  Colu <- c('score','Sample','mz','precursorMz','rt','target_compound_name')
  setcolorder(Iden, c(Colu,colnames(Iden)[!(colnames(Iden) %in% Colu)]))
  Iden$score[is.na(Iden$score)] <-0
  setorder(Iden, -score)
  write_xlsx(Iden,paste0(Path,Cell,"_MassbankResults.xlsx"))
  ## Unique the identification
  # IdenUni[score>0]
  # IdenUni <- Iden[!is.na(Iden$name)] # exclude the results missing names
  # names(Iden)
  # IdenUni <- unique(IdenUni, by=c('name')) %>% unique(.,by="Sample")
  # write.xlsx(IdenUni,paste0(Path,Cell,"_Results.xlsx"),sheetName="IdenUni",row.names=FALSE,append = T)
  ## Export the head to tail figure
  # for(id in c(1:nrow(F))){
  #   # id <- 1
  #   # print(id)
  #   lib_sel <- lib[lib$spectrum_id==F[id,]$target_spectrum_id]
  #   ex_sel <- ex_spectrum[ex_spectrum$feature_id==F[id,]$Sample]
  #   MirName <- paste(Cell, F[id,]$Sample, 'mz',round(F[id,]$precursorMz,digits = 2),
  #                    F[id,]$target_compound_name,'id',F[id,]$target_spectrum_id,
  #                    'dp=',round(F[id,]$score,digits = 2))
  #   MirName <- gsub(MirName, pattern = ":",replacement = "_",perl = TRUE)
  #   # MirName <- '391.28 3?-Hydroxy-12 Ketolithocholic Acid'
  #   MirName <- gsub(MirName, pattern = "\\?",replacement = "_",perl = TRUE)
  #   MirName <- gsub(MirName, pattern = "/",replacement = "_",perl = TRUE)
  #   png(paste0(MirName,'.png'),res=200,width = 2000,height = 2000)
  #   plotSpectraMirror(ex_sel[1], lib_sel ,main=MirName,
  #                     ppm = 20, labels = label_fun, labelPos = 2,
  #                     labelOffset = 0.2, labelSrt = -30)
  #   dev.off()
  # }
  return(Iden)
}
setwd(paste0(Path,'2020'))
Iden_2020 <- F_compareSpectra(dda = dda_2020, lib=mbank, Cell='2020')

# setwd(paste0(Path,'2022'))
# Iden_2022 <- F_compareSpectra(dda = dda_2022, lib=mbank, Cell='2022')


