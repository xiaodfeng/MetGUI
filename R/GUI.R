library(readr)
#detach("package:gWidgetstcltk")
#detach("package:gWidgets2tcltk")
#detach("package:gWidgets")
library(gWidgets2)
library(cairoDevice)
library(gWidgets2RGtk2)
require(gWidgets2RGtk2)
options(guiToolkit="RGtk2")
library(mzR)
library(Rcpp)
library(multtest)
library(MSnbase)
library(xcms) #version 3.3.3 
library(CAMERA)
library(RColorBrewer)
library(magrittr)
Sys.setenv(JAVA_HOME='c:/Program Files/Java/jre1.8.0_191')
library(rJava)
library(metfRag)
library(data.table)
library(stringr)
library(dplyr)
library (plyr)
library(ggplot2)
library(cowplot)
library(VennDiagram)
library(readr)
library(parallel)
library(MetaboAnalystR)
library(faahKO) # example dataset
library(later)
#library(png)
library(imager)
setwd("D:/github/GUI")
#register(SerialParam()) #disable paralell
BPPARAM  = SnowParam(detectCores()-3, progressbar = TRUE)
binSizeValue=fwhmValue=maxValue=snthreshValue=stepsValue=mzdiffValue=ppmValue=peakwidthMinValue=peakwidthMaxValue=snthreshCentWaveValue=noiseValue=binSizeRTValue=bwValue=binSizeGroupingValue=ppmIsotopeValue=ppmWindowValue=rtWindowValue=DatabaseValue=DatabaseSearchRelativeMassDeviationValue=FragmentPeakMatchAbsoluteMassDeviationValue=FragmentPeakMatchRelativeMassDeviationValue=ppmHomeBuiltValue=-1
options(warn=-1,show.error.messages=TRUE)
RangeRT=c(35.07705,1258.41516)
RangeMZ=c(250.1777,1746.5595)
Database<<-data.table(read.csv("d:/software/Database/lipidmaps2017Short.csv"))
win <- gwindow("Lipidomics identification workflow",expand=TRUE, width=1100,height=900,visible = T)

nb <- gnotebook(cont = win,visible=FALSE)
####ImportData####
gp_ImportData <- ggroup(cont = nb, label = "ImportData", horizontal = T,expand=TRUE)
#gp_ImportData$set_borderwidth(10L)
lo_ImportData <- glayout(cont = gp_ImportData, horizontal = F)
##Set Dataset frame in ImportData tab
lo_ImportData[1,1] <- gf_Dataset<-gframe("Dataset",horizontal = F,expand=TRUE) 
lo_Dataset <- glayout(cont = gf_Dataset, horizontal = F)
lo_Dataset[1,1] <- gbutton("Select dataset folder", handler = function(h,...) {
  my_path =choose.dir(caption = "Select folder")
  if(length(my_path)==0) {galert('Please Select the folder to be Processed',title = "File Selection Problems",delay = 1)}
  else {
    cdffiles <<- list.files(my_path, recursive = TRUE, full.names = TRUE)
    print(paste("You selected",cdffiles))
    name<<-list.files(my_path, recursive = TRUE, full.names = FALSE)
    #show the loaded files below
    ed_group<<-list()
    for (i in 1:length(cdffiles)){
      lo_Dataset[i+2,1:10]<-name[i]
      lo_Dataset[i+2,11]<-ed_group[i]<<-gedit(text = "wt", width = 4)
    }
  }
})
lo_Dataset[1,2] <- gbutton("?",icon="?", handler = function(h,...) {
  gmessage("Select the folder for all the raw data to be processed",title = "help")
})

lo_Dataset[2,1] <- gbutton("Set Group", handler = function(h,...) {
  smp_grp <- rep("wt", length(cdffiles))
  for (i in 1:length(cdffiles)){
    smp_grp[i] <- svalue(ed_group[[i]])
  }
  Group<<-smp_grp
  print("Group information already set")
})

lo_Dataset[2,2] <- gbutton("?",icon="?", handler = function(h,...) {
  gmessage("Set Group information of all the loaded files",title = "help")
})


###Set Experimental design frame in ImportData tab
lo_ImportData[1,2] <- gf_Design<-gframe("Experimental design", horizontal = F,expand=TRUE) 
lo_Design <- glayout(cont = gf_Design, horizontal = F)
lo_Design[1,1] <- gbutton("Import data", handler = function(h,...) {
  # Define a data.frame with sample descriptions
  pd <- data.frame(file = cdffiles, sample_group = Group)
  # Read the files
  data <- readMSData(cdffiles, pd = new("NAnnotatedDataFrame", pd),mode = "onDisk")
  Raw_data <<- data
  Raw_data  <<- Raw_data [grep(fData(Raw_data)$msLevel, pattern = "1")]
  group_colors <<- brewer.pal(12,"Set3") #set colours
  names(group_colors) <- name
  group_colors<<-group_colors
  print("Data imported successfully!")
})
lo_Design[1,2] <-gbutton("?",handler = function(h,...) {gmessage("Import the data into an object",title = "help")}) 
###Set TIC in Design frame
lo_Design[2,1] <-  gbutton("Plot TIC", handler = function(h,...) {
  F_Delete_Ggraphics()
  plot(chromatogram(Raw_data),col = group_colors[name])
  legend("topright",legend=name,col=group_colors[name], lty=1, cex=0.8) #,box.lty=0
  print("Plot Total ion current Done!")
})
lo_Design[2,2] <-gbutton("?",handler = function(h,...) {gmessage("Plot the total ion current of the raw data",title = "help")}) 
###Select export folder
lo_Design[3,1] <- gbutton("Select export folder", handler = function(h,...) {
  path =choose.dir(caption = "Select export folder")
  if(length(path)==0) {galert('Please Select the export folder',title = "File Selection Problems",delay = 1)}
  else {
    setwd(path)
    print("Export folder selected successfully!")
  }
})
lo_Design[3,2] <- gbutton("?",icon="?", handler = function(h,...) {gmessage("Select the folder to save all the exported results",title = "help")})
lo_Design[4,1] <- gbutton("export", handler = function(h,...) {
  ##extract MS1 according to function
  #fun_1 <- grep(fData(Raw_data)$msLevel, pattern = "1")
  #raw_data_all[fun_1]
  Fdata1<-data.table(fData(Raw_data))
  write.csv(Fdata1,"ImportData_Fdata1.csv")
  print("Results exported successfully!")
})
lo_Design[4,2] <- gbutton("?",icon="?", handler = function(h,...) {gmessage("Export base peak information from the raw dataset",title = "help")})


##Set Parameters frame in ImportData tab
lo_ImportData[2,1] <- gf_Parameters<-gframe("Parameters",horizontal = F,expand=TRUE) 
lo_Parameters <- glayout(cont = gf_Parameters, horizontal = F)
lo_Parameters[1,1] <- gbutton("Select Parameter File", handler = function(h,...) {
  Parameters =choose.files(caption = "Select Parameter File")
  if(length(Parameters)==0) {galert('Please Select the Parameter File to be Processed',title = "File Selection Problems",delay = 1)}
  else {
    Parameter<-read.csv(Parameters)
    Parameter<<-as.data.table(Parameter)
    print("Parameters loaded successfully!")
  }
})
lo_Parameters[1,2] <- gbutton("?",icon="?", handler = function(h,...) {
  gmessage("Select the parameter file",title = "help")
})
lo_Parameters[2,1] <- gbutton("ExportParameters", handler = function(h,...) {
  if(gconfirm(c("Export Parameters Now?", "Make sure you have already finished the identification before exporting parameters")))
  {ParametersExported<-data.table()
  ParametersExported<-Parameter
  UserValues=c(binSizeValue,fwhmValue,maxValue,snthreshValue,stepsValue,mzdiffValue,ppmValue,peakwidthMinValue,peakwidthMaxValue,snthreshCentWaveValue,noiseValue,binSizeRTValue,bwValue,binSizeGroupingValue,ppmIsotopeValue,ppmWindowValue,rtWindowValue,DatabaseValue,DatabaseSearchRelativeMassDeviationValue,FragmentPeakMatchAbsoluteMassDeviationValue,FragmentPeakMatchRelativeMassDeviationValue,ppmHomeBuiltValue)
  ParametersExported[,Values:=UserValues]
  write.csv(ParametersExported,"Parameters_ParametersExported.csv")
  print("Parameters exported successfully!")}
  else {print("Parameters are not exported!")}
})

####PeakDetection####
gp_PeakDetection <- ggroup(cont = nb, label = "PeakDetection", horizontal = T,expand=TRUE)
##Set method sub tab to select different methods
nb_Method<-gnotebook(cont = gp_PeakDetection,horizontal = F, tab.pos=2)
gp_MatchedFilter <- ggroup(cont = nb_Method, label = "MatchedFilter", horizontal = T,expand=TRUE)
gp_CentWave <- ggroup(cont = nb_Method, label = "CentWave", horizontal = T,expand=TRUE)
##set MatchedFilter parameters on gp_MatchedFilter
lo_MatchedFilter <- glayout(cont = gp_MatchedFilter, horizontal = T)
lo_MatchedFilter[1,1] <-"binSize"
lo_MatchedFilter[1,2]<- ed_binSize<-gedit(text = "0.04", width = 4)
lo_MatchedFilter[1,3]<-gbutton("?", handler = function(h,...) { gmessage("binSize. Default 0.1. Specifying the width of the bins/slices in m/z dimension.",title = "help")})

lo_MatchedFilter[2,1] <-glabel("fwhm")
lo_MatchedFilter[2,2]<-ed_fwhm<-gedit(text = "5", width = 4)
lo_MatchedFilter[2,3]<-gbutton("?", handler = function(h,...) { gmessage("fwhm. Default 30. specifying the full width at half maximum of matched filtration gaussian model peak, in seconds
                                                                         ",title = "help")})

lo_MatchedFilter[3,1] <-glabel("max")
lo_MatchedFilter[3,2]<-ed_max<-gedit(text = "10", width = 4)
lo_MatchedFilter[3,3]<-gbutton("?", handler = function(h,...) { gmessage("max. Default 5. representing the maximum number of peaks that are expected/will be identified per slice.
                                                                         ",title = "help")})

lo_MatchedFilter[4,1] <-glabel("snthresh")
lo_MatchedFilter[4,2]<-ed_snthresh<-gedit(text = "6", width = 4)
lo_MatchedFilter[4,3]<-gbutton("?", handler = function(h,...) { gmessage("snthresh. Default 10. defining the signal to noise cutoff to be used in the chromatographic peak detection step.
                                                                         ",title = "help")})

lo_MatchedFilter[5,1] <-glabel("steps")
lo_MatchedFilter[5,2]<-ed_steps<-gedit(text = "4", width = 4)
lo_MatchedFilter[5,3]<-gbutton("?", handler = function(h,...) { gmessage("steps. Default 2.  Defining the number of bins to be merged before filtration (i.e. the number of neighboring bins that will be joined to the slice in which filtration and peak detection will be performed).
                                                                         ",title = "help")})

lo_MatchedFilter[6,1] <-glabel("mzdiff")
lo_MatchedFilter[6,2]<-ed_mzdiff<-gedit(text = "0.008", width = 8)
lo_MatchedFilter[6,3]<-gbutton("?", handler = function(h,...) { gmessage("mzdiff. Default mzdiff = 0.8 - binSize * steps=0.6. Defining the minimum difference in m/z for peaks with overlapping retention times
                                                                         ",title = "help")})

lo_MatchedFilter[7,1]<-gbutton("Run", icon="help", handler = function(h,...) {
  binSizeValue<<-as.numeric(as.character(Parameter[Parameters=="binSize",Values]))
  fwhmValue<<-as.numeric(as.character(Parameter[Parameters=="fwhm",Values]))
  maxValue<<-as.numeric(as.character(Parameter[Parameters=="max",Values]))
  snthreshValue<<-as.numeric(as.character(Parameter[Parameters=="snthresh",Values]))
  stepsValue<<-as.numeric(as.character(Parameter[Parameters=="steps",Values]))
  mzdiffValue<<-as.numeric(as.character(Parameter[Parameters=="mzdiff",Values]))
  if (binSizeValue==-1) {binSizeValue<-as.numeric(svalue(ed_binSize))}
  if (fwhmValue==-1) {fwhmValue<-as.numeric(svalue(ed_fwhm))}
  if (maxValue==-1) {maxValue<-as.numeric(svalue(ed_max))}
  if (snthreshValue==-1) {snthreshValue<-as.numeric(svalue(ed_snthresh))}
  if (stepsValue==-1) {stepsValue<-as.numeric(svalue(ed_steps))}
  if (mzdiffValue==-1) {mzdiffValue<-as.numeric(svalue(ed_mzdiff))}  
  mfp <- MatchedFilterParam(
    binSize=as.numeric(binSizeValue),                #specifying the width of the bins/slices in m/z dimension.
    fwhm=as.numeric(fwhmValue),       #specifying the full width at half maximum of matched filtration gaussian model peak. Only used to calculate the actual sigma, see below.
    max=as.numeric(maxValue), # representing the maximum number of peaks that are expected/will be identified per slice.
    snthresh = as.numeric(snthreshValue), #defining the signal to noise cutoff to be used in the chromatographic peak detection step.
    steps=as.numeric(stepsValue),    # defining the number of bins to be merged before filtration (i.e. the number of neighboring bins that will be joined to the slice in which filtration and peak detection will be performed).
    mzdiff=as.numeric(mzdiffValue)
  )
  #xdata1<<- findChromPeaks(Raw_data, param = mfp, BPPARAM = BPPARAM)
  xdata1<<- findChromPeaks(Raw_data, param = MatchedFilterParam(), BPPARAM = BPPARAM)
  print("MatchedFilter Peak detection Done!")
})
lo_MatchedFilter[7,2]<- gbutton("export",handler = function(h,...) {
  write.csv(xdata1@msFeatureData[["chromPeaks"]],"PeakDetection_xdata1@msFeatureData_chromPeaks.csv")
  print("Results exported successfully!")
})
##set CentWave parameters on gp_CentWave
lo_CentWave <- glayout(cont = gp_CentWave, horizontal = T)
lo_CentWave[1,1] <-glabel("ppm")
lo_CentWave[1,2]<-ed_ppm<-gedit(text = "7", width = 4)
lo_CentWave[1,3]<-gbutton("?", handler = function(h,...) { gmessage("ppm. Default 25. maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition. Guidelines are 5–15 ppm for Orbitrap data, ~5 ppm for lock mass quadrupole time of flight (QTOF) data and 10–20 ppm for other QTOF instruments.
                                                                    ",title = "help")})
lo_CentWave[2,1] <-glabel("peakwidthMin")
lo_CentWave[2,2]<-ed_peakwidthMin<-gedit(text = "2.4", width = 4)
lo_CentWave[2,3]<-gbutton("?", handler = function(h,...) { gmessage("Peakwidth. Default 20-50. minimum and maximum chromatographic peak width detectable (seconds). This depends mainly on the type of chromatographic separation performed. For standard reverse-phase separations, a general guideline is 20–60 s, whereas for hydrophilic interaction liquid chromatography (HILIC), in which run times tend to be longer with broader peaks, we recommend 25–90 s. When running with UPLC, these values should drop markedly because of shorter run times and higher resolution, with suggested starting values between 2 and 5 s to a maximum of 30 s. If in doubt of values, check the raw chromatographic run for peak widths of some common compounds from each end of the trace. These values are not hard cutoffs and may be detected slightly out of this range, depending on the quality of the peak data
                                                                    ",title = "help")})

lo_CentWave[3,1] <-glabel("peakwidthMax")
lo_CentWave[3,2]<-ed_peakwidthMax<-gedit(text = "30", width = 4)

lo_CentWave[4,1] <-glabel("snthreshCentWave")
lo_CentWave[4,2]<-ed_snthreshCentWave<-gedit(text = "10", width = 4)
lo_CentWave[4,3]<-gbutton("?", handler = function(h,...) { gmessage("snthresh. Default 10,Signal/Noise ratio cutoff: ([maximum peak intensity] - [estimated baseline value]) / standard deviation of local chromatographic noise. Higher than default value to reduce false positive results.
                                                                    ",title = "help")})

lo_CentWave[5,1] <-glabel("noise")
lo_CentWave[5,2]<-ed_noise<-gedit(text = "100000", width = 4)
lo_CentWave[5,3]<-gbutton("?", handler = function(h,...) { gmessage("noise. Default 0. centroids with intensity < noise are omitted from ROI detection 
                                                                    ",title = "help")})
lo_CentWave[6,1]<- gbutton("Run",  handler = function(h,...) {
  
  ppmValue<<-as.numeric(as.character(Parameter[Parameters=="ppm",Values]))
  if (ppmValue==-1) {ppmValue<<-as.numeric(svalue(ed_ppm))}
  
  peakwidthMinValue<<-as.numeric(as.character(Parameter[Parameters=="peakwidthMin",Values]))
  if (peakwidthMinValue==-1) {peakwidthMinValue<-as.numeric(svalue(ed_peakwidthMin))}
  
  peakwidthMaxValue<<-as.numeric(as.character(Parameter[Parameters=="peakwidthMax",Values]))
  if (peakwidthMaxValue==-1) {peakwidthMaxValue<-as.numeric(svalue(ed_peakwidthMax))}
  
  snthreshCentWaveValue<<-as.numeric(as.character(Parameter[Parameters=="snthreshCentWave",Values]))
  if (snthreshCentWaveValue==-1) {snthreshCentWaveValue<<-as.numeric(svalue(ed_snthreshCentWave))}
  noiseValue<<-as.numeric(as.character(Parameter[Parameters=="noise",Values]))
  if (noiseValue==-1) {noiseValue<-as.numeric(svalue(ed_noise))}
  cwp <- CentWaveParam(
    ppm=as.numeric(ppmValue) ,                #camera 30, maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition.
    peakwidth=c(as.numeric(peakwidthMinValue),as.numeric(peakwidthMaxValue)),       #camera (5,12) minimum and maximum chromatographic peak width detectable (seconds), also use c(30, 80)
    snthresh = as.numeric(snthreshCentWaveValue),           # camera/default use 6,Signal/Noise ratio cutoff: ([maximum peak intensity] - [estimated baseline value]) / standard deviation of local chromatographic noise
    noise=as.numeric(noiseValue)  # centroids with intensity < noise are omitted from ROI detection
  )
  xdata1 <<- findChromPeaks(Raw_data, param = cwp, BPPARAM = BPPARAM)
  #Raw_data2@featureData@data$msLevel<-1
  #xdata2 <<- findChromPeaks(Raw_data2, param = cwp)
  print("CentWave Peak detection Done!")
})
lo_CentWave[6,2]<- gbutton("export",handler = function(h,...) {
  write.csv(xdata1@msFeatureData[["chromPeaks"]],"PeakDetection_xdata1@msFeatureData_chromPeaks.csv")
  print("Results exported successfully!")
})

####RTAlignment####
gp_RTAlignment <- ggroup(cont = nb, label = "RTAlignment", expand=T)
lo_RTAlignment <- glayout(cont = gp_RTAlignment)
###Set Obiwarp RT alignment frame in Retention Time Alignment tab
lo_RTAlignment[1,1] <- gf_Obiwarp<-gframe("Obiwarp RT Alignment",expand=T)
lo_Obiwarp <- glayout(cont = gf_Obiwarp)
lo_Obiwarp[1,1] <- glabel("binSize")   ### set binSize
lo_Obiwarp[1,2] <- ed_binSizeRT <- gedit(text = "0.04", width = 4)
lo_Obiwarp[2,1] <-  gbutton("Run", handler = function(h,...) {
  binSizeRTValue<<-as.numeric(as.character(Parameter[Parameters=="binSizeRT",Values]))
  if (binSizeRTValue==-1) {binSizeRTValue<<-as.numeric(svalue(ed_binSizeRT))}
  xdata1 <<- adjustRtime(xdata1, param = ObiwarpParam())# binSize defining the bin size (in mz dimension) to be used for the profile matrix generation. 
  #xdata1 <<- adjustRtime(xdata1, param = ObiwarpParam(binSize =as.numeric(binSizeRTValue))) 
  print("Obiwarp RT Alignment Done!")
})
lo_Obiwarp[2,2] <- gbutton("?",icon="?", handler = function(h,...) {
  gmessage("Obiwarp RT Alignment method performs retention time adjustment using the Obiwarp method. It supports alignment of multiple samples by aligning each against a center sample. The alignment is performed directly on the profile-matrix. binSize defining the bin size (in mz dimension) to be used for the profile matrix generation.",title = "help")
})

### Set Retention Time Alignment effect gframe in Retention Time Alignment tab
lo_RTAlignment[1,2] <- gf_Effect<-gframe("Retention Time Alignment effect",expand=T)
lo_Effect <- glayout(cont = gf_Effect, horizontal = T)
lo_Effect[1,1] <-  gbutton("TIC", handler = function(h,...) {
  F_Delete_Ggraphics()
  par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
  ##before time alignment
  plot(chromatogram(Raw_data), col = group_colors[name]) #,main="Before retention time alignment"
  legend("topright",legend=name,col=group_colors[name], lty=1, cex=0.8) #,box.lty=0
  ##after time alignment
  plot(chromatogram(xdata1), col = group_colors[name]) #,main="After retention time alignment"
  legend("topright",legend=name,col=group_colors[name], lty=1, cex=0.8) #,box.lty=0
  print("Plot RT Alignment effect Done!")
})
lo_Effect[1,2] <- gbutton("?",icon="?", handler = function(h,...) {
  gmessage("Plot the TIC before and after retention time alignment",title = "help")
})

lo_Effect[2,1] <- glabel(paste("rtMin (>",round(RangeRT[1])," s)",sep=""))  
lo_Effect[2,2] <- rtMinAlign <- gedit(text = "500", width = 6)
lo_Effect[3,1] <- glabel(paste("rtMax (<",round(RangeRT[2])," s)",sep=""))  
lo_Effect[3,2] <- rtMaxAlign <- gedit(text = "530", width = 6)
lo_Effect[4,1] <- glabel(paste("mzMin (>",round(RangeMZ[1]),")",sep=""))  
lo_Effect[4,2] <-  mzMinAlign <- gedit(text = "755.5", width = 6)
lo_Effect[5,1] <- glabel(paste("mzMax (<",round(RangeMZ[2]),")",sep=""))  
lo_Effect[5,2] <-  mzMaxAlign <- gedit(text = "755.6", width = 6)


lo_Effect[6,1] <-  gbutton("Extracted ion chromatogram", handler = function(h,...) {
  F_Delete_Ggraphics()
  par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
  ##before time alignment
  Before<-chromatogram(Raw_data, aggregationFun = "sum",mz = c(as.numeric(svalue(mzMinAlign)),as.numeric(svalue(mzMaxAlign))), rt = c(as.numeric(svalue(rtMinAlign)),as.numeric(svalue(rtMaxAlign))))
  plot(Before, col = group_colors[name]) #,main="Before retention time alignment"
  legend("topright",legend=name,col=group_colors[name], lty=1, cex=0.8) #,box.lty=0
  ##after time alignment
  After<-chromatogram(xdata1, aggregationFun = "sum",mz = c(as.numeric(svalue(mzMinAlign)),as.numeric(svalue(mzMaxAlign))), rt = c(as.numeric(svalue(rtMinAlign)),as.numeric(svalue(rtMaxAlign))))
  plot(After, col = group_colors[name]) #,main="After retention time alignment"
  legend("topright",legend=name,col=group_colors[name], lty=1, cex=0.8) #,box.lty=0
  print("Plot Extracted ion chromatogramt Done!")
})
lo_Effect[6,2] <- gbutton("?",icon="?", handler = function(h,...) {
  gmessage("Plot the EIC before and after retention time alignment",title = "help")
})






##Grouping tab
gp_Grouping <- ggroup(cont = nb, label = "Grouping", expand = T)
lo_Grouping <- glayout(cont = gp_Grouping, horizontal = T)
## Set Density grouping gframe in Grouping tab
lo_Grouping[1,1] <- gf_Density<-gframe("Density Grouping",expand=T)
lo_Density <- glayout(cont = gf_Density, horizontal = T)
lo_Density[1,1] <- glabel("BandWidth")   ### set binSize
lo_Density[1,2] <- ed_bw <- gedit(text = "10", width = 4)
lo_Density[2,1] <- glabel("binSizeGrouping")   ### set binSize
lo_Density[2,2] <-ed_binSizeGrouping<- gedit(text = "0.04", width = 8)
lo_Density[3,1] <- glabel("minFraction")   
lo_Density[3,2] <-ed_minFraction<- gedit(text = "0.5", width = 8)

lo_Density[4,1] <-  gbutton("Run", handler = function(h,...) {
  bwValue<<-as.numeric(as.character(Parameter[Parameters=="bw",Values]))
  if (bwValue==-1) {bwValue<<-as.numeric(svalue(ed_bw))}
  binSizeGroupingValue<<-as.numeric(as.character(Parameter[Parameters=="binSizeGrouping",Values]))
  if (binSizeGroupingValue==-1) {binSizeGroupingValue<<-as.numeric(svalue(ed_binSizeGrouping))}
  minFractionValue<<-as.numeric(svalue(ed_minFraction))
  pdp <- PeakDensityParam(
    minFraction=minFractionValue,
    sampleGroups = Group,
    bw=as.numeric(bwValue),       
    binSize=as.numeric(binSizeGroupingValue) #camera uses 0.015,default 0.25.width of overlapping m/z slices to use for creating peak density chromatograms and grouping peaks across samples
  )
  xdata1 <<- groupChromPeaks(xdata1, param = pdp)
  #xdata1 <<- groupChromPeaks(xdata1,param = PeakDensityParam(sampleGroups = Group, bw = 10))
  Feature <<-xdata1@msFeatureData[["featureDefinitions"]] #for plot XIC and grouping effect
  RangeRT<<-range(Feature$rtmed)
  RangeMZ<<-range(Feature$mzmed)
  #xdata2 <<- groupChromPeaks(xdata2, param = pdp)
  print("Grouping Done!")
})
lo_Density[4,2] <- gbutton("?",icon="?", handler = function(h,...) {
  gmessage("
           groupChromPeaks-density method performs correspondence (chromatographic peak grouping) based on the density (distribution) of identified peaks along the retention time axis within slices of overlapping mz ranges. All peaks (from the same or from different samples) being close on the retention time axis are grouped into a feature (peak group).bw: The maximum expected RT deviation across samples.  
           ")})
lo_Density[5,1] <- gbutton("export", handler = function(h,...) { 
  write.csv(xdata1@msFeatureData[["featureDefinitions"]],"Grouping_xdata1_msFeatureData_featureDefinitions.csv") 
  #write.csv(xdata2@msFeatureData[["featureDefinitions"]],"Grouping_xdata2_msFeatureData_featureDefinitions.csv") 
  print("Results exported successfully!")
})
## Set visualization gframe in Grouping tab
lo_Grouping[1,2] <- gf_XIC<-gframe("Visualization",expand=T)
lo_XIC <- glayout(cont = gf_XIC, horizontal = T)
lo_XIC[1,1] <- glabel(paste("rtMin (>",round(RangeRT[1])," s)",sep=""))  
lo_XIC[1,2] <- rtMinGroup <<- gedit(text = "500", width = 6)
lo_XIC[2,1] <- glabel(paste("rtMax (<",round(RangeRT[2])," s)",sep=""))  
lo_XIC[2,2] <- rtMaxGroup <<- gedit(text = "530", width = 6)
lo_XIC[3,1] <- glabel(paste("mzMin (>",round(RangeMZ[1]),")",sep=""))  
lo_XIC[3,2] <-  mzMinGroup <<- gedit(text = "755.5", width = 6)
lo_XIC[4,1] <- glabel(paste("mzMax (<",round(RangeMZ[2]),")",sep=""))  
lo_XIC[4,2] <-  mzMaxGroup <<- gedit(text = "755.6", width = 6)
lo_XIC[5,1] <-  gbutton("Mass trace", horizontal=FALSE, handler = function(h,...) {
  F_Delete_Ggraphics()
  xdata1 %>%
    filterRt(rt = c(as.numeric(svalue(rtMinGroup)),as.numeric(svalue(rtMaxGroup)))) %>%
    filterMz(mz = c(as.numeric(svalue(mzMinGroup)),as.numeric(svalue(mzMaxGroup)))) %>%
    plot(type = "XIC")
  print("Plot Mass trace Done!")
})
lo_XIC[5,2]<- gbutton("?",  handler = function(h,...) {gmessage("Plot mass trace of the raw data. The input MZ and RT value should within the dataset range.",title = "help")})  
lo_XIC[6,1] <-  gbutton("Extracted ion chromatogram", horizontal=FALSE, handler = function(h,...) {
  F_Delete_Ggraphics()
  EIC<-chromatogram(xdata1, aggregationFun = "sum",mz = c(as.numeric(svalue(mzMinGroup)),as.numeric(svalue(mzMaxGroup))), rt = c(as.numeric(svalue(rtMinGroup)),as.numeric(svalue(rtMaxGroup))))#
  plot(EIC, col = group_colors[name])
  legend("topright",legend=name,col=group_colors[name], lty=1, cex=0.8) #,box.lty=0
  print("Plot Extracted ion chromatogramt Done!")
})
lo_XIC[6,2]<- gbutton("?",  handler = function(h,...) {gmessage("Plot Extracted ion chromatogram of the raw data. The input MZ and RT value should within the dataset range.",title = "help")})  

## Set GroupingEffect gframe in Grouping tab
lo_Grouping[1,3] <- gf_GroupingEffect<-gframe("Grouping Effect",expand=T)
lo_GroupingEffect <- glayout(cont = gf_GroupingEffect, horizontal = T)
lo_GroupingEffect[1,1] <- glabel("bwr")  
lo_GroupingEffect[1,2] <-   bwr <- gedit(text = "6", width = 4)
lo_GroupingEffect[2,1] <-  glabel("binsizer:")  
lo_GroupingEffect[2,2] <- binsizer <- gedit(text = "0.04", width = 4)
lo_GroupingEffect[3,1] <-   gbutton("Plot", handler = function(h,...) {
  F_Delete_Ggraphics()
  plotChromPeakDensity(xdata1,col = group_colors[name],mz=c(as.numeric(svalue(mzMinGroup)),as.numeric(svalue(mzMaxGroup))),rt = c(as.numeric(svalue(rtMinGroup)),as.numeric(svalue(rtMaxGroup))),
                       param = PeakDensityParam(sampleGroups = Group,bw=as.numeric(svalue(bwr)),binSize = as.numeric(svalue((binsizer)))))
  legend("topright",legend=name,col=group_colors[name], lty=1, cex=0.8) #,box.lty=0
  print("Plot Grouping effect Done!")
})
lo_GroupingEffect[3,2] <-  gbutton("?", handler = function(h,...) {gmessage("Plot Grouping effect to see the parameters bw and Binsize's effect on grouping",title = "help")})  

##Fill missing peaks tab
gp_FillMissing <- ggroup(cont = nb, label = "FillMissing", expand=T)
lo_FillMissing <- glayout(cont = gp_FillMissing)
##Set Fill missing peaks default gframe in Fill missing peaks tab
lo_FillMissing[1,1] <- gf_FillMissingDefault<-gframe("Fill Missing peaks by Default",expand=T)
lo_FillMissingDefault <- glayout(cont = gf_FillMissingDefault)
lo_FillMissingDefault[1,1] <- gbutton("Run", handler = function(h,...) {
  xdata1 <<- fillChromPeaks(xdata1)
  print("Fill missing peaks Done!")
})
lo_FillMissingDefault[1,2] <-gbutton("?",icon="?",   handler = function(h,...) {
  gmessage("Integrate signal in the mz-rt area of a feature (chromatographic peak group) for samples in which no chromatographic peak for this feature was identified and add it to the chromPeaks. Such peaks will have a value of 1 in the is_filled column of the chromPeaks matrix of the object.
           ")})
lo_FillMissingDefault[2,1] <- gbutton("export", handler = function(h,...) { 
  write.csv(Feature,"FillMissing_Feature.csv") 
  print("Export Features Done!")
})

##Set IsotopeFilteration gframe in Fill missing peaks tab
lo_FillMissing[1,2] <- gf_IsotopeFilteration<-gframe("IsotopeFilteration",expand=T)
lo_IsotopeFilteration <- glayout(cont = gf_IsotopeFilteration)
lo_IsotopeFilteration[1,1] <- glabel("ppmIsotope")  
lo_IsotopeFilteration[1,2] <- ed_ppmIsotope <- gedit(text = "10", width = 4)
lo_IsotopeFilteration[2,1] <-  gbutton("Run", handler = function(h,...) {
  ppmIsotopeValue<<-as.numeric(as.character(Parameter[Parameters=="ppmIsotope",Values]))
  if (ppmIsotopeValue==-1) {ppmIsotopeValue<<-as.numeric(svalue(ed_ppmIsotope))}
  xset <- as(xdata1, "xcmsSet")
  xsa <- xsAnnotate(xset)# create an CAMERA object
  xsa_I<-findIsotopes(xsa,ppm =ppmIsotopeValue)#,mzabs=0.009,filter = FALSE
  
  peaklist_xsa_I<- getPeaklist(xsa_I)
  DT_peaklist<-data.table(peaklist_xsa_I)#transfer into datatable format
  intensity<-apply(DT_peaklist[,10:11],1,mean,na.rm = T)# get the mean of intensity
  DT_peaklist[,Intensity:=intensity]#add a new column Intensity to ms1
  setorder(DT_peaklist,-Intensity)
  IsoFilter<-str_replace_all(DT_peaklist$isotopes,pattern = ".*M\\]\\+$", replacement = "")
  DT_peaklist[,IsoFilter:=IsoFilter]
  DT_peaklist<-DT_peaklist[IsoFilter=="",]#
  Feature<<-unique(DT_peaklist,by=c("mz","rt"))
  print("Isotope Filteration Done!")
})
lo_IsotopeFilteration[2,2] <- gbutton("?",icon="?", handler = function(h,...) {
  gmessage("Annotate isotope peaks using CAMERA. ppm is the ppm error for the search ")})
##Set Remove Duplicates gframe in Fill missing peaks tab
lo_FillMissing[2,1] <- gf_RemoveDup <-gframe("Remove Duplicates",expand=T)
lo_RemoveDup <- glayout(cont = gf_RemoveDup)
lo_RemoveDup[1,1] <-glabel("ppmWindow")  
lo_RemoveDup[1,2] <- ed_ppmWindow <- gedit(text = "5", width = 4)                                                           
lo_RemoveDup[2,1] <-glabel("rtWindow")  
lo_RemoveDup[2,2] <- ed_rtWindow <- gedit(text = "5", width = 4)      
lo_RemoveDup[3,1] <- gbutton("Run", handler = function(h,...) {
  
  ppmWindowValue<<-as.numeric(as.character(Parameter[Parameters=="ppmWindow",Values]))
  if (ppmWindowValue==-1) {ppmWindowValue<<-as.numeric(svalue(ed_ppmWindow))}
  
  rtWindowValue<<-as.numeric(as.character(Parameter[Parameters=="rtWindow",Values]))
  if (rtWindowValue==-1) {rtWindowValue<<-as.numeric(svalue(ed_rtWindow))}
  
  Peaks<- data.table(Feature)
  Peaks<-unique(Peaks,by=c("mz","rt"))
  #Extract the duplicates in Feature_MZmine
  PeaksCopy<-Peaks
  Duplist<-list()
  Duplicates<-data.table()
  for (id in seq(1,nrow(Peaks),1)){
    #for each line of Peaks, get the duplicates within +/- 5PPM and +/- 20 rt
    Duplist[[id]]<-Peaks[(((mz-PeaksCopy$mz[id])/mz)*10^6>-ppmWindowValue)&(((mz-PeaksCopy$mz[id])/mz)*10^6<ppmWindowValue)&((rt-PeaksCopy$rt[id])<rtWindowValue)&((rt-PeaksCopy$rt[id])>-rtWindowValue),]
    Duplist[[id]]<-Duplist[[id]][-c(1)] #get rid of the first line
    #Bind all the candidates together
    Duplicates<-rbind(Duplicates,Duplist[[id]])	
    Duplicates<-unique(Duplicates,by=c("mz","rt"))
  }
  # get rid of duplicates from Feature_MZmine
  Feature<-anti_join(Peaks,Duplicates,by=c("mz","rt"))
  Feature<<-data.table(Feature)
  print(nrow(Peaks))
  print(nrow(Feature))
  print("Remove duplicates Done!")
})
lo_RemoveDup[3,2] <- gbutton("?",  handler = function(h,...) {gmessage("Remove duplicates in selected mz and rt window. ppm is the mz window, rt is the retention time window
                                                                       ",title = "help")})   
##Set Scatter Plot gframe in Fill missing peaks tab
lo_FillMissing[2,2] <- gf_Scatter <-gframe("Scatter Plot",expand=T)
lo_Scatter <- glayout(cont = gf_Scatter)
lo_Scatter[1,1]<-glabel("Log transformation.Log")
lo_Scatter[1,2:15]<-cb_log<- gcombobox(c("2","10"),expand=T)
lo_Scatter[2,1] <- gbutton("Plot", handler = function(h,...) {
  ###figure
  if(as.numeric(svalue(cb_log))==2){
    g=ggplot(Feature,aes(mz,rt))+ geom_point(aes(color=log2(Intensity)))+labs(title="Scatter Plot of detected features")#+scale_color_manual(values =brewer.pal(9,"Set1")) #+geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)
  }
  else if (as.numeric(svalue(cb_log))==10){
    g=ggplot(Feature,aes(mz,rt))+ geom_point(aes(color=log10(Intensity)))+labs(title="Scatter Plot of detected features")#+scale_color_manual(values =brewer.pal(9,"Set1")) #+geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)
  }
  F_Delete_Ggraphics()
  plot(g)
  ggsave("Scatter Plot.png",width=8,height=8)
  print("Scatter Plot Done!")
})
lo_Scatter[2,2] <- gbutton("?", handler = function(h,...) {gmessage("Scatter Plot of detected features, the colour indicates the intensiy of the features",title = "help")})  


##FilterNormalization tab
gp_FilterNormalization <- ggroup(cont = nb, label = "FilterNormalization", expand=T)
lo_FilterNormalization <- glayout(cont = gp_FilterNormalization)
##Set Initialize gframe in FilterNormalization tab
lo_FilterNormalization[1,1] <- gf_Initialize<-gframe("Initialize R object",expand=T)
lo_Initialize<-glayout(cont=gf_Initialize)
lo_Initialize[1,1]<-gbutton("Initialize", handler = function(h,...) {
  ### Export the feature table in the MetaboAnalyst format. Parameter 'label', defines the group assignment of the samples.
  exportMetaboAnalyst(xdata1, file = "met_test1.csv", label = xdata1$sample_group)
  # First step is to create the mSet Object, specifying that the data to be uploaded
  mSet <<- InitDataObjects("pktable", "stat", FALSE) # is a peak table ("pktable") and that statistical analysis will be performed ("stat").
  # The second step is to read in the processed data (created above)
  mSet <<- Read.TextData(mSet, "met_test1.csv", "colu", "disc")
  FilterationCheck<<-0
})
lo_Initialize[1,2]<-gbutton("?", handler = function(h,...) {
  gmessage("To create the mSet Object, specifying that the data to be uploaded",title = "help")
})
##Set Filteration gframe in FilterNormalization tab
lo_FilterNormalization[1,2] <- gf_Filteration<-gframe("Filteration",expand=T)
lo_Filteration<-glayout(cont=gf_Filteration)
lo_Filteration[1,1]<-bt_SanityCheck<-gbutton("SanityCheck", handler = function(h,...) {
  mSet <<- SanityCheckData(mSet)
  FilterationCheck<<-FilterationCheck+1
  print("Sanity Check done!")
})
lo_Filteration[1,2]<-gbutton("?", handler = function(h,...) {
  gmessage("SanityCheckData is used for data processing, and performs a basic sanity check of the uploaded content, ensuring that the data is suitable for further analysis. The function will return a message if the data has successfully passed the check and is deemed suitable for further analysis. If it fails, the function will return a 0. The function will perform the check directly onto the mSet$dataSet object, and must be performed immediately after reading in data. The sanity check function evaluates the accuracy of sample and class labels, data structure, deals with non-numeric values, removes columns that are constant across all samples (variance = 0), and by default replaces missing values with half of the original minimal positive value in your dataset.Before data analysis, a data integrity check is performed to make sure that all the necessary information
           has been collected. The class labels must be present and contain only two classes. If samples are paired,
           the class label must be from -n/2 to -1 for one group, and 1 to n/2 for the other group (n is the sample
           number and must be an even number). Class labels with same absolute value are assumed to be pairs.
           Compound concentration or peak intensity values should all be non-negative numbers. By default, all
           missing values, zeros and negative values will be replaced by the half of the minimum positive value
           found within the data
           ",title = "help")})
lo_Filteration[2,1]<-bt_ReplaceMissingValue<-gbutton("ReplaceMissingValue", handler = function(h,...) {
  FilterationCheck<<-FilterationCheck+1
  mSet<<-ReplaceMin(mSet)
  print("Replace Missing Value done!")
})
lo_Filteration[2,2]<-gbutton("?", handler = function(h,...) {
  gmessage("Too many zeroes or missing values will cause difficulties for downstream analysis.The default method replaces all the missing and zero values
           with a small values (the half of the minimum positive values in the original data) assuming to be
           the detection limit. The assumption of this approach is that most missing values are caused by low
           abundance metabolites (i.e.below the detection limit).",title = "help")})
lo_Filteration[3,1]<-bt_FilterVariable<-gbutton("FilterVariable", handler = function(h,...) {
  mSet<<-FilterVariable(mSet, "iqr", "F", 25)
  print("Filter Variable done!")
  FilterationCheck<<-FilterationCheck+1
  if (FilterationCheck>=3){
    enabled(gb_Normalization)<-TRUE
    enabled(gb_PlotNormalization)<-TRUE
  }
  else {print("SanityCheck, ReplaceMissingValue and FilterVariable should be performed before normalization")}
  
})
lo_Filteration[3,2]<-gbutton("?", handler = function(h,...) {
  gmessage("The purpose of the data filtering is to identify and remove variables that are unlikely to be of use
           when modeling the data. No phenotype information are used in the filtering process, so the result
           can be used with any downstream analysis. This step can usually improves the results. Data filter is
           strongly recommended for datasets with large number of variables (> 250) datasets contain much noise
           (i.e.chemometrics data).The function applies a filtering method, ranks the variables within the dataset, and removes variables based on its rank. The final dataset should contain no more than than 5000 variables for effective computing.Here we use iqr:interquantile for filter option; F: do not use QC samples; 25 is the relative standard deviation cut-off
           ",title = "help")})

##Set Normalization gframe in FilterNormalization tab
#delete(lo_FilterNormalization, gf_Normalization)
lo_FilterNormalization[2,1] <- gf_Normalization<-gframe("Normalization",expand=T)
lo_Normalization<-glayout(cont=gf_Normalization)
lo_Normalization[1,1]<-glabel("normalization method")
lo_Normalization[1,2:15]<-cb_rowNorm<- gcombobox(c("QuantileNorm","ProbNormT","ProbNormF","CompNorm","SumNorm","MedianNorm","SpecNorm"),expand=T)
lo_Normalization[1,16]<-gbutton("?", handler = function(h,...) {
  gmessage("Select the option for row-wise normalization, QuantileNorm for Quantile Normalization, ProbNormT for Probabilistic Quotient Normalization without using a reference sample, ProbNormF for Probabilistic Quotient Normalization based on a reference sample, CompNorm for Normalization by a reference feature, SumNorm for Normalization to constant sum, MedianNorm for Normalization to sample median, and SpecNorm for Normalization by a sample-specific factor.
           ",title = "help")})

lo_Normalization[2,1]<-glabel("transform method")
lo_Normalization[2,2:15]<-cb_transNorm<- gcombobox(c("LogNorm","CrNorm"))
lo_Normalization[2,16]<-gbutton("?", handler = function(h,...) {
  gmessage("Select option to transform the data, LogNorm for Log Normalization, and CrNorm for Cubic Root Transformation.",title = "help")})

lo_Normalization[3,1]<-glabel("scaling method")
lo_Normalization[3,2:15]<-cb_scaleNorm<- gcombobox(c("MeanCenter","AutoNorm","ParetoNorm","RangeNorm"))
lo_Normalization[3,16]<-gbutton("?", handler = function(h,...) {
  gmessage("Select option for scaling the data, MeanCenter for Mean Centering(mean-centered only), AutoNorm for Autoscaling(mean-centered and divided by standard deviation of each variable), ParetoNorm for Pareto Scaling(mean-centered and divided by the square root of standard deviation of each
           variable), amd RangeNorm for Range Scaling(mean-centered and divided by the value range of each variable).",title = "help")})



lo_Normalization[4,1]<- gb_Normalization<- gbutton("Run",do.functions=FALSE, handler = function(h,...) {
  mSet<<-PreparePrenormData(mSet)
  mSet<<-Normalization(mSet, svalue(cb_rowNorm), svalue(cb_transNorm), svalue(cb_scaleNorm), ratio=FALSE, ratioNum=20)
  print("Normalization done!")
})
enabled(gb_Normalization) <- FALSE
lo_Normalization[4,2]<-gbutton("?", handler = function(h,...) {
  gmessage("The data is stored as a table with one sample per row and one variable (bin/peak/metabolite) per
           column. The normalization procedures implemented below are grouped into four categories. Sample
           specific normalization allows users to manually adjust concentrations based on biological inputs (i.e.
           volume, mass); row-wise normalization allows general-purpose adjustment for differences among samples;
           data transformation and scaling are two different approaches to make features more comparable.This function performs row-wise normalization, transformation, and scaling of your metabolomic data.
           ",title = "help")})

lo_Normalization[5,1]<- gb_PlotNormalization<- gbutton("Plot Normalization Density", handler = function(h,...) {
  mSet <<- PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
  mSet <<- PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
  F_Delete_Ggraphics()
  plot.new()
  rasterImage(load.image(mSet[["imgSet"]][["summary_norm"]]),0,0,0.5,1)
  rasterImage(load.image(mSet[["imgSet"]][["norm"]]),0.5,0,1,1)
  print("Normalization Density done!")
})
enabled(gb_PlotNormalization) <- FALSE
lo_Normalization[5,2]<-gbutton("?", handler = function(h,...) {
  gmessage("To see the Normalization Density Effect.
           ",title = "help")})

#Statistics tab
gp_Statistics <- ggroup(cont = nb, label = "Statistics", expand=T)
lo_Statistics <- glayout(cont = gp_Statistics)

##Set FoldChange gframe in Statistics tab
lo_Statistics[1,2] <- gf_FoldChange<-gframe("Fold Change",expand=T)
lo_FoldChange<-glayout(cont=gf_FoldChange)
lo_FoldChange[1,1]<-"Fold change threshold"
lo_FoldChange[1,2]<-ed_FCThreshold<-gedit(text = "2", width = 4)  
lo_FoldChange[2,1]<- gbutton("Run", handler = function(h,...) {
  mSet<<-FC.Anal.unpaired(mSet, as.numeric(svalue(ed_FCThreshold)), 0)
  mSet<<-FC.Anal.unpaired(mSet, 2, 0)
  mSet<<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)
  F_Delete_Ggraphics()
  plot.new()
  #plot(1:10,main="Fold change Summary",ty="n",xlab="",ylab="")
  rasterImage(load.image(mSet[["imgSet"]][["fc"]]),0,0,1,1)
  print("Fold change done!")
})
lo_FoldChange[2,2]<-gbutton("?", handler = function(h,...) {
  gmessage("Calculate the fold change between WT and KO groups, each point represents a peak. The default Fold change threshold is set as 2
           ",title = "help")})
##Set Volcano button in FoldChange gframe
lo_FoldChange[3,1:2] <- gbutton("Run Volcano Plot", handler = function(h,...) {
  mSet<<-Volcano.Anal(mSet, paired=F, nonpar=as.logical(svalue(cb_nonpar)), threshp=as.numeric(svalue(ed_threshp)),fcthresh=as.numeric(svalue(ed_FCThreshold)),cmpType=1,equal.var=TRUE, pval.type="raw")
  mSet<<-PlotVolcano(mSet, "Volcano", 0,"png", dpi=72, width=NA)
  F_Delete_Ggraphics()
  plot.new()
  #plot(1:10,main="volcano",ty="n",xlab="",ylab="")
  rasterImage(load.image(mSet[["imgSet"]][["volcano"]]),0,0,1,1)
  print("Volcano Plot done!")
})



##Set Ttest gframe in Statistics tab
lo_Statistics[1,1] <- gf_Ttest<-gframe("Ttest",expand=T)
lo_Ttest<-glayout(cont=gf_Ttest)

lo_Ttest[1,1]<-"non-parametric test"
lo_Ttest[1,2:8]<-cb_nonpar<-gcombobox(c("TRUE","FALSE"))

lo_Ttest[2,1]<-"adjusted p-value (FDR) cutoff"
lo_Ttest[2,2]<-ed_threshp<-gedit(text = "0.05", width = 4) 

lo_Ttest[3,1]<-"paired"
lo_Ttest[3,2:8]<-cb_paired<-gcombobox(c("TRUE","FALSE"))

lo_Ttest[4,1]<- gbutton("Run", handler = function(h,...) {
  mSet<<-Ttests.Anal(mSet, as.logical(svalue(cb_nonpar)), as.numeric(svalue(ed_threshp)), as.logical(svalue(cb_paired)), TRUE) # perform t-test analysis
  mSet<<-PlotTT(mSet, "tt_0_", "png", 72, width=NA)
  F_Delete_Ggraphics()
  plot.new()
  #plot(1:10,main="T test",ty="n",xlab="",ylab="")
  rasterImage(load.image(mSet[["imgSet"]][["tt"]]),0,0,1,1)
  print("T test done!")
})
lo_Ttest[4,2]<-gbutton("?", handler = function(h,...) {
  gmessage("perform t-test analysis. This univariate analyses provide a preliminary overview about
           features that are potentially significant in discriminating the conditions under study.For paired fold change analysis, the algorithm rst counts the total number of pairs with fold changes
           that are consistently above/below the specified FC threshold for each variable. A variable will be reported
           as significant if this number is above a given count threshold (default > 75% of pairs/variable)",title = "help")})

##Set PCA gframe in Statistics tab
lo_Statistics[2,1] <- gf_PCA<-gframe("Principal Component Analysis",expand=T)
lo_PCA<-glayout(cont=gf_PCA)
lo_PCA[1,1]<- gbutton("RunPCA", handler = function(h,...) {
  mSet<<-PCA.Anal(mSet)
  print("Principle Component Analysis done!")
})

lo_PCA[1,2]<-gbutton("?", handler = function(h,...) {
  gmessage("PCA is an unsupervised method aiming to find the directions that best explain the variance in a data
           set (X) without referring to class labels (Y). The data are summarized into much fewer variables called
           scores which are weighted average of the original variables. The weighting profiles are called loadings.
           The PCA analysis is performed using the prcomp package. The calculation is based on singular value
           decomposition.",title = "help")})
lo_PCA[2,1]<- gbutton("PlotPCASummary", handler = function(h,...) {
  mSet<<-PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
  mSet<<-PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
  mSet<<-PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
  mSet<<-PlotPCALoading(mSet, "pca_loading_0_", format="png", dpi=72, width=NA, 1,2)
  mSet<<-PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1,2)
  mSet<<-PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
  F_Delete_Ggraphics()
  plot.new()
  #plot(1:10,main="PCA Summary",ty="n",xlab="",ylab="")
  rasterImage(load.image(mSet[["imgSet"]][["pca.pair"]]),0,0.5,0.33,1)
  rasterImage(load.image(mSet[["imgSet"]][["pca.scree"]]),0.33,0.5,0.66,1)
  rasterImage(load.image(mSet[["imgSet"]][["pca.score2d"]]),0.66,0.5,0.99,1)
  rasterImage(load.image(mSet[["imgSet"]][["pca.loading"]]),0,0,0.33,0.5)
  rasterImage(load.image(mSet[["imgSet"]][["pca.biplot"]]),0.33,0,0.66,0.5)
  rasterImage(load.image(mSet[["imgSet"]][["pca.score3d"]]),0.66,0,0.99,0.5)
  print("PCA Summary done!")
})
lo_PCA[2,2]<-gbutton("?", handler = function(h,...) {
  gmessage("pairwise score plots providing an overview of
           the various seperation patterns among the most significant PCs;",title = "help")})

##Set PLSDA gframe in Statistics tab
lo_Statistics[2,2] <- gf_PLSDA<-gframe("Partial Least Squares - Discriminant Analysis",expand=T)
lo_PLSDA<-glayout(cont=gf_PLSDA)
lo_PLSDA[1,1]<- gbutton("RunPLSDA", handler = function(h,...) {
  mSet<<-PLSR.Anal(mSet, reg=TRUE) # Perform PLS-DA
  mSet<<-PlotPLSPairSummary(mSet, "pls_pair_0_", "png", 72, width=NA, 5)
  mSet<<-PlotPLS2DScore(mSet, "pls_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
  mSet<<-PlotPLS3DScoreImg(mSet, "pls_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
  mSet<<-PlotPLSLoading(mSet, "pls_loading_0_", "png", 72, width=NA, 1, 2)
  mSet<<-PLSDA.CV(mSet, "L",3, "Q2")
  mSet<<-PlotPLS.Classification(mSet, "pls_cv_0_", "png", 72, width=NA)
  mSet<<-PlotPLS.Imp(mSet, "pls_imp_0_", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)
  print("Partial Least Squares - Discriminant Analysis done!")
})
lo_PLSDA[1,2]<-gbutton("?", handler = function(h,...) {
  gmessage("PLS is a supervised method that uses multivariate regression techniques to extract via linear combination
           of original variables (X) the information that can predict the class membership (Y).To assess the significance of class discrimination, a permutation test was performed. In each permu-
           tation, a PLS-DA model was built between the data (X) and the permuted class labels (Y) using the
           optimal number of components determined by cross validation for the model based on the original class
           assignment. MetaboAnalyst supports two types of test statistics for measuring the class discrimination.
           The first one is based on prediction accuracy during training. The second one is separation distance
           based on the ratio of the between group sum of the squares and the within group sum of squares (B/W-
           ratio). If the observed test statistic is part of the distribution based on the permuted class assignments,
           the class discrimination cannot be considered significant from a statistical point of view.6.",title = "help")})

lo_PLSDA[2,1]<- gbutton("Plot PLS Summary", handler = function(h,...) {
  F_Delete_Ggraphics()
  plot.new()
  rasterImage(load.image(mSet[["imgSet"]][["pls.pair"]]),0,0.5,0.33,1)
  rasterImage(load.image(mSet[["imgSet"]][["pls.class"]]),0.33,0.5,0.66,1)
  rasterImage(load.image(mSet[["imgSet"]][["pls.score2d"]]),0.66,0.5,0.99,1)
  rasterImage(load.image(mSet[["imgSet"]][["pls.loading"]]),0,0,0.33,0.5)
  rasterImage(load.image(mSet[["imgSet"]][["pls.imp"]]),0.33,0,0.66,0.5)
  rasterImage(load.image(mSet[["imgSet"]][["pls.score3d"]]),0.66,0,0.99,0.5)
  print("PLS Summary done!")
})
lo_PLSDA[2,2]<-gbutton("?", handler = function(h,...) {
  gmessage("Plot PLS Summary",title = "help")})

##Set PreparePDFReport gframe in Statistics tab
lo_Statistics[2,3] <-gf_Report<- gframe("PreparePDFReport",expand=T)
lo_Report<<-glayout(cont=gf_Report)
lo_Report[1,1]<- "Report Name"
lo_Report[1,2]<- ed_ReportName<-gedit(text = "Xiaodong", width = 8) 
lo_Report[2,1]<- gbutton("PreparePDFReport", handler = function(h,...) {
  PreparePDFReport(mSet, svalue(ed_ReportName)) #create a summary report of the statistical analysis 
  print("PreparePDFReport done!")
})

#Identification tab
gp_Identification <- ggroup(cont = nb, label = "Identification", expand=T)
lo_Identification <- glayout(cont = gp_Identification)
##Set Preprocessing gframe in Identification tab
lo_Identification[1,1] <- gf_Preprocessing<-gframe("Preprocessing",expand=T)
lo_Preprocessing<-glayout(cont=gf_Preprocessing)
lo_Preprocessing[1,1]<-gbutton("Extract Significant Features",  handler = function(h,...) {
  Plsda.coef.mat<-read.csv("plsda_coef.csv")   #Get the features from plsda_coef files
  Plsda.coef.mat<-data.table(Plsda.coef.mat)
  setnames(Plsda.coef.mat,"X","rownames") # change names
  significant.plsda<-Plsda.coef.mat[Overall>10,]
  #Get the features from XCMS object
  featureDefinitions<-data.table(rownames=xdata1@msFeatureData[["featureDefinitions"]]@rownames,mz=xdata1@msFeatureData[["featureDefinitions"]]@listData[["mzmed"]],rt=xdata1@msFeatureData[["featureDefinitions"]]@listData[["rtmed"]])
  #Only keep featureDefinitions features that exist in significant.plsda with column rownames
  significant.features<-inner_join(significant.plsda,featureDefinitions,by="rownames")
  write.csv(significant.features,"Identification_significant.features.csv")
  Feature<<-significant.features
  print("Extract Significant Features Done!")
  print("Significant features stored in：Identification_significant.features.csv")
})
lo_Preprocessing[2,1]<-gbutton("Calculate Neutral Mass", handler = function(h,...) {
  MH <- data.table(Feature,Adducts_S="M+H")
  MH[,NeutralMass:=(mz-1.0073)]
  MNa <- data.table(Feature,Adducts_S="M+Na")
  MNa[,NeutralMass:=(mz-22.9892)]
  MK <- data.table(Feature,Adducts_S="M+K")
  MK[,NeutralMass:=(mz-38.9632)]
  MNH4 <- data.table(Feature,Adducts_S="M+NH4")
  MNH4[,NeutralMass:=(mz-18.0338)]
  M2NaH <- data.table(Feature,Adducts_S="M+2Na-H")
  M2NaH[,NeutralMass:=(mz-44.9712)]
  MMH <- data.table(Feature,Adducts_S="2M+H")
  MMH[,NeutralMass:=(mz-1.0073)/2]
  MMNa <- data.table(Feature,Adducts_S="2M+Na")
  MMNa[,NeutralMass:=(mz-22.9892)/2]
  MMK <- data.table(Feature,Adducts_S="2M+K")
  MMK[,NeutralMass:=(mz-38.9632)/2]
  MMNH4 <- data.table(Feature,Adducts_S="2M+NH4")
  MMNH4[,NeutralMass:=(mz-18.0338)/2]
  MS1<<-rbind(MH,MNa,MK,MNH4) #,M2NaH,MMH,MMNa,MMK,MMNH4
  MS1$NeutralMass<-as.numeric(MS1$NeutralMass)
  MS1$rt<-as.numeric(MS1$rt)
  print("Calculate Neutral Mass Calculation Done!")
  write.csv(MS1,"Identification_MS1.csv")
  print("Results stored in：Identification_MS1.csv")
})


##Set HomeBuilt gframe in Identification tab
lo_Identification[1,2] <- gf_HomeBuilt<-gframe("Home built Database searching",expand=T)
lo_HomeBuilt<-glayout(cont=gf_HomeBuilt)
lo_HomeBuilt[2,1]<- glabel("DatabaseSearchRelativeMassDeviation (PPM):")  
lo_HomeBuilt[2,2]<- ed_ppmHomeBuilt <- gedit(text = "5", width = 6)
lo_HomeBuilt[3,1]<- gbutton("Run", handler = function(h,...) {
  
  print("Home built Identification is working...")
  ppmHomeBuiltValue<<-as.numeric(as.character(Parameter[Parameters=="DatabaseSearchRelativeMassDeviation",Values]))
  if (ppmHomeBuiltValue<0) {ppmHomeBuiltValue<<-as.numeric(svalue(ed_ppmHomeBuilt))}
  
  PPM<-as.numeric(ppmHomeBuiltValue)
  Candidate<-MS1
  Add<-list()
  Adducts<-data.table()
  for (id in seq(1,nrow(Candidate),1)){
    #for each line of Candidate, get the Database within PPM
    Add[[id]]<-Database[(((MonoisotopicMass-Candidate$NeutralMass[id])/Candidate$NeutralMass[id])*10^6>-PPM)&(((MonoisotopicMass-Candidate$NeutralMass[id])/Candidate$NeutralMass[id])*10^6<PPM),]
    #Add the compound name into Add
    Add[[id]][,mz:=Candidate$mz[id]]
    Add[[id]][,rt:=Candidate$rt[id]]
    Add[[id]][,Index:=Candidate$rownames[id]]
    Add[[id]][,Score:=Candidate$Overall[id]]
    Add[[id]][,Intensity:=Candidate$Intensity[id]]
    Add[[id]][,Adducts_S:=Candidate$Adducts_S[id]]
    Add[[id]][,NeutralMass:=Candidate$NeutralMass[id]]
    #Bind all the candidates together
    Adducts<-rbind(Adducts,Add[[id]])    
  }
  Adducts<-data.table(Adducts)
  identifications_lipidmaps<<-unique(Adducts,by=c("mz","rt","NeutralMass","Adducts_S"))
  print("Home built Identification Done!")
})
lo_HomeBuilt[3,2]<- gbutton("?", handler = function(h,...) {gmessage("Home Built Identification only use MS1 for the identificaiton",title = "help")})  
lo_HomeBuilt[4,1]<- gbutton("export", handler = function(h,...) { 
  write.csv(identifications_lipidmaps,"identifications_HomeBuilt.csv")
  print("Export identifications of HomeBuilt Done!")
})


##Set Metfrag gframe in Identification tab
lo_Identification[2,1] <- gf_Metfrag<-gframe("Metfrag Database searching",expand=T)
lo_Metfrag<-glayout(cont=gf_Metfrag)
lo_Metfrag[2,1] <-   gbutton("Extract MS2", handler = function(h,...) {
  raw_data_mem <- readMSData(files,pdata= new("NAnnotatedDataFrame")) #readMSData from packagee MSnbase, read the files using in memory mode
  results_dt <- ldply (raw_data_mem@assayData, data.frame)
  results_dt=data.table(results_dt)#Transfer results from list to data.table
  add.mz.rt<-function(x){
    id=as.character(x[1])
    rt=raw_data_mem@assayData[[id]]@rt
    precursorMz=raw_data_mem@assayData[[id]]@precursorMz
    x["rt"]=rt
    x["precursorMz"]=precursorMz
    return(x)
  }
  MS2<-apply(results_dt, 1, add.mz.rt)# for each line get the mz and rt according to .id, and add as a new column
  MS2<-data.table(t(MS2))
  MS2[,precursorMz:=(as.numeric(precursorMz))]
  MS2[,rt:=(as.numeric(rt))]
  MS2[,mz :=(as.numeric(mz))]
  MS2[,Intensity:=(as.numeric(i))] #intensity
  MS2<-MS2[Intensity>5,] #to make MS2 smaller 
  MS2<-unique(MS2,by=c("mz","rt","precursorMz"))
  MS3<<-list(MS1=MS1,MS2=MS2)#
  print("Extract MS2 Done!")
})
lo_Metfrag[3,1] <- glabel("DatabaseSearchRelativeMassDeviation")  
lo_Metfrag[3,2] <- ed_DatabaseSearchRelativeMassDeviation<-gedit(text = "5", width = 6)
lo_Metfrag[4,1] <- glabel("FragmentPeakMatchAbsoluteMassDeviation")  
lo_Metfrag[4,2] <-  ed_FragmentPeakMatchAbsoluteMassDeviation <- gedit(text = "0.005", width = 6)
lo_Metfrag[5,1] <-  glabel("FragmentPeakMatchRelativeMassDeviationv(PPM):")  
lo_Metfrag[5,2] <-  ed_FragmentPeakMatchRelativeMassDeviation <- gedit(text = "8", width = 6)
lo_Metfrag[6,1] <-  gbutton("Run", handler = function(h,...) {
  F_Identification<- function(x){ # first define the identification function
    candidates_null<-data.frame(Null="+")
    NeutralPrecursorMass<-as.numeric(x[18])#
    PrecuresorMZ<-as.numeric(x[1])
    RT<-as.numeric(x[4])
    rtmin<-as.numeric(x[5])
    rtmax<-as.numeric(x[6])
    #print(RT)
    settingsObject[["MetFragDatabaseType"]]<-c("LocalCSV") #can try ExtendedPubChem next time
    settingsObject[["LocalDatabasePath"]]<-c("d:/xcms/software/database/lipidmaps.csv") #
    settingsObject[["NeutralPrecursorMass"]]<-NeutralPrecursorMass #neutral monoisotopic precursor mass
    #peak<-subset(MS3$MS2,(rt-RT)<60&(rt-RT)>-60)
    peak<-subset(MS3$MS2,((rt<(rtmax+5)) & (rt>(rtmin-5)) & (precursorMz>(PrecuresorMZ-0.01))& (precursorMz<(PrecuresorMZ+0.01))))
    peak_mz_intensity<-peak[,c(2,6)]
    peak_mz_intensity<- as.matrix(peak_mz_intensity, ncol=2, byrow=TRUE)
    settingsObject[["PeakList"]]<-peak_mz_intensity # get the peak list from MS2
    candidates<-data.frame()
    scored.candidates<-data.frame()
    candidates<-run.metfrag(settingsObject)
    if(is.null(candidates$Score)) {scored.candidates<-candidates_null} #fill candidates with null created
    else {scored.candidates<-candidates# add the identified candidates to scored.candidates
    # add extra information to the identifications
    scored.candidates$RentionTime<-RT
    scored.candidates$rtmin<-rtmin
    scored.candidates$rtmax<-rtmax
    scored.candidates$NeutralPrecursorMass<-NeutralPrecursorMass
    scored.candidates$PrecuresorMZ<-PrecuresorMZ
    scored.candidates$mzmin<-x[2]
    scored.candidates$mzmax<-x[3]
    scored.candidates$npeaks <-x[7]
    scored.candidates$Intensity1 <-x[10]
    scored.candidates$Intensity2 <-x[11]
    scored.candidates$IntensityMean <-x[15]
    scored.candidates$isotops<-x[12]#
    scored.candidates$adduct<-x[17]
    }
    return(scored.candidates)
  }
  settingsObject<-list()
  
  DatabaseSearchRelativeMassDeviationValue<<-as.character(Parameter[Parameters=="DatabaseSearchRelativeMassDeviation",Values])
  if (DatabaseSearchRelativeMassDeviationValue==-1) {DatabaseSearchRelativeMassDeviationValue<<-svalue(ed_DatabaseSearchRelativeMassDeviation)}
  
  FragmentPeakMatchAbsoluteMassDeviationValue<<-as.character(Parameter[Parameters=="FragmentPeakMatchAbsoluteMassDeviation",Values])
  if (FragmentPeakMatchAbsoluteMassDeviationValue==-1) {FragmentPeakMatchAbsoluteMassDeviationValue<<-svalue(ed_FragmentPeakMatchAbsoluteMassDeviation)}
  
  FragmentPeakMatchRelativeMassDeviationValue<<-as.character(Parameter[Parameters=="FragmentPeakMatchRelativeMassDeviation",Values])
  if (FragmentPeakMatchRelativeMassDeviationValue==-1) {FragmentPeakMatchRelativeMassDeviationValue<<-svalue(ed_FragmentPeakMatchRelativeMassDeviation)}
  
  
  settingsObject[["DatabaseSearchRelativeMassDeviation"]]<-as.numeric(DatabaseSearchRelativeMassDeviationValue)
  settingsObject[["FragmentPeakMatchAbsoluteMassDeviation"]]<-as.numeric(FragmentPeakMatchAbsoluteMassDeviationValue)
  settingsObject[["FragmentPeakMatchRelativeMassDeviation"]]<-as.numeric(FragmentPeakMatchRelativeMassDeviationValue)
  results<-data.frame()
  results<-apply((MS3$MS1),1,F_Identification)
  results_dt <- ldply (results, data.frame) #Transfer results from list to data.frame
  identifications_lipidmaps<-data.table(results_dt)#
  identifications_lipidmaps<<-unique(identifications_lipidmaps,by=c("PrecuresorMZ","RentionTime","NeutralPrecursorMass","adduct"))
  identifications_lipidmaps<<-identifications_lipidmaps[!"+", on=.(Null)]# filterout the null candadates
  #identifications_lipidmaps$RentionTime<-identifications_lipidmaps$RentionTime/60
  setorder(identifications_lipidmaps,-RentionTime,-Score) #
  print("MetFrag Identification Done!")
})
lo_Metfrag[6,2] <-  gbutton("?", handler = function(h,...) {gmessage("MetFrag use both MS1 and MS2 for the identificaiton",title = "help")})  
lo_Metfrag[7,1]<- gbutton("export",  handler = function(h,...) { 
  write.csv(identifications_lipidmaps,"identifications_Metfrag.csv")
  print("Export identifications of Metfrag Done!")
})
##Set Parameters for MS-Pathway method



##About
gp_About <- ggroup(cont = nb, label = "About", horizontal = F)
addSpring(gp_About)
gtext("Lipidomics identification is a function set with Gui to process mass spectrometry data,developed by Xiaodong Feng (x.feng@umcg.nl).",font.attr=list(size="xx-large",family="monospace"),expand=TRUE,fill=TRUE,container = gp_About)  

#Plot tab
gg <<- ggraphics(cont=nb, label='Plot',visible=FALSE)
#plot.new() 
F_DefaultFigure<-function(){
  plot.new()
  rasterImage(load.image("d:/github/GUI/Logo/Logo.png"),0,0,1,1)
}


F_Delete_Ggraphics<-function(){
  delete(nb,gg)
  #dev.off()
  gg <<- ggraphics(cont=nb, label='Plot',visible=FALSE)
  #dev.new()
}

F_DefaultFigure()


