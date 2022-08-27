#' @name MetGUI
#' @title MetGUI for software
#' @description The MetGUI software.
#' @author Xiaodong Feng fengxdong@outlook.com
#' @references
#' https://doi.org/10.1101/2020.10.10.334342
#' @keywords GUI
#' @keywords Metabolomics
#' @keywords Statistical analysis
#' @return The output of GUI
#' @details  RGTK2 and GTK+ are required for MetGUI. We recommend the 64-bit
#' version of Windows 7, or newer for most users for RGTK2 installation.
#' The R 4.0.1 and RGtk2 2.20.36 sailed through the test.
#' @examples
#' if (interactive()) {
#'   statTargetGUI()
#' }
#' @export
Path <- 'd:/github/MetGUI'
setwd(Path)
register(SerialParam()) #disable paralell
MetGUI <- function() {
  #### Set GUI parameters####
  F_Capture <- function(Code) {
    capture.output(Code, file = "Output.log", type = "message", append = TRUE)
  }
  BPPARAM <- BiocParallel::SnowParam(detectCores() - 3, progressbar = TRUE)
  options(warn = -1, show.error.messages = TRUE)
  RangeRT <- c(35.07705, 1258.41516)
  RangeMZ <- c(250.1777, 1746.5595)
  Database <<- data.table(read.csv("input/lipidmaps2017Short.csv"))
  color <- grDevices::colors()
  linetype <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "F8", "431313", "22848222")
  linetype <- rep(linetype, 4)
  widgets <- list()
  val <- list()
  val$binSize = val$fwhm = val$max = val$snthresh = val$steps = val$mzdiff = val$ppm = val$peakwidthMin = val$peakwidthMax = val$snthreshCentWave = val$noise = val$binSizeRT = val$bw = val$binSizeGrouping = val$ppmIsotope = val$ppmWindow = val$rtWindow = val$Database = val$DatabaseSearchRelativeMassDeviation = val$FragmentPeakMatchAbsoluteMassDeviation = val$FragmentPeakMatchRelativeMassDeviation = val$ppmHomeBuilt = -1
    # win <- gwindow("MetGUI", visible=FALSE,
  #                handler = function(h, ...) {
  #                  gtkMainQuit() # stop main loop when windows is closed.
  #                }
  # )
  win <- gWidgets2::gwindow("MetGUI", expand = FALSE, fill = FALSE, height = 800, width = 800)
  nb <- gnotebook(cont = win, fill = TRUE, expand = TRUE)
  sb <- gstatusbar("Contact Us: fengxdong@outlook.com  University of Groningen", container = win)
  font(sb) <- list(size = 9, color = "blue")
  svalue(nb) <- 1

  #### ImportData####
  gp_ImportData <- ggroup(horizontal = FALSE, cont = nb, label = "Workflow")
  lo_ImportData <- glayout(cont = gp_ImportData, horizontal = F)
  ## Set Dataset frame in ImportData tab
  lo_ImportData[1, 1] <- gf_Dataset <- gframe("1 Data transformation", horizontal = F, expand = TRUE)
  lo_Dataset <- glayout(cont = gf_Dataset, horizontal = F)
  lo_Dataset[1, 1] <- widgets$bt_folder <- gbutton("Select dataset folder", handler = function(h, ...) {
    widgets$my_path <- choose.dir(caption = "Select folder")
    if (length(widgets$my_path) == 0) {
      galert("Please Select the folder to be Processed", title = "File Selection Problems", delay = 1)
    } else {
      cdffiles <<- list.files(widgets$my_path, recursive = TRUE, full.names = TRUE)
      print(paste("You selected", cdffiles))
      name <<- list.files(widgets$my_path, recursive = TRUE, full.names = FALSE)
      # show the loaded files below
      ed_group <<- list()
      for (i in 1:length(cdffiles)) {
        lo_Dataset[i + 4, 1] <- name[i]
        lo_Dataset[i + 4, 2] <- ed_group[i] <<- gedit(text = "Norm", width = 4)
      }
    }
  })
  tooltip(widgets$bt_folder) <- "Select the folder for all the raw data to be processed"
  ### Select export folder
  lo_Dataset[2, 1] <- widgets$bt_selectexportfolder <- gbutton("Select export folder", handler = function(h, ...) {
    path <- choose.dir(caption = "Select export folder")
    if (length(path) == 0) {
      galert("Please Select the export folder", title = "File Selection Problems", delay = 1)
    } else {
      setwd(path)
      print("Export folder selected successfully!")
    }
  })
  tooltip(widgets$bt_selectexportfolder) <- "Select the folder to save all the exported results"
  ## Set Parameters frame in ImportData tab
  lo_Dataset[3, 1] <- widgets$bt_selectparameterfile <- gbutton("Select parameter file", handler = function(h, ...) {
    Parameters <- choose.files(caption = "Select Parameter File")
    if (length(Parameters) == 0) {
      galert("Please Select the Parameter File to be Processed", title = "File Selection Problems", delay = 1)
    } else {
      Parameter <<- data.table(read.csv(Parameters))
      print("Parameters loaded successfully!")
      print(head(Parameter))
    }
  })
  tooltip(widgets$bt_selectparameterfile) <- "Select the parameter file"

  ### Set Experimental design frame in ImportData tab
  lo_Dataset[4, 1] <- widgets$bt_importdata <- gbutton("Import data", handler = function(h, ...) {
    # Define a data.frame with sample descriptions
    Group <<- c(rep('Hyp',3), rep('Norm',3))
    pd <- data.frame(file = cdffiles, sample_group = c(rep('Hyp',3), rep('Norm',3)))
    print(pd)
    # Read the files
    data <- readMSData(cdffiles, pd = new("NAnnotatedDataFrame", pd), mode = "onDisk")
    Raw_data <<- data
    Raw_data <<- Raw_data[grep(fData(Raw_data)$msLevel, pattern = "1")]
    group_colors <<- brewer.pal(12, "Set3") # set colours
    names(group_colors) <- name
    group_colors <<- group_colors
    print("Data imported successfully!")
  })
  tooltip(widgets$bt_importdata) <- "Import the data into an object"
  ### Set TIC in Design frame
  lo_Dataset[1, 2] <- widgets$bt_plotTIC <- gbutton("Plot TIC", handler = function(h, ...) {
    # F_Delete_Ggraphics()
    # pdf(paste0('Total ion current','.pdf'))
    svg(paste0('Total ion current','.svg'))
    bpis <- chromatogram(Raw_data, aggregationFun = "sum")
    plot(bpis, col = group_colors[name])
    legend("topright", legend = name, col = group_colors[name], lty = 1, cex = 0.8) # ,box.lty=0
    dev.off()
    print("Plot Total ion current Done!")
  })
  tooltip(widgets$bt_plotTIC) <- "Plot the total ion current of the raw data"
  lo_Dataset[2, 2] <- widgets$bt_plotBPC <- gbutton("Plot BPC", handler = function(h, ...) {
    # F_Delete_Ggraphics()
    # pdf(paste0('Total ion current','.pdf'))
    svg(paste0('Total ion current','.svg'))
    bpis <- chromatogram(Raw_data, aggregationFun = "max")
    plot(bpis, col = group_colors[name])
    legend("topright", legend = name, col = group_colors[name], lty = 1, cex = 0.8) # ,box.lty=0
    dev.off()
    print("Plot BPC Done!")
  })
  tooltip(widgets$bt_plotBPC) <- "Plot the BPC of the raw data"

  lo_Dataset[3, 2] <- widgets$bt_export <- gbutton("Export metadata", handler = function(h, ...) {
    ## extract MS1 according to function
    # fun_1 <- grep(fData(Raw_data)$msLevel, pattern = "1")
    # raw_data_all[fun_1]
    Fdata1 <- data.table(fData(Raw_data))
    write.csv(Fdata1, "ImportData_Fdata1.csv",row.names = FALSE)
    print("Results exported successfully!")
  })
  tooltip(widgets$bt_export) <- "Export base peak information from the raw dataset"
  
  lo_Dataset[4, 2] <- gbutton("Export parameter", handler = function(h, ...) {
    ParametersExported <- data.table()
    ParametersExported <- Parameter
    UserValues <- c(binSizeValue, fwhmValue, maxValue, snthreshValue, stepsValue, mzdiffValue, ppmValue, peakwidthMinValue, peakwidthMaxValue, snthreshCentWaveValue, noiseValue, binSizeRTValue, bwValue, binSizeGroupingValue, ppmIsotopeValue, ppmWindowValue, rtWindowValue, DatabaseValue, DatabaseSearchRelativeMassDeviationValue, FragmentPeakMatchAbsoluteMassDeviationValue, FragmentPeakMatchRelativeMassDeviationValue, ppmHomeBuiltValue)
    ParametersExported[, Values := UserValues]
    write.csv(ParametersExported, "Parameters_ParametersExported.csv")
    print("Parameters exported successfully!")
  })
  #### PeakDetection ####
  ## set CentWave parameters on gp_CentWave
  lo_ImportData[1, 2] <- gf_CentWave <- gframe("2 Peak detection", horizontal = F, expand = TRUE)
  lo_CentWave <- glayout(cont = gf_CentWave, horizontal = F)
  lo_CentWave[1, 1] <- glabel("Mass accuracy")
  lo_CentWave[1, 2] <- widgets$ed_ppm <- gedit(text = "1", width = 4)
  tooltip(widgets$ed_ppm) <- "ppm. Default 25. maximal tolerated mz deviation in consecutive scans in parts per million (ppm) for the initial ROI definition. Guidelines are 5–15 ppm for Orbitrap data, ~5 ppm for lock mass quadrupole time of flight (QTOF) data and 10–20 ppm for other QTOF instruments."

  lo_CentWave[2, 1] <- glabel("MRP MS1")
  lo_CentWave[2, 2] <- widgets$ed_MRPMS1<- gedit(text = "30000", width = 4)
  tooltip(widgets$ed_MRPMS1) <- "Mass resolving power"
  
  lo_CentWave[3, 1] <- glabel("Reference m/z")
  lo_CentWave[3, 2] <- widgets$ed_RefMZ<- gedit(text = "400", width = 4)
  tooltip(widgets$ed_RefMZ) <- "Reference m/z"
  
  lo_CentWave[4, 1] <- glabel("Instrument type")
  lo_CentWave[4, 2] <- widgets$ed_Instrument<- gedit(text = "2", width = 4)
  tooltip(widgets$ed_Instrument) <- "Instrument type"
  
  lo_CentWave[5, 1] <- glabel("Peakwidth min")
  lo_CentWave[5, 2] <- widgets$ed_peakwidthMin <- gedit(text = "2", width = 4)
  tooltip(widgets$ed_peakwidthMin) <- "Peakwidth. Default 20-50. minimum and maximum chromatographic peak width detectable (seconds). This depends mainly on the type of chromatographic separation performed. For standard reverse-phase separations, a general guideline is 20–60 s, whereas for hydrophilic interaction liquid chromatography (HILIC), in which run times tend to be longer with broader peaks, we recommend 25–90 s. When running with UPLC, these values should drop markedly because of shorter run times and higher resolution, with suggested starting values between 2 and 5 s to a maximum of 30 s. If in doubt of values, check the raw chromatographic run for peak widths of some common compounds from each end of the trace. These values are not hard cutoffs and may be detected slightly out of this range, depending on the quality of the peak data"

  lo_CentWave[6, 1] <- glabel("Peakwidth max")
  lo_CentWave[6, 2] <- ed_peakwidthMax <- gedit(text = "30", width = 4)

  lo_CentWave[7, 1] <- glabel("Signal/Noise ratio")
  lo_CentWave[7, 2] <- widgets$ed_snthreshCentWave <- gedit(text = "5", width = 4)
  tooltip(widgets$ed_snthreshCentWave) <- "snthresh. Default 10,Signal/Noise ratio cutoff: ([maximum peak intensity] - [estimated baseline value]) / standard deviation of local chromatographic noise. Higher than default value to reduce false positive results."

  lo_CentWave[8, 1] <- glabel("Noise level")
  lo_CentWave[8, 2] <- widgets$ed_noise <- gedit(text = "100", width = 4)
  tooltip(widgets$ed_noise) <- "noise. Default 0. centroids with intensity < noise are omitted from ROI detection"

  lo_CentWave[9, 1:2] <- gbutton("Run peak detection", handler = function(h, ...) {
    val$ppm <<- as.numeric(as.character(Parameter[Parameters == "ppm", Values]))
    if (val$ppm == -1) {
      val$ppm <<- as.numeric(svalue(widgets$ed_ppm))
    }
    val$MRPMS1 <<- as.numeric(svalue(widgets$ed_MRPMS1))
    val$RefMZ <<- as.numeric(svalue(widgets$ed_RefMZ))
    val$Instrument <<- as.numeric(svalue(widgets$ed_Instrument))
    
    val$peakwidthMin <<- as.numeric(as.character(Parameter[Parameters == "peakwidthMin", Values]))
    if (val$peakwidthMin == -1) {
      val$peakwidthMin <- as.numeric(svalue(widgets$ed_peakwidthMin))
    }

    val$peakwidthMax <<- as.numeric(as.character(Parameter[Parameters == "peakwidthMax", Values]))
    if (val$peakwidthMax == -1) {
      val$peakwidthMax <- as.numeric(svalue(ed_peakwidthMax))
    }

    val$snthreshCentWave <<- as.numeric(as.character(Parameter[Parameters == "snthreshCentWave", Values]))
    if (val$snthreshCentWave == -1) {
      val$snthreshCentWave <<- as.numeric(svalue(widgets$ed_snthreshCentWave))
    }
    val$noise <<- as.numeric(as.character(Parameter[Parameters == "noise", Values]))
    if (val$noise == -1) {
      val$noise <- as.numeric(svalue(widgets$ed_noise))
    }
    # MRP=70000 # Velos, 30000 for MS1 and 7500 for MS2 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3673017/
    # RefMZ=200 #  Velos is mz reference 400 while Q-Exactive is 200.
    # A <- 1 / (MRP * (RefMZ^0.5))
    # print(A)
    # B <- A / 2.35482
    # print(B)
    cwp <- xcms::CentWaveParam(
      A = 7.077682e-07, Instrument = 2, 
      ppm = as.numeric(val$ppm), # camera 30, maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition.
      peakwidth = c(as.numeric(val$peakwidthMin), as.numeric(val$peakwidthMax)), # camera (5,12) minimum and maximum chromatographic peak width detectable (seconds), also use c(30, 80)
      prefilter=c(3, 500),firstBaselineCheck = FALSE,
      snthresh = as.numeric(val$snthreshCentWave), # camera/default use 6,Signal/Noise ratio cutoff: ([maximum peak intensity] - [estimated baseline value]) / standard deviation of local chromatographic noise
      noise = as.numeric(val$noise) # centroids with intensity < noise are omitted from ROI detection
      )
    xdata1 <<- findChromPeaks(Raw_data, param = cwp, BPPARAM = BPPARAM)
    # Raw_data2@featureData@data$msLevel<-1
    # xdata2 <<- findChromPeaks(Raw_data2, param = cwp)
    print("CentWave Peak detection Done!")
  })
  lo_CentWave[10, 1:2] <- gbutton("Export peaks", handler = function(h, ...) {
    write.csv(xdata1@msFeatureData[["chromPeaks"]], "PeakDetection_xdata1@msFeatureData_chromPeaks.csv")
    print("Results exported successfully!")
  })
  #### Spectra pre-processing ####
  lo_ImportData[1, 3] <- gf_Processing <- gframe("3 Spectrum pre-processing", horizontal = F, expand = TRUE)
  lo_Processing <- glayout(cont = gf_Processing, horizontal = F)
  ### Set Obiwarp RT alignment frame in Retention Time Alignment tab
  lo_Processing[1, 1] <- glabel("binSize") ### set binSize
  lo_Processing[1, 2] <- widgets$ed_binSizeRT <- gedit(text = "1", width = 4)
  lo_Processing[2, 1:2] <- widgets$bt_Processing <- gbutton("Run alignment", handler = function(h, ...) {
    val$binSizeRT <<- as.numeric(as.character(Parameter[Parameters == "binSizeRT", Values]))
    if (val$binSizeRT == -1) {
      val$binSizeRT <<- as.numeric(svalue(widgets$ed_binSizeRT))
    }
    xdata1 <<- adjustRtime(xdata1, param = ObiwarpParam()) # binSize defining the bin size (in mz dimension) to be used for the profile matrix generation.
    # xdata1 <<- adjustRtime(xdata1, param = ObiwarpParam(binSize =as.numeric(binSizeRTValue)))
    print("Obiwarp RT Alignment Done!")
  })
  tooltip(widgets$bt_Processing) <- "Obiwarp RT Alignment method performs retention time adjustment using the Obiwarp method. It supports alignment of multiple samples by aligning each against a center sample. The alignment is performed directly on the profile-matrix. binSize defining the bin size (in mz dimension) to be used for the profile matrix generation."
  ### Set Retention Time Alignment effect gframe in Retention Time Alignment tab
  # lo_Processing[3, 1:2] <- widgets$bt_TIC <- gbutton("Plot alignment effect TIC", handler = function(h, ...) {
  #   # F_Delete_Ggraphics()
  #   svg(paste0('RT Alignment effect','.svg'))
  #   par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
  #   ## before time alignment
  #   plot(chromatogram(Raw_data), col = group_colors[name]) # ,main="Before retention time alignment"
  #   legend("topright", legend = name, col = group_colors[name], lty = 1, cex = 0.8) # ,box.lty=0
  #   ## after time alignment
  #   plot(chromatogram(xdata1), col = group_colors[name]) # ,main="After retention time alignment"
  #   legend("topright", legend = name, col = group_colors[name], lty = 1, cex = 0.8) # ,box.lty=0
  #   dev.off()
  #   print("Plot RT Alignment effect Done!")
  # })
  # tooltip(widgets$bt_TIC) <- "Plot the TIC before and after retention time alignment"
  lo_Processing[3, 1] <- glabel(paste("rtMin (s)", sep = ""))
  lo_Processing[3, 2] <- widgets$rtMin <- gedit(text = "760", width = 6)
  lo_Processing[4, 1] <- glabel(paste("rtMax (s)", sep = ""))
  lo_Processing[4, 2] <- widgets$rtMax <- gedit(text = "810", width = 6)
  lo_Processing[5, 1] <- glabel(paste("mzMin (Da)", sep = ""))
  lo_Processing[5, 2] <- widgets$mzMin <- gedit(text = "175.10", width = 6)
  lo_Processing[6, 1] <- glabel(paste("mzMax (Da)", sep = ""))
  lo_Processing[6, 2] <- widgets$mzMax <- gedit(text = "175.20", width = 6)
  lo_Processing[7, 1:2] <- ed_PlotAligned <- gradio(c('TIC','BPC','EIC'), horizontal=TRUE, width = 4)
  lo_Processing[8, 1:2] <- widgets$bt_ExtractedIonChromatogram <- gbutton("Plot raw and aligned EICs", handler = function(h, ...) {
    val$mzMin <<- as.numeric(svalue(widgets$mzMin))
    val$mzMax <<- as.numeric(svalue(widgets$mzMax))
    val$rtMin <<- as.numeric(svalue(widgets$rtMin))
    val$rtMax <<- as.numeric(svalue(widgets$rtMax))
    # F_Delete_Ggraphics()
    svg(paste0('Extracted ion chromatogramt','.svg'))
    par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
    ## before time alignment
    Before <- chromatogram(Raw_data, aggregationFun = "sum", mz = c(val$mzMin, val$mzMax),
                           rt = c(val$rtMin, val$rtMax))
    plot(Before, col = group_colors[name]) # ,main="Before retention time alignment"
    legend("topright", legend = name, col = group_colors[name], lty = 1, cex = 0.8) # ,box.lty=0
    ## after time alignment
    After <- chromatogram(xdata1, aggregationFun = "sum", mz = c(val$mzMin, val$mzMax), 
             rt = c(val$rtMin, val$rtMax))
    plot(After, col = group_colors[name],peakType = "none") # ,main="After retention time alignment"
    legend("topright", legend = name, col = group_colors[name], lty = 1, cex = 0.8) # ,box.lty=0
    print("Plot Extracted ion chromatogramt Done!")
    dev.off()
  })
  tooltip(widgets$bt_ExtractedIonChromatogram) <- "Plot the EIC before and after retention time alignment"


  ## Set Density Processing gframe in Processing tab
  lo_Processing[1, 3] <- glabel("Band width") ### set binSize
  lo_Processing[1, 4] <- ed_bw <- gedit(text = "10", width = 4)
  lo_Processing[2, 3] <- glabel("binSize") ### set binSize
  lo_Processing[2, 4] <- ed_binSizeProcessing <- gedit(text = "0.02", width = 8)
  lo_Processing[3, 3] <- glabel("minFraction")
  lo_Processing[3, 4] <- ed_minFraction <- gedit(text = "0.6", width = 8)

  lo_Processing[4, 3:4] <- widgets$bt_DensityProcessing <- gbutton("Peak grouping in dataset", handler = function(h, ...) {
    val$bw <<- as.numeric(as.character(Parameter[Parameters == "bw", Values]))
    if (val$bw == -1) {
      val$bw <<- as.numeric(svalue(ed_bw))
    }
    val$binSizeGrouping <<- as.numeric(as.character(Parameter[Parameters == "binSizeGrouping", Values]))
    if (val$binSizeGrouping == -1) {
      val$binSizeGrouping <<- as.numeric(svalue(ed_binSizeGrouping))
    }
    val$minFraction <<- as.numeric(svalue(ed_minFraction))
    pdp <- PeakDensityParam(
      minFraction = val$minFraction,
      minSamples = 2,
      sampleGroups = Group,
      bw = as.numeric(val$bw),
      binSize = as.numeric(val$binSizeGrouping) # camera uses 0.015,default 0.25.width of overlapping m/z slices to use for creating peak density chromatograms and grouping peaks across samples
    )
    xdata1 <<- groupChromPeaks(xdata1, param = pdp)
    # xdata1 <<- groupChromPeaks(xdata1,param = PeakDensityParam(sampleGroups = Group, bw = 10))
    Feature <<- xdata1@msFeatureData[["featureDefinitions"]] # for plot XIC and grouping effect
    RangeRT <<- range(Feature$rtmed)
    RangeMZ <<- range(Feature$mzmed)
    # xdata2 <<- groupChromPeaks(xdata2, param = pdp)
    print("Grouping Done!")
  })
  tooltip(widgets$bt_DensityGrouping) <- "groupChromPeaks-density method performs correspondence (chromatographic peak grouping) based on the density (distribution) of identified peaks along the retention time axis within slices of overlapping mz ranges. All peaks (from the same or from different samples) being close on the retention time axis are grouped into a feature (peak group).bw: The maximum expected RT deviation across samples."

  lo_Processing[5, 3:4] <- gbutton("Export grouped peaks", handler = function(h, ...) {
    write.csv(xdata1@msFeatureData[["featureDefinitions"]], "Grouping_xdata1_msFeatureData_featureDefinitions.csv")
    # write.csv(xdata2@msFeatureData[["featureDefinitions"]],"Grouping_xdata2_msFeatureData_featureDefinitions.csv")
    print("Results exported successfully!")
  })
  ## Set GroupingEffect gframe in Grouping tab
  lo_Processing[6, 3:4] <- widgets$bt_PlotGrouping <- gbutton("Plot peak grouping", handler = function(h, ...) {
    # F_Delete_Ggraphics()
    svg(paste0('Grouping effect','.svg'))
    plotChromPeakDensity(xdata1,
      col = group_colors[name], mz = c(val$mzMin, val$mzMax), rt = c(val$rtMin, val$rtMax),
      param = PeakDensityParam(sampleGroups = Group, bw = val$bw, binSize = val$binSizeGrouping)
    )
    legend("topright", legend = name, col = group_colors[name], lty = 1, cex = 0.8) # ,box.lty=0
    dev.off()
    print("Plot Grouping effect Done!")
  })
  tooltip(widgets$bt_PlotGrouping) <- "Plot Grouping effect to see the parameters bw and Binsize's effect on grouping"
  
  ## Set Fill missing peaks default gframe in Fill missing peaks tab
  lo_Processing[7, 3:4] <- widgets$bt_FillMissing <- gbutton("Fill missing peaks", handler = function(h, ...) {
    xdata1 <<- fillChromPeaks(xdata1)
    print("Fill missing peaks Done!")
  })
  tooltip(widgets$bt_FillMissing) <- "Integrate signal in the mz-rt area of a feature (chromatographic peak group) for samples in which no chromatographic peak for this feature was identified and add it to the chromPeaks. Such peaks will have a value of 1 in the is_filled column of the chromPeaks matrix of the object."


  #### Annotation ####
  lo_ImportData[2, 1] <- gf_Isotope <- gframe("4 Annotation", horizontal = F, expand = TRUE)
  lo_Isotope <- glayout(cont = gf_Isotope, horizontal = F)
  ## Set IsotopeFilteration gframe in Fill missing peaks tab
  # lo_Isotope[1, 1] <- glabel("ppm Isotope")
  # lo_Isotope[1, 2] <- ed_ppmIsotope <- gedit(text = "10", width = 4)
  lo_Isotope[2, 1:2] <- widgets$bt_Isotope <- gbutton("Isotope filteration", handler = function(h, ...) {
    val$ppmIsotope <<- as.numeric(as.character(Parameter[Parameters == "ppmIsotope", Values]))
    if (val$ppmIsotope == -1) {
      val$ppmIsotope <<- as.numeric(svalue(ed_ppmIsotope))
    }
    xset <- as(xdata1, "xcmsSet")
    xsa <- xsAnnotate(xset) # create an CAMERA object
    xsa_I <- findIsotopes(xsa, ppm = val$ppmIsotope) # ,mzabs=0.009,filter = FALSE
    peaklist_xsa_I <- getPeaklist(xsa_I)
    DT_peaklist <- data.table(peaklist_xsa_I) # transfer into datatable format
    intensity <- apply(DT_peaklist[, 11:16], 1, mean, na.rm = T) # get the mean of intensity
    DT_peaklist[, Intensity := intensity] # add a new column Intensity to ms1
    setorder(DT_peaklist, -Intensity)
    IsoFilter <- str_replace_all(DT_peaklist$isotopes, pattern = ".*M\\]\\+$", replacement = "")
    DT_peaklist[, IsoFilter := IsoFilter]
    DT_peaklist <- DT_peaklist[IsoFilter == "", ] #
    Feature <<- unique(DT_peaklist, by = c("mz", "rt"))
    print("Isotope Filteration Done!")
  })
  tooltip(widgets$bt_Isotope) <- "Annotate isotope peaks using CAMERA. ppm is the ppm error for the search "
  ## Set Remove Duplicates gframe in Fill missing peaks tab
  # lo_Isotope[3, 1] <- glabel("ppm Window")
  # lo_Isotope[3, 2] <- ed_ppmWindow <- gedit(text = "5", width = 4)
  lo_Isotope[4, 1] <- glabel("rt Window")
  lo_Isotope[4, 2] <- ed_rtWindow <- gedit(text = "5", width = 4)
  lo_Isotope[5, 1:2] <- widgets$bt_RemoveDup <- gbutton("Remove duplicates", handler = function(h, ...) {
    val$ppmWindow <<- as.numeric(as.character(Parameter[Parameters == "ppmWindow", Values]))
    if (val$ppmWindow == -1) {
      val$ppmWindow <<- as.numeric(svalue(ed_ppmWindow))
    }

    val$rtWindow <<- as.numeric(as.character(Parameter[Parameters == "rtWindow", Values]))
    if (val$rtWindow == -1) {
      val$rtWindow <<- as.numeric(svalue(ed_rtWindow))
    }
    Peaks <- as.data.table(Feature)
    Peaks <- unique(Peaks, by = c("mz", "rt"))
    # Extract the duplicates in Feature_MZmine
    PeaksCopy <- Peaks
    Duplist <- list()
    Duplicates <- data.table()
    MRP <- 30000
    RefMZ <- 400
    A <- 1 / (MRP * (RefMZ^0.5))
    B <- A / 2.35482
    for (id in seq(1, nrow(Peaks), 1)) {
      # for each line of Peaks, get the duplicates within +/- 5PPM and +/- 20 rt
      Duplist[[id]] <- Peaks[(abs(mz - PeaksCopy$mz[id]) < B * mz ^1.5 + mz / 1000000) &
        (abs(rt - PeaksCopy$rt[id]) < val$rtWindow), ]
      # Duplist[[id]] <- Peaks[(((mz - PeaksCopy$mz[id]) / mz) * 10^6 > -val$ppmWindow) & (((mz - PeaksCopy$mz[id]) / mz) * 10^6 < val$ppmWindow) & 
      #                          ((rt - PeaksCopy$rt[id]) < val$rtWindow) & ((rt - PeaksCopy$rt[id]) > -val$rtWindow), ]
      Duplist[[id]] <- Duplist[[id]][-c(1)] # get rid of the first line
      # Bind all the candidates together
      Duplicates <- rbind(Duplicates, Duplist[[id]])
      Duplicates <- unique(Duplicates, by = c("mz", "rt"))
    }
    # get rid of duplicates from Feature_MZmine
    Feature <- anti_join(Peaks, Duplicates, by = c("mz", "rt"))
    Feature <<- data.table(Feature)
    print(nrow(Peaks))
    print(nrow(Feature))
    print("Remove duplicates Done!")
  })
  tooltip(widgets$bt_RemoveDup) <- "Remove duplicates in selected mz and rt window. ppm is the mz window, rt is the retention time window"
  lo_Isotope[6, 1:2] <- gbutton("Export features", handler = function(h, ...) {
    write.csv(Feature, "Annotation_Feature.csv")
    print("Export Features Done!")
  })
  
  ## Set Scatter Plot gframe in Fill missing peaks tab
  lo_Isotope[7, 1] <- glabel("Log")
  lo_Isotope[7, 2] <- cb_log <- gcombobox(c("2", "10"), expand = T)
  lo_Isotope[8, 1:2] <- widgets$bt_Scatter <- gbutton("Plot annotation", handler = function(h, ...) {
    ### figure
    if (as.numeric(svalue(cb_log)) == 2) {
      g <- ggplot(Feature, aes(mz, rt)) +
        geom_point(aes(color = log2(Intensity))) +
        labs(title = "Scatter Plot of detected features") #+scale_color_manual(values =brewer.pal(9,"Set1")) #+geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)
    } else if (as.numeric(svalue(cb_log)) == 10) {
      g <- ggplot(Feature, aes(mz, rt)) +
        geom_point(aes(color = log10(Intensity))) +
        labs(title = "Scatter Plot of detected features") #+scale_color_manual(values =brewer.pal(9,"Set1")) #+geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)
    }
    # F_Delete_Ggraphics()
    svg(paste0('Scatter Plot','.svg'))
    plot(g)
    # ggsave("Scatter Plot.png", width = 8, height = 8)
    dev.off()
    print("Scatter Plot Done!")
  })
  tooltip(widgets$bt_Scatter) <- "Scatter Plot of detected features, the colour indicates the intensiy of the features"

  #### Metabolite identification####
  lo_ImportData[2, 2] <- gf_Identification <- gframe("5 Metabolite identification", horizontal = F, expand = TRUE)
  lo_Identification <- glayout(cont = gf_Identification, horizontal = F)
  ## Set Preprocessing gframe in Identification tab
  lo_Identification[1, 1:2] <- gbutton("Create query library", handler = function(h, ...) {
    print("Link MS1 with MS2 is working...")
    cdfs <<- dir("d:/github/MetGUI/input/NormalVSHYP", full.names = TRUE,recursive = TRUE)
    pa <<- purityA(cdfs)
    # if (length(unique(msLevel(xdata1))) != 1) {
    #   xdata1 <- filterMsLevel(xdata1, msLevel = 1)
    # } # to use MS1
    pa <<- frag4feature(pa, xdata1) # link the MS2 with MS1
    pa <<- filterFragSpectra(pa) # filter the fragmentation spectra , allfrag = TRUE
    pa <<- averageAllFragSpectra(pa) # treat the inter and intra fragmentation scans the same
    # pa@puritydf
    print("Link MS1 with MS2 done!")
    print("Create query database is working...")
    q_dbPthValue <<- msPurity::createDatabase(pa, xdata1, dbName = "q_dbPth.sqlite")
    q_dbPthValue <<- 'd:/github/MetGUI/output/q_dbPth.sqlite'
    qon <- DBI::dbConnect(RSQLite::SQLite(), q_dbPthValue)
    ## Extract from the query database the "c_peak_groups", which will be used for subset.
    # c_peak_groups <- qon %>%
    #   dplyr::tbl("c_peak_groups") %>%
    #   dplyr::collect() %>%
    #   as.data.table(.)
    # names(c_peak_groups)
    # q_xcmsGroupsValue <- c_peak_groups$grpid # 455 # use this index for mass query
    print("Create query database done!")
  })
  lo_Identification[2, 1:2] <- gbutton("Create consensus library", handler = function(h, ...) {
    print("Create library database is working...")
    l_dbPthValue <<- "d:/github/MetGUI/input/MoNA-export-LC-MS-MS_Spectra.sqlite"
    con <- DBI::dbConnect(RSQLite::SQLite(), l_dbPthValue)
    print("Create library database done!")
  })
  
  ## Set HomeBuilt gframe in Identification tab
  lo_Identification[3, 1] <- glabel("MRP MS2")
  lo_Identification[3, 2] <- ed_MRPMS2 <- gedit(text = "7500", width = 4)
  lo_Identification[4, 1:2] <- ed_MS2 <- gradio(c('Exp. MS2','In-silico MS2'), horizontal=TRUE, width = 4)
  lo_Identification[5, 1:2] <- widgets$bt_HomeBuilt <- gbutton("Run identification", handler = function(h, ...) {
    print("Experimental Identification is working...")
    result <- spectralMatching(
      q_dbPth = q_dbPthValue, l_dbPth = l_dbPthValue
      # l_accessions = l_accessionsValue,
      # q_xcmsGroups = q_xcmsGroupsValue
    )
    # names(result)
    # result
  })
  tooltip(widgets$bt_HomeBuilt) <- "Experimental Identification"
  
  ## Set Metfrag gframe in Identification tab
  
  # lo_Identification[7, 1] <- glabel("DatabaseSearchRelativeMassDeviation")
  # lo_Identification[7, 2] <- ed_DatabaseSearchRelativeMassDeviation <- gedit(text = "5", width = 6)
  # lo_Identification[8, 1] <- glabel("FragmentPeakMatchAbsoluteMassDeviation")
  # lo_Identification[8, 2] <- ed_FragmentPeakMatchAbsoluteMassDeviation <- gedit(text = "0.005", width = 6)
  # lo_Identification[9, 1] <- glabel("FragmentPeakMatchRelativeMassDeviationv(PPM):")
  # lo_Identification[9, 2] <- ed_FragmentPeakMatchRelativeMassDeviation <- gedit(text = "8", width = 6)
  # lo_Identification[5, 1:2] <- widgets$bt_Metfrag <- gbutton("Run metfrag", handler = function(h, ...) {
  #   F_Identification <- function(x) { # first define the identification function
  #     candidates_null <- data.frame(Null = "+")
  #     NeutralPrecursorMass <- as.numeric(x[18]) #
  #     PrecuresorMZ <- as.numeric(x[1])
  #     RT <- as.numeric(x[4])
  #     rtmin <- as.numeric(x[5])
  #     rtmax <- as.numeric(x[6])
  #     # print(RT)
  #     settingsObject[["MetFragDatabaseType"]] <- c("LocalCSV") # can try ExtendedPubChem next time
  #     settingsObject[["LocalDatabasePath"]] <- c("inst/extdata/Database/lipidmaps2017Short.csv") #
  #     settingsObject[["NeutralPrecursorMass"]] <- NeutralPrecursorMass # neutral monoisotopic precursor mass
  #     # peak<-subset(MS3$MS2,(rt-RT)<60&(rt-RT)>-60)
  #     peak <- subset(MS3$MS2, ((rt < (rtmax + 5)) & (rt > (rtmin - 5)) & (precursorMz > (PrecuresorMZ - 0.01)) & (precursorMz < (PrecuresorMZ + 0.01))))
  #     peak_mz_intensity <- peak[, c(2, 6)]
  #     peak_mz_intensity <- as.matrix(peak_mz_intensity, ncol = 2, byrow = TRUE)
  #     settingsObject[["PeakList"]] <- peak_mz_intensity # get the peak list from MS2
  #     candidates <- data.frame()
  #     scored.candidates <- data.frame()
  #     candidates <- run.metfrag(settingsObject)
  #     if (is.null(candidates$Score)) {
  #       scored.candidates <- candidates_null
  #     } # fill candidates with null created
  #     else {
  #       scored.candidates <- candidates # add the identified candidates to scored.candidates
  #       # add extra information to the identifications
  #       scored.candidates$RentionTime <- RT
  #       scored.candidates$rtmin <- rtmin
  #       scored.candidates$rtmax <- rtmax
  #       scored.candidates$NeutralPrecursorMass <- NeutralPrecursorMass
  #       scored.candidates$PrecuresorMZ <- PrecuresorMZ
  #       scored.candidates$mzmin <- x[2]
  #       scored.candidates$mzmax <- x[3]
  #       scored.candidates$npeaks <- x[7]
  #       scored.candidates$Intensity1 <- x[10]
  #       scored.candidates$Intensity2 <- x[11]
  #       scored.candidates$IntensityMean <- x[15]
  #       scored.candidates$isotops <- x[12] #
  #       scored.candidates$adduct <- x[17]
  #     }
  #     return(scored.candidates)
  #   }
  #   settingsObject <- list()
  #   
  #   val$DatabaseSearchRelativeMassDeviation <<- as.character(Parameter[Parameters == "DatabaseSearchRelativeMassDeviation", Values])
  #   if (val$DatabaseSearchRelativeMassDeviation == -1) {
  #     val$DatabaseSearchRelativeMassDeviation <<- 10
  #   }
  #   
  #   val$FragmentPeakMatchAbsoluteMassDeviation <<- as.character(Parameter[Parameters == "FragmentPeakMatchAbsoluteMassDeviation", Values])
  #   if (val$FragmentPeakMatchAbsoluteMassDeviation == -1) {
  #     val$FragmentPeakMatchAbsoluteMassDeviation <<- 0.005
  #   }
  #   
  #   val$FragmentPeakMatchRelativeMassDeviation <<- as.character(Parameter[Parameters == "FragmentPeakMatchRelativeMassDeviation", Values])
  #   if (val$FragmentPeakMatchRelativeMassDeviation == -1) {
  #     val$FragmentPeakMatchRelativeMassDeviation <<- 10
  #   }
  #   
  #   
  #   settingsObject[["DatabaseSearchRelativeMassDeviation"]] <- as.numeric(val$DatabaseSearchRelativeMassDeviation)
  #   settingsObject[["FragmentPeakMatchAbsoluteMassDeviation"]] <- as.numeric(val$FragmentPeakMatchAbsoluteMassDeviation)
  #   settingsObject[["FragmentPeakMatchRelativeMassDeviation"]] <- as.numeric(val$FragmentPeakMatchRelativeMassDeviation)
  #   results <- data.frame()
  #   results <- apply((MS3$MS1), 1, F_Identification)
  #   results_dt <- ldply(results, data.frame) # Transfer results from list to data.frame
  #   identifications_lipidmaps <- data.table(results_dt) #
  #   identifications_lipidmaps <<- unique(identifications_lipidmaps, by = c("PrecuresorMZ", "RentionTime", "NeutralPrecursorMass", "adduct"))
  #   identifications_lipidmaps <<- identifications_lipidmaps[!"+", on = .(Null)] # filterout the null candadates
  #   # identifications_lipidmaps$RentionTime<-identifications_lipidmaps$RentionTime/60
  #   setorder(identifications_lipidmaps, -RentionTime, -Score) #
  #   print("MetFrag Identification Done!")
  # })
  lo_Identification[6, 1] <- glabel("FDR")
  lo_Identification[6, 2] <- ed_FDR <- gedit(text = "0.05", width = 4)
  tooltip(ed_FDR) <- "The false discovery rate used for quality control, the value should be positive larger than 0.001 
  and are usually set as 0.05 (5% FDR) or 0.01 (1% FDR)"
  # lo_Identification[6, 1:2] <- ed_FDR <- gradio(c('0.1%','1%','5% FDR'), horizontal=TRUE, width = 4)
  # outputRadio <- gradio(c("HTML", "LATEX", "TEXT"), horizontal=TRUE, container=reportDialogOutFrame)
  lo_Identification[7, 1:2] <- gbutton("Visualize target decoy", handler = function(h, ...) {
    print("Visualize target decoy Done!")
  })
  
  lo_Identification[8, 1:2] <- gbutton("Export results", handler = function(h, ...) {
    write.csv(result,'Experimental Identification.csv')
    # write.csv(identifications_lipidmaps, "identifications_Metfrag.csv")
    print("Export results Done!")
  })
  ## Set Parameters for MS-Pathway method
  #### Statistics ####
  lo_ImportData[2, 3] <- gf_Statistics <- gframe("6 Statistics", horizontal = F, expand = TRUE)
  lo_Statistics <- glayout(cont = gf_Statistics, horizontal = F)
  ## Set Filteration gframe in FilterNormalization tab
  lo_Statistics[1, 1:2] <- widgets$bt_SanityCheck <- gbutton("Sanity check", handler = function(h, ...) {
    ### Export the feature table in the MetaboAnalyst format. Parameter 'label', defines the group assignment of the samples.
    exportMetaboAnalyst(xdata1, file = "met_test1.csv", label = xdata1$sample_group)
    # First step is to create the mSet Object, specifying that the data to be uploaded
    mSet <<- InitDataObjects("pktable", "stat", FALSE) # is a peak table ("pktable") and that statistical analysis will be performed ("stat").
    # The second step is to read in the processed data (created above)
    mSet <<- Read.TextData(mSet, "met_test1.csv", "colu", "disc")
    mSet <<- SanityCheckData(mSet)
    print("Sanity Check done!")
  })
  tooltip(widgets$bt_SanityCheck) <- "SanityCheckData is used for data processing, and performs a basic sanity check of the uploaded content, ensuring that the data is suitable for further analysis. The function will return a message if the data has successfully passed the check and is deemed suitable for further analysis. If it fails, the function will return a 0. The function will perform the check directly onto the mSet$dataSet object, and must be performed immediately after reading in data. The sanity check function evaluates the accuracy of sample and class labels, data structure, deals with non-numeric values, removes columns that are constant across all samples (variance = 0), and by default replaces missing values with half of the original minimal positive value in your dataset.Before data analysis, a data integrity check is performed to make sure that all the necessary information
           has been collected. The class labels must be present and contain only two classes. If samples are paired,
           the class label must be from -n/2 to -1 for one group, and 1 to n/2 for the other group (n is the sample
           number and must be an even number). Class labels with same absolute value are assumed to be pairs.
           Compound concentration or peak intensity values should all be non-negative numbers. By default, all
           missing values, zeros and negative values will be replaced by the half of the minimum positive value
           found within the data
           "
  lo_Statistics[2, 1:2] <- widgets$bt_ReplaceMissingValue <- gbutton("Replace missing value", handler = function(h, ...) {
    mSet <<- ReplaceMin(mSet)
    print("Replace Missing Value done!")
  })
  tooltip(widgets$bt_ReplaceMissingValue) <- "Too many zeroes or missing values will cause difficulties for downstream analysis.The default method replaces all the missing and zero values
           with a small values (the half of the minimum positive values in the original data) assuming to be
           the detection limit. The assumption of this approach is that most missing values are caused by low
           abundance metabolites (i.e.below the detection limit)."

  lo_Statistics[3, 1:2] <- widgets$bt_FilterVariable <- gbutton("Filter variable", handler = function(h, ...) {
    mSet <<- FilterVariable(mSet, "iqr", "F", 25)
    print("Filter Variable done!")
  })
  tooltip(widgets$bt_FilterVariable) <- "The purpose of the data filtering is to identify and remove variables that are unlikely to be of use
           when modeling the data. No phenotype information are used in the filtering process, so the result
           can be used with any downstream analysis. This step can usually improves the results. Data filter is
           strongly recommended for datasets with large number of variables (> 250) datasets contain much noise
           (i.e.chemometrics data).The function applies a filtering method, ranks the variables within the dataset, and removes variables based on its rank. The final dataset should contain no more than than 5000 variables for effective computing.Here we use iqr:interquantile for filter option; F: do not use QC samples; 25 is the relative standard deviation cut-off
           "
  ## Set Normalization gframe in FilterNormalization tab
  lo_Statistics[4, 1] <- glabel("Normalization")
  lo_Statistics[4, 2] <- widgets$cb_rowNorm <- gcombobox(c("QuantileNorm", "ProbNormT", "ProbNormF", "CompNorm", "SumNorm", "MedianNorm", "SpecNorm"), expand = T)
  tooltip(widgets$cb_rowNorm) <- "Select the option for row-wise normalization, QuantileNorm for Quantile Normalization, ProbNormT for Probabilistic Quotient Normalization without using a reference sample, ProbNormF for Probabilistic Quotient Normalization based on a reference sample, CompNorm for Normalization by a reference feature, SumNorm for Normalization to constant sum, MedianNorm for Normalization to sample median, and SpecNorm for Normalization by a sample-specific factor."

  lo_Statistics[5, 1] <- glabel("Transform")
  lo_Statistics[5, 2] <- widgets$cb_transNorm <- gcombobox(c("LogNorm", "CrNorm"),)
  tooltip(widgets$cb_transNorm) <- "Select option to transform the data, LogNorm for Log Normalization, and CrNorm for Cubic Root Transformation."

  lo_Statistics[6, 1] <- glabel("Scaling")
  lo_Statistics[6, 2] <- widgets$cb_scaleNorm <- gcombobox(c("MeanCenter", "AutoNorm", "ParetoNorm", "RangeNorm"))
  tooltip(widgets$cb_scaleNorm) <- "Select option for scaling the data, MeanCenter for Mean Centering(mean-centered only), AutoNorm for Autoscaling(mean-centered and divided by standard deviation of each variable), ParetoNorm for Pareto Scaling(mean-centered and divided by the square root of standard deviation of each
           variable), amd RangeNorm for Range Scaling(mean-centered and divided by the value range of each variable)."

  lo_Statistics[7, 1:2] <- widgets$bt_Normalization <- gbutton("Run normalization",  handler = function(h, ...) {
    mSet <<- PreparePrenormData(mSet)
    mSet <<- Normalization(mSet, svalue(widgets$cb_rowNorm), svalue(widgets$cb_transNorm), svalue(widgets$cb_scaleNorm), ratio = FALSE, ratioNum = 20)
    print("Normalization done!")
  })
  tooltip(widgets$bt_Normalization) <- "The data is stored as a table with one sample per row and one variable (bin/peak/metabolite) per
           column. The normalization procedures implemented below are grouped into four categories. Sample
           specific normalization allows users to manually adjust concentrations based on biological inputs (i.e.
           volume, mass); row-wise normalization allows general-purpose adjustment for differences among samples;
           data transformation and scaling are two different approaches to make features more comparable.This function performs row-wise normalization, transformation, and scaling of your metabolomic data.
           "

  lo_Statistics[8, 1:2] <- widgets$bt_PlotNormalization <- gbutton("Show normalization effect", handler = function(h, ...) {
    mSet <<- PlotNormSummary(mSet, "norm_0_", "png", 72, width = NA)
    mSet <<- PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width = NA)
    # F_Delete_Ggraphics()
    svg(paste0('Normalization Density','.svg'))
    par(mar = c(0, 0, 0, 0))
    plot.new()
    rasterImage(load.image(mSet[["imgSet"]][["summary_norm"]]), 0, 0, 0.5, 1)
    rasterImage(load.image(mSet[["imgSet"]][["norm"]]), 0.5, 0, 1, 1)
    dev.off()
    print("Normalization Density done!")
  })
  tooltip(widgets$bt_PlotNormalization) <- "To see the Normalization Density Effect."
  lo_Statistics[9, 1:2] <- widgets$bt_Mummichog <- gbutton("Mummichog analysis", handler = function(h, ...) {
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
  })
  tooltip(widgets$bt_PlotNormalization) <- "mummichog pathway enrichment"
  
  lo_Statistics[10, 1:2] <- widgets$bt_PathwayNetwork <- gbutton("Pathway network", handler = function(h, ...) {
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
  })
  tooltip(widgets$bt_PathwayNetwork) <- "Network view of meta-analysis"
  
  

  ##Set Volcano button in Statistics tab
  lo_Statistics[1,3]<-"Fold change"
  lo_Statistics[1,4]<-ed_FCThreshold<-gedit(text = "2", width = 4)  
  tooltip(ed_FCThreshold) <- "Calculate the fold change between WT and KO groups, each point represents a peak. 
                                    The default Fold change threshold is set as 2"
  lo_Statistics[2,3]<-"Non-parametric"
  lo_Statistics[2,4:5]<-cb_nonpar<-gcombobox(c("FALSE","TRUE"))
  
  lo_Statistics[3,3]<-"p-value cutoff"
  lo_Statistics[3,4]<-ed_threshp<-gedit(text = "0.05", width = 4) 
  
  lo_Statistics[4,3]<-"paired"
  lo_Statistics[4,4:5]<-cb_paired<-gcombobox(c("FALSE","TRUE"))
  tooltip(ed_threshp) <- "perform t-test analysis. This univariate analyses provide a preliminary overview about
           features that are potentially significant in discriminating the conditions under study.For paired fold change analysis, the algorithm rst counts the total number of pairs with fold changes
           that are consistently above/below the specified FC threshold for each variable. A variable will be reported
           as significant if this number is above a given count threshold (default > 75% of pairs/variable)"
  lo_Statistics[5,3:4] <- gbutton("Run volcano plot", handler = function(h,...) {
    mSet<<-Volcano.Anal(mSet, fcthresh=as.numeric(svalue(ed_FCThreshold)),
                        nonpar=as.logical(svalue(cb_nonpar)), 
                        threshp=as.numeric(svalue(ed_threshp)),
                        paired=as.logical(svalue(cb_paired)), 
                        cmpType=1,equal.var=TRUE, pval.type="raw")
    mSet<<-PlotVolcano(mSet, "Volcano", 0,"png", dpi=72, width=NA)
    # F_Delete_Ggraphics()
    svg(paste0('Volcano Plot','.svg'))
    par(mar=c(0,0,0,0)) 
    plot.new()
    #plot(1:10,main="volcano",ty="n",xlab="",ylab="")
    rasterImage(load.image(mSet[["imgSet"]][["volcano"]]),0,0,1,1)
    dev.off()
    print("Volcano Plot done!")
  })
  ##Set PCA gframe in Statistics tab
  lo_Statistics[6,3:4]<- widgets$bt_PCA <- gbutton("Run PCA", handler = function(h,...) {
    mSet<<-PCA.Anal(mSet)
    print("Principle Component Analysis done!")
  })
  tooltip(widgets$bt_PCA) <- "PCA is an unsupervised method aiming to find the directions that best explain the variance in a data
           set (X) without referring to class labels (Y). The data are summarized into much fewer variables called
           scores which are weighted average of the original variables. The weighting profiles are called loadings.
           The PCA analysis is performed using the prcomp package. The calculation is based on singular value
           decomposition."
  
  lo_Statistics[7,3:4]<- widgets$bt_PlotPCASummary <- gbutton("Plot PCA summary", handler = function(h,...) {
    mSet<<-PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
    mSet<<-PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
    mSet<<-PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
    mSet<<-PlotPCALoading(mSet, "pca_loading_0_", format="png", dpi=72, width=NA, 1,2)
    mSet<<-PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1,2)
    mSet<<-PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
    # F_Delete_Ggraphics()
    svg(paste0('PCA Summary','.svg'))
    par(mar=c(0,0,0,0)) 
    plot.new()
    #plot(1:10,main="PCA Summary",ty="n",xlab="",ylab="")
    rasterImage(load.image(mSet[["imgSet"]][["pca.pair"]]),0,0.5,0.33,1)
    rasterImage(load.image(mSet[["imgSet"]][["pca.scree"]]),0.33,0.5,0.66,1)
    rasterImage(load.image(mSet[["imgSet"]][["pca.score2d"]]),0.66,0.5,0.99,1)
    rasterImage(load.image(mSet[["imgSet"]][["pca.loading"]]),0,0,0.33,0.5)
    rasterImage(load.image(mSet[["imgSet"]][["pca.biplot"]]),0.33,0,0.66,0.5)
    rasterImage(load.image(mSet[["imgSet"]][["pca.score3d"]]),0.66,0,0.99,0.5)
    dev.off()
    print("PCA Summary done!")
  })
  tooltip(widgets$bt_PlotPCASummary) <- "pairwise score plots providing an overview of
           the various seperation patterns among the most significant PCs;"
  
  ##Set PLSDA gframe in Statistics tab
  lo_Statistics[8,3:4]<- widgets$bt_PLSDA <- gbutton("Run PLS-DA", handler = function(h,...) {
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
  tooltip(widgets$bt_PLSDA) <- "PLS is a supervised method that uses multivariate regression techniques to extract via linear combination
           of original variables (X) the information that can predict the class membership (Y).To assess the significance of class discrimination, a permutation test was performed. In each permu-
           tation, a PLS-DA model was built between the data (X) and the permuted class labels (Y) using the
           optimal number of components determined by cross validation for the model based on the original class
           assignment. MetaboAnalyst supports two types of test statistics for measuring the class discrimination.
           The first one is based on prediction accuracy during training. The second one is separation distance
           based on the ratio of the between group sum of the squares and the within group sum of squares (B/W-
           ratio). If the observed test statistic is part of the distribution based on the permuted class assignments,
           the class discrimination cannot be considered significant from a statistical point of view.6."
  
  lo_Statistics[9,3:4]<- widgets$bt_PlotPLSSummary <-gbutton("Plot PLS-DA Summary", handler = function(h,...) {
    # F_Delete_Ggraphics()
    svg(paste0('PLS Summary','.svg'))
    par(mar=c(0,0,0,0)) 
    plot.new()
    rasterImage(load.image(mSet[["imgSet"]][["pls.pair"]]),0,0.5,0.33,1)
    rasterImage(load.image(mSet[["imgSet"]][["pls.class"]]),0.33,0.5,0.66,1)
    rasterImage(load.image(mSet[["imgSet"]][["pls.score2d"]]),0.66,0.5,0.99,1)
    rasterImage(load.image(mSet[["imgSet"]][["pls.loading"]]),0,0,0.33,0.5)
    rasterImage(load.image(mSet[["imgSet"]][["pls.imp"]]),0.33,0,0.66,0.5)
    rasterImage(load.image(mSet[["imgSet"]][["pls.score3d"]]),0.66,0,0.99,0.5)
    dev.off()
    print("PLS Summary done!")
  })
  tooltip(widgets$bt_PlotPLSSummary) <- "Plot PLS Summary"
  
  ##Set PreparePDFReport gframe in Statistics tab
  # lo_Statistics[10,1]<- "Report name"
  # lo_Statistics[10,2]<- ed_ReportName<-gedit(text = "MetGUI", width = 4) 
  lo_Statistics[10,3:4]<- gbutton("Export report", handler = function(h,...) {
    PreparePDFReport(mSet, "MetGUI") #create a summary report of the statistical analysis 
    print("PreparePDFReport done!")
  })

  #### About####
  gp_About <- ggroup(horizontal = FALSE, cont = nb, label = "About")
  # addSpring(gp_About)
  gtext("MetGUI is a function set with Gui to process mass spectrometry data", font.attr = list(size = "xx-large", family = "monospace"), expand = TRUE, fill = TRUE, container = gp_About)
  #### Plot tab####
  # gg <<- ggraphics(cont=nb, label='Plot',visible=FALSE)
  # gp_Plot= ggroup(horizontal = FALSE, cont = nb, label = "Plot")
  # F_Delete_Ggraphics<-function(){
  #   delete(nb,gg)
  #   #dev.off()
  #   gg <<- ggraphics(cont=nb, label='Plot',visible=FALSE)
  # }
  # par(mar=c(0,0,0,0)) # To solve the Error in plot.new() : figure margins too large
  # plot.new()
  # rasterImage(imager::load.image("extdata/Logo/Logo.png"),0,0,1,1)
}
