#' Checking quality of internal control probes
#'
#' @description Used to check the quality of internal control
#' probes and return a data frame which includes whether
#' samples have failed on different categories of internal controls
#'
#' @inheritParams UniD.load
#' @param sampleType Indicate the input data are "FFPE" or "other"
#' @param loading List generated from the \code{UniD.load()}
#' @return a data frame with quality checking results for internal
#' control probes
#' @examples
#' \dontrun{
#' samQC <- UniD.intctl(loading = loading, dataDir = "~/Desktop/IDAT/",
#' outDir = "~/Desktop/output/", arrayType = "450k", sampleType = "FFPE",
#' write = T)
#'
#' samQC <- UniD.intctl(loading = loading, dataDir = "~/Desktop/IDAT/",
#' outDir = NULL, arrayType = "EPIC", sampleType = "other", write = F)
#' }
#' @importFrom minfi getRed getGreen
#' @importFrom utils write.table
#' @export
UniD.intctl <- function(loading = loading,
                        dataDir,
                        outDir,
                        arrayType,
                        sampleType,
                        write)
{

  message("")
  message("====Internal Control Checking Start====")

  #check if loading exist
  if(exists("loading") == FALSE){

    if(file.exists(paste0(outDir,"/loading.RData")) == TRUE){
      #if loading not exist in current working directory, load it from outDir
      message("#Load RData from \"", outDir, "\"")
      load(paste0(outDir,"/loading.RData"))
      message("#Load RData finished")
    }else{
      #no loading in current environment nor outDir, read raw data from dataDir
      message("Read raw data from \"", dataDir, "\" using UniD.load.R")
      loading <- UniD.load(dataDir, outDir, arrayType, write)
    }
  }


  r <- getRed(loading$rgSet)
  g <- getGreen(loading$rgSet)
  sampleSet <- GetsampleSet(loading$rgSet)
  sentrix_barcode <- Getsentrix_barcode(loading$rgSet)
  sentrix_position <- Getsentrix_position(loading$rgSet)

  ####Generate Internal controls####
  #####Internal controls (not distinguish arrayType)
  #restoration
  res.ctl <- g["41636384",]

  #staining
  biotin_high <- g["41666334",]
  biotin_bkg <- g["34648333",]
  staining_green <- biotin_high/biotin_bkg   #>=5

  DNP_high <- r["27630314",]
  DNP_bkg <- r["43603326",]
  staining_red <- DNP_high/DNP_bkg    #>=5

  ####Internal controls (distinguish arrayType) generation
  if(arrayType == "EPIC"){

    #extension green or red as bkg
    bkg_g_AT_max <- apply(g[c("63642461","21752326"),], 2, function(x) max(x))
    bkg_r_CG_max <- apply(r[c("12719506","74666473"),] , 2, function(x) max(x))

    #restoration, if FFPE samples, threshold = 1; others as 0
    Restoration_green <- (res.ctl/(bkg_g_AT_max + 3000))

    #extension
    extension_green <- apply(g[c("12719506","74666473"),] , 2,
                             function(x) min(x))/bkg_g_AT_max  #>=5
    extension_red   <- apply(r[c("63642461","21752326"),], 2,
                             function(x) min(x))/bkg_r_CG_max   #>=5

    #hybridization
    hybridization_green_HM <- g["28684356",]/g["39782321",]  #>=1
    hybridization_green_ML <- g["39782321",]/g["21771417",]  #>=1

    #target removal
    targetRemoval_green1 <- (bkg_g_AT_max+3000)/g["39773404",]  #>=1
    targetRemoval_green2 <- (bkg_g_AT_max+3000)/g["42790394",]  #>=1

    #bisulfite conversion
    bisulfite_conversion_I_green1 <- apply( g[c("22795447","56682500"),], 2,
                                            function(x) min(x))/
      apply( g[c("24637490","33665449"),], 2, function(x) max(x))  #>=1
    bisulfite_conversion_I_green2 <- (bkg_g_AT_max + 3000)/
      apply( g[c("24637490","33665449"),], 2, function(x) max(x)) #>=1

    bisulfite_conversion_I_red1 <- apply( r[c("54705438","49720470","26725400")
                                            ,], 2, function(x) min(x))/
      apply( r[c("57693375","15700381","33635504"),], 2, function(x) max(x))  #>=1
    bisulfite_conversion_I_red2 <- (bkg_r_CG_max + 3000)/
      apply( r[c("57693375","15700381","33635504"),], 2, function(x) max(x))  #>=1

    bisulfite_conversion_II_1 <- apply( r[c("43720395","70664314","71718498","12722428"),],
                                        2, function(x) min(x))/
      apply( g[c("43720395","70664314","71718498","12722428"),], 2, function(x) max(x)) #>=1
    bisulfite_conversion_II_2 <- (bkg_g_AT_max + 3000)/
      apply( g[c("43720395","70664314","71718498","12722428"),], 2, function(x) max(x)) #>=1

    #specificity
    specificity_I_green <- apply( g[c("65735497","51804467","61624401"),], 2,
                                  function(x) min(x))/
      apply( g[c("50611399","30684480","28649458"),], 2, function(x) max(x))  #>=1
    specificity_I_red   <- apply( r[c("46779338","59783305","28618334"),], 2,
                                  function(x) min(x))/
      apply( r[c("51745378","65797428","52712334"),], 2, function(x) max(x))   #>=1
    specificity_II_1    <- apply( r[c("29662396","17661470","34730329"),], 2,
                                  function(x) min(x))/
      apply( g[c("29662396","17661470","34730329"),], 2, function(x) max(x))   #>=1
    specificity_II_2    <- (bkg_g_AT_max + 3000)/
      apply( g[c("29662396","17661470","34730329"),], 2, function(x) max(x))  #>=1

    #nonploymorphic
    nonpolymorphic_green <- apply( g[c("23663352","38796356"),], 2,
                                   function(x) min(x))/
      apply( g[c("24701411","18773482"),], 2, function(x) max(x))   #>=5
    nonpolymorphic_red   <- apply( r[c("24701411","18773482"),], 2,
                                   function(x) min(x))/
      apply( r[c("23663352","38796356"),], 2, function(x) max(x))   #>=5

  }else if(arrayType == "450k"){

    #extension green or red as bkg
    bkg_g_AT_max <- apply( g[c("63642461","47640365"),], 2, function(x) max(x))
    bkg_r_CG_max <- apply( r[c("31698466","74666473"),], 2, function(x) max(x))

    #restoration, distinguish on sampleType
    Restoration_green <- res.ctl/(bkg_g_AT_max + 3000)

    #extension
    extension_green <- apply( g[c("31698466","74666473"),], 2, function(x) min(x))/bkg_g_AT_max  #>=5
    extension_red   <- apply( r[c("63642461","47640365"),], 2, function(x) min(x))/bkg_r_CG_max   #>=5

    #hybridization
    hybridization_green_HM <- g["28684356",] / g["26772442",]  #>=1
    hybridization_green_ML <- g["26772442",] / g["21771417",]  #>=1

    #target removal
    targetRemoval_green1 <- (bkg_g_AT_max+3000) / g["13643320",]  #>=1
    targetRemoval_green2 <- (bkg_g_AT_max+3000) / g["42790394",]  #>=1

    #bisulfite conversion
    bisulfite_conversion_I_green1 <- apply( g[c("22711390","56682500","22795447"),],
                                            2, function(x) min(x))/
      apply( g[c("24637490","33665449","46651360"),], 2, function(x) max(x))  #>=1
    bisulfite_conversion_I_green2 <- (bkg_g_AT_max + 3000)/
      apply( g[c("24637490","33665449","46651360"),], 2, function(x) max(x)) #>=1

    bisulfite_conversion_I_red1 <- apply( r[c("54705438","49720470","26725400"),],
                                          2, function(x) min(x))/
      apply(r[c("57693375","15700381","33635504"),], 2, function(x) max(x))   #>=1
    bisulfite_conversion_I_red2 <- (bkg_r_CG_max + 3000)/
      apply( r[c("57693375","15700381","33635504"),], 2, function(x) max(x))  #>=1

    bisulfite_conversion_II_1 <- apply( r[c("43720395","70664314","71718498","30724412"),],
                                        2, function(x) min(x))/
      apply(g[c("43720395","70664314","71718498","12722428"),], 2, function(x) max(x))  #>=1
    bisulfite_conversion_II_2 <- (bkg_g_AT_max + 3000)/
      apply( g[c("43720395","70664314","71718498","30724412"),], 2, function(x) max(x))  #>=1

    #specificity
    specificity_I_green <- apply(g[c("23777311","51804467","10673427"),],2,
                                 function(x) min(x))/
      apply(g[c("67672371","30684480","58661465"),],2, function(x) max(x))  #>=1
    specificity_I_red   <- apply(r[c("46779338","59783305","53740460"),],2,
                                 function(x) min(x))/
      apply(r[c("51745378","65797428","12808347"),],2, function(x) max(x))   #>=1
    specificity_II_1    <- apply(r[c("29662396","17661470","34730329"),],2,
                                 function(x) min(x))/
      apply(g[c("29662396","17661470","34730329"),],2, function(x) max(x))   #>=1
    specificity_II_2    <- (bkg_g_AT_max + 3000)/
      apply( g[c("29662396","17661470","34730329"),], 2, function(x) max(x))  #>=1

    #nonploymorphic
    nonpolymorphic_green <- apply(g[c("23663352","70645401"),],2, function(x) min(x))/
      apply(g[c("24701411","18773482"),],2, function(x) max(x))   #>=5
    nonpolymorphic_red   <- apply(r[c("24701411","18773482"),],2, function(x) min(x))/
      apply(r[c("23663352","70645401"),],2, function(x) max(x))   #>=5

  }

  #generate sample internal control data summary
  samQ <- data.frame("Sentrix_Barcode" = sentrix_barcode,
                     "Sentrix_Position" = sentrix_position,
                     Restoration_green,
                     staining_green,
                     staining_red,
                     extension_green,
                     extension_red,
                     hybridization_green_HM,
                     hybridization_green_ML,
                     targetRemoval_green1,
                     targetRemoval_green2,
                     bisulfite_conversion_I_green1,
                     bisulfite_conversion_I_green2,
                     bisulfite_conversion_I_red1,
                     bisulfite_conversion_I_red2,
                     bisulfite_conversion_II_1,
                     bisulfite_conversion_II_2,
                     specificity_I_green,
                     specificity_I_red,
                     specificity_II_1,
                     specificity_II_2,
                     nonpolymorphic_green,
                     nonpolymorphic_red)

  rownames(samQ) <- sampleSet

  ####Evaluate Internal controls####
  samQC <- samQ[, -c(1,2)]

  #restoration
  if(sampleType == "FFPE"){
    samQC$Restoration_green[samQC$Restoration_green <= 1] <- "Failed"
  }else if(sampleType == "other"){
    samQC$Restoration_green[samQC$Restoration_green <= 0] <- "Failed"
  }

  #staining
  samQC$staining_green[samQC$staining_green <=5 ] <- "Failed"
  samQC$staining_red[samQC$staining_red     <=5 ] <- "Failed"


  #extension
  samQC$extension_green[samQC$extension_green <=5 ] <- "Failed"
  samQC$extension_red[samQC$extension_red     <=5 ] <- "Failed"

  #hybridization
  samQC$hybridization_green_HM[samQC$hybridization_green_HM <= 1] <- "Failed"
  samQC$hybridization_green_ML[samQC$hybridization_green_ML <= 1] <- "Failed"

  #target removal
  samQC$targetRemoval_green1[samQC$targetRemoval_green1 <=1 ] <- "Failed"
  samQC$targetRemoval_green2[samQC$targetRemoval_green2 <=1 ] <- "Failed"

  #bisulfite conversion
  samQC$bisulfite_conversion_I_green1[samQC$bisulfite_conversion_I_green1 <= 1] <- "Failed"
  samQC$bisulfite_conversion_I_green2[samQC$bisulfite_conversion_I_green2 <= 1] <- "Failed"
  samQC$bisulfite_conversion_I_red1[samQC$bisulfite_conversion_I_red1     <= 1] <- "Failed"
  samQC$bisulfite_conversion_I_red2[samQC$bisulfite_conversion_I_red2     <= 1] <- "Failed"
  samQC$bisulfite_conversion_II_1[samQC$bisulfite_conversion_II_1         <= 1] <- "Failed"
  samQC$bisulfite_conversion_II_2[samQC$bisulfite_conversion_II_2         <= 1] <- "Failed"

  #specificity
  samQC$specificity_I_green[samQC$specificity_I_green <= 1] <- "Failed"
  samQC$specificity_I_red[samQC$specificity_I_red     <= 1] <- "Failed"
  samQC$specificity_II_1[samQC$specificity_II_1       <= 1] <- "Failed"
  samQC$specificity_II_2[samQC$specificity_II_2       <= 1] <- "Failed"

  #non-polymorphic
  samQC$nonpolymorphic_green[samQC$nonpolymorphic_green <= 5] <- "Failed"
  samQC$nonpolymorphic_red[samQC$nonpolymorphic_red     <= 5] <- "Failed"

  samQC[samQC != "Failed" ] <- ""

  #print out samples with suspect Internal control issue
  message("Failed sample in Internal Control:")
  message("")
  for(i in 1:ncol(samQC)){
    if( length( samQC[,i][samQC[,i] == "NA"]) != nrow(samQC) ){
      message(paste0( colnames(samQC)[i] , " : ", paste(rownames(samQC)[
        which(samQC[,i] == "Failed")], collapse = "; ")))
    }
  }
  message("")

  #write data out
  if(write){
    write.table(samQC, file = paste0(outDir, "/samQC_Internal.Control.csv"),
                sep = ",", col.names = NA, quote = F)
    message("Internal Controls checking results saved as: ", outDir,
            "/samQC_Internal.Control.csv")
  }

  message("====Internal Control Checking Finished====")
  message("")
  return(samQC)
}
