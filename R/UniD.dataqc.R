#' Probe quality evaluation
#'
#' @description Used to evaluate the probe quality in term of
#' probe detection p-value and beadcount. For probes which failed the
#' cut point, they will be set as missing value. The data with missing
#' value will be returned/saved as beta value.
#'
#' @param loading List generated from the \code{UniD.load()}
#' @param detP.cut The cut point for detection P-value, default is 0.05.
#' If the probe's detection p-value > 0.05, the probe will be set as missing
#' value in a sample-wise fashion.
#' @param bc.cut The cut point for beadcount per probe, default is 3.
#' If beadcount < 3, the probe will be set as missing value sample-wise.
#' control probes
#' @inheritParams UniD.load
#' @return A data frame with beta values for all samples and probes.
#' @examples
#' \dontrun{
#' Beta.raw <- UniD.dataqc(loading = loading, outDir = "~/Desktop/output/",
#' detP.cut = 0.05, bc.cut = 3, arrayType = "450k", write = T)
#'
#' Beta.raw <- UniD.dataqc(loading = loading,  outDir = NULL, detP.cut = 0.05,
#' bc.cut = 3, arrayType = "450k", write = F)
#' }
#' @importFrom minfi getBeta
#' @import wateRmelon
#' @importFrom utils write.table


#' @export
UniD.dataqc <- function(loading = loading,
                        outDir,
                        detP.cut = 0.05,
                        bc.cut = 3,
                        arrayType,
                        write
)
{

  message("")
  message("====Sample-wise Probe Quality Assessment Start====")

  #check if loading exist
  if(exists("loading") == FALSE){
    stop("loading object is missing")
  }

  detP <- loading$detP
  sampleSet <- GetsampleSet(loading$rgSet)

  Beta.raw <- getBeta(loading$Mset)   #Here, Beta value only
  Beta.raw <- as.data.frame(Beta.raw)

  bc <- beadcount(loading$rgSet) #get beadcount data from rgSet
  bc[bc < bc.cut ] <- NA  #set NA

  Beta.raw[detP > detP.cut] <- NA
  Beta.raw[is.na(bc)] <- NA

  #assess probe quality in detection P-value and beadcount number
  p.fail <- data.frame(Fail.Frac.detP = apply(detP, 2, function(x)
    length(which(x > detP.cut)))/nrow(detP),
                       Fail.Frac.beadcount  = colSums(is.na(bc))/nrow(bc),
                       Fail.Frac.NA = colSums(is.na(Beta.raw))/nrow(Beta.raw))

  print("The failed fraction per sample (failed detP and bc may overlap): ")
  print(p.fail)

  if(write){
    save(Beta.raw, file = paste0(outDir, "/UniD_Beta.raw.RData"))
    message("The raw Beta value were saved as: ", outDir,
            "/UniD_Beta.raw.RData")
    write.table(p.fail, file = paste0(outDir, "/Failed_probe_Fraction.csv"),
                sep = ",", col.names = NA, quote = F)
    message(paste0("The fraction of failed probes per sample saved as: ",
                   outDir, "/Failed_probe_Fraction.csv"))
  }

  message("====Sample-wise Probe Quality Assessment Finished====")
  message("")

  return(Beta.raw)
}
