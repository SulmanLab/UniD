#' Load raw data
#'
#' @description Used to load the raw data and return a list.
#'
#' @param dataDir directory where raw data stored
#' @param outDir directory where output data should be saved if write = T
#' @param arrayType the platform which raw data generated, can be "450k" or "EPIC"
#' @param write whether the output should be saved, highly recommended
#' @return A list which includes rgSet, Mset, and detP

#' @examples
#' \dontrun{
#' loading <- UniD.load(dataDir = "~/Desktop/IDAT/", outDir = "~/Desktop/output/",
#' arrayType = "450k", write = T)
#'
#' loading <- UniD.load(dataDir = "~/Desktop/IDAT/", outDir = NULL,
#' arrayType = "EPIC", write = F)
#' }
#' @importFrom minfi read.metharray.sheet read.metharray.exp detectionP preprocessRaw
#' @import minfi
#' @import IlluminaHumanMethylation450kmanifest
#' @import IlluminaHumanMethylationEPICmanifest
#' @importFrom utils packageVersion

#' @export
UniD.load <- function(dataDir,
                      outDir,
                      arrayType,
                      write)
{

  message("")
  message("====Data Loading Start====")
  message("Package version of UniD is: ", packageVersion("UniD"))
  message("Loading data from \"", dataDir, "\"")

  targets <- read.metharray.sheet(dataDir)  #the sample sheet has specific format
  rgSet <- read.metharray.exp(targets = targets, extended = T)

  # set annotation based on array type
  if (arrayType == "450k")
    rgSet@annotation <- c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")
  if (arrayType == "EPIC")
    rgSet@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b3.hg19")

  detP <- detectionP(rgSet)
  Mset <- preprocessRaw(rgSet)

  loading <- list()
  loading[["rgSet"]] <- rgSet
  loading[["Mset"]] <- Mset
  loading[["detP"]] <- detP

  # write loadings out
  if (write) {
    save(loading, file = paste0(outDir, "/loading.RData"))
    message(paste0("loading(rgSet, Mset, detP) saved as: ", outDir, "/loading.RData"))
  }

  message("====Data Loading Finsihed====")
  message("")
  return(loading)
}
