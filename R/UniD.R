#' UniD: A package for unified diagnostic platform for gliomas.
#'
#' The UniD package use the DNA methylation microarray data to predict other
#' genomic information, including IDH mutation, ATRX mutation, TERT promoter
#' mutation, chromosome 1p/19q co-deletion, TCGA gene expression subtypes (
#' Classical, Mesenchymal, and Proneural) with their probability, and MGMT
#' gene expression level.
#'
#' UniD package providing two categories of functions:
#'
#' (1) Data processing and quality control:
#'
#'       \code{UniD.intctl}, \code{UniD.dataqc}, \code{UniD.probefilter}
#'
#' (2) Genomic information prediction:
#'
#'       \code{UniD.pred}
#'
#' @section Functions:
#' \code{UniD.load}: load the raw data
#'
#' \code{UniD.intctl}: check the internal control probes which can reflect
#' different aspects of experiments quality.
#'
#' \code{UniD.dataqc}: check data quality in three aspects: detection P-value,
#'  beadcount per probes, and missing values.
#'
#' \code{UniD.probefilter}: filter out probes belong to different categories
#' which may have quality issues, and probes with high proportion of missing
#' values.
#'
#' \code{UniD.BMIQ}: normalize type I and type II probes based on BMIQ method.
#'
#' \code{UniD.pred}: predict othergenomic information using Illumina DNA methylation
#' microarray data, currently support 450k and EPIC platform
#'
#' @docType package
#' @name UniD
#' @seealso The BMIQ normalization method is adapted from \code{\link[wateRmelon]{BMIQ}}.
#' @inheritParams UniD.load
#' @inheritParams UniD.intctl
#' @return A data frame with all predicted results
#' @examples
#' \dontrun{
#'  pred.result <- UniD(dataDir = "~/Desktop/input/", outDir = "~/Desktop/output/",
#'  arrayType = "450k", sampleType = "other", write = T)
#' }

#' @export
UniD <- function(dataDir,
                 outDir,
                 arrayType,
                 sampleType,
                 write)
{
  loading <- UniD.load(dataDir = dataDir, outDir = outDir, arrayType = arrayType,
                       write = write)
  samQC <- UniD.intctl(loading = loading, dataDir = dataDir,outDir = outDir,
                       arrayType = arrayType, sampleType = sampleType,
                       write = write)
  Beta.raw <- UniD.dataqc(loading = loading, outDir = outDir,
                          arrayType = arrayType, write = write)
  Beta.BMIQ <- UniD.BMIQ(Beta.raw, outDir = outDir,
                         arrayType = arrayType, write = write)
  Beta.clean <- UniD.probefilter(Beta.raw, outDir = outDir,
                                 arrayType = arrayType, write = T)
  Pred <- UniD.pred(inputdata = Beta.raw, inputvalueType = "B",
                    inputdata.BMIQ = Beta.BMIQ,  outDir = outDir, write = write)
  return(Pred)
}
