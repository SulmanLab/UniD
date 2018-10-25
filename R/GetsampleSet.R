#' Get sample Set from rgSet
#'
#' @description \code{GetsampleSet} is used to extract the samples set from
#' the rgSet. \code{Getsentrix_barcode} is used to extract sentrix barcode from
#' rgSet. \code{Getsentrix_position} is used to extract sentrix position from rgSet
#'
#' @param rgSet The RGChannelSetExtended with in the list \code{loading}
#'
#' @return \code{sampleSet} sample set of the raw data
#' @return \code{sentrix_barcode} the sentrix barcode of samples
#' @return \code{sentrix_position} the sentrix position of samples

#'@export
GetsampleSet <- function(rgSet){
  sampleSet <- rgSet@colData@listData$Sample_Name
  return(sampleSet)
}

#' @rdname GetsampleSet
Getsentrix_barcode <- function(rgSet){
  sentrix_barcode <- rgSet@colData@listData$Slide
  return(sentrix_barcode)
}

#' @rdname GetsampleSet
Getsentrix_position <- function(rgSet){
  sentrix_position <- rgSet@colData@listData$Array
  return(sentrix_position)
}
