#' Imputation of missing value
#'
#' @description K-nearest neighbor averaging (KNN) is used to impute missing vlaues
#' of the input data.
#'
#' @param inputdata.small input data set which needs imputation of
#' missing values
#'
#' @return \code{inputedata2} the data set after imputation, without missing
#' values now
#' @importFrom impute impute.knn
#' @export
imputeData <- function(inputdata.small){
  if(sum(is.na(inputdata.small)) > 0){
    message("KNN used to impute ", sum(is.na(inputdata.small)), " missing value")
    inputdata2 <- impute::impute.knn(inputdata.small, rowmax = 1, colmax = 1)$data
    message("Imputation Done (Check missing value fraction for each sample)")
  }else{
    inputdata2 <- inputdata.small
    message("No missing value for predictor.")
  }
  return(inputdata2)
}
