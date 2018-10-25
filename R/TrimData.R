#' Trim data set
#'
#' @description Trim the data set with probes for specific prediction models
#'
#' @param inputdata input data
#' @param probeNames the probe list needed for prediction models
#'
#' @return \code{inputdata.small} The smaller data set with only specific
#' probes included
#'

#' @export
TrimData <- function(probeNames,  inputdata){
  inputdata <- as.matrix(inputdata)
  inputdata.small <- matrix(ncol = ncol(inputdata), nrow = length(probeNames))
  rownames(inputdata.small) <- probeNames
  colnames(inputdata.small) <- colnames(inputdata)
  for(i in 1:length(probeNames)){
    if(rownames(inputdata.small)[i] %in% rownames(inputdata)){
      inputdata.small[i,] <- inputdata[rownames(inputdata.small)[i],]
    }else if(!rownames(inputdata.small)[i] %in% rownames(inputdata)){
      inputdata.small[i,] <- rep(NA, ncol(inputdata))
    }
  }
  return(inputdata.small)
}
