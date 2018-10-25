#' Change Data Type
#'
#' @description Change data value type between "Beta" and "M". As we have:
#' B = 2^M/(2^M + 1) and M = log2(B/(1-B))
#'
#' @param InputValueType data value type of input data
#' @param TargetValueType data value type the output data should be
#' @param data the input data set
#'
#' @return \code{data} the output data set with wanted data value type
#
#' @export
ChangeDataType <- function(InputValueType, TargetValueType, data){
  if(InputValueType == "B"){
    if(TargetValueType == "M"){
      message("Change B value to M value")
      data[data[,] == 0] <- 0.0001
      data[data[,] == 1] <- 0.9999
      data <- log2(data/(1-data))
    }
  }else if(InputValueType == "M"){
    if(TargetValueType == "B"){
      message("Change M value to B value")
      data <- 2^data/(2^data + 1)
    }
  }
  return(data)
}
