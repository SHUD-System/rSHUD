#' Generate the default Boundary Condition
#' \code{BC} 
#' @param x matrix of BC. c(days, value)
#' @param start Start time, POSIXct.
#' @return Time-series data of BC
#' @export
BC <- function(x=cbind(0:100, 0), start=as.POSIXct('2000-01-01')){
  tsd=xts::as.xts(x[,-1], 
                  order.by=start+x[,1]*86400)
  return(tsd)
}