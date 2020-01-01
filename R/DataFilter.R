
#' Filter the river data with a filter
#' \code{datafilter.riv} 
#' @param x  Input data
#' @param filter Data filter
#' @param plot Whether plot the data
#' @importFrom grDevices dev.off graphics.off png rgb topo.colors 
#' @importFrom graphics grid hist lines par plot points 
#' @importFrom methods as 
#' @importFrom stats dist rnorm time 
#' @importFrom utils read.table 
#' @return Matrix information, c('ID','Vmin','Vmax', 'Filter')
#' @export
datafilter.riv <-function(x, filter=NULL, plot=TRUE){
  msg=paste0('datafilter.riv::')
  # y=x[['YRivstage']]
  y = x
  # plot(y)
  pr=readriv()
  cb = readcalib()
  tid = pr@river[,'Type']
  uid = sort(unique(tid))
  st=pr@rivertype[tid, 'Depth'] + cb['RIV_DPTH']
  if( is.null(filter) ){
    filter = st
  }
  ymax = apply(y, 2, max, na.rm=T)
  ymin = apply(y, 2, min, na.rm=T)
  id = which(ymax > filter) 
  ret = data.frame(id, ymin[id], ymax[id], filter[id])
  colnames(ret) = c('ID','Vmin','Vmax', 'Filter')
  rownames(ret) = id
  ylim=range(c(filter, y))
  if(plot ){
    if(length(id) > 0){
      id = id
      message(msg, length(id), ' rivers are filtered.')
    }else{
      id = 1:ncol(x)
    }
    yv = sort(( unique(filter) ))
    ny = length(yv)
    col = uid
    zoo::plot.zoo(y[,id], col=col[tid[id]], ylim=ylim, screen=1)
    graphics::abline( h=yv, col=col, lwd=3, lty=2)
  }
  ret
}
