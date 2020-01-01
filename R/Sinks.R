
#' Get the Sink cells in SHUD mesh
#' \code{meshSinks} 
#' @param pm SHUD mesh
#' @param details Whether return details of sinks infomration.
#' @return Vector of TRUE/FALSE or data.frame
#' @export
meshSinks <- function(pm = readmesh(), details=FALSE){
  ze = getElevation(pm=pm)
  nab = pm@mesh[4 + 1:3]
  zi=cbind(ze, ze, ze)
  zn=zi*NA
  i=1
  for(i in 1:3){
    nid=nab[,i]
    ltid = which(nid > 0)
    zn[ltid, i] = ze[ nid[ltid] ]
  }
  dz = zn - zi
  tf=apply(dz, 1, FUN=function(x) {all( x >= 0, na.rm = TRUE)} )
  if(details){
    ret = data.frame(ze, dz, apply(dz, 1, min, na.rm=T), (tf * 1) )
    colnames(ret) = c('Zsurf', paste0('dz', 1:3), 'minDZ', 'Sink')
  }else{
    ret = tf
  }
  ret
}
