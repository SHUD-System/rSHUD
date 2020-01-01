
#' Simplify Polylines by length
#' \code{sp.simplifyLen} 
#' @param sp a SpatialLines or SpatialPoygons object
#' @param tol Tolerance
#' @param threshold Minimum length for a segment
#' @export
sp.simplifyLen <-function(sp, tol, threshold=tol/2){
  getdiag<-function(x, pos){
    d=row(x) - col(x)
    r = x[d == pos]
    r
  }
  getDist <- function(xy, acc=T){
    dm=as.matrix(dist(xy[,], diag=T))
    if(acc){
      dv = cumsum(  c(0, getdiag(dm, 1)) )
    }else{
      dv = dm[,1]
    }
    dv
  }
  
  simplify.line<- function(sll, tol=NULL, threshold=tol/2){
    coord = sp::coordinates(sll)
    xd = getDist(coord)
    nx = length(xd)
    if(is.null(tol)){
      tol = max(xd)/10
    }
    
    seqdist = unique(c(seq(0, max(xd), by=tol), max(xd) ) )
    nq = length(seqdist)
    
    if(nq<=2){
      rt = coord[c(1,nx),]
    }else{
      if(diff(seqdist[-1:0 + nq]) < threshold ){
        #if the P(n-1) too close to P(n), remove the P(n-1) point.
        seqdist = seqdist[-1 * ( nq-1)]
      }
      np = length(seqdist)
      rt = matrix(0, nrow=np, ncol=2)
      
      rt[1,] = coord[1, ]
      rt[np,] = coord[nx, ]
      for(i in 2:(np-1)){
        id = order( abs(seqdist[i] - xd))[1]
        rt[i, ] = coord[id, ]
      }
      
    }
    rt
  }
  shp = methods::as(sp, 'SpatialLines')
  nshp=length(shp)
  shp.new=shp;
  
  for(i in 1:nshp){
    # i=1
    sl = shp@lines[[i]]@Lines
    # sl = pls@Lines
    nsl =length(sl)
    for(j in 1:nsl){
      sll = sl[[j]]
      sll.simp = simplify.line(sll, tol=tol)
      shp.new@lines[[i]]@Lines[[j]]@coords=sll.simp
    }
  }
  shp.new 
}
