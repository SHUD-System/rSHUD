#' Generate triangle mesh domain
#' \code{shud.domain} 
#' @param wb SpatialPolygon or SpatialLines which define the watershed boundary
#' @param riv SpatialLines of river network, optional
#' @param dem Elevation data.
#' @param lk SpatialPolygon of lake.
#' @param q minimum angle of triangle
#' @param pts Extra pts to build triangular mesh.
#' @param ... more options in RTriangle::triangulate()
#' @return A object with class triangulation.
#' @export
shud.triangle <- function(wb, dem = NULL, 
                                  riv=NULL, lk=NULL, 
                                  pts=NULL,
                                  q=30, ...){
  #x is a spatialpolygons
  # wb=wb$wb.simp
  # riv=riv$Riv.simp
  # wb=wb.simp
  # lk=lk.sim
  # plot(wb);plot(lk,add=T)
  ps1 = sp2PSLG(wb)
  if(!is.null(dem)){
    raster::movingFun()
  }
  # if(!is.null(lk)){
  #   ps3 = sp2PSLG(lk)
  #   n1 = nrow(ps1$P)
  #   ps=RTriangle::pslg(P = rbind(ps1$P, ps3$P),
  #           S = rbind(ps1$S, ps3$S + n1) )
  # }else{
  #   
  # }
  if(!is.null(riv)){
    ps2 = sp2PSLG(riv)
    n1 = nrow(ps1$P)
    ps=list('P' = rbind(ps1$P, ps2$P),
            'S' = rbind(ps1$S, ps2$S + n1) )
  }else{
    ps = ps1
  }
  if(!is.null(pts) ){
    ps$P = rbind(ps$P, pts)
  }
  p = RTriangle::pslg(P=ps$P,
           S = ps$S)
  # tri <- RTriangle::triangulate(p, a=500000, q=20)
  # plot(tri, asp=1, type='n')
  # dim(tri$T)
  
  if(q >35){
    q = 35;
  }
  tri <- RTriangle::triangulate(p, q=q,...)
  # plot(tri, asp=1)
  # points(ps1$P, col=2)
  # points(ps2$P, col=3)
  # plot(riv, add=TRUE, col=4)
  tri
}

#' Generate triangle mesh domain
#' \code{sp.mesh2Shape} 
#' @param pm SHUD mesh object
#' @param dbf attribute table of the mesh triangles.
#' @param crs  Projection parameters.
#' @return SpatialPolygons object
#' @export
sp.mesh2Shape <- function(pm=readmesh(), dbf=NULL, crs=NULL){
  tt = pm@mesh[,2:4]
  pp=pm@point[,2:3]
  dd = pm@point[,4]
  zz = pm@point[,5]
  tri = list('T'=tt, 'P'=pp)
  if(is.null(dbf)){
    aqd = (dd[tt[,1]] + dd[tt[,2]] + dd[tt[,3]] ) /3
    zs = (zz[tt[,1]] + zz[tt[,2]] + zz[tt[,3]] ) /3
    dbf = cbind(pm@mesh, 'AqDepth'= aqd, 'Zsurf'=zs)
  }
  ret <- sp.Tri2Shape(tri, dbf=dbf)
  if(!is.null(crs)){
    raster::crs(ret) = crs
  }
  return(ret)
}
#' Generate triangle mesh domain
#' \code{sp.Tri2Shape} 
#' @param tri triangulate 
#' @param dbf attribute table, data.frame or matrix
#' @param crs Projection 
#' @return Coordinates of triangles centroids
#' @export
sp.Tri2Shape <- function(tri, dbf=NULL, crs=NA){
  ta = tri$T
  pt = tri$P
  ncell = nrow(ta)
  p.x = pt[,1]
  p.y = pt[,2]
  ipt = t (ta[,c(1:3,1)] )
  xp = matrix( p.x[ipt], nrow=4)
  yp = matrix( p.y[ipt], nrow=4)
  s.tri=matrix(paste(as.numeric(xp), as.numeric(yp)), nrow=4)
  str=paste('GEOMETRYCOLLECTION(', 
            paste(
              paste('POLYGON((',apply(s.tri, 2, paste, collapse=','),'))'),
              collapse = ','),
            ')')
  SRL = rgeos::readWKT(str)
  ia=rgeos::gArea(SRL, byid = TRUE)
  dbf=cbind(dbf, 'Area'=ia)
  ret = sp::SpatialPolygonsDataFrame(Sr=SRL, data=as.data.frame(dbf), match.ID = FALSE)
  if(!is.na(crs) ){
    raster::crs(ret) = crs
  }
  ret
}
