#' Calculate the elevation of mesh cells
#' \code{getElevation} 
#' @param pm \code{shud.mesh}
#' @return elevation
#' @export
getElevation <- function(pm = readmesh()){
  tt <- pm@mesh[,2:4];
  zz = pm@point[,5]
  zs = (zz[tt[,1]] + zz[tt[,2]] + zz[tt[,3]] ) /3
  ret <- zs
  return(ret)
}
#' Calculate the AquiferDepth of mesh cells
#' \code{getAquiferDepth} 
#' @param pm \code{shud.mesh}
#' @return aquifer thickness
#' @export
getAquiferDepth <- function(pm = readmesh()){
  tt <- pm@mesh[,2:4]
  dd = pm@point[,4]
  aqd = (dd[tt[,1]] + dd[tt[,2]] + dd[tt[,3]] ) /3
  ret <- aqd
  return(ret)
}
#============
#' Get the Lake Element Index
#' \code{getLakeEleID} 
#' @param pa \code{shud.att}
#' @return index of lake element index
#' @export
getLakeEleID <-function(pa = readatt() ){
  id = which(pa$LAKE > 0)
  return(id)
}
#============
#' Calculate the area of the cells
#' \code{getArea} 
#' @param pm \code{shud.mesh}
#' @return Area of cells
#' @export
getArea <-function(pm=readmesh() ){
  spm=sp.mesh2Shape(pm=pm)
  ia=rgeos::gArea(spm, byid = TRUE)
  return(ia)
}

#============
#' Calculate the Vertex of the cells
#' \code{getVertex} 
#' @param pm \code{shud.mesh}
#' @return Vertex of \code{shud.mesh}, dim = c(ncell, vertex, 4), 3-dimension = c('x','y','AqD', 'zmax')
#' @export
getVertex <- function(pm = readmesh() ){
  msh <- pm@mesh;
  pts <- pm@point;
  nabr<-pm@mesh[,5:7]   #nabor or each cell.
  node<-pm@mesh[,2:4]    #vetex of each cell.
  
  xv       =cbind(pts[node[,1],2],pts[node[,2],2],pts[node[,3],2]); #end with v, means vertex of triangles.
  yv       =cbind(pts[node[,1],3],pts[node[,2],3],pts[node[,3],3]);
  aqd    =cbind(pts[node[,1],4],pts[node[,2],4],pts[node[,3],4]);
  zsurfv   =cbind(pts[node[,1],5],pts[node[,2],5],pts[node[,3],5]);
  
  ret = abind::abind(xv, yv, aqd, zsurfv,
                     along=3)
  dimnames(ret) = list(paste0('Cell.', 1:nrow(msh) ),
                       paste0('Vertex.', 1:3),
                       c('X','Y', 'AqD', 'ZMAX') )
  ret
  return(ret)
}

#' Get the From/To nodes of the river.
#' \code{getRiverNodes} 
#' @param spr SpatialLine* of river streams.
#' @return a list, c(points, FT_ID)
#' @export
getRiverNodes <- function(spr=readriv.sp() ){
  crs.pcs = raster::crs(spr)
  pts=extractCoords(spr)
  ft0=FromToNode(spr, coord=pts)
  ft.pts = pts[sort(unique(as.numeric(ft0))),]; 
  colnames(ft.pts)=c('X', 'Y')
  
  ft=FromToNode(spr, coord = ft.pts)
  
  ll = ProjectCoordinate(ft.pts, proj4string = crs.pcs, P2G=TRUE)
  ft.pts = cbind(ft.pts, ll)
  riv.pts = data.frame('ID'=1:nrow(ft.pts), ft.pts)
  riv.ft = data.frame('ID'=1:nrow(ft), ft)
  return(list(points=riv.pts, 
         FT_ID=riv.ft) )
}
#' Get the centroids of the cells
#' \code{getCentroid} 
#' @param pm \code{shud.mesh}
#' @return centroids of \code{shud.mesh}, ncell x 4 c('x','y','AqD', 'zmax')
#' @export
getCentroid <- function(pm = readmesh() ){
  x = getVertex(pm=pm);
  xc      =rowMeans(x[,,1]);      #end of c, means centroid of the triangles.
  yc      =rowMeans(x[,,2]);
  if(ncol(pm@mesh) ==  9){
    aqd   =pm@mesh[,8]
    zsurfc  =pm@mesh[,9]
  }else{
    aqd   =rowMeans(x[,,3]);
    zsurfc  =rowMeans(x[,,4]);
  }
  ret = cbind(1:length(xc), xc, yc, aqd, zsurfc)
  colnames(ret) = c('ID', 'X','Y', 'aqd', 'ZMAX')
  ret
  return(ret)
}
