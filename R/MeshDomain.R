#' Generate the mesh data from the triangulation
#' \code{shud.mesh}
#' @param tri Triangles defination
#' @param dem Elevation. Projection of the DEM raster must be same as watershed boundary.
#' @param AqDepth Aquifer depth, numeric.
#' @param r.aq  Aquifer Thickness. Raster object
#' @return Triangle mesh of model domain.
#' @export
#' @examples
#' library(raster)
#' data(sac)
#' wbd=sac[['wbd']]
#' dem=sac[['dem']]
#' a.max = 1e6 * 1;
#' q.min = 33;
#' tol.riv = 200
#' tol.wb = 200
#' tol.len = 500
#' AqDepth = 10
#'
#' wbbuf = rgeos::gBuffer(wbd, width = 2000)
#' dem = raster::crop(dem, wbbuf)
#'
#' wb.dis = rgeos::gUnionCascaded(wbd)
#' wb.simp = rgeos::gSimplify(wb.dis, tol=tol.wb, topologyPreserve = TRUE)
#'
#'
#' tri = shud.triangle(wb=wb.simp,q=q.min, a=a.max)
#' plot(tri, asp=1)
#'
#' # generate SHUD .mesh
#' pm=shud.mesh(tri,dem=dem, AqDepth = AqDepth)
#' sm = sp.mesh2Shape(pm)
#' raster::plot(sm)
shud.mesh <- function(tri, dem, AqDepth = 10, r.aq = dem * 0 + AqDepth){
  topo=tri$NB
  topo[topo<0]=0
  pt = tri$P;
  pid=tri$T;
  
  x = (pt[pid[,1],1] + pt[pid[,2],1] + pt[pid[,3], 1]) / 3
  y = (pt[pid[,1],2] + pt[pid[,2],2] + pt[pid[,3], 2]) / 3
  
  cxy = cbind(x, y)
  zc=raster::extract(dem, cxy)
  
  m = data.frame(1:nrow(topo), tri$T, topo[,1:3], zc)
  colnames(m) = c('ID', paste0('Node', 1:3), paste0('Nabr',1:3),'Zmax' )
  
  zs=raster::extract(dem, pt)
  aq=raster::extract(r.aq, pt)
  pt = data.frame(1:nrow(pt), pt, aq, zs)
  colnames(pt) = c('ID', 'X','Y', 'AqDepth', 'Elevation');
  mm=SHUD.MESH(mesh=m, point=pt)
}

#' Centroids of the triangulation
#' \code{Tri2Centroid}
#' @param tri Triangles defination
#' @return Centroids of the triangles, m x 2;
#' @export
Tri2Centroid <- function(tri){
  tt=tri$T;
  pt=tri$P;
  xc= (pt[tt[,1],1] + pt[tt[,2],1] + pt[tt[,3],1]) / 3;
  yc= (pt[tt[,1],2] + pt[tt[,2],2] + pt[tt[,3],2]) / 3;
  ret <- cbind(xc,yc)
}
#' extract Coordinates of SpatialLines or  SpatialPolygons
#' \code{shud.att}
#' @param tri Triangles defination
#' @param r.soil raster of soil classification
#' @param r.geol raster of geology layer
#' @param r.lc raster of land cover, LAI and Roughness Length
#' @param r.forc raster of forcing data
#' @param r.mf raster of melt factor
#' @param r.BC raster of boundary condition
#' @param r.SS raster of Source/Sink
#' @param sp.lake SpatialPolygon of lake layers.
#' @return data.frame of SHUD .att
#' @export
shud.att <- function(tri, r.soil =NULL, r.geol=NULL, r.lc=NULL, r.forc=NULL,
                     r.mf = NULL, r.BC = NULL, r.SS =NULL, sp.lake=NULL){
  p.centroids = Tri2Centroid(tri)
  ncell = nrow(p.centroids)
  atthead=c( "INDEX",  "SOIL", "GEOL", "LC", 
             'FORC', 'MF', 'BC', 'SS', 'LAKE')
  nh = length(atthead)
  att = data.frame(cbind(1:ncell, 1, 1, 1, 
              1, 1, 0, 0, 0) )
  extract.id <- function(r, p.centroids){
    id = raster::extract(r, p.centroids)
    if(is.matrix(id) | is.data.frame(id)){
      ret=id[,2]
    }else{
      ret =id
    }
    ret
  }
  apply.raster <- function(rr, pxy, v0){
    xv = rep(v0, nrow(pxy))
    if (!is.null(rr)){
      if( is.numeric(rr)){
        xv = xv * rr
      }else{
        xv = extract.id(rr, pxy)
      }
    }
    xv[is.na(xv)] = v0
    return(xv)
  }
  colnames(att) = atthead;
  nx = nrow(p.centroids)
  att$SOIL = apply.raster(r.soil, p.centroids, v0=1)
  att$GEOL = apply.raster(r.geol, p.centroids, v0=1)
  att$LC = apply.raster(r.lc, p.centroids, v0=1)
  att$FORC = apply.raster(r.forc, p.centroids, v0=1)
  att$MF = apply.raster(r.mf, p.centroids, v0=1)
  att$BC = apply.raster(r.BC, p.centroids, v0=0)
  att$SS = apply.raster(r.SS, p.centroids, v0=0)
  att$LAKE = apply.raster(sp.lake, p.centroids, v0=0)
  return(att)
}
