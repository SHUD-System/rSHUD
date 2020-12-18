#' Write ESRI shapefile out
#' \code{writeshape}
#' @param shp Spatial file
#' @param crs projection
#' @param file file path, without '.shp'.
#' @export
#' @examples
#' library(sp)
#' library(rgeos)
#' library(rgdal)
#' sp1 = readWKT("POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))")
#' raster::crs(sp1) =sp::CRS("+init=epsg:4326")
#' writeshape(sp1, file=file.path(tempdir(), 'sp1'))
#' sp2=readOGR(file.path(tempdir(), 'sp1.shp'))
#' plot(sp2)
writeshape <- function(shp, file=NULL, crs = raster::crs(shp)){
  msg='writeshape::'
  if(grepl(class(shp)[1],'SpatialPolygons' ) ){
    # shp = methods::as(shp, "SpatialPolygonsDataFrame")
    shp=sp::SpatialPolygonsDataFrame(shp, 
                                     data=data.frame('ID'=1:length(shp)),
                                     match.ID = FALSE)
  }else if ( grepl(class(shp)[1],'SpatialLines' )   ){
    # shp = methods::as(shp, "SpatialLinesDataFrame")
    shp=sp::SpatialLinesDataFrame(shp, data=data.frame('ID'=1:length(shp)),match.ID = FALSE)
  }
  if( is.null(file) ){
    # message(msg, 'No file exported')
  }else{
    path = dirname(file)
    fn = basename(file)
    if(!dir.exists(path)){
      dir.create(path, showWarnings = T, recursive = T)
    }
    raster::crs(shp) = crs;
    prj = sp::proj4string(shp)
    rgdal::writeOGR(obj=shp, driver = 'ESRI Shapefile',
                    layer=fn,
                    dsn=path, overwrite_layer = T)
    if( !is.na(crs) ){
      fn.prj = file;
      raster::extension(fn.prj) = '.prj'
      invisible(rgdal::showWKT(prj, file = fn.prj))
    }
    message(msg, file, ' is saved')
  }
}
#' Re-project coordinates betwen GCS and PCS
#' \code{ProjectCoordinate} 
#' @param  x 2-column matrix of coordinates.
#' @param  proj4string proj4string
#' @param  P2G if TRUE, a cartographic projection into lat/long, otherwise projects from lat/long into a cartographic projection.
#' @return Basic model infomation, figures and tables
#' @export
ProjectCoordinate <- function(x, proj4string, P2G=TRUE){
  # Transformed data
  x=as.matrix(x)
  y <- proj4::project(x, proj4string, inverse=P2G)
  if(P2G){
    colnames(y)=c('Lon','Lat')
  }else{
    colnames(y)=c('X','Y')
  }
  y
}
#' SpatialData to Raster
#' \code{sp2raster}
#' @param sp SpatialPolygon
#' @param mask Raster mask of mesh domain.
#' @param ngrids Number of grid along x direction.
#' @param resolution Resolution, defaul = NULL, resolution = extent / ngrids
#' @param field Index of field
#' @return Raster map
#' @export
sp2raster <- function (sp, mask = get('MASK', envir = .shud),
                       ngrids=200, 
                       resolution=NULL, field=1) {
  if( is.null(mask) ){
    ext <-  raster::extent (sp)
    xlim=ext[1:2]
    ylim=ext[3:4]
    if ( resolution<=0 || is.null(resolution)){
      dx=diff(xlim) / ngrids;
    }else{
      dx=resolution
    }
    r <- raster::raster(ext, res=dx)
  }else{
    r = mask
  }
  ## Rasterize the shapefile
  rr <-raster::rasterize(sp, r, field=field)
  return(rr)
}

#' Generate the raster mask of Mesh domain
#' \code{shud.mask}
#' @param pm \code{shud.mesh}
#' @param n Number of grid
#' @param rr Default mask in .shud environment
#' @param cellsize Resolution, defaul = NULL, resolution = extent / ngrids
#' @param proj Projection parameter
#' @return Raster map
#' @export
shud.mask  <- function (pm = readmesh(), proj=NULL,
                        rr = get('MASK', envir=.shud),
                        n=10000, cellsize=NULL){
  # mesh=readmesh(shp=TRUE); ngrids=100; resolution=0
  if(is.null(rr)){
    spm =sp.mesh2Shape(pm)
    sp0=rgeos::gUnaryUnion(spm)
    if(is.null(cellsize)){
      # grd <- as.data.frame(sp::spsample(spm, "regular", n=n))
      # grd <- as.data.frame(sp::spsample(spm, "regular", nsig=2, n=n))
      grd <- sp::makegrid(sp0, n = n)
    }else{
      # grd <- as.data.frame(sp::spsample(spm, "regular", cellsize = cellsize))
      grd <- sp::makegrid(sp0, cellsize = cellsize)
    }
    names(grd)       <- c("X", "Y")
    sp::coordinates(grd) <- c("X", "Y")
    sp::gridded(grd)     <- TRUE  # Create SpatialPixel object
    sp::fullgrid(grd)    <- TRUE  # Create SpatialGrid object
    rr=raster::raster(grd); rr[]=1
    rr=raster::mask(rr, sp0)
    if(!is.null(proj)){
      raster::crs(rr) = proj
    }
    assign('MASK', rr, envir=.shud)
  }else{
    rr = rr
  }
  rr
}

#' SpatialData to Raster
#' \code{MeshData2Raster}
#' @param x vector or matrix, length/nrow is number of cells.
#' @param rmask mask of the mesh file
#' @param stack Whether export the stack, only when the x is a matrix, i.e. (Ntime x Ncell).
#' @param proj Projejction parameter
#' @param pm shud mesh
#' @param plot Whether plot the result.
#' @return Raster map
#' @export
MeshData2Raster <- function(x=getElevation(),
                            rmask=shud.mask(proj=proj), 
                            pm=readmesh(), proj=NULL,
                            stack=FALSE, method='ide',
                            plot =FALSE){
  
  if(stack){
    ret <- raster::stack(apply(x, 1, FUN = MeshData2Raster) )
  }else{
    if( is.matrix(x) | is.data.frame(x)){
      x = as.numeric(x[nrow(x),])
    }
    if(any(is.na(x)) ){
      x[is.na(x)] = 0
    }
    if (any(is.infinite(x))){
      x[is.infinite(x)] = 0
    }
    xy=getCentroid(pm=pm)[,2:3]
    
    if(grepl('idw', tolower(method))){
      val= data.frame(xy, x)
      colnames(val) = c('X', 'Y', 'Z')
      sp::coordinates(val) = c('X', 'Y')
      grd=methods::as(rmask, 'SpatialGrid')
      # if(grepl(method, 'idw')){
      # Interpolate the grid cells using a power value of 2 (idp=2.0)
      dat <- gstat::idw(Z ~ 1, val, newdata=grd, idp=2.0)
      r = raster::raster(dat)
    }
    if(grepl('linear', tolower(method))){
      xr = raster::rasterToPoints(rmask)
      ext=raster::extent(rmask);res=raster::res(rmask); hr = res/2
      r0=rmask; r0[]=1
      xyo=raster::rasterToPoints(r0)
      xx=interp::interp(x=xy[,1], y=xy[,2], z=x, xo=xyo[,1], yo=xyo[,2] )
      r = raster::setValues(rmask, as.numeric(xx$z))
    }
    if(grepl('ide', tolower(method))){
      tps <- fields::Tps(xy, x)
      # use model to predict values at all locations
      r <- raster::interpolate(rmask, tps)
    }
    ret <- raster::mask(r,rmask)
  }
  if(plot){
    raster::plot(ret)
  }
  if(!is.null(proj)){ raster::crs(ret) <- proj }
  return(ret)
}


#' Remove the holes in polygons
#' \code{removeholes}
#' @param sp SpatialPolygons or SpatialPolygonDataFrame
#' @return Raster map
#' @export
#' @examples
#' library(sp)
#' p.out = Polygon(cbind(c(4,4,6,7,4),c(5,3,2,5,5))  )
#' p.hole = Polygon(cbind(c(5,6,6,5,5),c(4,4,3,3,4) ), hole = TRUE)
#' sp <- SpatialPolygons(list(Polygons(list(p.out, p.hole), "1")))
#' s = removeholes(sp)
#' par(mfrow=c(1,2))
#' plot(sp)
#' plot(s)
removeholes <- function(sp){
  x = sp
  nx = length(x)
  # vn = rownames(x@data)
  rl = list()
  for(i in 1:nx){
    spg = x@polygons[[i]]@Polygons
    npg = length(spg)
    ypg = list()
    k=1
    for(j in 1:npg){
      if(spg[[j]]@hole){
      }else{
        ypg[[k]] = spg[[j]]
        k = k+1
      }
    }
    rl[[i]] = sp::Polygons(ypg, ID=i)
  }
  ysp = sp::SpatialPolygons(rl)
  # ret = SpatialPolygonsDataFrame(ysp, data=x@data)
  ret = ysp
  if( !is.na(raster::crs(x)) & ! is.null(raster::crs(x)) ){
    raster::crs(ret) = raster::crs(x)
  }
  return(ret)
}
#' Generatue fishnet
#' \code{fishnet}
#' @param ext Extension of the fishnet. c(xmin, xmax, ymin, ymax)
#' @param n Number of dx
#' @param dx Interval of x direction
#' @param dy Interval of y direction
#' @param crs Projection
#' @param polygons Whether to export SpatialPolygons
#' @param points Whether to export SpatialPoints
#' @param lines Whether to export SpatialLines
#' @export
#' @examples
#' library(raster)
#' pg=fishnet(ext=c(-80,80, -50,50), dx=5)
#' plot(pg)
fishnet <- function(ext, crs=sp::CRS("+init=epsg:4326"),
                    n=10,
                    dx=diff(ext[1:2])/n, dy=dx,
                    lines=FALSE, polygons=TRUE, points=FALSE){
  xmin = ext[1]
  xmax = ext[2]
  ymin = ext[3]
  ymax = ext[4]
  dx = min(dx, diff(ext[1:2]))
  dy = min(dy, diff(ext[3:4]))
  xx=seq(ext[1], ext[2], by=dx)
  yy=seq(ext[3], ext[4], by=dy)
  nx = length(xx)
  ny = length(yy)
  if(lines){
    vline=cbind(xx, ymin, xx, ymax)
    hline =  cbind(xmin, yy, xmax, yy)
    mat = rbind(vline, hline)
    str=paste(paste('(', mat[,1], mat[,2],',', mat[,3], mat[,4], ')'), collapse = ',')
    spl=rgeos::readWKT(paste('MULTILINESTRING(', str, ')'))
    df = as.data.frame(mat)
    colnames(df) = c('x1', 'y1', 'x2','y2')
    spdf=sp::SpatialLinesDataFrame(spl, data = df)
    ret = spdf
    raster::crs(ret) = crs
    return(ret)
  }
  if(polygons){
    xy=expand.grid(xx,yy)
    xm = matrix(xy[,1], nx,ny)
    ym = matrix(xy[,2], nx, ny)
    xloc=abind::abind(as.matrix(xm[-nx, -ny]), as.matrix(xm[-nx, -1]), as.matrix(xm[-1, -1]),
                      as.matrix(xm[-1, -ny]), as.matrix(xm[-nx, -ny]), along=3)
    yloc=abind::abind(as.matrix(ym[-nx, -ny]), as.matrix(ym[-nx, -1]), as.matrix(ym[-1, -1]),
                      as.matrix(ym[-1, -ny]), as.matrix(ym[-nx, -ny]), along=3)
    
    # df=as.data.frame(matrix(0, nrow=(nx-1)*(ny-1), 6))
    df=data.frame(as.numeric(apply(xloc, 1:2, min)),
                  as.numeric(apply(xloc, 1:2, max)),
                  as.numeric(apply(yloc, 1:2, min)),
                  as.numeric(apply(yloc, 1:2, max)))
    df = data.frame(df, rowMeans(df[,1:2]), rowMeans(df[,1:2+2]) )
    colnames(df) = c('xmin','xmax','ymin', 'ymax','xcenter','ycenter')
    str=paste('GEOMETRYCOLLECTION(',
              paste(paste('POLYGON((',
                          paste(xm[-nx, -ny], ym[-nx, -ny], ',' ),
                          paste(xm[-nx, -1],  ym[-nx, -1], ','),
                          paste(xm[-1, -1],   ym[-1, -1], ','),
                          paste(xm[-1, -ny],  ym[-1, -ny], ','),
                          paste(xm[-nx, -ny], ym[-nx, -ny], '' ), '))' )
                    , collapse =','),
              ')' )
    # str=paste('MULTIPOLYGON(', paste(xt, collapse = ', '), ')')
    SRL = rgeos::readWKT(str)
    # x1 = x0@polygons[[1]]@Polygons
    # SRL =lapply(1:length(x1),  function(x, i) {Polygons(list(x[[i]]), ID=i)},  x=x1 )
    ret = sp::SpatialPolygonsDataFrame( Sr=SRL, data=df, match.ID = TRUE)
    raster::crs(ret) = crs
    return(ret)
  }
  return(NA)
}
#' Add holes into Polygons
#' \code{AddHoleToPolygon}
#' @param poly SpatialPolygons
#' @param hole Hole Polygon
#' @export
AddHoleToPolygon <-function(poly,hole){
  # https://stackoverflow.com/questions/29624895/how-to-add-a-hole-to-a-polygon-within-a-spatialpolygonsdataframe
  # invert the coordinates for Polygons to flag it as a hole
  coordsHole <-  hole@polygons[[1]]@Polygons[[1]]@coords
  newHole <- sp::Polygon(coordsHole,hole=TRUE)
  
  # punch the hole in the main poly
  listPol <- poly@polygons[[1]]@Polygons
  listPol[[length(listPol)+1]] <- newHole
  punch <-sp::Polygons(listPol,poly@polygons[[1]]@ID)
  
  # make the polygon a SpatialPolygonsDataFrame as the entry
  new <- sp::SpatialPolygons(list(punch),proj4string=poly@proj4string)
  new <- sp::SpatialPolygonsDataFrame(new,data=as(poly,"data.frame"))
  return(new)
}
#' Cut sptialLines with threshold.
#' \code{sp.CutSptialLines}
#' @param sl SpatialLines or SpatialLineDataFrame
#' @param tol Tolerence. If the length of segment is larger than tolerance, cut the segment until the maximum segment is shorter than tolerance.
#' @export
#' @examples
#' library(rSHUD)
#' library(sp)
#' x=1:1000/100
#' l1 = Lines(Line(cbind(x, sin(x)) ), ID='a' )
#' sl = SpatialLines( list(l1) )
#' tol1=5;
#' tol2 =2
#' sl1 = sp.CutSptialLines(sl, tol1)
#' sl2 = sp.CutSptialLines(sl, tol2)
#' par(mfrow=c(1,2))
#' plot(sl1, col=1:length(sl1));title(paste0('Tol=', tol1))
#' plot(sl2, col=1:length(sl2));title(paste0('Tol=', tol2))
#'
#' data(sh)
#' riv=sh$riv
#' x = sp.CutSptialLines(riv, tol=5)
#' par(mfrow=c(2,1))
#' plot(riv, col=1:length(riv), lwd=3);
#'
#' plot(riv, col='gray', lwd=3);
#' plot(add=TRUE, x, col=1:length(x))
sp.CutSptialLines <- function(sl, tol){
  msg='sp.CutSptialLines::'
  ll = rgeos::gLength(sl, byid = TRUE)
  if(all(ll < tol) ){
    ret = sl
  }else{
    nsp = length(sl)
    xsl = list(); ik=1
    for(i in 1:nsp){
      sx = sl[i, ]
      pxy = extractCoords(sx,unique = TRUE)
      np = nrow(pxy)
      dacc = cumsum( sp::LineLength(pxy, sum = FALSE))
      # dacc =getDist(pxy)
      tol = max(c(tol, min(dacc) ) )
      len= rgeos::gLength(sx)
      if(len > tol){
        nsplit = ceiling(len / tol)
      }else{
        nsplit = 1
      }
      dd = len / nsplit
      v0 = 1  # Vetex 0, Vetex 1
      message(msg, i, '/', nsp, '\t', nsplit, '\t', round(dd, 2) )
      for(k in 1:nsplit){
        if(v0 >=np){
          break
        }
        # message(msg, '\t', k, '/', nsplit)
        dk = dd * k
        v1 = order(abs(dacc - dk), decreasing = FALSE)[1] + 1
        if(v1 + 1>np){
          v1 = np
        }
        message(msg, v0,'\t', v1)
        if(v0 == v1){
          next;
        }
        # plot(sl[i, ]);points(pxy); points(pxy[c(v0, v1), ], pch=2, col=2)
        xsl[[ik]]= sp::Lines(sp::Line( pxy[c(v0:v1), ]), ID=ik)
        ik=ik+1
        # points(pxy[v0:v1,], col=k)
        v0=v1
      }
    }
    nsl = length(xsl)
    tmp = sp::SpatialLines(xsl, proj4string = raster::crs(sl))
    ilen = rgeos::gLength(tmp, byid=TRUE)
    att=data.frame('INDEX'=1:length(tmp), 'Length'=ilen)
    ret = sp::SpatialLinesDataFrame(tmp, data = att)
  }
  return(ret)
}

#' Extract values on Raster map. The line is a straight line between (0,1). 
#' \code{extractRaster}
#' @param r Raster
#' @param xy coordinates of the line, dim=(Npoints, 2); x and y must be in [0, 1]
#' @param ext extension of value xy.
#' @param plot Whether plot result.
#' @importFrom grDevices dev.off graphics.off png rgb topo.colors
#' @importFrom graphics grid hist lines par plot points
#' @importFrom methods as
#' @importFrom stats dist rnorm time
#' @importFrom utils read.table
#' @export
#' @examples
#' library(raster)
# r <- raster(ncol=36, nrow=18)
# r[] <- 1:ncell(r)
# extractRaster(r)
extractRaster<-function(r, xy=NULL, ext = raster::extent(r), plot=T){
  if(is.null(xy)){
    ndim = dim(r)
    x=0:ndim[2] / ndim[2]
    y = rep(0.5, length(x))
    xy = cbind(x,y)
  }
  x = ext[1] + xy[,1] * (ext[2]- ext[1] )
  y = ext[3] + xy[,2] * (ext[4]- ext[3] )
  if(plot){
    raster::plot(r);
    points(x, y, col=2)
    nx=length(x)
    points(x,y)
    graphics::arrows(x[1], y[1], x[nx], y[nx], lty=3, lwd=1.5, col=2)
    # lines(x,y, lwd=1.5, col=2, lty=2)
  }
  v = raster::extract(r, cbind(x,y))
  ret = cbind('x'=x,'y'=y,'z'=v)
  return(ret)
}

#' Simplify SpatialData.
#' @param x SpatialData
#' @return Simplified SpatialData
#' @export
SimpleSpatial <-function(x){
  # n1=length(x@polygons)
  # nj=unlist(lapply(1:n1, function(i){ length(x@polygons[[i]]@Polygons) } ))
  # x@polygons[[1]]@Polygons[[1]]@coords
  msg='SimpleSpatial'
  ni = length(x@polygons)
  k=1
  sl=list()
  for(i in 1:ni){
    nj = length(x@polygons[[1]]@Polygons)
    for(j in 1:nj){
      cd = x@polygons[[i]]@Polygons[[j]]@coords
      np=nrow(cd)
      message(msg, i,'-',j ,'\t', np)
      sl[[k]] = paste('POLYGON((', paste( paste(cd[,1], cd[,2]), collapse = ',' ), '))')
      if(k==1){
        str = sl[[k]]
      }else{
        str = paste(str, ',', sl[[k]] )
      }
      k=k+1
    }
  }
  r=rgeos::readWKT(paste('GEOMETRYCOLLECTION(', str, ')'))
}

#' Find the points in distance less than tol.
#' \code{PointInDistance}
#' @param pt 2-column coordinates (x,y).
#' @param tol Tolerance
#' @return Index of points that within tol
#' @export
PointInDistance <- function(pt, tol){
  msg='PointInDistance'
  dm = as.matrix(stats::dist(pt, diag  = TRUE))
  dm[dm==0]=NA
  # View(dm)
  dmin=apply(dm, 1, min, na.rm=T)
  id=which(dmin < 100)
  id1=id2=NULL
  tmp1=tmp=dmin[id]
  i=2
  for(i in 1:length(id)){
    if(id[i] %in% id2){next }
    id1 = c(id1, i); id1
    tmp1[id1] = NA; tmp1
    id[which(tmp1 %in% tmp[id1])]
    id2=unique(c(id2,  id[which(tmp1 %in% tmp[id1])]))
    id2
    tmp1=tmp
  }
  cbind('P1'=id[id1], P2=id2)
}

#' Conver the MULTIPOLYGONS to SINGLEPOLYGONS.
#' @param x the spatialpolygon*
#' @param id Index of the sorted (decreasing) polygons to return. default = 0;
#' @export
SinglePolygon <- function(x, id=0){
  n1 = length(x)
  y1 = list()
  y2 = list()
  k=1
  for(i in 1:n1){
    message('level 1:', i, '/', n1)
    x1 = x@polygons[[i]]
    
    n2=length(x1@Polygons)
    for(j in 1:n2){
      message('level 2:', j, '/', n2)
      x2 = x1@Polygons[[j]]
      y1[[k]] = Polygons( list(x2), ID=k)
      k=k+1
    }
  }
  
  y=sp::SpatialPolygonsDataFrame(SpatialPolygons(y1), data=data.frame('ID'=2:k-1))
  
  if(id < 1){
    return(y)
  }else{
    ia=rgeos::gArea(y, byid = TRUE)
    id=order(ia, decreasing = TRUE)[1]
    return(y[id,])
    
  }
}

#' Remove the duplicated lines, which means the FROM and TO points are identical.
#' \code{rmDuplicatedLines}
#' @param x ShapeLine*
#' @param ... More options in duplicated()
#' @return ShapeLine* without duplicated lines.
#' @export
rmDuplicatedLines <- function(x, ...){
  # x=spi.riv
  cd = extractCoords(x)
  # dim(cd)
  ft = FromToNode(x, cd)
  id=which(duplicated(ft, MARGIN = 1, ...))
  if(length(id)>0){
    r =x[-id,]
  }else{
    r = x
  }
  return(r)
}

#' Generate Thiesson Polygons from a point data
#' \code{voronoipolygons}
#' @param x ShapePoints* in PCS
#' @param pts Coordinates (x,y) of the points
#' @param rw extent
#' @param crs projection parameters
#' @return ShapePolygon*
#' @export
#' @examples 
#' library(rgeos)
#' library(rSHUD)
#' n=10
#' xx = rnorm(n)
#' yy = rnorm(n)
#' str = paste('MULTIPOINT(', paste(paste('(', xx, yy, ')'), collapse = ','), ')')
#' x=readWKT(str)
#' vx = voronoipolygons(x)
#' raster::plot(vx, axes=TRUE)
#' raster::plot(add=TRUE, x, col=2)
#' #' ====END====
#' 
#' x=1:5
#' y=1:5
#' xy=expand.grid(x,y)
#' vx=voronoipolygons(pts=xy)
#' plot(vx, axes=TRUE)
#' points(xy)
#' #' ====END====
#' 
#' library(rgdal)
#' library(raster)
#' library(rgeos)
#' n=10
#' xx = rnorm(n)
#' yy = rnorm(n)
#' str = paste('MULTIPOINT(', paste(paste('(', xx, yy, ')'), collapse = ','), ')')
#' x=readWKT(str)
#' y=readWKT(paste('MULTIPOINT(', paste(paste('(', xx+2, yy+2, ')'), collapse = ','), ')'))
#' e1 = extent(y)
#' e2 = extent(x)
#' rw = c(min(e1[1], e2[1]),
#'        max(e1[2], e2[2]),
#'        min(e1[3], e2[3]),
#'        max(e1[4], e2[4]) ) + c(-1, 1, -1, 1)
#' vx=voronoipolygons(x=x, rw=rw)
#' plot(vx); plot(add=TRUE, x, col=2); plot(add=TRUE, y, col=3)
voronoipolygons = function(x, pts = x@coords, rw=NULL, crs=NULL) {
  z = deldir::deldir(pts[,1], pts[,2], rw=rw)
  w = deldir::tile.list(z)
  polys = vector(mode='list', length=length(w))
  for (i in seq(along=polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = sp::Polygons(list(sp::Polygon(pcrds)), ID=as.character(i))
  }
  SP = sp::SpatialPolygons(polys)
  voronoi = sp::SpatialPolygonsDataFrame(SP, 
                                         data=data.frame(x=pts[,1], 
                                                         y=pts[,2], row.names=sapply(slot(SP, 'polygons'), 
                                                                                     function(x) slot(x, 'ID'))))
  if(!is.null(crs)){
    raster::crs(voronoi) = crs;
  }
  return(voronoi)
}



#' Generate the coverage map for forcing sites.
#' \code{ForcingCoverage}
#' @param sp.meteoSite ShapePoints* in PCS
#' @param pcs Projected Coordinate System
#' @param gcs Geographic Coordinate System
#' @param dem DEM raster
#' @param wbd watershed boundary
#' @param enlarge enlarge factor for the boundary.
#' @return ShapePolygon*
#' @export
ForcingCoverage <- function(sp.meteoSite=NULL, pcs, gcs=sp::CRS('+init=epsg:4326'), 
                            dem, wbd, enlarge = 10000){
  if( is.null(sp.meteoSite) ){
    sp.meteoSite = rgeos::gCentroid(wbd)
  }
  x.pcs = spTransform(sp.meteoSite, pcs)
  x.gcs = spTransform(sp.meteoSite, gcs)
  
  ll = coordinates(x.gcs)
  xy = coordinates(x.pcs)
  z = raster::extract(dem, x.pcs)
  
  att = data.frame(1:length(sp.meteoSite), ll, xy, z, paste0(sp.meteoSite@data$ID, '.csv'))
  colnames(att) = c('ID', 'Lon', 'Lat', 'X', 'Y','Z', 'Filename')
  
  e1 = extent(wbd)
  e2 = extent(x.pcs)
  rw = c(min(e1[1], e2[1]), 
         max(e1[2], e2[2]),
         min(e1[3], e2[3]),
         max(e1[4], e2[4]) ) + c(-1, 1, -1, 1) * enlarge
  sp.forc=voronoipolygons(x=x.pcs, rw=rw, crs=pcs)
  att[is.na(att)] = -9999
  sp.forc@data = att
  return(sp.forc)
}

