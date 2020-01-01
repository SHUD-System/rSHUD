#' Plot map of time-series data
#' \code{ts2map} 
#' @param x Time-series data
#' @param fun Functions to summarize the TS data.
#' @param raster Whether to make output as Raster
#' @param ... More options in fun
#' @export
ts2map <- function(x, fun = mean, raster=TRUE, ...){
  y = apply(x, 2, fun, ...)
  if(raster){
    r=MeshData2Raster(y)
    raster::plot(r)
  }else{
    dbf = cbind(y); colnames(dbf) = 'tsvalue'
    r = sp.mesh2Shape(dbf=y)
    plot_sp(r, zcol='tsvalue')
  }
  return(r)
}
#' Plot multiple maps
#' \code{compareMaps} 
#' @param r List of raster or SpatialData
#' @param mfrow mfrow in par()
#' @param contour Whether plot the contour.
#' @param ... More options in par()
#' @export
#' @examples 
#' library(raster)
#' data(volcano)
#' r <- raster(volcano)
#' extent(r) <- c(0, 610, 0, 870)
#' r1= sin(r/100)
#' r2= cos(r/100)
#' compareMaps(list(r,r1), mfrow=c(1,2))
#' compareMaps(list(r,r1,r2, r1+r2), mfrow=c(2,2))
#' compareMaps(list(r,r1,r2, r1+r2), mfrow=c(2,2), contour = TRUE)
compareMaps <- function(r, mfrow, contour=FALSE, ...){
  nr = length(r)
  is.Raster <- function(x)  {
    return((class(x)=="RasterLayer" || class(x)=="RasterBrick" || class(x)=="RasterStack"))
  }
  par(mfrow = mfrow, ...)
  for(i in 1:nr){
    raster::plot(r[[i]], axes=FALSE, box=FALSE)
    if(contour){
      if(is.Raster(r[[i]])){
        raster::contour(r[[i]], add=TRUE)
      }
    }
  }
  grid()
  par(mfrow=c(1,1))
}

#' Plot SpatialPolygons maps
#' \code{plot.sp} 
#' @param x Spatial*
#' @param zcol The column to plot
#' @param col.fun Functions of coloring.
#' @param ncol Number of colors in color function.
#' @export
plot_sp <-function(x, zcol, col.fun = topo.colors, ncol = 20){
  z = x@data[, zcol]
  col = col.fun(length(z))
  ord = order(z)
  raster::plot(x[ord, ], col=col)
}

#' Plot animation maps
#' \code{plot.animate} 
#' @param x Time-series data.
#' @param id Which row(s) to plot in the animation.
#' @param nmap Number of maps to be plot.
#' @param rlist RasterStack of the maps.
#' @export
plot_animate <- function(x, id = NULL, nmap = 10, rlist = NULL ){
  if(is.null(rlist)){
    nr = nrow(x)
    if(is.null(id)){
      if(nr < 20){
        id = 1:nr
      }else{
        id = sort(unique( round(seq(1, nrow(x), length.out = nmap) ) ) )
      }
    }
    y = x[id,]
    tx = paste(time(y) )
    rlist = raster::stack(apply(y, 1, FUN=function(xr){ MeshData2Raster(xr)  } ))
    names(rlist) = tx
  }else{
    rlist = raster::stack(rlist)
  }
  par(mfrow=c(1,1))
  raster::animate(rlist)
  ret = rlist
}


#' Highlight elements with id.
#' \code{highlight_id} 
#' @param EleID Index of element
#' @param RivID Index of river reaches
#' @export
highlight_id <- function(EleID=NULL, RivID=NULL){
  spm=sp.mesh2Shape()
  raster::plot(spm)
  if(!is.null(EleID)){
    raster::plot(spm[EleID, ], add=T, col=2)
  }
  inpath=shud.filein()['inpath']
  fn = file.path(inpath, 'gis', 'riv.shp')
  if(file.exists(fn)){
    spr=rgdal::readOGR()
    raster::plot(spr, add=T, col=3)
    if(!is.null(RivID)){
      raster::plot(spr[RivID, ], add=T, col=2, lwd=2)
    }
  }
}
