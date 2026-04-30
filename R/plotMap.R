#' Plot map of time-series data
#' \code{ts2map} 
#' @param x Time-series data
#' @param fun Functions to summarize the TS data.
#' @param raster Whether to make output as Raster
#' @param ... More options in fun
#' @export
ts2map <- function(x, fun = mean, raster = TRUE, ...){
  y = apply(x, 2, fun, ...)
  if(raster){
    r = MeshData2Raster(y)
    terra::plot(r)
  }else{
    dbf = cbind(y); colnames(dbf) = 'tsvalue'
    r = mesh_to_sf(dbf = y)
    plot_polygons(r, field = 'tsvalue')
  }
  return(r)
}

#' Plot multiple maps
#' \code{compareMaps} 
#' @param r List of raster or SpatialData
#' @param mfrow mfrow in par()
#' @param contour Whether plot the contour.
#' @param ... More options in par()
#' @noRd
#' @examples 
#' library(raster)
#' data(volcano)
#' r <- raster(volcano)
#' extent(r) <- c(0, 610, 0, 870)
#' r1 = sin(r / 100)
#' r2 = cos(r / 100)
#' compareMaps(list(r, r1), mfrow = c(1, 2))
#' compareMaps(list(r, r1, r2, r1 + r2), mfrow = c(2, 2))
#' compareMaps(list(r, r1, r2, r1 + r2), mfrow = c(2, 2), contour = TRUE)
compareMaps_legacy <- function(r, mfrow, contour = FALSE, ...){
  nr = length(r)
  is_raster_like <- function(x) inherits(x, c("SpatRaster", "RasterLayer", "RasterBrick", "RasterStack"))
  par(mfrow = mfrow, ...)
  for(i in 1:nr){
    ri <- if (inherits(r[[i]], "SpatRaster")) r[[i]] else if (is_raster_like(r[[i]])) terra::rast(r[[i]]) else r[[i]]
    if (inherits(ri, "SpatRaster")) {
      terra::plot(ri, axes = FALSE, box = FALSE)
    } else {
      plot(ri, axes = FALSE, box = FALSE)
    }
    if(contour){
      if(inherits(ri, "SpatRaster")){
        terra::contour(ri, add = TRUE)
      }
    }
  }
  grid()
  par(mfrow = c(1, 1))
}



#' Plot animation maps
#' \code{plot.animate} 
#' @param x Time-series data.
#' @param id Which row(s) to plot in the animation.
#' @param nmap Number of maps to be plot.
#' @param rlist RasterStack of the maps.
#' @export
plot_animate <- function(x, id = NULL, nmap = 10, rlist = NULL){
  if(is.null(rlist)){
    nr = nrow(x)
    if(is.null(id)){
      if(nr < 20){
        id = 1:nr
      }else{
        id = sort(unique(round(seq(1, nrow(x), length.out = nmap))))
      }
    }
    y = x[id, ]
    tx = paste(time(y))
    rlist = terra::rast(lapply(seq_len(nrow(y)), function(i) MeshData2Raster(y[i, ])))
    names(rlist) = tx
  }else{
    rlist = if (inherits(rlist, "SpatRaster")) rlist else terra::rast(rlist)
  }
  par(mfrow = c(1, 1))
  terra::animate(rlist)
  ret = rlist
}


#' Highlight elements with id (DEPRECATED - INCOMPLETE)
#' 
#' @description
#' \strong{This function is deprecated and incomplete.} It uses legacy spatial
#' libraries (raster/sp) and has incomplete implementation. This function will
#' be removed in a future version.
#' 
#' Users should create custom highlighting plots using modern spatial packages:
#' - Use \code{terra::plot()} for raster visualization
#' - Use \code{sf::plot()} or \code{ggplot2} for vector visualization
#' 
#' @param EleID Index of element to highlight
#' @param RivID Index of river reaches to highlight
#' @return NULL (plots to current device)
#' @export
#' @section Deprecated:
#' This function is deprecated due to incomplete implementation and use of
#' legacy spatial libraries. It will be removed in version 2.4.0.
#' 
#' @examples
#' \dontrun{
#' # This function is deprecated
#' # Use modern plotting instead:
#' library(sf)
#' library(terra)
#' 
#' # Read and plot mesh
#' mesh_sf <- st_read("mesh.shp")
#' plot(mesh_sf)
#' 
#' # Highlight specific elements
#' plot(mesh_sf[EleID, ], col = "red", add = TRUE)
#' }
highlight_id <- function(EleID = NULL, RivID = NULL){
  .Deprecated(msg = paste(
    "highlight_id() is deprecated and incomplete.",
    "Use terra::plot() or sf::plot() for custom visualization instead.",
    "This function will be removed in version 2.4.0."
  ))
  
  warning("This function uses legacy spatial libraries and may not work correctly.")
  
  spm <- mesh_to_sf()
  plot(sf::st_geometry(spm))
  if(!is.null(EleID)){
    plot(sf::st_geometry(spm[EleID, ]), add = TRUE, col = 2)
  }
  inpath <- shud.filein()["inpath"]
  fn <- file.path(inpath, "gis", "riv.shp")
  if(file.exists(fn)){
    spr <- sf::st_read(fn, quiet = TRUE)
    plot(sf::st_geometry(spr), add = TRUE, col = 3)
    if(!is.null(RivID)){
      plot(sf::st_geometry(spr[RivID, ]), add = TRUE, col = 2, lwd = 2)
    }
  }
}
