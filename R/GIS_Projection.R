#' Build an Albers Equal Area projection based on the spatial data or extent.
#' \code{crs.Albers}
#' @param spx Spatial Data
#' @param ext  extent of the data.
#' @return CRS in Albers Equal Area projection.
#' @export
#' @examples
#' crs.Albers(ext=c(90, 100, 35, 40))
crs.Albers <- function(spx, ext=NULL){
  if(is.null(ext)){
    crs0 = sp::CRS('+init=epsg:4326')
    sp.gcs = sp::spTransform(spx, crs0)
    ext= raster::extent(sp.gcs)
  }
  x0 = round(mean(ext[1:2]), 1)
  y0 = round(mean(ext[1:2+2]), 1)
  dx = round(diff(ext[1:2]), 2)
  dy = round(diff(ext[1:2+2]), 2)
  if(dx < 1){
    my = 0.25
  }else{
    my=round(dy/4, 2)
  }
  lat1 = y0 + my
  lat2 = y0 - my
  str = paste0(' +proj=', 'aea',
               ' +lat_1=', lat1,
               ' +lat_2=', lat2,
               ' +lon_0=', x0,
               ' +datum=', 'WGS84',
               " +units=", "m"
  )
  print(str)
  ret = sp::CRS(str)
  return(ret)
}
#' Build a Lambert Equal Area projection based on the spatial data or extent.
#' \code{crs.Lambert}
#' @param spx Spatial Data
#' @param ext  extent of the data.
#' @return CRS in Lambert
#' @export
#' @examples
#' crs.Lambert(ext=c(90, 100, 35, 40))
crs.Lambert <- function(spx, ext=NULL){
  if(is.null(ext)){
    crs0 = sp::CRS('+init=epsg:4326')
    sp.gcs = sp::spTransform(spx, crs0)
    ext= raster::extent(sp.gcs)
  }
  x0 = round(mean(ext[1:2]), 1)
  y0 = round(mean(ext[1:2+2]), 1)
  dx = round(diff(ext[1:2]), 2)
  dy = round(diff(ext[1:2+2]), 2)
  
  ret = sp::CRS(paste0('+proj=', 'leac',
                       ' +lat_1=', y0,
                       ' +lon_0=', x0,
                       ' +datum=', 'WGS84',
                       " +units=", "m"))
  return(ret)
}

#' Get UTM zone frome longitutde
#' \code{crs.long2utmZone} 
#' @param  lon Longitude in degree
#' @return Zone number
#' @export
#' @examples 
#' long = seq(-180, 179, 6)+1
#' crs.long2utmZone(long)
crs.long2utmZone <- function(lon){
  zid =floor((lon + 180) / 6) + 1  
  return(zid)
}
#' Get UTM projection parameters in CRS, from longtitude and latitude.
#' \code{crs.long2utm} 
#' @param  lon Longitude in degree
#' @param lat Latitude in degree, default = 30, North.
#' @return UTM projection parameters in CRS
#' @export
#' @examples 
#' long = seq(-180, 179, 6)+1
#' x =crs.long2utm(long)
crs.long2utm <- function(lon, lat=30){
  if(any(lon > 180 | lon < -180)){
    message('\nRange of data:', paste(range(lon), collapse = ', '))
    message('Acceptable range of longitude is (-180, 180) \n')
    stop()
  }
  if(length(lon)>1){
    x = lapply(lon, crs.long2utm)
  }else{
    zid =crs.long2utmZone(lon)
    if(lat>=0){
      str = paste0('+proj=utm +zone=', zid, ' +datum=WGS84 +units=m +no_defs')
    }else{
      str = paste0('+proj=utm +zone=', zid, ' +south +datum=WGS84 +units=m +no_defs')
    }
    x = sp::CRS(str)
  }
  return(x)
}



