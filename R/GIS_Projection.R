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



