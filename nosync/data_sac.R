library(sf)
library(terra)
library(rSHUD)
data(sac)

convert_raster <- function(r, fallback_val=1) {
  if (inherits(r, 'RasterLayer')) {
    if (!raster::inMemory(r)) {
      tryCatch({
        r <- raster::readAll(r)
      }, error = function(e) {
        warning(paste('Could not read raster data for', names(r), 'creating a fallback constant raster.'))
        r_new <- raster::raster(raster::extent(r), res=raster::res(r), crs=raster::crs(r))
        r_new[] <- fallback_val
        r <<- r_new
      })
    }
    return(terra::rast(r))
  }
  return(r)
}

sac$dem = convert_raster(sac$dem)
sac$rsoil = convert_raster(sac$rsoil)
sac$rgeol = convert_raster(sac$rgeol)
sac$rlc = convert_raster(sac$rlc)

if(inherits(sac$riv, 'Spatial')) sac$riv = sf::st_as_sf(sac$riv)
if(inherits(sac$wbd, 'Spatial')) sac$wbd = sf::st_as_sf(sac$wbd)
if(inherits(sac$forc, 'Spatial')) sac$forc = sf::st_as_sf(sac$forc)

unlink('data/sac.rda')
usethis::use_data(sac, overwrite = TRUE)
