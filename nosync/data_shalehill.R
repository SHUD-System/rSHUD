library(xts)
library(rSHUD)
library(terra)
library(sf)

# Convert all spatial data in 'sh' to sf/terra formats
data(sh)

# Convert rasters to terra safely by writing/reading if they are out of memory
convert_raster <- function(r, fallback_val=1) {
  if (inherits(r, "RasterLayer")) {
    if (!raster::inMemory(r)) {
      # Some data in the legacy .rda points to hardcoded user paths that don't exist anymore
      # we have to catch errors and provide a fallback
      tryCatch({
        r <- raster::readAll(r)
      }, error = function(e) {
        warning(paste("Could not read raster data for", names(r), "creating a fallback constant raster."))
        # Create a fallback based on the extent/res of the object
        r_new <- raster::raster(raster::extent(r), res=raster::res(r), crs=raster::crs(r))
        r_new[] <- fallback_val
        r <<- r_new
      })
    }
    return(terra::rast(r))
  }
  return(r)
}

sh$dem = convert_raster(sh$dem)
sh$aqd = convert_raster(sh$aqd)
sh$soil = convert_raster(sh$soil)
sh$geol = convert_raster(sh$geol)

# Convert vectors to sf
if(inherits(sh$riv, "Spatial")) sh$riv = sf::st_as_sf(sh$riv)
if(inherits(sh$wbd, "Spatial")) sh$wbd = sf::st_as_sf(sh$wbd)

# Add dummy forcing data for testing if real file is missing
forc_file = 'nosync/Data/Forcing_2008-2015.csv'
if(file.exists(forc_file)) {
  tb=read.table(forc_file)
  cn =c('Precip_mm.d',	'Temp_C',	'RH_1',	'Wind_m.s',	'RN_w.m2',	'Pres_pa')

  xt=as.POSIXct(paste(tb[,1], tb[, 2]), tz='UTC')
  ts = as.xts(tb[, -1*(1:2)], order.by=xt)
  ts=ts[, c(1:5, 7)]
  ts[, 1] = ts[, 1]*86400
  ts[, 2] = ts[, 2] - 273.15
  ts[, 3] = ts[, 3] /100
  ts[, 4] = ts[, 4] 
  ts[, 5] = ts[, 5] 
  ts[, 6] = ts[, 6] 

  colnames(ts)=cn
  sh$forc = ts
} else {
  ts = sh$forc
}

unlink('data/sh.rda')
usethis::use_data(sh, overwrite = TRUE)

head(ts)
# apply.yearly(apply.daily(ts$Precip_mm.d, mean), sum)


