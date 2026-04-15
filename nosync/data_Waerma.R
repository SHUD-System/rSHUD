library(raster)
library(rSHUD)
library(sp)
library(rgeos)
library(sf)
library(terra)
library(sf)

# Load existing waerma data to convert to sf/terra if original files missing
data(waerma)

convert_raster <- function(r, fallback_val=1) {
  if (inherits(r, "RasterLayer")) {
    if (!raster::inMemory(r)) {
      tryCatch({
        r <- raster::readAll(r)
      }, error = function(e) {
        warning(paste("Could not read raster data for", names(r), "creating a fallback constant raster."))
        r_new <- raster::raster(raster::extent(r), res=raster::res(r), crs=raster::crs(r))
        r_new[] <- fallback_val
        r <<- r_new
      })
    }
    return(terra::rast(r))
  }
  return(r)
}

waerma$dem = convert_raster(waerma$dem)
waerma$soil = convert_raster(waerma$soil)
waerma$geol = convert_raster(waerma$geol)
waerma$lc = convert_raster(waerma$lc)

if(inherits(waerma$riv, "Spatial")) waerma$riv = sf::st_as_sf(waerma$riv)
if(inherits(waerma$wbd, "Spatial")) waerma$wbd = sf::st_as_sf(waerma$wbd)
if(inherits(waerma$meteosite, "Spatial")) waerma$meteosite = sf::st_as_sf(waerma$meteosite)

# Try to load raw files if they exist (local developer environment)
r.dem.file = '/Users/leleshu/CloudDrive/Commondata/China/黄河源/WaErMa/box.tif'
if(file.exists(r.dem.file)) {
  r.dem = terra::rast(r.dem.file)
  sp.wbd = sf::st_read('/Users/leleshu/CloudDrive/Commondata/China/黄河源/WaErMa/delineation/wbd_dem.shp', quiet=TRUE)
  sp.riv = sf::st_read('/Users/leleshu/CloudDrive/Commondata/China/黄河源/WaErMa/delineation/stm_dem.shp', quiet=TRUE)
  r.geol=terra::rast('/Users/leleshu/CloudDrive/文章天下事/Drafts/2020_rSHUD/Waerma/ETV/soil.tif')
  r.soil=terra::rast('/Users/leleshu/CloudDrive/文章天下事/Drafts/2020_rSHUD/Waerma/ETV/geol.tif')
  r.lc=terra::rast('/Users/leleshu/CloudDrive/文章天下事/Drafts/2020_rSHUD/Waerma/ETV/landuse.tif')
  sp.meteo = sf::st_read('/Users/leleshu/CloudDrive/文章天下事/Drafts/2020_rSHUD/Waerma/Modeling/GCS/meteo.shp', quiet=TRUE)

  at.lc = read.df('/Users/leleshu/CloudDrive/文章天下事/Drafts/2020_rSHUD/Waerma/ETV/landuse.csv')[[1]]
  at.soil = read.df('/Users/leleshu/CloudDrive/文章天下事/Drafts/2020_rSHUD/Waerma/ETV/soil.csv')[[1]]
  at.geol = read.df('/Users/leleshu/CloudDrive/文章天下事/Drafts/2020_rSHUD/Waerma/ETV/geol.csv')[[1]]

  tsd.lai = read.tsd('/Users/leleshu/CloudDrive/文章天下事/Drafts/2020_rSHUD/Waerma/ETV/lai.csv')[[1]]

  fns=list.files('/Users/leleshu/CloudDrive/文章天下事/Drafts/2020_rSHUD/Waerma/ETV/TSD/', 
                 pattern = glob2rx('*.csv'), full.names = TRUE)
  tsd.forc=list(); nf=length(fns)
  for(i in 1:nf){
    tsd.forc[[i]] = read.tsd(fns[i])[[1]][1:(8*731), ]
  }
  names(tsd.forc) = basename(fns)

  waerma = list(dem = r.dem, 
            wbd = sp.wbd, 
            riv = sp.riv,
            meteosite= sp.meteo,
            soil=r.soil,
            geol=r.geol,
            lc = r.lc,
            att=list(soil=at.soil, geol=at.geol, lc=at.lc),
            tsd=list(forc=tsd.forc, lai=tsd.lai)
            )
}

unlink('data/waerma.rda')
usethis::use_data(waerma, overwrite = TRUE)
