#' Do Watershed Delineation
#' \code{watershedDelineation} 
#' @param fn.wbd  file of watershed SpatialPolygons
#' @param fnr.dem  file of DEM in raster
#' @param fsp.outlets  file of watershed outlet
#' @param dirout  dir for out data(dem_filled.tif, dem_stm.shp, dem_wbd.shp, dem_outlets.shp) saved. Default = tempdir()
#' @param dirtemp  dir for temporary data save. Default = tempdir()
#' @param FlowAccCell.min  Minimum number of cells for river generation. Default  = Number of DEM / 100
#' @return list of path to dem-filled, wbd, riv, outlet
#' @export
watershedDelineation <- function(fn.wbd, 
                                 fnr.dem, 
                                 fsp.outlets= NULL,
                                 dirout = tempdir(),
                                 dirtemp = tempdir(),
                                 FlowAccCell.min = NULL, 
                                 plot=FALSE
){
  # library(raster)
  # library(whitebox)
  # library(terra)
  dir.create(dirtemp, showWarnings = FALSE, recursive = TRUE)
  whitebox::wbt_init()
  fnr.filled = file.path(dirout, "dem.filled.tif")
  fsp.stm = file.path(dirout, 'dem_stm.shp')
  fsp.wbd = file.path(dirout, 'dem_wbd.shp')
  
  fnr.smoothing = file.path(dirtemp, "dem.smoothed.tif")
  fnr.breached = file.path(dirtemp, "dem.breached.tif")
  fnr.d8fa = file.path(dirtemp, 'd8fa.tif')
  fnr.d8ptr = file.path(dirtemp, 'd8ptr.tif')
  fnr.stm =  file.path(dirtemp, 'stream.tif')
  fnr.stmclip =  file.path(dirtemp, 'stream_clip.tif')
  fnr.subs =  file.path(dirtemp, 'subbasins.tif')
  fnr.flood = file.path(dirtemp, 'flood.tif')
  fnr.wbd = file.path(dirtemp, 'wbd.tif')
  
  plotr <- function(x, ...){raster::plot(raster(x), ...)}
  plotv <- function(x, ...){raster::plot(rgdal::readOGR(x), ...)}
  
  # # 1. Fill Pits
  # writelog(paste0('1. Fill Pits'), caller=caller)
  whitebox::wbt_feature_preserving_smoothing(dem = fnr.dem, output = fnr.smoothing, filter = 9)
  whitebox::wbt_breach_depressions_least_cost(dem = fnr.smoothing,  output = fnr.breached,  dist = 5,  fill = TRUE)
  whitebox::wbt_fill_depressions_wang_and_liu(dem = fnr.breached, output = fnr.filled)
  plot(raster(fnr.filled))
  
  
  # # 2. Flow Accumulation and Pointer Grids
  # writelog(paste0('2. Flow Accumulation and Pointer Grids'), caller=caller)
  whitebox::wbt_d8_flow_accumulation(input = fnr.filled, output = fnr.d8fa)
  whitebox::wbt_d8_pointer(dem = fnr.filled, output = fnr.d8ptr)
  # plot(raster(fnr.d8fa))
  
  # 3. Watershed.
  # writelog(paste0('3. Watershed.'), caller=caller)
  do_outlets <- function(){
    r = raster(fnr.d8fa)
    sp.wbd = rgdal::readOGR(fn.wbd)
    spp = sp::spTransform(sp.wbd, CRSobj = crs(r))
    r = raster::mask(r, spp)
    # plot(r)
    maxval = cellStats(r, max, na.rm=TRUE)
    idx = which.max(Which(r >= maxval))
    # idx = which.max(r, na.rm=TRUE)
    ll.outlets = xyFromCell(r,idx)
    sp.outlets = rSHUD::xy2shp(ll.outlets)
    rSHUD::writeshape(sp.outlets, file=fsp.outlets)
  }
  if(is.null(fsp.outlets) ){
    fsp.outlets = file.path(dirout, 'dem_outlets.shp')
    do_outlets()
  }else{
    if(file.exists(fsp.outlets) ){
      #void
    }else{
      fsp.outlets = file.path(dirout, 'dem_outlets.shp')
      do_outlets()
    }
  }
  # plotr(fnr.wbd)
  # plotv(fsp.wbd, add=T)
  
  whitebox::wbt_watershed(d8_pntr = fnr.d8ptr, pour_pts = fsp.outlets, output = fnr.wbd)
  whitebox::wbt_raster_to_vector_polygons(fnr.wbd, output = fsp.wbd)
  
  # writelog(paste0('4. Extract Streams'), caller=caller)
  # # 4. Extract Streams
  if(is.null(FlowAccCell.min)){
    r=raster::raster(fnr.dem)
    FlowAccCell.min = max(10, round(length(r)/500) )
  }else{
    FlowAccCell.min = FlowAccCell.min
  }
  whitebox::wbt_extract_streams(flow_accum = fnr.d8fa, output = fnr.stm, 
                                threshold = FlowAccCell.min, command_only=F)
  whitebox::wbt_raster_streams_to_vector(streams=fnr.stm, d8_pntr = fnr.d8ptr, output = fsp.stm)
  sp.wbd = rgdal::readOGR(fsp.wbd)
  sp.stm = rgdal::readOGR(fsp.stm)
  spx = raster::crop(sp.stm, sp.wbd)
  writeshape(spx, fsp.stm)
  
  # writelog(paste0('Plot watershed_delineation.png'), caller=caller)
  if(plot){
    plotr(fnr.dem)
    plotv(fsp.wbd, add=T)
    plotv(fsp.outlets, add=T, col='red', cex=3)
    plotv(fsp.stm, add=T, col=4)
  }
  
  # writelog(paste0('Finished.'), caller=caller)
  ret  = list(dem = fnr.filled, 
              stm = fsp.stm, 
              wbd = fsp.wbd)
  return(ret)
}

