
#' Automatically run SHUDgis based on the input data and parameters.
#' \code{shud_auto_build}
#' @param indata  Input data, list of data.
#' @param forcfiles Filenames of .csv forcing data
#' @param prjname Projectname
#' @param outdir Output directory
#' @param a.max  Maximum area of triangles. Unit=m2
#' @param q.min Minimum of angles for each triangle
#' @param tol.riv Tolerance to simplify river lines, in meter
#' @param tol.wb Tolerance to simplify watershed boundary, in meter
#' @param tol.len Tolerance to simplyfy river by longituanal length, meter.
#' @param AqDepth Depth of aquifer bottom, which is surface elemvation minus bedrock elevation. meter.
#' @param years Years to generate the LAI, Meltfactor, Roughness Length, which must be same as the period of forcing data.
#' @param clean Whether clean the existing model files in output directory.
#' @param cfg.para  model configuration, parameter
#' @param cfg.calib model calibration
#' @param rm.outlier Whether to remove the outlier in soil/geol;
#' @param mf Meltfactor
#' @param quiet Whether to ask confirmation when multiple outlets exist.
#' @importFrom grDevices dev.off graphics.off png rgb topo.colors
#' @importFrom graphics grid hist lines par plot points
#' @importFrom methods as
#' @importFrom stats dist rnorm time
#' @importFrom utils read.table
#' @return \code{SHUD.mesh}
#' @export
#' @examples
#' data(sac)
#' indata = sac
#' sp.forc = indata[['forc']]
#' forc.fns = paste0(sp.forc@data[, 'NLDAS_ID'], '.csv')
#' forc.fns
# shud_auto_build(sac, forcfiles = forc.fns, outdir=tempdir())
shud_auto_build <- function(
  indata,
  forcfiles,
  prjname = 'sac',
  outdir = getwd(),
  a.max = 1e6 * .2,
  q.min = 33,
  tol.riv = 0,
  tol.wb = 200,
  tol.len = 0,
  AqDepth = 10,
  years = 2000:2010,
  clean = TRUE,
  cfg.para = shud.para(nday = 365 * length(years) + round(length(years) / 4)),
  cfg.calib = shud.calib(),
  mf = MeltFactor(years = years),
  rm.outlier = TRUE, 
  quiet = FALSE){
  msg = paste0('autoBuildModel(', prjname, '):: ')
  message(msg, 'Automatic build the SHUD model.')
  dir.out = file.path(outdir)
  dir.create(dir.out, showWarnings = FALSE, recursive = TRUE)
  ny = length(years)
  
  fin <- shud.filein(prjname, inpath = dir.out, outpath = dir.out)
  # x=list.files(dir.out, pattern = utils::glob2rx(paste0(prjname, '.*.*')), full.names = T)
  
  pngout = file.path(dir.out, 'fig')
  gisout = file.path(dir.out, 'gis')
  dir.create(dir.out, showWarnings = FALSE, recursive = TRUE)
  dir.create(pngout, showWarnings = FALSE, recursive = TRUE)
  dir.create(gisout, showWarnings = FALSE, recursive = TRUE)
  
  message(msg, 'Get data ...')
  # data(indata)
  wbd = indata[['wbd']]
  riv = indata[['riv']]
  dem = indata[['dem']]
  rsoil = indata[['rsoil']]
  rgeol = indata[['rsoil']]
  rlc = indata[['rlc']]
  alc = indata[['alc']]
  sp.forc = indata[['forc']]

  # Check for wrapped objects and unwrap them (handles data loading issue)
  if (inherits(dem, "PackedSpatRaster")) dem <- terra::unwrap(dem)
  if (inherits(rsoil, "PackedSpatRaster")) rsoil <- terra::unwrap(rsoil)
  if (inherits(rgeol, "PackedSpatRaster")) rgeol <- terra::unwrap(rgeol)
  if (inherits(rlc, "PackedSpatRaster")) rlc <- terra::unwrap(rlc)
  if (inherits(wbd, "PackedSpatVector")) wbd <- terra::unwrap(wbd)
  if (inherits(riv, "PackedSpatVector")) riv <- terra::unwrap(riv)
  if (inherits(sp.forc, "PackedSpatVector")) sp.forc <- terra::unwrap(sp.forc)

  if(!quiet) terra::plot(sp.forc)
  
  message(msg, 'Data buffer')
  # Use sf buffer
  wbbuf = sf::st_buffer(sf::st_as_sf(wbd), dist = max(c(2000, tol.wb)))
  # crop raster dem
  dem = terra::crop(terra::rast(dem), terra::vect(wbbuf))
  
  png(filename = file.path(pngout, 'data_0.png'), height = 11, width = 11, res = 100, units = 'in')
  terra::plot(dem); plot(sf::st_geometry(sf::st_as_sf(wbd)), add = TRUE, border = 2, lwd = 2); plot(sf::st_geometry(sf::st_as_sf(riv)), add = TRUE, lwd = 2, col = 4)
  dev.off()
  
  message(msg, 'River network')
  if(tol.riv > 0){
    riv.s1 = sf::st_simplify(sf::st_as_sf(riv), dTolerance = tol.riv, preserveTopology = TRUE)
  }else{
    riv.s1 = sf::st_as_sf(riv)
  }
  if(tol.len > 0){
    riv.s2 = sp.simplifyLen(riv.s1, tol.len)
  }else{
    riv.s2 = riv.s1
  }
  # plot(st_geometry(riv.s1)); plot(st_geometry(riv.s2), add=T, col=3)
  message(msg, 'Simplying boundary and river network...')
  wb.dis = sf::st_union(sf::st_as_sf(wbd))
  if(tol.wb > 0){
    wb.s1 = sf::st_simplify(wb.dis, dTolerance = tol.wb, preserveTopology = TRUE)
  }else{
    wb.s1 = wb.dis
  }
  # wb.s2 = sp.simplifyLen(wb.s1, tol.len)
  wb.s2 = wb.s1
  
  png(filename = file.path(pngout, 'data_1.png'), height = 11, width = 11, res = 100, units = 'in')
  terra::plot(dem); plot(sf::st_geometry(wb.s2), add = TRUE, border = 2, lwd = 2);
  plot(sf::st_geometry(riv.s2), add = TRUE, lwd = 2, col = 4)
  dev.off()
  
  
  # shp.riv =raster::crop(riv.simp, wb.simp)
  # shp.wb = raster::intersect( wb.simp, riv.simp)
  wb.simp = wb.s2
  riv.simp = riv.s2
  message(msg, 'Triangulating ... ')
  tri = shud.triangle(wb = wb.simp, q = q.min, a = a.max)
  ncells = nrow(tri$T)
  # graphics.off()
  # plot(tri, asp=1)
  message(msg, 'Number of Triangles = ', nrow(tri$T))
  # cin = readline(prompt = "See the plot. Go on?\n")
  # if( cin %in% c('n','N')){
  #   return(NULL);
  # }
  # generate SHUD .mesh
  pm = shud.mesh(tri, dem = dem, AqDepth = AqDepth)
  spm = mesh_to_sf(pm, crs = sf::st_crs(wb.s2))
  writeshape(spm, sf::st_crs(sf::st_as_sf(wbd)), file = file.path(gisout, 'domain'))
  if(!quiet) plot(sf::st_geometry(spm))
  png(filename = file.path(pngout, 'data_Mesh.png'), height = 11, width = 11, res = 100, units = 'in')
  plot(sf::st_geometry(spm))
  dev.off()
  
  # generate SHUD .att
  if(is.list(rsoil)){ 
    xsoil = 1:ncells
  }else{
    xsoil = rsoil
  }
  if(is.list(rgeol)){ 
    xgeol = 1:ncells
  }else{
    xgeol = rgeol
  }
  pa = shud.att(tri, r.soil = xsoil, r.geol = xgeol, r.lc = rlc, r.forc = sp.forc )
  
  write_forc(forcfiles,  file = fin['md.forc'], startdate = paste0(min(years), '0101'))
  
  # generate SHUD .riv
  AA = as.numeric(sf::st_area(sf::st_as_sf(wbd))) * 1e-6
  pr = shud.river(riv.simp, dem, AREA = AA)
  oid = getOutlets(pr)
  message(msg, 'Number of Rivers = ', nrow(pr@river))
  # Correct river slope to avoid negative slope
  
  # shud.riv to Shapefile
  spr = riv.simp
  if(is(spr, 'SpatialLines')){
    spr = sp::SpatialLinesDataFrame(spr, 
                                   data.frame('Riv' = pr@river, 'Tp' = pr@rivertype[pr@river$Type,]),
                                   match.ID = FALSE)
  }else if (inherits(spr, "sf")){
    spr <- cbind(spr, data.frame('Riv' = pr@river, 'Tp' = pr@rivertype[pr@river$Type,]))
  }else{
    spr@data = data.frame(spr@data, 'Riv' = pr@river, 'Tp' = pr@rivertype[pr@river$Type,], match.ID = FALSE)
  }
  writeshape(spr, sf::st_crs(sf::st_as_sf(wbd)), file = file.path(gisout, 'river'))
  
  if(length(oid) > 1){
    message(msg, "There are ", length(oid), ' outlets in streams')
    dev.off();
    plot(sf::st_geometry(sf::st_as_sf(spr))); plot(sf::st_geometry(sf::st_as_sf(spr)[oid, ]), add = TRUE, col = 2, lwd = 3)
    if(!quiet){
      flag = readline('Continue?')
      if(flag =='N' | flag =='n'){
        message(msg, 'Exit the autoBuildModel')
        stop('Exit')
      }else{
        message(msg, 'Continure the ModelBuilder')
      }
    }
  }
  # Cut the rivers with triangles
  message(msg, 'Cut the rivers with triangles.')
  sp.seg = shud.rivseg(spm, sf::st_as_sf(spr))
  writeshape(sp.seg, sf::st_crs(sf::st_as_sf(wbd)), file = file.path(gisout, 'seg'))
  
  prs = sp.seg
  
  # Generate initial condition
  message(msg, 'Initial conditions')
  pic = shud.ic(nrow(pm@mesh), nrow(pr@river), AqD = AqDepth)
  
  # Generate shapefile of river
  # spp.riv = sp.riv2shp(pr);
  
  
  # Generate shapefile of mesh domain
  message(msg, 'Generate shapefile of mesh domain')
  sp.dm = mesh_to_sf(pm)
  xcentroid = getCentroid(pm)
  png(filename = file.path(pngout, 'data_2.png'), height = 11, width = 11, res = 100, units = 'in')
  zz = sp.dm$Zsurf
  ord = order(zz)
  col = grDevices::terrain.colors(length(sp.dm$Zsurf))
  plot(sf::st_geometry(sp.dm[ord, ]), col = col)
  plot(sf::st_geometry(sf::st_as_sf(spr)), col = pr@river$Type + 1 , add = TRUE, lwd = 3)
  dev.off()
  
  #soil/geol/landcover
  message(msg, 'Generate land cover')
  lc = unlist(alc)
  para.lc = lc.NLCD(lc)
  
  if(is.list(rsoil)){
    asoil = cbind(terra::extract(terra::rast(rsoil[[1]]), xcentroid[, 2:3])[, 2], 
                  terra::extract(terra::rast(rsoil[[2]]), xcentroid[, 2:3])[, 2], 
                  terra::extract(terra::rast(rsoil[[3]]), xcentroid[, 2:3])[, 2], 
                  terra::extract(terra::rast(rsoil[[4]]), xcentroid[, 2:3])[, 2])
  }else{
    asoil = indata[['asoil']]
  }
  if(is.list(rgeol)){
    ageol = cbind(terra::extract(terra::rast(rgeol[[1]]), xcentroid[, 2:3])[, 2], 
                  terra::extract(terra::rast(rgeol[[2]]), xcentroid[, 2:3])[, 2], 
                  terra::extract(terra::rast(rgeol[[3]]), xcentroid[, 2:3])[, 2], 
                  terra::extract(terra::rast(rgeol[[4]]), xcentroid[, 2:3])[, 2])
  }else{
    ageol = indata[['asoil']]
  }
  para.soil = PTF.soil(asoil, rm.outlier = rm.outlier)
  para.geol = PTF.geol(ageol, rm.outlier = rm.outlier)
  
  
  # 43-mixed forest in NLCD classification
  # 23-developed, medium
  # 81-crop land
  # 11-water
  lr = fun.lairl(lc, years = years)
  graphics.off()
  png(filename = file.path(pngout, 'data_lairl.png'), height = 11, width = 11, res = 100, units = 'in')
  graphics::par(mfrow = c(2,1))
  col = 1:length(lc)
  
  zoo::plot.zoo(lr$LAI, col = col, main = 'LAI');
  graphics::legend('top', paste0(lc), col = col, lwd = 1)
  zoo::plot.zoo(lr$RL, col = col, main = 'Roughness Length');
  graphics::legend('top', paste0(lc), col = col, lwd = 1)
  dev.off()
  message(msg, 'Write SHUD model input files.')
  write_tsd(lr$LAI, file = fin['md.lai'])
  write_tsd(lr$RL, file = fin['md.rl'])
  
  write_tsd(mf, file = fin['md.mf'])
  
  # write SHUD input files.
  write_mesh(pm, file = fin['md.mesh'])
  write_river(pr, file = fin['md.riv'])
  write_ic(pic, file = fin['md.ic'])
  
  write_df(pa, file = fin['md.att'])
  write_df(prs, file = fin['md.rivseg'])
  write_df(para.lc, file = fin['md.lc'])
  write_df(para.soil, file = fin['md.soil'])
  write_df(para.geol, file = fin['md.geol'])
  
  write_config(cfg.para, fin['md.para'])
  write_config(cfg.calib, fin['md.calib'])
  message(msg, 
          'Ncell=', nrow(pm@mesh), '\t',
          'Nriv=', nrow(pr@river), '\t',
          'Nseg=', nrow(prs), '\n'
  )
  pm <- pm
}
# 
# 
# pm=autoBuildModel(indata,
#                   forcfiles = forc.fns,
#                   prjname = prjname,
#                   outdir= dirs$model,
#                   a.max = a.max,
#                   AqDepth = 20,
#                   tol.len = 0,
#                   tol.riv = 0,
#                   tol.wb = sqrt(a.max)/3,
#                   quiet = TRUE)
