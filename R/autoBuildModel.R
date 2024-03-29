
#' Automatically run SHUDgis based on the input data and parameters.
#' \code{autoBuildModel}
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
#' sp.forc =indata[['forc']]
#' forc.fns = paste0(sp.forc@data[, 'NLDAS_ID'], '.csv')
#' forc.fns
# autoBuildModel(sac, forcfiles = forc.fns, outdir=tempdir())
autoBuildModel <- function(
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
  cfg.para = shud.para(nday = 365*length(years) +round(length(years)/4) ),
  cfg.calib = shud.calib(),
  mf = MeltFactor(years = years),
  rm.outlier = TRUE, 
  quiet = FALSE
){
  msg = paste0('autoBuildModel(', prjname, '):: ')
  message(msg, 'Automatic build the SHUD model.')
  dir.out = file.path(outdir)
  dir.create(dir.out, showWarnings = FALSE, recursive = TRUE)
  ny=length(years)
  
  fin <- shud.filein(prjname, inpath = dir.out, outpath=dir.out)
  # x=list.files(dir.out, pattern = utils::glob2rx(paste0(prjname, '.*.*')), full.names = T)
  
  pngout = file.path(dir.out, 'fig')
  gisout = file.path(dir.out, 'gis')
  dir.create(dir.out, showWarnings = F, recursive = T)
  dir.create(pngout, showWarnings = F, recursive = T)
  dir.create(gisout, showWarnings = F, recursive = T)
  
  message(msg, 'Get data ...')
  # data(indata)
  wbd=indata[['wbd']]
  riv=indata[['riv']]
  dem=indata[['dem']]
  rsoil=indata[['rsoil']]
  rgeol=indata[['rsoil']]
  rlc=indata[['rlc']]
  alc =indata[['alc']]
  sp.forc =indata[['forc']]
  raster::plot(sp.forc)
  
  message(msg, 'Data buffer')
  wbbuf = rgeos::gBuffer(wbd, width = max(c(2000, tol.wb) ) )
  dem = raster::crop(dem, wbbuf)
  
  png(filename = file.path(pngout, 'data_0.png'), height=11, width=11, res=100, units='in')
  raster::plot(dem); raster::plot(wbd, add=T, border=2, lwd=2); raster::plot(riv, add=T, lwd=2, col=4)
  dev.off()
  
  message(msg, 'River network')
  if(tol.riv > 0){
    riv.s1 = rgeos::gSimplify(riv, tol=tol.riv, topologyPreserve = TRUE)
  }else{
    riv.s1 = riv
  }
  if(tol.len > 0){
    riv.s2 = sp.simplifyLen(riv.s1, tol.len)
  }else{
    riv.s2 = riv.s1
  }
  # raster::plot(riv.s1); raster::plot(riv.s2, add=T, col=3)
  message(msg, 'Simplying boundary and river network...')
  wb.dis = rgeos::gUnionCascaded(wbd)
  if(tol.wb>0){
    wb.s1 = rgeos::gSimplify(wb.dis, tol=tol.wb, topologyPreserve = TRUE)
  }else{
    wb.s1 = wb.dis
  }
  # wb.s2 = sp.simplifyLen(wb.s1, tol.len)
  wb.s2= wb.s1
  
  png(filename = file.path(pngout, 'data_1.png'), height=11, width=11, res=100, units='in')
  raster::plot(dem); raster::plot(wb.s2, add=TRUE, border=2, lwd=2);
  raster::plot(riv.s2, add=T, lwd=2, col=4)
  dev.off()
  
  
  # shp.riv =raster::crop(riv.simp, wb.simp)
  # shp.wb = raster::intersect( wb.simp, riv.simp)
  wb.simp = wb.s2
  riv.simp = riv.s2
  message(msg, 'Triangulating ... ')
  tri = shud.triangle(wb=wb.simp,q=q.min, a=a.max)
  ncells = nrow(tri$T)
  # graphics.off()
  # plot(tri, asp=1)
  message(msg, 'Number of Triangles = ', nrow(tri$T))
  # cin = readline(prompt = "See the plot. Go on?\n")
  # if( cin %in% c('n','N')){
  #   return(NULL);
  # }
  # generate SHUD .mesh
  pm=shud.mesh(tri,dem=dem, AqDepth = AqDepth)
  spm = sp.mesh2Shape(pm, crs = raster::crs(wb.s2))
  writeshape(spm, raster::crs(wbd), file=file.path(gisout, 'domain'))
  raster::plot(spm)
  png(filename = file.path(pngout, 'data_Mesh.png'), height=11, width=11, res=100, units='in')
  raster::plot(spm)
  dev.off()
  
  # generate SHUD .att
  if(is.list(rsoil)){ 
    xsoil = 1:ncells
  }else{
    xsoil=rsoil
  }
  if(is.list(rgeol)){ 
    xgeol = 1:ncells
  }else{
    xgeol=rgeol
  }
  pa=shud.att(tri, r.soil = xsoil, r.geol = xgeol, r.lc = rlc, r.forc = sp.forc )
  
  write.forc(forcfiles,  file=fin['md.forc'], startdate = paste0(min(years), '0101'))
  
  # generate SHUD .riv
  AA = rgeos::gArea(wbd) * 1e-6
  pr=shud.river(riv.simp, dem, AREA = AA)
  oid = getOutlets(pr)
  message(msg, 'Number of Rivers = ', nrow(pr@river))
  # Correct river slope to avoid negative slope
  
  # shud.riv to Shapefile
  spr = riv.simp
  if(is(spr, 'SpatialLines')){
    spr =sp::SpatialLinesDataFrame(spr, 
                                   data.frame('Riv'=pr@river, 'Tp'=pr@rivertype[pr@river$Type,]),
                                   match.ID = FALSE)
  }else{
    spr@data = data.frame(spr@data, 'Riv'=pr@river, 'Tp'=pr@rivertype[pr@river$Type,], match.ID = FALSE)
  }
  writeshape(spr, raster::crs(wbd), file=file.path(gisout, 'river'))
  
  if(length(oid)>1){
    message(msg, "There are ", length(oid), ' outlets in streams')
    dev.off();
    raster::plot(spr); raster::plot(spr[oid, ], add=TRUE, col=2, lwd=3)
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
  sp.seg = sp.RiverSeg(spm, spr)
  writeshape(sp.seg, raster::crs(wbd), file=file.path(gisout, 'seg'))
  
  # Generate the River segments table
  prs = shud.rivseg(sp.seg)
  
  # Generate initial condition
  message(msg, 'Initial conditions')
  pic = shud.ic(nrow(pm@mesh), nrow(pr@river), AqD = AqDepth)
  
  # Generate shapefile of river
  # spp.riv = sp.riv2shp(pr);
  
  
  # Generate shapefile of mesh domain
  message(msg, 'Generate shapefile of mesh domain')
  sp.dm = sp.mesh2Shape(pm)
  xcentroid = getCentroid(pm)
  png(filename = file.path(pngout, 'data_2.png'), height=11, width=11, res=100, units='in')
  zz = sp.dm@data[,'Zsurf']
  ord=order(zz)
  col=grDevices::terrain.colors(length(sp.dm))
  raster::plot(sp.dm[ord, ], col = col)
  raster::plot(spr, col=pr@river$Type+1 , add=TRUE, lwd=3)
  dev.off()
  
  #soil/geol/landcover
  message(msg, 'Generate land cover')
  lc = unlist(alc)
  para.lc = PTF.lc(lc)
  
  if(is.list(r.soil)){
    asoil=cbind(extract(rsoil[[1]], xcentroid[, 2:3]), 
                extract(rsoil[[2]], xcentroid[, 2:3]), 
                extract(rsoil[[3]], xcentroid[, 2:3]), 
                extract(rsoil[[4]], xcentroid[, 2:3]))
  }else{
    asoil=indata[['asoil']]
  }
  if(is.list(r.geol)){
    ageol=cbind(extract(rgeol[[1]], xcentroid[, 2:3]), 
                extract(rgeol[[2]], xcentroid[, 2:3]), 
                extract(rgeol[[3]], xcentroid[, 2:3]), 
                extract(rgeol[[4]], xcentroid[, 2:3]))
  }else{
    ageol=indata[['asoil']]
  }
  para.soil = PTF.soil(asoil, rm.outlier = rm.outlier)
  para.geol = PTF.geol(ageol, rm.outlier = rm.outlier)
  
  
  # 43-mixed forest in NLCD classification
  # 23-developed, medium
  # 81-crop land
  # 11-water
  lr=fun.lairl(lc, years=years)
  graphics.off()
  png(filename = file.path(pngout, 'data_lairl.png'), height=11, width=11, res=100, units='in')
  graphics::par(mfrow=c(2,1))
  col=1:length(lc)
  
  zoo::plot.zoo(lr$LAI, col=col, main='LAI');
  graphics::legend('top', paste0(lc), col=col, lwd=1)
  zoo::plot.zoo(lr$RL, col=col, main='Roughness Length');
  graphics::legend('top', paste0(lc), col=col, lwd=1)
  dev.off()
  message(msg, 'Write SHUD model input files.')
  write.tsd(lr$LAI, file = fin['md.lai'])
  write.tsd(lr$RL, file = fin['md.rl'])
  
  write.tsd(mf, file=fin['md.mf'])
  
  # write SHUD input files.
  write.mesh(pm, file = fin['md.mesh'])
  write.riv(pr, file=fin['md.riv'])
  write.ic(pic, file=fin['md.ic'])
  
  write.df(pa, file=fin['md.att'])
  write.df(prs, file=fin['md.rivseg'])
  write.df(para.lc, file=fin['md.lc'])
  write.df(para.soil, file=fin['md.soil'])
  write.df(para.geol, file=fin['md.geol'])
  
  write.config(cfg.para, fin['md.para'])
  write.config(cfg.calib, fin['md.calib'])
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
