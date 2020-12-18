#' Prepare Input file of SHUD model
#' \code{shud.filein} 
#' @param projname Character, project name, default= PRJNAME which is a global variable.
#' @param inpath SHUD model input directory, default = inpath which is global variable
#' @param outpath SHUD model output directory, default = outpath which is global variable
#' @return Character of full path of input files for SHUD model
#' @export
shud.filein <- function(projname=get('PRJNAME', envir=.shud),
                        inpath = get('inpath', envir=.shud), 
                        outpath = get('outpath', envir=.shud)
                        ){
  dir.create(inpath, showWarnings = FALSE, recursive = TRUE)
  fn.mesh = file.path(inpath, paste0(projname, '.sp.mesh' ) );
  fn.att = file.path(inpath, paste0(projname, '.sp.att' ) );
  fn.pt = file.path(inpath, paste0(projname, '.sp.points' ) )
  fn.edge = file.path(inpath, paste0(projname, '.sp.edges' ) )
  fn.mdriv = file.path(inpath, paste0(projname, '.sp.riv' ) )
  fn.mdseg = file.path(inpath, paste0(projname, '.sp.rivseg' ) )
  
  fn.ic = file.path(inpath, paste0(projname, '.cfg.ic' ) )
  fn.para = file.path(inpath, paste0(projname, '.cfg.para' ) )
  fn.calib = file.path(inpath, paste0(projname, '.cfg.calib' ) )
  
  fn.soil = file.path(inpath, paste0(projname, '.para.soil' ) )
  fn.geol = file.path(inpath, paste0(projname, '.para.geol' ) )
  fn.lc = file.path(inpath, paste0(projname, '.para.lc' ) )
  
  fn.forc = file.path(inpath, paste0(projname, '.tsd.forc' ) )
  fn.bc = file.path(inpath, paste0(projname, '.tsd.bc' ) )
  fn.lai = file.path(inpath, paste0(projname, '.tsd.lai' ) )
  fn.rl = file.path(inpath, paste0(projname, '.tsd.rl' ) )
  fn.mf = file.path(inpath, paste0(projname, '.tsd.mf' ) )
  gisdir = file.path(inpath, 'gis')
  fns = c(
    inpath, outpath, gisdir, 
    fn.mesh, fn.pt, fn.edge, fn.mdseg, fn.mdriv, fn.att,
    fn.ic,fn.para, fn.calib,
    fn.soil, fn.geol, fn.lc,
    fn.forc, fn.bc, fn.lai, fn.rl, fn.mf)
  names(fns) = c(
    'inpath', 'outpath', 'gispath',
    "md.mesh","md.pt", "md.edge", "md.rivseg", "md.riv", "md.att",
                 'md.ic', 'md.para', 'md.calib',
                 'md.soil', 'md.geol', 'md.lc',
                 'md.forc', 'md.bc', 'md.lai', 'md.rl', 'md.mf' )
  fns
}
#' Prepare configure file for autoSHUD
#' \code{read.auto} 
#' @return Character of full path of input files for SHUD model
#' @export
read.auto<-function(fn){
  msg = 'read.auto: '
  message(msg,'Reading rawdata information from file: ', fn)
  if(file.exists(fn)){
    
    tmp=as.matrix(read.table(fn, header = F, row.names = 1))
    cdir = tmp[,1]
    # print(cdir)
  }else{
    message(msg, 'File missing.')
    stop('File missing: ', fn)
  }
  # dir.rawdata=cdir['dir.rawdata']
  # dir.soil = cdir['dir.soil']
  # dir.ldas = cdir['dir.ldas']
  # dir.out = cdir['dir.out']
  # 
  # prjname=cdir['prjname']
  # years=as.numeric(cdir['startyear']): as.numeric(cdir['endyear'])
  # 
  # fsp.wbd = cdir['fsp.wbd']
  # fsp.stm = cdir['fsp.stm']
  # 
  # fr.dem=cdir['fr.dem']
  # fr.landuse = cdir['fr.landuse']
  # NumCells = as.numeric(cdir['NumCells'])
  # if(is.null(NumCells) | is.na(NumCells)){
  #   NumCells = 1000;
  # }
  # AqDepth = as.numeric(cdir['AqDepth'])
  # if(is.null(AqDepth) | is.na(AqDepth)){
  #   AqDepth = 10;
  # }
  # # years=2017:2018
  # dir.png =file.path(dir.out, 'Image')
  # dir.pihmgis = file.path(dir.out, 'PIHMgis' )
  # dir.pihmin <- file.path(dir.out, 'input', prjname)
  # dir.pihmout <- file.path(dir.out, 'output', paste0(prjname, '.out') )
  # # dir.forc <- file.path(dir.out, 'forcing')
  # dir.forc=cdir['dir.forc']
  # 
  # tmp=lapply(list(dir.out,dir.png, dir.pihmgis, dir.pihmin, dir.pihmout, dir.forc), dir.create, showWarnings=F, recursive=T)
  as.list(cdir)
}
