#============ 
#' Setup .shud environment.
#' \code{shud.env} 
#' @param inpath path of SHUD model input directory
#' @param outpath path of SHUD model output directory
#' @param prjname charactor of SHUD model project name.
#' @param anapath charactor of SHUD model project name.
#' @export
shud.env <- function(prjname, inpath = file.path('input', prjname), 
                 outpath = file.path('output', paste0(prjname, '.out')), 
                 anapath = file.path(outpath, 'SHUDtb')){
  assign('PRJNAME', prjname, envir=.shud)
  assign('inpath', inpath, envir=.shud)
  assign('outpath', outpath, envir=.shud)
  
  dir.create(anapath, showWarnings = F, recursive = T)
  assign('anapath', anapath, envir = .shud)
  assign('MASK', NULL, envir=.shud)
  return(list(prjname=prjname, inpath=inpath, outpath=outpath, anapath=anapath))
}
