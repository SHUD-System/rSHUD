# Suppress R CMD check notes for ggplot2 NSE in plot functions
utils::globalVariables(c("Time", "Value", "Period", "Precipitation", "value", "variable"))

#============ 
#' Setup .shud environment.
#' \code{shud.env} 
#' @param prjname Project name
#' @param inpath path of SHUD model input directory
#' @param outpath path of SHUD model output directory
#' @param prjname character of SHUD model project name.
#' @param anapath character of SHUD model project name.
#' @param version version of SHUD model.
#' @export
shud.env <- function(prjname, inpath = file.path('input', prjname), 
                     outpath = file.path('output', paste0(prjname, '.out')), 
                     anapath = file.path(outpath, 'postprocessing'),
                     version = 2.0){
  prjname = as.character(prjname)
  inpath = as.character(inpath)
  outpath = as.character(outpath)
  version = as.numeric(version)
  assign('PRJNAME', prjname, envir = .shud)
  assign('inpath', inpath, envir = .shud)
  assign('outpath', outpath, envir = .shud)
  assign('version', version, envir = .shud)
  
  dir.create(anapath, showWarnings = FALSE, recursive = TRUE)
  assign('anapath', anapath, envir = .shud)
  assign('MASK', NULL, envir = .shud)
  
  # Use custom environment
  .shud$PRJNAME <- prjname
  .shud$inpath <- inpath
  .shud$outpath <- outpath
  .shud$anapath <- anapath
  .shud$version <- version
  
  return(list(prjname = prjname, inpath = inpath, 
              outpath = outpath, anapath = anapath, version = version))
}

