# #' Write a backup file if the file exist
# #' \code{filebackup}
# #' @param file file name
# #' @param backup if backup, TRUE/FALSE
# #' @export
# filebackup <- function(file, backup = FALSE){
#   if(file.exists(file) & backup){
#     bakFile <- paste0(file,format(Sys.time(), '%Y%m%d.%H%M%S'),'-',round(rnorm(1),2));
#     file.copy(file,bakFile)
#   }
# }
#' Write xts data out into file
#' \code{write.xts}
#' @param x xts data
#' @param file file name
#' @param append whether append
#' @export
write.xts <- function(x, file, append = F){
    msg='write.ts::'
    tt = stats::time(x)
    md=data.frame('TIME'=tt, x) 
    write.table(md, file = file, append = append,
                sep = '\t', col.names = TRUE, row.names = FALSE) 
    
}
#' Write xts data out into file
#' \code{write.tsd}
#' @param x xts data
#' @param file file name
#' @param dt_units Time interval of the 1 unit change in time column.
#' @param append whether append
#' @param header Header of the file. Default is the dimension of data.frame.
#' @param quite TRUE/FALSE, if quiet mode
#' @export
write.tsd <- function(x, file, dt_units='days', append = F, quite=F, header = NULL){
    msg='write.tsd::'
  
   # x=lr$LAI
  mat = as.matrix(rbind(x))
  nr = nrow(x)
  nc = ncol(x)
  if(!quite){
    message(msg, 'Writing to file ', file)
  }
  tt = stats::time(x)
  if(grepl('minute', dt_units)){
    dt_sec = 60 #  60 sec
  }else if(grepl('second', dt_units)){
    dt_sec = 1
  }else{
    dt_sec = 3600*24
  }
  # time_tag = as.numeric( difftime(tt, tt[1], units = 'days') )
  time_tag = as.numeric(tt - tt[1]) / dt_sec
  if(is.null(header)){
    t0 = format(time(x)[1], '%Y%m%d')
    t1 = format(time(x)[nrow(x)], '%Y%m%d')
    header = c(nr, nc+1, t0, t1, dt_sec)
  }
  dd = data.frame('Time_interval' = time_tag, mat)
  write(header,file = file, ncolumns = length(header), append = append, sep = '\t')
  write(colnames(dd), file = file, ncolumns = nc+1, append = TRUE, sep = '\t')
  suppressWarnings(write.table(dd, file = file, append = TRUE, sep = '\t', col.names = FALSE, row.names = FALSE) )
}
#' Write data.frame out into file
#' \code{write.df}
#' @param x data.frame
#' @param file file name
#' @param append whether append
#' @param header Header of the file. Default is the dimension of data.frame.
#' @param quite TRUE/FALSE, if quiet mode
#' 
#' @export
write.df <- function(x, file, append = F, quite=F, header = NULL){
    msg='write.df::'
  
  x=as.matrix(rbind(x))
  nr = nrow(x)
  nc = ncol(x)
  if(is.null(header)){
    header = c(nr,nc)
  }
  if(!quite){
    message(msg, 'Writing to file ', file)
  }
  write(header, file = file, append = append, sep = '\t',  ncolumns = length(header))
  write(colnames(x), file = file, ncolumns = nc,append = T, sep = '\t')
  write(t(x), file = file, ncolumns = nc, append = T, sep = '\t')
}
#' Write  .mesh file
#' \code{write.mesh}
#' @param pm SHUD.MESH class
#' @param file file name
#' 
#' @export
write.mesh <- function(pm, file){
    msg='writemesh::'
  
  ncell= nrow(pm@mesh)
  np = nrow(pm@point)
  message(msg, 'Writing to file ', file)
  write.df(pm@mesh, file=file, append=F, quite = T)
  write.df(pm@point, file=file, append=T, quite = T)
  # write(c(ncell, np),file = file, append = F)
  # write.table(pm@mesh, file = file, append = T, quote = F, row.names = F)
  # write.table(pm@point, file = file, append = T, quote = F, row.names = F)
}
#' Write SHUD .mesh file
#' \code{writemesh}
#' @param pr SHUDmesh class
#' @param file file name
#' 
#' @export
write.riv <- function(pr, file){
    msg='writeriv::'
  
  message(msg, 'Writing to file ', file)
  write.df(pr@river, file=file, append=F, quite = T)
  write.df(pr@rivertype, file=file, append=T, quite = T)
  if(length(pr@point) >0 ){
    write.df(pr@point, file=file, append=T, quite = T)
  }
}
#' Write SHUD .ic file
#' \code{writeinit}
#' @param x Initial condition, list()
#' @param file file name
#' 
#' @export
write.ic <- function(x, file){
    msg='writeinit::'
  
  message(msg, 'Writing to file', file)
  write.df(x[[1]], file=file, append=F, quite = T)
  write.df(x[[2]], file=file, append=T, quite = T)
}
#' Write SHUD configuration files (.para, .calib, etc.)
#' @param x SHUD model configure parameter or calibration
#' @param file file name
#' @importFrom utils type.convert write.table
#' @export
write.config <-function(x, file){
    msg='write.config::'
  
  message(msg, 'Writing to file ', file)
  out=cbind(names(x), t(x))
  write.table(out, file, append = F, sep = '\t', quote = F, 
              row.names = FALSE, col.names = FALSE)
}

#' Write SHUD .forc file
#' \code{write.forc}
#' @param x  data.frame of the forcing sites information. c(id, lon, lat, x, y, z, filename)
#' @param path Common path of the files.
#' @param startdate Start Date. Character. e.g. 20000101
#' @param file file name
#' 
#' @export
write.forc <- function(x, file, path='', startdate='20000101'){
    msg='writeforc::'
  
  nf=nrow(x)
  nc=ncol(x)
  message(msg, 'Writing to file ', file)
  write( paste(nf, startdate), file=file, append=F)
  write( path, file=file, append=T,  ncolumns = 1)
  write(colnames(x), file = file, ncolumns = nc,append = T, sep = '\t')
  write(t(x), file = file, ncolumns = nc, append = T, sep = '\t')
  # write( fns, file=file, append=T,  ncolumns = 1)
}


