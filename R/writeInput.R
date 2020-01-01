#' Write a backup file if the file exist
#' \code{filebackup}
#' @param file file name
#' @param backup if backup, TRUE/FALSE
#' @export
filebackup <- function(file, backup = TRUE){
  if(file.exists(file) & backup){
    bakFile <- paste0(file,format(Sys.time(), '%Y%m%d.%H%M%S'),'-',round(rnorm(1),2));
    file.copy(file,bakFile)
  }
}
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
#' @param append whether append
#' @param header Header of the file. Default is the dimension of data.frame.
#' @param quite TRUE/FALSE, if quiet mode
#' @param backup TRUE/FALSE, if backup the existing file
#' @export
write.tsd <- function(x, file, append = F, quite=F, header = NULL, backup=TRUE){
    msg='write.tsd::'
  filebackup(file, backup = backup)
   # x=lr$LAI
  mat = as.matrix(rbind(x))
  nr = nrow(x)
  nc = ncol(x)
  if(!quite){
    message(msg, 'Writing to file ', file)
  }
  tt = stats::time(x)
  tday = as.numeric( difftime(tt, tt[1], units = 'days') )
  if(is.null(header)){
    t0 = format(time(x)[1], '%Y%m%d')
    header = c(nr, nc+1, t0)
  }
  dd = data.frame('Time_Day' = tday, mat)
  write(header,file = file, ncolumns = length(header), append = append, sep = '\t')
  write(colnames(dd), file = file, ncolumns = nc+1, append = TRUE, sep = '\t')
  suppressWarnings(
    write.table(dd, file = file, append = TRUE, sep = '\t', col.names = FALSE, row.names = FALSE) 
    )
}
#' Write data.frame out into file
#' \code{write.df}
#' @param x data.frame
#' @param file file name
#' @param append whether append
#' @param header Header of the file. Default is the dimension of data.frame.
#' @param quite TRUE/FALSE, if quiet mode
#' @param backup TRUE/FALSE, if backup the existing file
#' @export
write.df <- function(x, file, append = F, quite=F, header = NULL, backup=TRUE){
    msg='write.df::'
  filebackup(file, backup = backup)
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
#' @param backup TRUE/FALSE, if backup the existing file
#' @export
write.mesh <- function(pm, file,backup = TRUE){
    msg='writemesh::'
  filebackup(file, backup = backup)
  ncell= nrow(pm@mesh)
  np = nrow(pm@point)
  message(msg, 'Writing to file ', file)
  write.df(pm@mesh, file=file, append=F, quite = T, backup = F)
  write.df(pm@point, file=file, append=T, quite = T, backup = F)
  # write(c(ncell, np),file = file, append = F)
  # write.table(pm@mesh, file = file, append = T, quote = F, row.names = F)
  # write.table(pm@point, file = file, append = T, quote = F, row.names = F)
}
#' Write SHUD .mesh file
#' \code{writemesh}
#' @param pr SHUDmesh class
#' @param file file name
#' @param backup TRUE/FALSE, if backup the existing file
#' @export
write.riv <- function(pr, file,backup = TRUE){
    msg='writeriv::'
  filebackup(file, backup = backup)
  message(msg, 'Writing to file ', file)
  write.df(pr@river, file=file, append=F, quite = T, backup = F)
  write.df(pr@rivertype, file=file, append=T, quite = T, backup = F)
  if(length(pr@point) >0 ){
    write.df(pr@point, file=file, append=T, quite = T, backup = F)
  }
}
#' Write SHUD .ic file
#' \code{writeinit}
#' @param x Initial condition, list()
#' @param file file name
#' @param backup TRUE/FALSE, if backup the existing file
#' @export
write.ic <- function(x, file, backup = TRUE){
    msg='writeinit::'
  filebackup(file, backup = backup)
  message(msg, 'Writing to file', file)
  write.df(x[[1]], file=file, append=F, quite = T, backup = F)
  write.df(x[[2]], file=file, append=T, quite = T, backup = F)
}
#' Write SHUD configuration files (.para, .calib, etc.)
#' @param x SHUD model configure parameter or calibration
#' @param file file name
#' @param backup TRUE/FALSE
#' @importFrom utils type.convert write.table
#' @export
write.config <-function(x, file, backup=TRUE){
    msg='write.config::'
  filebackup(file, backup = backup)
  message(msg, 'Writing to file ', file)
  out=cbind(names(x), t(x))
  write.table(out, file, append = F, sep = '\t', quote = F, 
              row.names = FALSE, col.names = FALSE)
}

#' Write SHUD .forc file
#' \code{write.forc}
#' @param fns  filenames; Vector of character,
#' @param path Common path of the files.
#' @param startdate Start Date. Character. e.g. 20000101
#' @param file file name
#' @param backup TRUE/FALSE, if backup the existing file
#' @export
write.forc <- function(fns, path='', startdate='20000101', file, backup=TRUE){
    msg='writeforc::'
  filebackup(file, backup = backup)
  nf=length(fns)
  message(msg, 'Writing to file ', file)
  write( paste(nf, startdate), file=file, append=F)
  write( path, file=file, append=T,  ncolumns = 1)
  write( fns, file=file, append=T,  ncolumns = 1)
}
