#' Read the matrix or data.frame file which has (nrow, ncol) in header
#' \code{read.df} 
#' @param file full path of file
#' @param text text return from readLines(file)
#' @return a list of matrix
#' @export
read.df <-function(file, text = readLines(file)){
  # fn=fin['md.mesh']
  # text=readLines(fn)
  r0 = 1
  nrow = length(text)
  xl = list()
  for(i in 1: 100){
    ndim = as.numeric(utils::read.table(text = text[r0] ))
    nr = ndim[1]
    nc = ndim[2]
    if(nr <= 0){
      nr = nrow-2
    }else{
    }
    xl[[i]] = utils::read.table(text = text[0:nr + 1 + r0], header = T)
    # xl[[i]] = as.matrix(utils::read.table(text = text[0:nr + 1 + r0], header = T))
    r0 = r0 + nr + 2;
    if(r0 + 1 > nrow){
      break
    }
  }
  xl
}
#' Read the .mesh file
#' \code{readmesh} 
#' @param file full path of .mesh file
#' @return Mesh domain
#' @export
readmesh <-function(file = shud.filein()['md.mesh'] ){
  d = read.df(file = file);
  ret <- SHUD.MESH(mesh = data.frame(d[[1]]),
                   point = data.frame(d[[2]]) )
}

#' Read the .att file
#' \code{readmesh} 
#' @param file full path of file
#' @return matrix
#' @export
readatt <-function( file = shud.filein()['md.att']){
  x = read.df(file)
  ret = x[[1]]
}

#' Read the .riv file
#' \code{readriv} 
#' @param file full path of .riv file
#' @return SHUD.RIVER
#' @export
readriv <-function(file = shud.filein()['md.riv'] ){
  d = read.df(file = file);
  ret <- SHUD.RIVER(river = data.frame(d[[1]]), 
                    rivertype = data.frame(d[[2]]),
                    point = data.frame() #data.frame(d[[3]]) 
                    )
}
#' Read the river shapefile
#' \code{readriv.sp} 
#' @param file full path of spatialLines file
#' @return spatialLines file
#' @export
readriv.sp <-function(file = file.path(shud.filein()['inpath'], 'gis', 'river.shp' ) ){
  spr = rgdal::readOGR(file)  
}
#' Read the .rivchn file
#' \code{readchannel} 
#' @param file full path of .rivchn file
#' @return matrix
#' @export
readchannel <-function(file = shud.filein()['md.rivseg'] ){
  d = read.df(file = file);
  ret <- d[[1]]
}

#============
#' Read the .para file
#' \code{readpara} 
#' @param file full path of file
#' @return .para
#' @export
readpara <- function(file = shud.filein()['md.para'] ){
  return(readconfig(file))
}

#============
#' Read the .calib file
#' \code{readcalib} 
#' @param file full path of file
#' @return .calib
#' @export
readcalib <- function(file = shud.filein()['md.calib'] ){
  return(readconfig(file))
}

#============
#' Read the configuration (.para, .calib) file
#' \code{readconfig} 
#' @param file full path of file
#' @return .para or .calib
#' @importFrom utils type.convert write.table
#' @export
readconfig <- function(file = shud.filein()['md.para']){
  tline = readLines(file, skipNul = TRUE)
  tline=tline[!grepl('^#', tline)]
  tmp = which(grepl('[:Alpha:]', tline) | grepl('[:alpha:]', tline))
  x =  utils::read.table(text = tline, header=FALSE, stringsAsFactors = FALSE)
  xdf = data.frame(rbind(t(x[, -1]), NULL), stringsAsFactors = FALSE)
  colnames(xdf) =toupper(as.character(as.matrix(x[, 1])) )
  # ret = type.convert(xdf)
  ret = xdf;
  return(ret)
}
#============
#' Read the .ic file
#' \code{readic} 
#' @param file full path of file
#' @return .ic
#' @export
readic <- function(
  file = shud.filein()['md.ic']){
  # file = fin['md.ic']
  x=read.df(file)
  cn = c('minit', 'rinit', 'linit')
  names(x) = cn[1:length(x)]
  x
}

#' Read the .soil file
#' \code{readsoil} 
#' @param file full path of file
#' @return .soil
#' @export
readsoil <-function(file =  shud.filein()['md.soil']){
  x=read.df(file)
  ret = x[[1]]
}

#============
#' Read the .geol file
#' \code{readgeol} 
#' @param file full path of file
#' @return .geol
#' @export
readgeol <-function(file =  shud.filein()['md.geol']){
  x=read.df(file)
  ret = x[[1]]
}

#============ 
#' Read the .lc file
#' \code{readlc} 
#' @param file full path of file
#' @return .lc
#' @export
readlc <-function( file = shud.filein()['md.lc']){
  x=read.df(file)
  ret = x[[1]]
}
#============ 
#' Read the .forc file
#' \code{readlc} 
#' @param file full path of file
#' @return files of forcing .csv.
#' @export
readforc.fn <-function( file = shud.filein()['md.forc']){
  text = readLines(file)
  hd=read.table(text = text[1])
  path=text[2]
  fns = text[-1 * 1:2]
  ret = as.matrix(file.path(path, fns), ncol=1)
  colnames(ret) = hd[2]
  ret
}
#============ 
#' Read the .csv files in .forc file
#' \code{readlc} 
#' @param file full path of file
#' @param id  Index of the forcing sites.
#' @return forcing data, list.
#' @export
readforc.csv <-function(file = shud.filein()['md.forc'], id=NULL){
  msg='readforc.csv::'
  fns=readforc.fn(file=file)
  tstr = colnames(fns)
  t0 = as.POSIXct(tstr, format = '%Y%m%d')
  if( is.null(id)){
    fn = fns
    id=1:length(fns)
  }else{
    fn = fns[id]
  }
  nf = length(fn)
  ret = list()
  for(i in 1:nf){
    message(msg, i, '/', nf,'\t', fn[i])
    x=read.df(fns[id[i]]) 
    y=x[[1]]
    xt = t0+y[,1]*86400
    tsd=zoo::zoo(y[,-1], xt)
    ret[[i]] = tsd
  }
  names(ret) = basename(fns)
  ret
}
