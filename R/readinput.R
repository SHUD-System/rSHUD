#' Read the matrix or data.frame file which has (nrow, ncol) in header
#' \code{read.df} 
#' @param file full path of file
#' @param text text return from readLines(file)
#' @param sep seperator of the table.
#' @return a list of matrix
#' @noRd
read.df_legacy <-function(file, text = readLines(file), sep='\t'){
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
    # print(text[0:nr + 1 + r0])
    xl[[i]] = rbind(utils::read.table(text = text[0:nr + 1 + r0], header = T, sep = sep))
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
#' @noRd
readmesh_legacy <-function(file = shud.filein()['md.mesh'] ){
  d = read.df_legacy(file = file);
  ret <- SHUD.MESH(mesh = data.frame(d[[1]]),
                   point = data.frame(d[[2]]) )
}

#' Read the .att file
#' \code{readmesh} 
#' @param file full path of file
#' @return matrix
#' @noRd
readatt_legacy <-function( file = shud.filein()['md.att']){
  x = read.df_legacy(file)
  ret = x[[1]]
}

#' Read the .riv file
#' \code{readriv} 
#' @param file full path of .riv file
#' @return SHUD.RIVER
#' @noRd
readriv_legacy <-function(file = shud.filein()['md.riv'] ){
  d = read.df_legacy(file = file);
  ret <- SHUD.RIVER(river = data.frame(d[[1]]), 
                    rivertype = data.frame(d[[2]]),
                    point = data.frame() #data.frame(d[[3]]) 
                    )
}
#' Read the river shapefile
#' \code{readriv.sp} 
#' @param file full path of spatialLines file
#' @return spatialLines file
#' @noRd
readriv.sp_legacy <-function(file = file.path(shud.filein()['inpath'], 'gis', 'river.shp' ) ){
  # Check if file exists
  if(!file.exists(file)){
    stop('River shapefile does not exist: ', file)
  }
  
  spr = as(sf::st_read(file, quiet = TRUE), "Spatial")
}
#' Read the .rivseg file
#' \code{readrivseg} 
#' @param file full path of .rivseg file
#' @return matrix
#' @noRd
readrivseg_legacy <-function(file = shud.filein()['md.rivseg'] ){
  d = read.df_legacy(file = file);
  ret <- d[[1]]
}

#============
#' Read the .para file
#' \code{readpara} 
#' @param file full path of file
#' @return .para
#' @noRd
readpara_legacy <- function(file = shud.filein()['md.para'] ){
  return(readconfig_legacy(file))
}

#============
#' Read the .calib file
#' \code{readcalib} 
#' @param file full path of file
#' @return .calib
#' @noRd
readcalib_legacy <- function(file = shud.filein()['md.calib'] ){
  return(readconfig_legacy(file))
}

#============
#' Read the configuration (.para, .calib) file
#' \code{readconfig} 
#' @param file full path of file
#' @return .para or .calib
#' @importFrom utils type.convert write.table
#' @noRd
readconfig_legacy <- function(file = shud.filein()['md.para']){
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
#' @noRd
readic_legacy <- function(
  file = shud.filein()['md.ic']){
  # file = fin['md.ic']
  x=read.df_legacy(file)
  cn = c('minit', 'rinit', 'linit')
  names(x) = cn[1:length(x)]
  x
}

#' Read the .soil file
#' \code{readsoil} 
#' @param file full path of file
#' @return .soil
#' @noRd
readsoil_legacy <-function(file =  shud.filein()['md.soil']){
  x=read.df_legacy(file)
  ret = x[[1]]
}

#============
#' Read the .geol file
#' \code{readgeol} 
#' @param file full path of file
#' @return .geol
#' @noRd
readgeol_legacy <-function(file =  shud.filein()['md.geol']){
  x=read.df_legacy(file)
  ret = x[[1]]
}

#============ 
#' Read the .lc file
#' \code{readlc} 
#' @param file full path of file
#' @return .lc
#' @noRd
readlc_legacy <-function( file = shud.filein()['md.lc']){
  x=read.df_legacy(file)
  ret = x[[1]]
}
#============ 
#' Read the .forc file
#' \code{readlc} 
#' @param file full path of file
#' @return data.frame of forcing sites.
#' @noRd
readforc.fn_legacy <-function( file = shud.filein()['md.forc']){
  txt = readLines(file)
  hd=read.table(text = txt[1])
  path=txt[2]
  df=read.table(text=txt[-1 * 1:2 ], header = TRUE)
  nc=ncol(df)
  df[, nc] = file.path(path, df[, nc])
  ret = list('StartTime' = as.character(hd[2]) ,
             Sites = df)
  return(ret)
}
#============ 
#' Read the .csv files in .forc file
#' \code{readforc.csv} 
#' @param file full path of file
#' @param id  Index of the forcing sites.  default = NULL, which return average values of all sites.
#' @return forcing data, list.
#' @noRd
readforc.csv_legacy <-function(file = shud.filein()['md.forc'], id=NULL){
  msg='readforc.csv::'
  xf = readforc.fn_legacy(file=file)
  fns = xf$Sites[, ncol(xf$Sites)]
  tstr = xf$StartTime
  t0 = as.POSIXct(tstr, format = '%Y%m%d')
  if( is.null(id)){
    fn = fns
    RID=1:length(fns)
  }else{
    fn = fns[id]
    RID = id
  }
  nf = length(fn)
  xl = list()
  for(i in 1:nf){
    message(msg, i, '/', nf,'\t', fn[i])
    # x=read.df(fns[RID[i]]) 
    # y=x[[1]]
    # xt = t0+y[,1]*86400
    # tsd=zoo::zoo(y[,-1], xt)
    tsd = read_tsd(fns[RID[i]])[[1]]
    xl[[i]] = tsd
  }
  names(xl) = basename(fns[RID])
  
  if( is.null(id) ){
    mx = xl[[1]] * 0
    for(i in 1:nf){
      mx = mx + xl[[i]]
    }
    return(mx/nf)
  }else{
    return(xl)
  }
}
