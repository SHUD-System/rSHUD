#============ 
#' Read the time-series file, which has c(nrow, ncol, timestring) in header
#' \code{read.tsd} 
#' @param file full path of file
#' @return xts
#' @export
read.tsd <- function(file){
  # file=fin['md.lai']
  text=readLines(file)
  r0 = 1
  nrow = length(text)
  xl = list()
  for(i in 1: 100){
    header = utils::read.table(text = text[r0] )
    nr = as.numeric(header[1])
    nc = as.numeric(header[2])
    t0 = header[3]
    tmp =  as.matrix(utils::read.table(text = text[0:nr + 1 + r0], header = T))
    tsd = xts::as.xts(tmp[,-1], 
                      order.by = as.POSIXct(paste(t0), format="%Y%m%d") + tmp[,1] * 86400)
    xl[[i]] =tsd
    r0 = r0 + nr + 2;
    if(r0 + 1 > nrow){
      break
    }
  }
  xl
}

#============ 
#' Read the .ts.lai Time-Series file
#' \code{readlai} 
#' @param file full path of file
#' @return list of Time-Series data
#' @export
readlai <- function( file = shud.filein()['ts.lai'] ){
  x=read.tsd(file)
  names(x) = c('LAI','RL')
  x
}

