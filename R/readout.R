#' Read output file from SHUD model
#' @param keyword keyword of the output file. File name format is: projectname.keyword.dat
#' @param file  Full path of output file. 
#' @param path path of the outputfile
#' @param ASCII If TRUE, convert the binary file to ASCII format (.csv)
#' @keywords read output. Could be used for reading mesh and river format.
#' @importFrom grDevices dev.off graphics.off png rgb topo.colors 
#' @importFrom graphics grid hist lines par plot points 
#' @importFrom methods as 
#' @importFrom stats dist rnorm time 
#' @importFrom utils read.table 
#' @return A TimeSeries data. This list require support of xts packages.
#' @export  
readout <- function(keyword,
                    path=get('outpath', envir=.shud) , ASCII=FALSE,
                    file = file.path(path, paste0(get('PRJNAME', envir=.shud),'.', keyword,'.dat') ) 
){
  msg='readout::'
  fid=file(file, 'rb');
  # nc=readBin(fid, what=integer(), n=1, size = 4)
  # st=readBin(fid, what=integer(), n=1, size = 8) #Long integer
  tmp=readBin(fid, what=numeric(), n=1e9, size=8)
  dat=tmp[-1 * (1:2)]
  nc=tmp[1]
  st=tmp[2]
  close(fid)
  
  nrr = length(dat)/(nc+1)
  nr = round(nrr)
  if(nr < nrr){
    message(msg, 'File may not completed. ', nrr, "X", nc+1)
  }
  mat=t(matrix(dat[1:( nr*(nc+1) )], nrow=nc+1))
  if(ASCII){
    fn.asc=paste0(substr(file, 1, nchar(file)-3), 'csv')
    write(paste(nc, st), file = fn.asc, append = FALSE)
    tmp=mat
    colnames(tmp)=c('T_MIN', paste0('X', 1:nc))
    suppressWarnings(
      write.table(tmp, file = fn.asc, append = TRUE, col.names = TRUE, row.names = FALSE)
    )
  }
  tmove = diff(mat[,1])
  tmove = c(tmove, tmove[length(tmove)])
  tsec =   ( mat[,1]) * 60 
  xt = as.POSIXct(as.character(st), format='%Y%m%d') + tsec
  if(nc<=1){
    tsd = xts::as.xts(cbind(mat[,-1]), order.by = xt)
  }else{
    tsd = xts::as.xts(rbind(mat[,-1]), order.by = xt)
  }
  tsd
}

#' Read multiply SHUD model output files and do time-series plot
#' @param varname vector of output keywords 
#' @param plot Whether do the time-series plot
#' @param imap Whether do the raster plot for Element data. Only works for element data
#' @param return Whether return the data. Some the results are too huge to load in memoery at once.
#' @param iRDS Whether save RDS file.
#' @param sp.riv River SpatialLine*
#' @param rdsfile Save RDS file
#' @param w.focal forcal matrix
#' @keywords read output.
#' @return A list of TimeSeries data. 
#' @export  
BasicPlot <- function(
  varname=c(paste0('eley',c( 'surf','unsat', 'gw', 'snow') ), 
            paste0('elev',c('prcp','infil', 'rech') ),
            paste0('elev',c('etp', 'eta', 'etev', 'ettr', 'etic') ),
            paste0('rivq',c('down', 'sub', 'surf')),
            paste0('rivy','stage')
  ) ,
  sp.riv=NULL,
  rdsfile = file.path(get('outpath', envir = .shud), 'BasicPlot.RDS'),
  plot=TRUE, imap=FALSE, iRDS = TRUE,
  return=T, w.focal=matrix(1/9, 3, 3)){
  msg='BasicPlot::'
  graphics.off()
  varname = tolower(varname)
  print(varname)
  nv=length(varname)
  prjname = get('PRJNAME', envir = .shud)
  if(plot || imap){
    path = file.path(get('anapath', envir = .shud), 'BasicPlot')
    if(!dir.exists(path)){
      dir.create(path, recursive = T, showWarnings = F)
    }
  }
  
  vtype = substr(varname,4,4)
  att=cbind(varname, vtype, unit='')
  att[vtype %in% 'y', 3] = 'm'
  att[vtype %in% 'v', 3] = 'm/d'
  att[vtype %in% 'q', 3] = 'm3/d'
  
  ret=list()
  for(i in 1:nv){
    vn=varname[i]
    fn= paste0(prjname,'_', vn, '.png')
    message(msg, i, '/', nv, '\t', vn, '\t', fn)
    x=readout(vn);
    xm = xts::as.xts(apply(x, 1, mean), order.by=time(x))
    if(return){
      ret[[i]] = x;
    }
    if(plot){
      time =time(x)
      t0=min(time)
      t1=max(time)
      tr=paste( strftime(t0) ,'-', strftime(t1))
      png(file.path(path, fn), width=11, height=9, res=100, units = 'in')
      pp=xts::plot.xts(x, main=paste0(att[i, 1], ' (', att[i,3], ')') )
      pp = lines(xm, col=2, lwd=2, lty=2, type='b', pch=1)
      print(pp)
      dev.off()
    }
    if(imap ){
      # if( grepl('^eley', vn) | grepl('^elev', vn)  | grepl('^eleq', vn) ){
      if( grepl('^ele', vn)){
        fn= paste0('Map.', prjname,'_', vn, '.png')
        y = colMeans(x)
        png(file.path(path, fn), width=11, height=9, res=100, units = 'in')
        # map2d(y)
        r = MeshData2Raster(y, stack=FALSE)
        if(!is.null(w.focal)){
          r = raster::focal(r, w=w.focal)
        }
        raster::plot(r)
        if(!is.null(sp.riv)){
          raster::plot(sp.riv, col=rgb(1,0,0,0.7), lwd=2, add=T)
        }
        dev.off()
      }
    }
  }
  if(return){
    names(ret) = varname;
    if(iRDS){
    saveRDS(ret, file = rdsfile)
    }
  }
  ret <- ret
}
