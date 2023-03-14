#'  Read the time tage from NC file.
#' @param ncid  the ncid from nc_open()
#' @return a POSIXct date-time objects
#' @importFrom lubridate seconds minutes hours days month year
#' @export
readnc.time <- function(ncid) {
  # modified after https://stackoverflow.com/questions/46001573/convert-a-netcdf-time-variable-to-an-r-date-object
  # require(lubridate)
  ncdims <- names(ncid$dim) #get netcdf dimensions
  timevar <- ncdims[which(ncdims %in% c("time", "Time", "datetime", "Datetime", "date", "Date"))[1]] #find time variable
  times <- ncdf4::ncvar_get(ncid, timevar)
  if (length(timevar)==0) stop("ERROR! Could not identify the correct time variable")
  timeatt <- ncdf4::ncatt_get(ncid, timevar) #get attributes
  timedef <- strsplit(gsub('^ *', '', timeatt$units), " ")[[1]]
  timeunit <- timedef[1]
  if(length(timedef) < 5){
    # cat("Warning:", tz, "not a valid timezone. Assuming UTC\n")
    tz <- "UTC"
  }else{
    tz <- timedef[5]
  }
  timestart <- strsplit(timedef[4], ":")[[1]]
  if (length(timestart) != 3 || timestart[1] > 24 || timestart[2] > 60 || timestart[3] > 60 || any(timestart < 0)) {
    cat("Warning:", timestart, "not a valid start time. Assuming 00:00:00\n")
    warning(paste("Warning:", timestart, "not a valid start time. Assuming 00:00:00\n"))
    timedef[4] <- "00:00:00"
  }
  timestart <- lubridate::ymd_hms(paste(timedef[3], timedef[4]), tz=tz)
  f <- switch(tolower(timeunit), #Find the correct lubridate time function based on the unit
              seconds=lubridate::seconds, second=lubridate::seconds, sec=lubridate::seconds,
              minutes=lubridate::minutes, minute=lubridate::minutes, min=lubridate::minutes,
              hours=lubridate::hours,     hour=lubridate::hours,     h=lubridate::hours,
              days=lubridate::days,       day=lubridate::days,       d=lubridate::days,
              months=lubridate::months,   month=lubridate::months,   m=lubridate::months,
              years=lubridate::years,     year=lubridate::years,     yr=lubridate::years,
              NA
  )
  suppressWarnings(if (is.na(f)) stop("Could not understand the time unit format"))
  timestart + f(times)
}
#'  Read the time tage from NC file.
#' @param ncid  the ncid from nc_open()
#' @param varid The varid.
#' @param ext box coordinates of the subset.
#' @return a POSIXct date-time objects
#' @importFrom graphics  grid lines points
#' @export
#'
readnc<-function(ncid, varid=NULL,  ext = NULL){
  msg= 'readnc:: '
  vars = names(ncid$var)
  nvars = length(vars)
  if(is.null(varid)){ # read all
    varid = varid[!(vars %in% 'time_bnds')] # don't need the time_bnds
  }else if(is.character(varid)){ # read VARID (character) by user
    if(!all(varid %in% vars)){ # validate the input chars.
      message(msg, 'ERROR:: some varid is missing in the dataset.\n')
      print(varid[! (varid %in% vars)])
      stop('Stop with error.')
    }
    varid = varid
  }else if(is.numeric(varid)) {  # read VARID (index) by user
    if(max(varid)>nvars || min(varid) < 1){
      message(msg, 'ERROR:: Wrong value in varid.\n')
      stop('Stop with error.')
    }
    message('Reading VARID = ', vars[varid])
    varid = vars[varid]
  }else{  # ERROR
    message(msg, 'ERROR:: Wrong format of varid.\n')
    stop('Stop with error.')
  }
  
  ncdims = names(ncid$dim)
  var.lon <- ncdims[which(grepl('lon', tolower(ncdims)))]
  var.lat <- ncdims[which(grepl('lat', tolower(ncdims)))]
  
  lon <- ncdf4::ncvar_get(ncid, varid = var.lon)
  lat <- ncdf4::ncvar_get(ncid, varid = var.lat)
  dx = mean(diff(lon)); dy = mean(diff(lat))
  xmin = min(lon - dx/2); xmax = max(lon + dx/2)
  ymin = min(lat - dy/2); ymax = max(lat + dy/2)
  if(is.null(ext)){
    ext= c(min(lon), max(lon), min(lat), max(lat))
  }
  if(ext[1] < xmin | ext[2] > xmax | ext[3] < ymin | ext[4] > ymax){
    warning(paste('Extent required is larger than the boundbox of dataset'))
    message(paste(ext, collaps=','))
    message(paste(c(xmin,xmax,ymin, ymax), collaps=','))
  }
  xid = min(which(abs(lon  - ext[1]) <= dx/2)):max(which(abs(lon  - ext[2]) <= dx/2))
  yid = min(which(abs(lat  - ext[3]) <= dy/2)):max(which(abs(lat  - ext[4]) <= dy/2))
  nx = length(xid); ny = length(yid)
  x.cord = lon[xid]; y.cord = lat[yid]
  # plot.check <- function(){
  #   # plot(0, type='n',
  #   #      xlim= range(c(x.cord, ext[1:2] )),
  #   #      ylim=range(c(y.cord, ext[3:4] )), xlab='Lon', ylab='Lat'); grid()
  #   sp.fn = fishnet(xx=x.cord, yy=y.cord, type='point')
  #   raster::plot(sp.fn, axes=TRUE); grid()
  #   lines(cbind(c(ext[1], ext[1], ext[2], ext[2], ext[1]),
  #               c(ext[3], ext[4], ext[4], ext[3], ext[3])), col=4)
  #   points(x.cord, rep(ext[3], length(xid)), col=2)
  #   points(rep(ext[1], length(yid)), y.cord, col=3)
  # }
  # arr = array(0, dim=c(ny, nx, length(varid)), dimnames= list(y.cord, x.cord, c(varid)))
  arr = array(0, dim=c(nx, ny, length(varid)), dimnames= list(x.cord, y.cord, c(varid)))
  ndims = ncid$ndims
  start = c(min(xid), min(yid), 1)
  count = c(nx, ny, 1)
  for(i in 1:length(varid)){  #reading file
    vn=varid[i]
    arr[, , i]=ncdf4::ncvar_get(ncid, vn, start = start, count = count)
  }
  tx = readnc.time(ncid = ncid)
  rt = list('x' = x.cord, 'y' = y.cord, 'arr' = arr, 'time' = tx)
  return(rt)
}
#'  Read the time tage from (x, y, matrix).
#' @param x  list(x, y, arr) or x coordinates
#' @param y y coordinates
#' @param arr array. dim = (nx, ny)
#' @param Dxy Resolution
#' @param flip Whether flip the matrix.
#' @param plot Whether plot the raster.
#' @return a Raster/RasterStack
#' @export
xyz2Raster <- function(x, y=NULL, arr=NULL,Dxy=NULL,
                       flip=TRUE, plot=TRUE){
  if(is.null(y) | is.null(arr)){
    nc = x;
  }else{
    nc = list(x=x, y=y, arr=arr)
  }
  dims = dim(nc$arr)
  ndims = length(dims)
  x = nc$x; y = nc$y
  
  if(is.null(Dxy) & (length(x) == 1 | length(y) == 1) ){
    message('Dxy cannot be NULL. Dxy MUST be given.')
    stop()
  }
  
  if( ndims > 2){
    # multiple layers
    rl = list()
    for(i in 1:dims[3]){
      rl[[i]] = xyz2Raster(x = x, y = y, Dxy=Dxy,
                           arr=matrix(nc$arr[, , i], nrow=dim(nc$arr)[1], ncol=dim(nc$arr)[2]) )
    }
    rs = raster::stack(rl)
  }else{
    # single layer
    # val = matrix(arr, nrow=nrow(arr), ncol=ncol(arr))
    val = arr
    if(!is.null(Dxy)){
      dx = Dxy[1]
      if(length(Dxy)==1){
        dy = Dxy[1]
      }else{
        dy = Dxy[2]
      }
    }else{
      dx = abs(mean(diff(x)));
      dy = abs(mean(diff(y)))
    }
    nx = length(x);   ny = length(y)
    r = raster::raster(ncols=nx, nrows=ny)
    raster::extent(r) = c(min(x), max(x), min(y), max(y)) + c(-dx, dx, -dy, dy)/2
    # raster::res(r) = c(dx, dy)
    # r = raster::setValues(r, t(val[, ny:1]))
    if(flip){    idx = ny:1
    }else{    idx = 1:ny  }
    rs = raster::setValues(r, t(val[, idx ]) )
  }
  rs
}