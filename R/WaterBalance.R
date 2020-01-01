#' Calculate the water balance
#' \code{wb.all}
#' @param xl List of data. Five variables are included: prcp (eleqprcp), etic (eleqetic), ettr (eleqettr), etev (eleqetev) and discharge (rivqflx)
#' @param fun function to process the time-series data. Default = apply.daily.
#' @param ic Initial Condition.
#' @param plot Whether plot the result
#' @return A matrix, contains the colums of water balance factors
#' @export
wb.all <-function(
  xl=BasicPlot(varname=c(paste0('elev', c('prcp', 'etic', 'ettr', 'etev', 'etp') )
                         , 'rivqdown', 
                         paste0('eley', c('surf', 'unsat', 'gw') ),
                         paste0('rivy', 'stage') ),plot = FALSE, return = TRUE), 
               ic=readic(),
  fun = xts::apply.monthly, plot=TRUE
){
  func <-function(x, w){
    aa= sum(ia)
    y = sweep(x, 2, w, '*')
  }
  ia=getArea()
  aa=sum(ia)
  w = ia/aa
  pr = readriv()
  oid=getOutlets(pr)
  P = fun( func(xl$elevprcp, w ), FUN=sum)
  Q = fun(xl$rivqdown[,oid], FUN=sum) / aa
  IC = fun( func(xl$elevetic, w ), FUN=sum)
  ET = fun( func(xl$elevettr, w ), FUN=sum)
  EV = fun( func(xl$elevetev, w ), FUN=sum)
  PET = fun( func(xl$elevetp, w ), FUN=sum)
  AET=IC+EV+ET
  ds = wb.DS(xl=xl, ic=ic)
  dh = P-Q-IC-EV-ET
  x=cbind(dh, P, Q, AET, PET,IC, ET,EV)
  # colnames(x)=c('P', 'PET','Q','ET_IC', 'ET_TR','ET_EV')
  # y = cbind(dh, x)
  colnames(x)=c('DH', 'P', 'Q','AET','PET','ET_IC', 'ET_TR','ET_EV')
  if(plot){
      hydrograph(x)
  }
  dse = sum( apply(ds$Ele, 2, sum) * w)
  dsr = sum(ds$Riv)
  dss = dse + dsr
  y=c(apply(x, 2, sum, na.rm=TRUE), DS=dss, DS.Ele=dse, DS.Riv = dsr)
  mat = rbind('H'=y, '%'=y/y[2]*100)
  print(mat)
  return(x)
}

#' Convert Time-Series data to data.frame
#' \code{ts2df}
#' @param x Time-Series data.
#' @return data.frame. First column is Time.
#' @export
ts2df <- function(x){
  d1 = as.data.frame(x)
  colnames(d1) = colnames(x)
  d2 = data.frame('Time'=time(x), d1)
}

#' Calculate the water balance
#' \code{wb.riv}
#' @param xl List of data. Five variables are included: prcp (eleqprcp),and discharge ()
#' @param fun function to process the time-series data. Default = apply.daily.
#' @param plot Whether plot the result
#' @return A matrix, contains the colums of water balance factors
#' @export
wb.riv <-function(
  xl=BasicPlot(varname=paste0('rivq', c('sub', 'surf', 'down') ),
               plot = FALSE, return = TRUE),
  fun = xts::apply.yearly, plot=TRUE
){
  fun.read <- function(xx, val){
    cn=names(xx)
    if( val %in% cn ){
      return(xx[[val]])
    }else{
      return(readout(val))
    }
  }
  # xl=BasicPlot(varname=paste0('rivq', c('sub', 'surf', 'flx') ) )
  # fun = xts::apply.daily

  func <-function(x, w){
    aa= sum(ia)
    y = sweep(x, 2, w, '*')
  }
  ia=getArea()
  aa=sum(ia)
  pr = readriv()
  oid=getOutlets(pr)
  qr = fun.read(xl, 'rivqdown')[,oid]
  qsf = fun.read(xl, 'rivqsurf')[,]
  qsb = fun.read(xl, 'rivqsub')[,]

  Qo = fun(qr, FUN=sum)
  Qsf = fun(qsf, FUN=sum)
  Qsub = fun(qsb, FUN=sum)
  q3=cbind(Qo, -Qsf, -Qsub)  / aa
  # plot(q3)

  dh = -( Qo + Qsf + Qsub ) /aa
  # plot(dh)

  x=cbind(dh, q3)
  colnames(x)=c('DH', 'Qout','Qin_sf','Qin_gw')
  if(plot){
    hydrograph(x)
  }
  x
}

#' Calculate the water balance
#' \code{wb.ele}
#' @param xl List of data. Five variables are included: 
#' @param fun function to process the time-series data. Default = apply.daily.
#' @param period Period of the waterbalance
#' @param plot Whether plot the result
#' @return A matrix, contains the colums of water balance factors
#' @export
wb.ele <-function(
  xl=BasicPlot(varname=c(paste0('elev', c('prcp', 'etic', 'ettr', 'etev') ),
                         paste0('eleq', c('surf', 'sub') ),
                         paste0('eley', c('surf', 'unsat', 'gw')  ) ) ),
  fun = xts::apply.yearly, period = 'years', plot=TRUE ){
  # xl=BasicPlot(varname=c(paste0('elev', c('prcp', 'etic', 'ettr', 'etev')),
  #                        paste0('eleq', c('surf', 'sub') ),
  #                        paste0('eley', c('surf', 'unsat', 'gw') )  ) )
  # fun = xts::apply.daily
  func <-function(x, w){
    aa= sum(ia)
    y = sweep(x, 2, w, '*')
  }
  fun.first<- function(x, period='month'){
    r <- do.call(rbind, lapply(split(x, period), xts::first))
  }
  fun.last<- function(x, period='month'){
    r <- do.call(rbind, lapply(split(x, period), xts::last))
  }
  ia=getArea()
  aa=sum(ia)
  w=1/ia
  fun.read <- function(xx, val){
    cn=names(xx)
    if( val %in% cn ){
      return(xx[[val]])
    }else{
      return(readout(val))
    }
  }
  P = fun(fun.read(xl, 'elevprcp'), FUN=mean)
  IC = fun( fun.read(xl, 'elevetic'), FUN=mean)
  ET = fun( fun.read(xl, 'elevettr'), FUN=mean)
  EV = fun( fun.read(xl, 'elevetev'), FUN=mean)
  QS = fun( fun.read(xl, 'eleqsurf'), FUN=mean)
  Qg = fun( fun.read(xl, 'eleqsub'), FUN=mean)
  dys = fun(xts::diff.xts(xl$eleysurf), FUN=colSums, na.rm=T)
  dyg = fun(xts::diff.xts(xl$eleygw), FUN=colSums, na.rm=T)
  dyu = fun(xts::diff.xts(xl$eleyunsat), FUN=colSums, na.rm=T)
  arr=abind::abind(P, IC, ET, EV, QS/aa, Qg/aa, dys, dyu, dyg, along=3)
  dimnames(arr)[[3]] = c('P','ET_IC', 'ET_TR','ET_EV', 'qs', 'qg','dys','dyu', 'dyg')
  arr
}

#' Calculate the Change of Storage.
#' \code{DeltaS}
#' @param x Time-Serres Matrix
#' @param x0 Intial condition or the values of first time-step
#' @param t1 Time-step one, default = 1
#' @param t2 Time-step two, default = nrow(x)
#' @return A vector
#' @export
DeltaS<-function(x, x0=x[t1, ], t1=1, t2=nrow(x)){
  ds = (as.numeric(x[t2,]) - as.numeric(x0))
  return(ds)
}

#' Calculate the Change of Storage.
#' \code{wb.DS}
#' @param xl List of data. Five variables are included: surface (eleysurf), unsat zone (eleyunsat), groundwater (eleygw) and river stage (rivystage)
#' @param ic Initial Condition.
#' @return A list. list(Ele, Riv)
#' @export
wb.DS<-function(xl=BasicPlot(varname = c(paste0('eley', c('surf', 'unsat', 'gw')),
                                         paste0('rivy', 'stage')),plot = FALSE, return = TRUE),
                ic=readic() ){
  ic=readic()
  g=readgeol()
  att=readatt()
  cfg.calib=readcalib()
  pr=readriv()
  rtype = pr@river$Type
  ra = pr@rivertype$Width[rtype] * pr@river$Length
  AA = sum(getArea())
  por=g$ThAETS.m3_m3.[att$GEOL] * cfg.calib$GEOL_THAETS
  ds.sf=DeltaS(xl$eleysurf, x0=ic$minit$Surface)
  ds.us=DeltaS(xl$eleyunsat, x0=ic$minit$Unsat) * por
  ds.gw=DeltaS(xl$eleygw, x0=ic$minit$GW) * por
  ds.riv=DeltaS(xl$rivystage, x0=ic$rinit$Stage) * ra / AA
  r = list(Ele = rbind(ds.sf, ds.us, ds.gw),
           Riv = rbind(ds.riv)
  )
  return(r)
}
