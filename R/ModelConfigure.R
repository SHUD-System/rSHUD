#' extract Coordinates of SpatialLines or  SpatialPolygons
#' \code{shud.ic} 
#' @param ncell Number of triangles in SHUD.MESH domain
#' @param nriv Number of river lines
#' @param AqD Default aquifer depth 
#' @param stage Default river stage
#' @param p1 Fraction of unsat zone, [0,1]
#' @param p2 Fraction of Groundwater, [0,1]
#' @return List of initial condition
#' @export
shud.ic <- function(ncell, nriv, AqD = 10, stage = 0.1, p1=0.4, p2 = p1){
  mi = data.frame(cbind(1:ncell, rep(0, ncell), 0, 0,
             AqD * p1, AqD * p2))
  ri = data.frame(cbind(1:nriv, rep(stage, nriv)))
  colnames(mi) = c('Index', 'Canopy', 'Snow', 'Surface', 'Unsat', 'GW')
  colnames(ri) = c('Index', 'Stage')
  ret = list('mic' = mi, 'ric' = ri)
  return(ret)
}
#' Generate the default model configuration
#' \code{shud.para} 
#' @param nday Simulation periods (days)
#' @return model configureation, Vector
#' @export
shud.para <- function( nday = 10){
  dts = c(
    paste0('dt_',
           c(paste0('ye_', c('snow', 'surf', 'unsat', 'gw') ), 
             paste0('Qe_', c('surf', 'sub') ),
             paste0('qe_et'), 
             paste0('qe_', c('prcp', 'infil', 'rech') ),
             paste0('yr_', c('stage')),
             paste0('Qr_', c('down', 'surf', 'sub', 'up'))
           ))
  )
  
  vdt = rep(1440, length(dts))
  
  vn = c('VERBOSE', 'INIT_MODE', 'CloseBoundary', 
         'ASCII_OUTPUT', 'Binary_OUTPUT', 
         'SpinupDay', 'NUM_OPENMP', 'SCR_INTV',
         'ABSTOL', 'RELTOL', 
         'INIT_SOLVER_STEP', 'MAX_SOLVER_STEP', 'LSM_STEP', 
         'START', 'END', 
         dts)
  val = c(0,  3, 1,
        0, 1,
        0, 8, 1440,
        1e-4, 1e-4,
        1e-2, 10, 60,
        0, nday,
        vdt
  )
  val=data.frame(rbind(val))
  names(val) = toupper( vn )
  val
}
#' Generate the default model calibration
#' \code{shud.calib} 
#' @return Default calibration values for model
#' @export
shud.calib <- function(){
  cn = toupper( c(
    'GEOL_KSATH', 'GEOL_KSATV', 'GEOL_KMACSATH','GEOL_MACVF', 'GEOL_THETAS','GEOL_THETAR',  'GEOL_DMAC',
    'SOIL_KINF',  'SOIL_KMACSATV', 'SOIL_DINF', 'SOIL_ALPHA', 'SOIL_BETA', 'SOIL_MACHF', 
    'LC_VEGFRAC', 'LC_ALBEDO', 'LC_ROUGH', 'LC_DROOT','LC_ISMAX', 'LC_ImpAF', 'LC_SoilDgd', 
    'TS_PRCP', 'TS_LAI',
    'TS_SFCTMP+', 
    'ET_ETP', 'ET_IC', 'ET_TR', 'ET_Soil',
    'RIV_ROUGH', 'RIV_KH', 
    'RIV_SINU', 'RIV_CWR', 'RIV_BedThick',
    'RIV_BSLOPE+','RIV_DPTH+','RIV_WDTH+', 
    'IC_GW+', 'IC_RIV+',
    'AQ_DEPTH+' 
    ) )
  v=data.frame(rbind(rep(1, length(cn))))
  names(v) = toupper(cn)
  id=which(grepl('\\+', cn))
  # id
  v[id] = 0
  return(v)
}
