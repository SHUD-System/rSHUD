#' Pedotransfer function to generate soil/geol parameter from soil texture
#' \code{PTF} 
#' @param  x  data.frame or matrix, column = c(Silt_perc, Clay_perc, OrganicMatter _perc, BulkDensity g/cm3)
#' @param topsoil TRUE/FALSE, default = TRUE
#' @param xmin Lower limitation of the value of c(Silt_perc, Clay_perc, OrganicMatter _perc, BulkDensity g/cm3)
#' @return Hydraulic parameters, matrix
#' @export
#' @examples 
#' x=cbind(10,40, 1:20, 1.5)
#' y=PTF(x)
#' apply(y[,-1], 2, summary)
#' plot(x[,3], y[,2])
PTF <- function (x=t(matrix(c(33, 33, 2, 1.4), ncol=5, nrow=4) ), topsoil=TRUE, 
                 xmin=c(.1, .1, 1.3, 1.3) ){
  msg='PTF:: '
  #Wösten, J. H. M., Pachepsky, Y. a., & Rawls, W. J. (2001). Pedotransfer functions: Bridging the gap between available basic soil data and missing soil hydraulic characteristics. Journal of Hydrology, 251(3–4), 123–150. https://doi.org/10.1016/S0022-1694(01)00464-4
  if(is.matrix(x) || is.data.frame(x)){
    if(ncol(x) < 4){
      stop('Input X must be N*4 matrix or data.frame')
    }else if(ncol(x) == 4){
      y=x
    }else{
      y = x[, ncol(x) - 3:0 ]
    }
    
  }else{
    y = as.matrix(x, ncol=4);
  }
  nsoil  = nrow(y)
  ly =  matrix(ifelse(topsoil,1,0), nrow= nsoil, ncol=1)
  for (i in  1:nsoil){
    id= i;
    applymin <- function(x, xb, name){
      if (x <=xb) { #must be positive. minimum value 0.1_perc.
        message(msg, name, ' value ', i, ' is less than limit ', xb )
        x = xb
      } 
      return (x)
    }
    S = applymin(y[i, 1], xmin[1], 'SILT')
    C = applymin(y[i, 2], xmin[2], 'CLAY')
    OM = applymin(y[i, 3], xmin[3], 'Organic Mater')
    D = applymin(y[i, 4], xmin[4], 'Bulk Density')
    
    # S = y[i,1] #Silt. percentage to ratio
    # if (S <=xmin[1]) { #must be positive. minimum value 0.1_perc.
    #   warning('Non-positive SILT percentage, type ', i)
    #   S = 1/10
    # }
    # C = y[i,2] # Clay ratio. percentage to ratio
    # if (C <=xmin[2]){ #must be positive. minimum value 0.1_perc.
    #   warning('Non-positive CLAY percentage, type ', i)
    #   C = 1/10
    # }
    # OM = y[i,3] #Percentage Organic matter
    # if(OM <=xmin[3]) { #must be positive or zero
    #   warning('Non-positive OM percentage, type ', i)
    #   OM = 1.3
    # }
    # D = y[i,4]  # Bulk Density
    # if (D <= xmin[4]) { # cannot be negative. default = 1.3 g/cm3
    #   warning('Non-positive bulk density, type ', i)
    #   D = 1.3
    # }
    topsoil = ly[i];
    #KsatV
    #outData[i][1]= (7.755+0.03252*S+0.93*topsoil-0.967*D*D-0.000484*C*C-0.000322*S*S+0.001/S-0.0748/OM-0.643*log(S)-0.01398*D*C-0.1673*D*OM+0.02986*topsoil*C-0.03305*topsoil*S);
    KsatV=exp(7.755+0.03252*S+0.93*topsoil-0.967*D*D-0.000484*C*C-0.000322*S*S+0.001/S-0.0748/OM-0.643*log(S)-0.01398*D*C-0.1673*D*OM+0.02986*topsoil*C-0.03305*topsoil*S);
    #outData[i][1]=outData[i][1]/100;
    KsatV = KsatV / 100; #cm/day -> m/day
    #ThetaS
    #outData[i][2]=(0.7919+0.001691*C-0.29619*D-0.000001491*S*S+0.0000821*OM*OM+0.02427/C+0.01113/S+0.01472*log(S)-0.0000733*OM*C-0.000619*D*C-0.001183*D*OM-0.0001664*topsoil*S);
    ThetaS=  (0.7919+0.001691*C-0.29619*D-0.000001491*S*S+0.0000821*OM*OM+0.02427/C+0.01113/S+0.01472*log(S)-0.0000733*OM*C-0.000619*D*C-0.001183*D*OM-0.0001664*topsoil*S);
    #ThetaR
    #ThetaR=0.01;
    #InfD
    #InfD=0.10;
    #Alpha
    #outData[i][5]=log(-14.96+0.03135*C+0.0351*S+0.646*OM+15.29*D-0.192*topsoil-4.671*D*D-0.000781*C*C-0.00687*OM*OM+0.0449/OM+0.0663*log(S)+0.1482*log(OM)-0.04546*D*S-0.4852*D*OM+0.00673*topsoil*C);
    Alpha=100*exp(-14.96+0.03135*C+0.0351*S+0.646*OM+15.29*D-0.192*topsoil-4.671*D*D-0.000781*C*C-0.00687*OM*OM+0.0449/OM+0.0663*log(S)+0.1482*log(OM)-0.04546*D*S-0.4852*D*OM+0.00673*topsoil*C);
    #outData[i][5]=   exp(-14.96+0.03135*C+0.0351*S+0.646*OM+15.29*D-0.192*topsoil-4.671*D*D-0.000781*C*C-0.00687*OM*OM+0.0449/OM+0.0663*log(S)+0.1482*log(OM)-0.04546*D*S-0.4852*D*OM+0.00673*topsoil*C);
    # exp in Alpha function is what is used in this package .
    #Beta
    # if(Alpha <0){
    #   print(Alpha);
    # }
    Beta=1+exp(-25.23-0.02195*C+0.0074*S-0.1940*OM+45.5*D-7.24*D*D+0.0003658*C*C+0.002885*OM*OM-12.81/D-0.1524/S-0.01958/OM-0.2876*log(S)-0.0709*log(OM)-44.6*log(D)-0.02264*D*C+0.0896*D*OM+0.00718*topsoil*C);
    #hAreaF
    #hAreaF=0.01;
    #macKsatV
    #macKsatV=100*outData[i][1];
    #macKsatV = KsatV * 100
    #val = c(id, KsatV, ThetaS,ThetaR,InfD, Alpha, Beta, hAreaF, macKsatV)
    val = c(id, KsatV, ThetaS, Alpha, Beta)
    if( i==1 ){
      mat =matrix(val, ncol = 5);
    }else{
      mat= rbind(mat, val);
    }
  }
  #colnames(mat) = c('INDEX', 'KsatV(m_d)', 'ThetaS(m3_m3)', 'ThetaR(m3_m3)', 'InfD(m)', 'Alpha(1_m)', 'Beta', 'hAreaF(m2_m2)', 'macKsatV(m_d)')
  colnames(mat) = c('INDEX', 'KsatV(m_d)', 'ThetaS(m3_m3)',  'Alpha(1_m)', 'Beta')
  rownames(mat) = paste(1:nsoil)
  return(mat)
}

#' Pedotransfer function to generate soil parameter
#' \code{PTF.soil} 
#' @param topsoil TRUE/FALSE, default = TRUE
#' @param  x  data.frame or matrix, column = c(Silt_perc, Clay_perc, OrganicMatter _perc, BulkDensity g/cm3)
#' @param rm.outlier Whether replace the outliers in each column with the mean value of the column.
#' @return Hydraulic parameters, matrix
#' @export
#' @examples 
#' x=cbind(10,40, 1:20, 1.5)
#' y=PTF.soil(x)
#' apply(y[,-1], 2, summary)
#' plot(x[,3], y[,2])
PTF.soil <- function(x=t(matrix(c(33, 33, 2, 1.4), ncol=5, nrow=4) ),
                     topsoil=TRUE,
                     rm.outlier=FALSE
                     ){
  y = as.matrix(x, ncol = 4)
  nsoil  = nrow(y)
  ly = topsoil * matrix(1, nrow= nsoil)
  ptf = PTF(x, topsoil=topsoil)
  
  ret = matrix(0, ncol=9, nrow=nsoil)
  colnames(ret) = c('INDEX', 'KsatV(m_d)', 'ThetaS(m3_m3)', 
                    'ThetaR(m3_m3)', 'InfD(m)', 'Alpha(1_m)',
                    'Beta', 'hAreaF(m2_m2)', 'macKsatV(m_d)')
  ret[,1:3]= ptf[,1:3]
  ret[,2]= ptf[,2]    #
  ret[,4]=0.01        #ThetaR - residual
  ret[,5]= 0.1    #infiltration depth 10cm
  ret[,6:7]= ptf[,4:5] #alpha, beta in van genuchten.
  ret[,8] = 0.01; #macropore fraction, 1_perc
  ret[,9] = ptf[,2] * 100 # Kmp = kmx * 100;
  if(rm.outlier){
    for(i in 2:ncol(ret)){
      oid = which_outliers(ret[,i])
      if(i==2){
        oid = c(oid, which(ret[,i] < 1e-3))
      }
      ret[oid,i] = mean(ret[-oid, i])
    }
  }
  return(ret)
}


#' Pedotransfer function to generate soil parameter
#' \code{PTF.geol} 
#' @param  x  data.frame or matrix, column = c(Silt_perc, Clay_perc, OrganicMatter _perc, BulkDensity g/cm3)
#' @param topsoil TRUE/FALSE, default = FALSE
#' @param rm.outlier Whether replace the outliers in each column with the mean value of the column.
#' @return Hydraulic parameters, matrix
#' @export
#' @examples 
#' x=cbind(10,40, 1:20, 1.5)
#' y=PTF.geol(x)
#' apply(y[,-1], 2, summary)
#' plot(x[,3], y[,2])
PTF.geol <- function(x=t(matrix(c(33, 33, 2, 1.4), ncol=5, nrow=4) ),
                     topsoil=FALSE, rm.outlier = FALSE){
  y = as.matrix(x, ncol = 4)
  nsoil  = nrow(y)
  ly = topsoil * matrix(1, nrow= nsoil)
  ptf = PTF(x, topsoil=topsoil)
  
  ret = matrix(0, ncol=8, nrow=nsoil)
  colnames(ret) = c('INDEX', 'KsatH(m_d)','KsatV(m_d)',
                    'ThetaS(m3_m3)', 'ThetaR(m3_m3)', 
                    # 'Alpha(1_m)','Beta', 
                    'vAreaF(m2_m2)', 'macKsatH(m_d)', 'Dmac(m)')
  ret[,1] = ptf[,1]    #INDEX
  ret[,2] = ptf[,2] * 10    #Horizontal K Kh = Kv * 10
  ret[,3] = ptf[,2]      #Vertical K    /100 -- cm/day -> m/day
  ret[,4] = ptf[,3]       #porosity
  ret[,5] = 0.01        #ThetaR - residual
  # ret[,6:7] = ptf[,4:5] #alpha, beta in van genuchten.
  ret[,6] = 0.01;         #vertical macropore fraction, 1_perc
  ret[,7] = ptf[,2] * 1e3 # Kmp = kmx * 100,000;
  ret[,8] = 1            # Depth of Macropore;
  if(rm.outlier){
    for(i in 2:ncol(ret)){
      oid = which_outliers(ret[,i])
      if(i %in% c(2,3) ){
        oid = c(oid, which(ret[,i] < 1e-3))
      }
      ret[oid,i] = mean(ret[-oid, i])
    }
  }
  return(ret)
}
#' The default land cover parameters from University of Marryland Classification
#' \code{lc.default} 
#' @return Default land cover parameters
#' @export
lc.default <- function(){
  cn = toupper(c('INDEX', 'LAIMAX', 'RMIN', 'RSREF', 'ALBEDO',
                 'VEGFRAC', 'ROUGH', 'RZD', 'SOILDGRD', 'IMPAF') )
  v = rbind( 
    c(1,0.8540044,0.00202374,12613990.87,0.14262859,0.0799863,4.34E-07,0.6,0,0),
    c(2,6.8320352,0.00201167,6565959.705,0.19602872,0.6398904,6.37E-07,0.6,0,0),
    c(3,4.80592291,0.00189796,8615796.962,0.28386765,0.727328,4.28E-07,0.6,0,0),
    c(4,4.86573019,0.00186856,8555289.368,0.27348679,0.731848,3.41E-07,0.6,0,0),
    c(5,4.93749891,0.00183329,8482680.254,0.26102975,0.737272,2.37E-07,0.6,0,0),
    c(6,4.99730619,0.00180389,8422172.66,0.25064889,0.741792,1.50E-07,0.6,0,0),
    c(7,0.09427062,0.00202904,13382620.83,0.23489923,0.0863528,4.14E-07,0.6,0,0),
    c(8,6.2166,0.00197916,7188601.405,0.256808,0.770208,6.71E-07,0.6,0,0),
    c(9,4.983,0.00180555,8436646.37,0.243008,0.830208,6.71E-07,0.6,0,0),
    c(10,7.2126,0.00197916,6180938.642,0.236708,0.767208,6.02E-07,0.6,0,0),
    c(11,1.5321585,0.00121528,6536697.724,0.14975106,0.375003,3.12E-07,0.6,0,0),
    c(12,4.3038,0.00171875,7776000,0.259218,0.652968,4.17E-07,0.6,0,0),
    c(13,2.8692,0.00081325,5184000,0.14994972,0.5012199,2.78E-07,0.6,0,0),
    c(14,2.8692,0.00081325,5184000,0.14994972,0.5012199,2.78E-07,0.6,0,0),
    c(15,1.5325585,0.00202546,11927491.17,0.24261089,0.404959,4.75E-07,0.6,0,0),
    c(16,0.70052962,0.00205233,12769262.77,0.25277193,0.160861,4.75E-07,0.6,0,0) 
  )
  colnames(v) = cn
  return(v)
}

#' Generate the default land cover parameters of NLCD classes.
#' \code{PTF.lc} 
#' @param lc NLCD (2001 and later version) land use code. 
#' @return Default land cover parameters of NLCD classes.
#' @export
PTF.lc <- function(lc){
  dtab = rbind(
  # http://glcf.umiacs.umd.edu/data/landcover/
    # Index		LAIMAX		RMIN		RSREF		ALBEDO		VEGFRAC		ROUGH		Drz		SOILDGRD		IMPAF	
    c(  0	,	0	,	0.00202546	,	13477995.32	,	0.135	,	0	,	4.05E-07	,	0.00	,	0.00	,	0.00	), #Water
    c(  1	,	10.76	,	0.00202546	,	2592000	,	0.182	,	0.8	,	8.10E-07	,	0.60	,	0.00	,	0.00	), #Evergreen Needleleaf
    c(  2	,	5.117	,	0.00173611	,	8301077.283	,	0.213	,	0.9	,	8.10E-07	,	0.60	,	0.00	,	0.00	), #Evergreen Broadleaf
    c(  3	,	10.76	,	0.00202546	,	2592000	,	0.182	,	0.8	,	8.10E-07	,	0.60	,	0.00	,	0.00	), #Deciduous Needleleaf
    c(  4	,	7.173	,	0.00202546	,	6221002.342	,	0.236	,	0.8	,	8.10E-07	,	0.60	,	0.00	,	0.00	), #Deciduous Broadleaf
    c(  5	,	8.833	,	0.00202546	,	4541564.403	,	0.2025	,	0.795	,	6.94E-07	,	0.60	,	0.00	,	0.00	), #Mixed Forest
    c(  6	,	8.540044	,	0.00200822	,	4837950.801	,	0.2112859	,	0.799863	,	6.94E-07	,	0.00	,	0.00	,	0.00	), #Woodland
    c(  7	,	6.2395131	,	0.00195775	,	7165420.002	,	0.25245035	,	0.801837	,	5.79E-07	,	0.40	,	0.00	,	0.00	), #Wooded Grassland
    c(  8	,	2.5535975	,	0.00202546	,	10894496.21	,	0.2495851	,	0.625005	,	5.21E-07	,	0.40	,	0.00	,	0.00	), #Closed Shrubland
    c(  9	,	1.1668827	,	0.00207025	,	12297448.88	,	0.26652016	,	0.218175	,	5.21E-07	,	0.40	,	0.00	,	0.00	), #Open Shrubland
    c( 10	,	4.782	,	0.00190972	,	8640000	,	0.28802	,	0.72552	,	4.63E-07	,	0.40	,	0.1	,	0.00	), #Grassland
    c( 11	,	4.782	,	0.00135542	,	8640000	,	0.2499162	,	0.8353665	,	4.63E-07	,	0.40	,	0.50	,	0.00	), #Cropland
    c( 12	,	0.001	,	0.00202546	,	13476983.61	,	0.23214958	,	0.07489	,	4.05E-07	,	0.05	,	0.60	,	0.00	), #Bare Ground
    c( 13	,	5.0212291	,	0.00179213	,	8397969.622	,	0.24649654	,	0.7436	,	1.16E-07	,	0.05	,	0.90	,	0.90	) #Urban and Build
    )
  
  y= t(lc_EQ( t(dtab), lc))
  y[,1] = 1:nrow(y)
  colnames(y) = toupper(c('INDEX', 'LAIMAX', 'RMIN', 'RSREF', 'ALBEDO',
                   'VEGFRAC', 'ROUGH', 'RZD', 'SOILDGRD', 'IMPAF') )
  rownames(y)=lc
  return(y)
}
#' Generate the default melt factor
#' \code{MeltFactor} 
#' @param years Years.
#' @return Default default melt factor.
#' @export
#' @examples 
#' years = 2000:2001
#' x = MeltFactor(years=years)
#' plot(x)
MeltFactor <- function(years){
  mf=c(0.001308019, 0.001633298,  0.002131198, 0.002632776, 0.003031171,  0.003197325, 0.003095839, 0.00274524,     0.002260213, 0.001759481, 0.001373646,  0.001202083);
  years=sort(c(years,max(years)+1))
  yrlim=range(years);
  ny = length(years)
  t1=as.Date(paste(yrlim[1],'-01-01',sep=''))
  t2=as.Date(paste(yrlim[2]+1,'-01-01',sep=''))
  tdaily = seq.Date(t1,t2,by=1)
  DataDaily= xts::as.xts(numeric(length(tdaily)),order.by=tdaily)
  DataMon =  xts::apply.monthly(DataDaily,FUN=sum)
  #tmon = time(DataMon)- days_in_month(time(DataMon))+1
  tmon = as.Date( format(stats::time(DataMon), "%Y-%m-01"))
  ret = xts::as.xts(c(rep(mf, ny), mf[1]), order.by=tmon)
  colnames(ret) = 'MF'
  return(ret)
}