
#' Generate the default LAI and Roughness Length
#' \code{LaiRf.NLCD}
#' @param lc land classe codes 
#' @param years numeric years.
#' @param class The classification of the Landuse of NLCD
#' @return Default LaiRf.NLCD parameters, a list $LAI, $RL
#' @export
#' @examples
#' lc = c(43, 23, 81, 11)
#' lr=LaiRf.NLCD(lc, years=2000:2001)
#' par(mfrow=c(2,1))
#' col=1:length(lc)
#' plot(lr$LAI, col=col, main='LAI');
#' legend('top', paste0(lc), col=col, lwd=1)
#' plot(lr$RL, col=col, main='Roughness Length');
#' legend('top', paste0(lc), col=col, lwd=1)
LaiRf.NLCD <- function(lc, years=2000){
  # http://glcf.umiacs.umd.edu/data/landcover/
  rltbl =matrix(c(
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
    1.112,	1.103,	1.088,	1.082,	1.076,	1.068,	1.073,	1.079,	1.082,	1.088,	1.103,	1.112,
    2.653,	2.653,	2.653,	2.653,	2.653,	2.653,	2.653,	2.653,	2.653,	2.653,	2.653,	2.653,
    1.112,	1.103,	1.088,	1.082,	1.076,	1.068,	1.073,	1.079,	1.082,	1.088,	1.103,	1.112,
    0.52,	0.52,	0.666,	0.91,	1.031,	1.044,	1.042,	1.037,	1.036,	0.917,	0.666,	0.52,
    0.816,	0.8115,	0.877,	0.996,	1.0535,	1.056,	1.0575,	1.058,	1.059,	1.0025,	0.8845,	0.816,
    0.7602524,	0.7551426,	0.7772204,	0.8250124,	0.846955,	0.8449668,	0.8471342,	0.8496604,	0.8514252,	0.8299022,	0.7857734,	0.7602744,
    0.35090494,	0.34920916,	0.36891486,	0.40567288,	0.42336056,	0.42338372,	0.42328378,	0.42485112,	0.42631836,	0.40881268,	0.37218526,	0.35096866,
    0.05641527,	0.05645892,	0.05557872,	0.05430207,	0.05425842,	0.05399002,	0.05361482,	0.0572041,	0.05892068,	0.05821407,	0.05709462,	0.05645892,
    0.03699235,	0.03699634,	0.03528634,	0.03272533,	0.03272134,	0.03270066,	0.03268178,	0.03907616,	0.04149324,	0.04032533,	0.03823134,	0.03699634,
    0.0777,	0.0778,	0.0778,	0.0779,	0.0778,	0.0771,	0.0759,	0.0766,	0.0778,	0.0779,	0.0778,	0.0778,
    0.0777,	0.0778,	0.0778,	0.0779,	0.0778,	0.0771,	0.0759,	0.0766,	0.0778,	0.0779,	0.0778,	0.0778,
    0.0112,	0.0112,	0.0112,	0.0112,	0.0112,	0.0112,	0.0112,	0.0112,	0.0112,	0.0112,	0.0112,	0.0112,
    0.1947138,	0.19413424,	0.20831414,	0.23348558,	0.24574614,	0.24605016, 0.24538258,	0.24630454,	0.247455,	0.23527388,	0.20963734,	0.19478494
  ), ncol =14, nrow=12)
  laitbl = matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    8.76,  9.16,  9.827,  10.093,  10.36,  10.76,  10.493,  10.227,  10.093,  9.827,  9.16,  8.76,
    5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117,
    8.76,  9.16,  9.827,  10.093,  10.36,  10.76,  10.493,  10.227,  10.093,  9.827,  9.16,  8.76,
    0.52,  0.52,  0.867,  2.107,  4.507,  6.773,  7.173,  6.507,  5.04,  2.173,  0.867,  0.52,
    4.64,  4.84,  5.347,  6.1,  7.4335,  8.7665,  8.833,  8.367,  7.5665,  6,  5.0135,  4.64,
    5.276088,  5.528588,  6.006132,  6.4425972,  7.2448806,  8.3639474,  8.540044,  8.126544,  7.2533006,  6.3291908,  5.6258086,  5.300508,
    2.3331824,  2.4821116,  2.7266101,  3.0330155,  3.8849492,  5.5212224,  6.2395131,  5.7733017,  4.1556703,  3.1274641,  2.6180116,  2.4039116 ,
    0.580555,  0.6290065,  0.628558,  0.628546,  0.919255,  1.7685454,  2.5506969,  2.5535975,  1.7286418,  0.9703975,  0.726358,  0.6290065 ,
    0.3999679,  0.4043968,  0.3138257,  0.2232945,  0.2498679,  0.3300675,  0.4323964,  0.7999234,  1.1668827,  0.7977234,  0.5038257,  0.4043968,
    0.782,  0.893,  1.004,  1.116,  1.782,  3.671,  4.782,  4.227,  2.004,  1.227,  1.004,  0.893,
    0.782,  0.893,  1.004,  1.116,  1.782,  3.671,  4.782,  4.227,  2.004,  1.227,  1.004,  0.893,
    0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001 ,
    1.2867143,  1.3945997,  1.5506977,  1.7727263,  2.5190228,  4.1367678,  5.0212291,  4.5795799,  2.8484358,  1.8856229,  1.5178736,  1.3656797
  ),
  ncol = 14, nrow = 12)
  
  years=sort(c(years, max(years, na.rm=TRUE)+1), decreasing = FALSE)
  nlc = length(lc)
  ny  = length(years)
  ym = rbind(expand.grid(1:12, years), c(1, max(years)+1) )
  tstr = paste(ym[,2], ym[,1], '01', sep='-')
  tday = zoo::as.Date(tstr)
  mlai = NULL
  mrl = NULL
    for (i in 1:nlc){
      ilc= lc[i]
      lai = lc_EQ(laitbl, ilc)
      rl = lc_EQ (rltbl, ilc)
      mlai = cbind(mlai, lai)
      mrl = cbind(mrl, rl)
    }
    colnames(mlai) = lc
    colnames(mrl) = lc

  rep.row<-function(x,n){
    ret=NULL
    for(i in 1:n){
      ret=rbind(ret, x)
    }
    return(ret)
  }

  ts.lai = xts::as.xts(rbind(rep.row(mlai, ny), mlai[1,]), order.by =tday)
  ts.rl = xts::as.xts(rbind(rep.row(mrl, ny), mrl[1,]), order.by =tday)
  # ts.lai = zoo::zoo(rbind(rep.row(mlai, ny), mlai[1,]), order.by =tday)
  # ts.rl = zoo::zoo(rbind(rep.row(mrl, ny), mrl[1,]), order.by =tday)
  colnames(ts.lai) = paste(lc)
  colnames(ts.rl) = paste(lc)
  ret = list(LAI = ts.lai, RL=ts.rl, nLAI=nlc)
  ret
}


#' Generate the default LAI and Roughness Length for USGS GLC
#' \code{LaiRf.GLC}
#' @param lc land classe codes 
#' @param years numeric years.
#' @param if.daily whether daily output
#' @return Default LaiRf.GLC parameters, a list $LAI, $RL
#' @export
#' @examples 
#' lc = 1:14
#' lr=LaiRf.GLC(lc, years=2000:2001)
#' par(mfrow=c(2,1))
#' col=1:length(lc)
#' plot(lr$LAI, col=col, main='LAI');
#' legend('top', paste0(lc), col=col, lwd=1)
#' plot(lr$RL, col=col, main='Roughness Length');
#' legend('top', paste0(lc), col=col, lwd=1)
LaiRf.GLC <- function(lc=NULL,  years=2000+1:2,  if.daily=FALSE){
  # from UMD to GLC classification
  mapid = c(1:11, 1, 12, 14, 12, 1, 13)
  rep.row<-function(x,n){
    for(i in 1:n){ if(i==1){        ret = x;      }else{        ret=rbind(ret, x)      } }
    return(ret)  }
  rep.col<-function(x,n){
    for(i in 1:n){ if(i==1){        ret = x;      }else{        ret=cbind(ret, x)      } }
    return(ret)  }
  years=sort(c(years,max(years)+1))
  yrlim = range(years);
  ny = length(years)
  t1 = as.Date(paste(yrlim[1],'-01-01',sep=''))
  t2 = as.Date(paste(yrlim[2],'-12-31',sep=''))
  tdaily = seq.Date(t1,t2,by=1)
  DataDaily = xts::as.xts(numeric(length(tdaily)),order.by=tdaily)
  DataMon = xts::apply.monthly(DataDaily,FUN=sum)
  tmon =as.Date( format(time(DataMon), "%Y-%m-01"))
  #tmon = time(DataMon)- days_in_month(time(DataMon))+1
  nlc = length(lc)
  l = matrix(0, nrow=12, ncol=nlc)
  r = matrix(0, nrow=12, ncol=nlc)
  lai0 = t(tsLaiRf(type=1)[mapid, ])
  rl0 = t(tsLaiRf(type=2)[mapid, ])
  # for (i in 1:nlc){
  #   l[,i] = tsLaiRf(type=1)[lc[i], ]
  #   r[,i] = tsLaiRf(type=2)[lc[i], ]
  # }
  if(is.null(lc)){ lc = 1:17}
  lmat = xts::as.xts(rep.row(lai0[, lc], ny), order.by=tmon)
  rmat = xts::as.xts(rep.row(rl0[, lc], ny), order.by=tmon)
  colnames(lmat) = lc-1
  colnames(rmat) = lc-1
  ret = list('LAI'=lmat, 'RL'=rmat)
  if(if.daily){
    ld = NA*rep.col(DataDaily, nlc);
    rd = NA*rep.col(DataDaily, nlc);
    ld[time(lmat),]=lmat
    rd[time(rmat),]=rmat
    ld=zoo::na.approx(ld)
    rd=zoo::na.approx(ld)
    colnames(ld)=lc
    colnames(rd)=lc
    ret=list('LAI'=ld, 'RL'=rd)
  }
  return(ret)
}


#' Generate the default LAI and Roughness Length for UMD classification system
#' \code{tsLaiRf}
#' @param type 1: LAI, 2:RL
#' @return Default LAI and RL
#' @export
#' @examples 
#' tsLaiRf(type=1)
#' tsLaiRf(type=2)
tsLaiRf <- function (type=1) {
  dlai=rbind(c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             c(   8.76,  9.16,  9.827,  10.093,  10.36,  10.76,  10.493,  10.227,  10.093,  9.827,  9.16,  8.76),
             c(   5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117,  5.117),
             c(    8.76,  9.16,  9.827,  10.093,  10.36,  10.76,  10.493,  10.227,  10.093,  9.827,  9.16,  8.76),
             c(    0.52,  0.52,  0.867,  2.107,  4.507,  6.773,  7.173,  6.507,  5.04,  2.173,  0.867,  0.52),
             c(    4.64,  4.84,  5.347,  6.1,  7.4335,  8.7665,  8.833,  8.367,  7.5665,  6,  5.0135,  4.64),
             c(    5.276088,  5.528588,  6.006132,  6.4425972,  7.2448806,  8.3639474,  8.540044,  8.126544,  7.2533006,  6.3291908,  5.6258086,  5.300508),
             c(   2.3331824,  2.4821116,  2.7266101,  3.0330155,  3.8849492,  5.5212224,  6.2395131,  5.7733017,  4.1556703,  3.1274641,  2.6180116,  2.4039116 ),
             c(   0.580555,  0.6290065,  0.628558,  0.628546,  0.919255,  1.7685454,  2.5506969,  2.5535975,  1.7286418,  0.9703975,  0.726358,  0.6290065 ),
             c(    0.3999679,  0.4043968,  0.3138257,  0.2232945,  0.2498679,  0.3300675,  0.4323964,  0.7999234,  1.1668827,  0.7977234,  0.5038257,  0.4043968),
             c(    0.782,  0.893,  1.004,  1.116,  1.782,  3.671,  4.782,  4.227,  2.004,  1.227,  1.004,  0.893),
             c(    0.782,  0.893,  1.004,  1.116,  1.782,  3.671,  4.782,  4.227,  2.004,  1.227,  1.004,  0.893),
             c(   0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001 ),
             c(    1.2867143,  1.3945997,  1.5506977,  1.7727263,  2.5190228,  4.1367678,  5.0212291,  4.5795799,  2.8484358,  1.8856229,  1.5178736,  1.3656797)
  );
  drl=rbind(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            c(	1.112, 1.103, 1.088, 1.082, 1.076, 1.068, 1.073, 1.079, 1.082, 1.088, 1.103, 1.112),
            c(	2.653, 2.653, 2.653, 2.653, 2.653, 2.653, 2.653, 2.653, 2.653, 2.653, 2.653, 2.653),
            c(	1.112, 1.103, 1.088, 1.082, 1.076, 1.068, 1.073, 1.079, 1.082, 1.088, 1.103, 1.112),
            c(	0.52, 0.52, 0.666, 0.91, 1.031, 1.044, 1.042, 1.037, 1.036, 0.917, 0.666, 0.52),
            c(	0.816, 0.8115, 0.877, 0.996, 1.0535, 1.056, 1.0575, 1.058, 1.059, 1.0025, 0.8845, 0.816),
            c(	0.7602524, 0.7551426, 0.7772204, 0.8250124, 0.846955, 0.8449668, 0.8471342, 0.8496604, 0.8514252, 0.8299022, 0.7857734, 0.7602744),
            c(	0.35090494, 0.34920916, 0.36891486, 0.40567288, 0.42336056, 0.42338372, 0.42328378, 0.42485112, 0.42631836, 0.40881268, 0.37218526, 0.35096866),
            c(	0.05641527, 0.05645892, 0.05557872, 0.05430207, 0.05425842, 0.05399002, 0.05361482, 0.0572041, 0.05892068, 0.05821407, 0.05709462, 0.05645892),
            c(	0.03699235, 0.03699634, 0.03528634, 0.03272533, 0.03272134, 0.03270066, 0.03268178, 0.03907616, 0.04149324, 0.04032533, 0.03823134, 0.03699634),
            c(	0.0777, 0.0778, 0.0778, 0.0779, 0.0778, 0.0771, 0.0759, 0.0766, 0.0778, 0.0779, 0.0778, 0.0778),
            c(	0.0777, 0.0778, 0.0778, 0.0779, 0.0778, 0.0771, 0.0759, 0.0766, 0.0778, 0.0779, 0.0778, 0.0778),
            c(	0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112),
            c(	0.1947138, 0.19413424, 0.20831414, 0.23348558, 0.24574614, 0.24605016, 0.24538258, 0.24630454, 0.247455, 0.23527388, 0.20963734, 0.19478494)
  );
  tab=switch(type,'lai'=dlai, 'rl'=drl)
  return(tab)
}

#' convert equation from Univeristy of Marryland classes to NLCD classes
#' \code{lc_EQ}
#' @param dtab conversion table.
#' @param lc land classe codes in NLCD classes.
#' @return Converted values.
#' @references Lele Shu(2017), Gopal Bhatt(2012)
#' @export
lc_EQ <-function(dtab, lc){
  #Source: http://www.pihm.psu.edu/EstimationofVegitationParameters.htm
  lc = as.numeric(lc)
  tab = matrix(c( #0, 0.00, 0, 0.00, #00
    0	,	0	,	0	,	0		,	#	1	0
    0	,	0	,	0	,	0		,	#	2	0
    0	,	0	,	0	,	0		,	#	3	0
    0	,	0	,	0	,	0		,	#	4	0
    0	,	0	,	0	,	0		,	#	5	0
    0	,	0	,	0	,	0		,	#	6	0
    0	,	0	,	0	,	0		,	#	7	0
    0	,	0	,	0	,	0		,	#	8	0
    0	,	0	,	0	,	0		,	#	9	0
    0	,	1	,	0	,	0		,	#	10	1
    0	,	0.9	,	6	,	0.1		,	#	11	1	Open Water
    0	,	0.9	,	6	,	0.1		,	#	12	1	Perennial Ice/Snow
    0	,	0	,	0	,	0		,	#	13	0
    0	,	0	,	0	,	0		,	#	14	0
    0	,	0	,	0	,	0		,	#	15	0
    0	,	0	,	0	,	0		,	#	16	0
    0	,	0	,	0	,	0		,	#	17	0
    0	,	0	,	0	,	0		,	#	18	0
    0	,	0	,	0	,	0		,	#	19	0
    13	,	0.65	,	0	,	0		,	#	20	1
    13	,	0.2	,	10	,	0.8		,	#	21	1	Developed, Open Space
    13	,	0.35	,	10	,	0.65		,	#	22	1	Developed, Low Intensity
    13	,	0.65	,	10	,	0.35		,	#	23	1	Developed, Medium Intensity
    13	,	0.9	,	10	,	0.1		,	#	24	1	Developed High Intensity
    0	,	0	,	0	,	0		,	#	25	0
    0	,	0	,	0	,	0		,	#	26	0
    0	,	0	,	0	,	0		,	#	27	0
    0	,	0	,	0	,	0		,	#	28	0
    0	,	0	,	0	,	0		,	#	29	0
    12	,	1	,	0	,	0		,	#	30	1
    12	,	0.92	,	9	,	0.08		,	#	31	1	Barren Land (Rock/Sand/Clay)
    12	,	1	,	9	,	0		,	#	32	1
    0	,	0	,	0	,	0		,	#	33	0
    0	,	0	,	0	,	0		,	#	34	0
    0	,	0	,	0	,	0		,	#	35	0
    0	,	0	,	0	,	0		,	#	36	0
    0	,	0	,	0	,	0		,	#	37	0
    0	,	0	,	0	,	0		,	#	38	0
    0	,	0	,	0	,	0		,	#	39	0
    6	,	0.65	,	0	,	0		,	#	40	1
    4	,	0.6	,	10	,	0.4		,	#	41	1	Deciduous Forest
    2	,	0.6	,	10	,	0.4		,	#	42	1	Evergreen Forest
    5	,	0.6	,	10	,	0.4		,	#	43	1	Mixed Forest
    0	,	0	,	0	,	0		,	#	44	0
    0	,	0	,	0	,	0		,	#	45	0
    0	,	0	,	0	,	0		,	#	46	0
    0	,	0	,	0	,	0		,	#	47	0
    0	,	0	,	0	,	0		,	#	48	0
    0	,	0	,	0	,	0		,	#	49	0
    8	,	1	,	0	,	0		,	#	50	1
    9	,	0.4	,	10	,	0.6		,	#	51	1	Dwarf Scrub
    8	,	0.6	,	10	,	0.4		,	#	52	1	Shrub/Scrub
    0	,	0	,	0	,	0		,	#	53	0
    0	,	0	,	0	,	0		,	#	54	0
    0	,	0	,	0	,	0		,	#	55	0
    0	,	0	,	0	,	0		,	#	56	0
    0	,	0	,	0	,	0		,	#	57	0
    0	,	0	,	0	,	0		,	#	58	0
    0	,	0	,	0	,	0		,	#	59	0
    0	,	0	,	0	,	0		,	#	60	0
    0	,	0	,	0	,	0		,	#	61	0
    0	,	0	,	0	,	0		,	#	62	0
    0	,	0	,	0	,	0		,	#	63	0
    0	,	0	,	0	,	0		,	#	64	0
    0	,	0	,	0	,	0		,	#	65	0
    0	,	0	,	0	,	0		,	#	66	0
    0	,	0	,	0	,	0		,	#	67	0
    0	,	0	,	0	,	0		,	#	68	0
    0	,	0	,	0	,	0		,	#	69	0
    10	,	0.85	,	0	,	0		,	#	70	1
    10	,	0.9	,	12	,	0.1		,	#	71	1	Grassland/Herbaceous
    10	,	0.9	,	12	,	0.1		,	#	72	1	Sedge/Herbaceous
    10	,	0.9	,	12	,	0.1		,	#	73	1	Lichens
    10	,	0.9	,	12	,	0.1		,	#	74	1	Moss
    0	,	0	,	0	,	0		,	#	75	0
    0	,	0	,	0	,	0		,	#	76	0
    0	,	0	,	0	,	0		,	#	77	0
    0	,	0	,	0	,	0		,	#	78	0
    0	,	0	,	0	,	0		,	#	79	0
    10	,	0.9	,	0	,	0		,	#	80	1
    11	,	0.6	,	10	,	0.4		,	#	81	1	Pasture/Hay
    11	,	0.9	,	12	,	0.1		,	#	82	1	Cultivated Crops
    0	,	0	,	0	,	0		,	#	83	0
    0	,	0	,	0	,	0		,	#	84	0
    0	,	0	,	0	,	0		,	#	85	0
    0	,	0	,	0	,	0		,	#	86	0
    0	,	0	,	0	,	0		,	#	87	0
    0	,	0	,	0	,	0		,	#	88	0
    0	,	0	,	0	,	0		,	#	89	0
    8	,	0.6	,	0	,	0.4		,	#	90	1	Woody Wetlands
    6	,	0.6	,	0	,	0		,	#	91	1
    7	,	0.6	,	0	,	0		,	#	92	1
    6	,	0.6	,	0	,	0		,	#	93	1
    7	,	0.6	,	0	,	0		,	#	94	1
    10	,	0.8	,	0	,	0.2		,	#	95	1	Emergent Herbaceous Wetlands
    10	,	0.8	,	2	,	0.2		,	#	96	1
    10	,	0.8	,	4	,	0.2		,	#	97	1
    10	,	0.2	,	0	,	0.8		,	#	98	1
    10	,	0.2	,	0	,	0.8		)	#	99	1
    , ncol=99, nrow=4)
  c1= tab[2, lc]
  v1 = dtab[, tab[1, lc]+1 ]
  c2 = tab[4, lc]
  v2 =  dtab[, tab[3,lc]+1  ]
  ret <- c1 * v1 + c2 * v2;
}


#' The default land cover parameters from UMD
#' \code{lc.UMD} 
#' @return Default land cover parameters
#' @export
lc.UMD <- function(){
  dtab = rbind(
    # http://glcf.umiacs.umd.edu/data/landcover/
    # Index		LAIMAX		RMIN		RSREF		ALBEDO		VEGFRAC		ROUGH(S/M^(1/3))		Drz		SOILDGRD		IMPAF	
    c(  1	,	0	,	0.00202546	,	13477995.32	,	0.135	,	0	,	3.5E-02	,	0.00	,	0.00	,	0.00	), #Water
    c(  2	,	10.76	,	0.00202546	,	2592000	,	0.182	,	0.8	,	7.0E-02	,	0.60	,	0.00	,	0.00	), #Evergreen Needleleaf
    c(  3	,	5.117	,	0.00173611	,	8301077.283	,	0.213	,	0.9	,	7.0E-02	,	0.60	,	0.00	,	0.00	), #Evergreen Broadleaf
    c(  4	,	10.76	,	0.00202546	,	2592000	,	0.182	,	0.8	,	7.0E-02	,	0.60	,	0.00	,	0.00	), #Deciduous Needleleaf
    c(  5	,	7.173	,	0.00202546	,	6221002.342	,	0.236	,	0.8	,	7.0E-02	,	0.60	,	0.00	,	0.00	), #Deciduous Broadleaf
    c(  6	,	8.833	,	0.00202546	,	4541564.403	,	0.2025	,	0.795	,	6.0E-02	,	0.60	,	0.00	,	0.00	), #Mixed Forest
    c(  7	,	8.540044	,	0.00200822	,	4837950.801	,	0.2112859	,	0.799863	,	6.0E-02	,	0.00	,	0.00	,	0.00	), #Woodland
    c(  8	,	6.2395131	,	0.00195775	,	7165420.002	,	0.25245035	,	0.801837	,	5.0E-02	,	0.40	,	0.00	,	0.00	), #Wooded Grassland
    c(  9	,	2.5535975	,	0.00202546	,	10894496.21	,	0.2495851	,	0.625005	,	4.5E-02	,	0.40	,	0.00	,	0.00	), #Closed Shrubland
    c(  10	,	1.1668827	,	0.00207025	,	12297448.88	,	0.26652016	,	0.218175	,	4.5E-02	,	0.40	,	0.00	,	0.00	), #Open Shrubland
    c( 11	,	4.782	,	0.00190972	,	8640000	,	0.28802	,	0.72552	,	4.0E-02	,	0.40	,	0.1	,	0.00	), #Grassland
    c( 12	,	4.782	,	0.00135542	,	8640000	,	0.2499162	,	0.8353665	,	4.0E-02	,	0.40	,	0.50	,	0.00	), #Cropland
    c( 13	,	0.001	,	0.00202546	,	13476983.61	,	0.23214958	,	0.07489	,	3.5E-02	,	0.05	,	0.60	,	0.00	), #Bare Ground
    c( 14	,	5.0212291	,	0.00179213	,	8397969.622	,	0.24649654	,	0.7436	,	1.0E-02	,	0.05	,	0.90	,	0.90	) #Urban and Build
  )
  cn = toupper(c('INDEX', 'LAIMAX', 'RMIN', 'RSREF', 'ALBEDO',
                 'VEGFRAC', 'ROUGH', 'RZD', 'SOILDGRD', 'IMPAF') )
  colnames(dtab) = cn
  dtab = dtab[, -1 * (2:4)]
  return(dtab)
}

#' Generate the default land cover parameters of NLCD classes.
#' \code{PTF.NLCD} 
#' @param lc NLCD (2001 and later version) land use code. 
#' @return Default land cover parameters of NLCD classes.
#' @export
lc.NLCD <- function(lc){
  dtab = lc.UMD()
  y= t(lc_EQ( t(dtab), lc))
  y[,1] = 1:nrow(y)
  colnames(y) = toupper(c('INDEX', 'ALBEDO',
                          'VEGFRAC', 'ROUGH', 'RZD', 'SOILDGRD', 'IMPAF') )
  rownames(y)=lc
  return(y)
}

#' The default land cover parameters from USGS GLC
#' \code{lc.GLC} 
#' @return Default land cover parameters
#' @export
lc.GLC <- function(){
  cn = toupper(c('INDEX', 'ALBEDO', 'VEGFRAC', 'ROUGH', 'RZD', 'SOILDGRD', 'IMPAF') )
  mapid = c(1:11, 1, 12, 14, 12, 1, 13)
  v = lc.UMD()[mapid, ]
  v[, 1] = 1:nrow(v)-1
  return(v)
}

