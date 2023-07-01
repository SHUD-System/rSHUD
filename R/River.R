#' calculate river order, downstream, slope and length
#' \code{shud.river}
#' @param sl SpatialLines*
#' @param dem Raster of elevation
#' @param rivord Order of each river reach
#' @param rivdown Downstream Index of each river reach.
#' @param AREA Area of the watershed, for estimating the width/depth of river.
#' @return SHUD.RIVER
#' @export
shud.river <- function(sl, dem, rivord  = NULL, rivdown=NULL, AREA=NULL){
  msg='shud.river::'
  # sp.slt = sp.Polylines2Lines(sp)
  sp.slt = sl
  nsp = length(sp.slt)
  xy = data.frame(extractCoords(sp.slt))
  # z = raster::extract(dem, xy)
  # pt = cbind(1:nrow(xy), xy, z)
  # colnames(pt) = c('Index','X', 'Y', 'Zmax')
  if(is.null(rivord)){
    message(msg, '\nCalculate river order ...')
    rivord = sp.RiverOrder(sp.slt)
  }
  if(is.null(rivdown)){
    message(msg, '\nIdentify the downstream ...')
    rivdown = sp.RiverDown(sp.slt)
  }
  message(msg, '\nFrom/To nodes ...') 
  ft = rbind(FromToNode(sp.slt, simplify=TRUE)[, 2:3])
  message(msg, '\nSlope and length of river ...')
  zf = raster::extract(dem, xy[ft[,1], ])
  zt = raster::extract(dem, xy[ft[,2], ])
  len = rgeos::gLength(sp.slt, byid = TRUE)
  slp = (zf - zt) / len;
  row.names(sp.slt) = paste(1:nsp)
  
  df = data.frame('Index'=1:nsp,
                  'Down' = rivdown,
                  'Type' = rivord,
                  'Slope' = slp,
                  'Length' = len,
                  'BC' = len *0)
  rownames(df) = row.names(sp.slt)
  
  p.df = data.frame(xy[ft[, 1], ], zf, xy[ft[,2], ], zt)
  colnames(p.df) = c('From.x','From.y','From.z', 'To.x', 'To.y', 'To.z')
  rownames(p.df) = row.names(sp.slt)
  
  sp.df = sp::SpatialLinesDataFrame(sp.slt, data=df)
  ntype = max(rivord)
  if(is.null(AREA)){
    rtype = RiverType(ntype)
  }else{
    fx  = function(a, n = 10){
      dd =  (1 / (1:n) ) ^ 0.8
      rev(8 * log10(a + 1)  * a  ^ 0.25 * dd)
    }
    wd = round(fx(AREA, ntype), 2)
    rtype = RiverType(ntype, width=wd)
  }
  p.z=raster::extract(dem, xy)
  SHUD.RIVER(river = df,
             rivertype = data.frame(rtype), 
             point=p.df
  )
}
#' Calculate the river length
#' \code{RiverLength}
#' @param pr SHUD.RIVER class
#' @return River length, vector
#' @export
RiverLength<- function(pr){
  rv = pr@river
  pt = pr@point
  n = nrow(rv)
  pfr = rv[,2]
  pto = rv[,3]
  p0 = pt[pfr, 2:3]
  p1 = pt[pto, 2:3]
  d = rep(0, n)
  for(i in 1:n){
    d[i] <- Eudist(p0[i,],p1[i,])
  }
  d
}
#' Calculate the river slope
#' \code{RiverSlope}
#' @param pr SHUD.RIVER class
#' @return River slopes, vector
#' @export
RiverSlope <- function(pr){
  rv = pr@river
  rt = pr@rivertype
  pt = pr@point

  pfr = rv[,'FrNode']
  pto = rv[,'ToNode']
  z0 = pt[pfr, 'Zmax']
  z1 = pt[pto, 'Zmax']
  len = RiverLength(pr)
  s0 = (z0 - z1) / len

  z = 0.5 * (z0 + z1)
  #=============
  idown = rv[, 'Down']
  zz = pt[pto, 'Zmax']
  oid=which(idown <= 0)
  idown[oid] = oid

  z0 = z
  z1 = z[idown]

  downlen = 0.5 * (len + len[idown] )
  s1 = (z0 - z1) / downlen
  s1[oid] = s0[oid]
  ret <- cbind(s0, s1)
}
#' Convert the bed slope of river, to avoid negative/zero slope
#' \code{correctRiverSlope}
#' @param pr SHUD.RIVER class
#' @param minSlope Minimum slope
#' @param maxrun Maximum number to run the loops
#' @return SpatialLinesDataFrame
#' @export
correctRiverSlope <- function(pr, minSlope = 1e-5, maxrun = 500){
    msg='correctRiverSlope::'
  EPS = 1e-9
  # pr = ppr
  # minSlope = 1e-4
  nflag = 1
  pfr = pr@river[,'FrNode']
  pto = pr@river[,'ToNode']
  idown = pr@river[,'Down']
  len = RiverLength(pr)
  oid = getOutlets(pr)
  nloop=0
  pz0 = pr@point[, 'Zmax']
  while(nflag > 0 & maxrun > nloop){
    nloop=nloop+1
    pz = pr@point[, 'Zmax']
    s=RiverSlope(pr)
    rid = which(s[,1] < minSlope )
    nflag = length(rid)
    if(nflag >= 1){
      message(msg, nflag, ' rivers in type 1 (reverse)')
      idn = idown[rid]
      od = which(rid %in% oid)
      if(length(od) > 0){
        otid = rid[od]

        p1 = pfr[otid]
        p2 = pto[otid]
        ll = len[otid]

        pz[p2] = pz[p1] - minSlope*ll - EPS
        #pz[pto[rid]] = pz[pfr[rid]] - len[rid] * minSlope - EPS
        pr@point[,'Zmax'] = pz
      }else{
        # for(ip in 1:length(rid)){
        # ip=1
        #   p1 = pfr[rid[ip]]
        #   p2 = pto[rid[ip]]
        #   p3 = pto[idn[ip]]
        #   l12 = len[rid[ip]]
        #   l23 = len[idn[ip]]
        #   ll = l12 + l23
          # pz[p2] = (pz[p1] * l23 + pz[p3] * l12) / ll
          pz[pto[rid]] = pz[pfr[rid]] - len[rid] * minSlope - EPS
        # }
        pr@point[,'Zmax'] = pz
        # print(cbind(rid,idn,p1,p2,p3, pz[p1], pz[p2], pz[p3]))
        # readline("go on?")
      }
    }
  }

  flag = 1
  idown = pr@river[,'Down']
  nloop=0
  while(nflag > 0 & maxrun > nloop){
    nloop=nloop+1
    pz = pr@point[, 'Zmax']
    s=RiverSlope(pr)
    rid = which(s[,2] < minSlope )
    nflag = length(rid)
    if(nflag > 0){
      message(msg, nflag, ' rivers in type 2 (revers down)')
      p1 = pfr[rid]
      p2 = pto[rid]
      p3 = pto[idown[rid]]
      print(rid)
      pz[p2] = 0.5 * (pz[p1] + pz[p2])
      pr@point[,'Zmax'] = pz
    }
  }
  dz = pz0 - pz
  dz0 = dz[ dz!= 0]
  message(msg, length(dz0), ' points was changed.')
  print(summary(dz0))
  pr
}

#' Convert the SHUD.RIVER class to SpatialLines
#' \code{sp.riv2shp}
#' @param pr SHUD.RIVER class
#' @param dbf Attribute data for exported SpatialLines
#' @return SpatialLinesDataFrame
#' @export
sp.riv2shp <- function(pr = readriv(), dbf=NULL){
  # pt = pr@point[,2:3]
  # rt = pr@rivertype
  # riv = pr@river
  # 
  # nr = nrow(riv)
  # nt = nrow(rt)
  # np = nrow(pt)
  # 
  # ft = riv[,2:3]
  # sl=list()
  # sls = list()
  # for(i in 1:nr){
  #   coord = rbind(pt[ft[i,1], ], pt[ft[i,2], ])
  #   sl[[i]] = sp::Line( coord )
  #   sls[[i]] = sp::Lines(sl[[i]], i)
  # }
  # if(is.null(dbf)){
  #   df = data.frame(riv, rt[riv[, 'Type'], ])
  # }else{
  #   df = dbf
  # }
  # spl = sp::SpatialLines(sls)
  # rownames(df) = names(spl)
  # sp = sp::SpatialLinesDataFrame(spl, data=df)
  fn = file.path(get('inpath', envir = .shud), 'gis', 'river.shp')
  sp = rgdal::readOGR(fn)
  sp
}

#' Generate the river segments table
#' \code{shud.rivseg}
#' @param sl SpatialLinesDataFrame
#' @return river segments table
#' @export
shud.rivseg <- function(sl){
  df = data.frame('Index' = 1:length(sl),
                  sl@data,
                  'Length' = rgeos::gLength(sl, byid = T))
  df
}
#' Get outlet ID(s)
#' \code{getOutlets}
#' @param pr SHUD.RIVER class
#' @return Index of outlets, numeric
#' @export
getOutlets <- function(pr=readriv()){
  idown = pr@river[,'Down']
  ids = which(idown < 0)
  ids
}
