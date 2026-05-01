#' Calculate river order, downstream, slope, and length
#' \code{shud.river}
#' @param sl \code{sf} LINESTRING object; legacy \code{SpatialLines*} is
#'   accepted for compatibility.
#' @param dem Elevation as a \code{terra::SpatRaster}
#' @param rivord Order of each river reach
#' @param rivdown Downstream Index of each river reach.
#' @param AREA Area of the watershed, for estimating the width/depth of river.
#' @return SHUD.RIVER
#' @export
shud.river <- function(sl, dem, rivord = NULL, rivdown = NULL, AREA = NULL){
  rivers_sf <- if (inherits(sl, "sf")) sl else sf::st_as_sf(sl)
  dem_rast <- if (inherits(dem, "SpatRaster")) dem else terra::rast(dem)

  network <- build_river_network(
    rivers = rivers_sf,
    dem = dem_rast,
    area = AREA,
    river_order = rivord,
    downstream = rivdown
  )

  as_shud_river(network)
}

#' Calculate the river length
#' \code{RiverLength}
#' @param pr SHUD.RIVER class
#' @return River length, vector
#' @export
RiverLength <- function(pr){
  rv = pr@river
  pt = pr@point
  n = nrow(rv)
  pfr = rv[, 2]
  pto = rv[, 3]
  p0 = pt[pfr, 2:3]
  p1 = pt[pto, 2:3]
  d = rep(0, n)
  for(i in 1:n){
    d[i] <- Eudist(p0[i, ], p1[i, ])
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

  pfr = rv[, 'FrNode']
  pto = rv[, 'ToNode']
  z0 = pt[pfr, 'Zmax']
  z1 = pt[pto, 'Zmax']
  len = RiverLength(pr)
  s0 = (z0 - z1) / len

  z = 0.5 * (z0 + z1)
  #=============
  idown = rv[, 'Down']
  zz = pt[pto, 'Zmax']
  oid = which(idown <= 0)
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
#' @return Corrected \code{SHUD.RIVER} object
#' @export
correctRiverSlope <- function(pr, minSlope = 1e-5, maxrun = 500){
    msg = 'correctRiverSlope::'
  EPS = 1e-9
  # pr = ppr
  # minSlope = 1e-4
  nflag = 1
  pfr = pr@river[,'FrNode']
  pto = pr@river[,'ToNode']
  idown = pr@river[,'Down']
  len = RiverLength(pr)
  oid = getOutlets(pr)
  nloop = 0
  pz0 = pr@point[, 'Zmax']
  while(nflag > 0 & maxrun > nloop){
    nloop = nloop + 1
    pz = pr@point[, 'Zmax']
    s = RiverSlope(pr)
    rid = which(s[, 1] < minSlope )
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

        pz[p2] = pz[p1] - minSlope * ll - EPS
        #pz[pto[rid]] = pz[pfr[rid]] - len[rid] * minSlope - EPS
        pr@point[,'Zmax'] = pz
      }else{
        # for(ip in 1:length(rid)){
        # ip = 1
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
  nloop = 0
  while(nflag > 0 & maxrun > nloop){
    nloop = nloop + 1
    pz = pr@point[, 'Zmax']
    s = RiverSlope(pr)
    rid = which(s[, 2] < minSlope )
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
  dz0 = dz[ dz != 0]
  message(msg, length(dz0), ' points was changed.')
  print(summary(dz0))
  pr
}

#' Deprecated: read the SHUD.RIVER shapefile as legacy SpatialLines
#' \code{sp.riv2shp}
#' @param pr SHUD.RIVER class
#' @param dbf Attribute data for exported line features; retained for
#'   compatibility.
#' @return Legacy \code{SpatialLinesDataFrame}
#' @keywords deprecated
#' @export
sp.riv2shp <- function(pr = read_river(), dbf = NULL){
  # pt = pr@point[,2:3]
  # rt = pr@rivertype
  # riv = pr@river
  # 
  # nr = nrow(riv)
  # nt = nrow(rt)
  # np = nrow(pt)
  # 
  # ft = riv[,2:3]
  # sl = list()
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
  sp = as(sf::st_read(fn, quiet = TRUE), "Spatial")
  sp
}

#' Get outlet ID(s)
#' \code{getOutlets}
#' @param pr SHUD.RIVER class
#' @return Index of outlets, numeric
#' @export
getOutlets <- function(pr = read_river()){
  idown = pr@river[,'Down']
  ids = which(idown < 0)
  ids
}
