#' Write ESRI shapefile out
#' \code{writeshape}
#' @param shp \code{sf} object, or legacy \code{sp} \code{Spatial*} object (converted
#'   with \code{sf::st_as_sf()} before writing).
#' @param crs CRS to assign before write (default: CRS taken from \code{shp} via
#'   \code{sf::st_crs()} or \code{raster::crs()} for \code{sp} objects). Use
#'   \code{NA} to skip \code{st_set_crs} (not recommended for new data).
#' @param file file path, with or without \code{.shp} extension.
#' @details
#' All outputs are written with \code{sf::st_write()}. \code{SpatialPolygons},
#' \code{SpatialLines}, and \code{SpatialPoints} without a data slot are given
#' a synthetic \code{ID} column so conversion to \code{sf} is well-defined.
#' @export
#' @examples
#' \dontrun{
#' f <- file.path(tempdir(), "out")
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   nc <- sf::read_sf(system.file("shape/nc.shp", package = "sf"))
#'   writeshape(nc, file = f)
#' }
#' }
writeshape <- function(shp, file = NULL, crs = NULL) {
  msg <- "writeshape::"
  if (is.null(crs)) {
    crs <- if (inherits(shp, "sf")) {
      sf::st_crs(shp)
    } else {
      raster::crs(shp)
    }
  }
  if (!is.null(file) && tolower(tools::file_ext(file)) == "shp") {
    file <- sub("\\.[Ss][Hh][Pp]$", "", file)
  }
  if (is.null(file)) {
    return(invisible(shp))
  }

  path <- dirname(file)
  fn <- basename(file)
  if (!dir.exists(path)) {
    dir.create(path, showWarnings = TRUE, recursive = TRUE)
  }

  if (!inherits(shp, "sf")) {
    if (inherits(shp, "SpatialPolygons") &&
        !inherits(shp, "SpatialPolygonsDataFrame")) {
      shp <- sp::SpatialPolygonsDataFrame(
        shp,
        data = data.frame(ID = seq_len(length(shp)), stringsAsFactors = FALSE),
        match.ID = FALSE
      )
    } else if (inherits(shp, "SpatialLines") &&
               !inherits(shp, "SpatialLinesDataFrame")) {
      shp <- sp::SpatialLinesDataFrame(
        shp,
        data = data.frame(ID = seq_len(length(shp)), stringsAsFactors = FALSE),
        match.ID = FALSE
      )
    } else if (inherits(shp, "SpatialPoints") &&
               !inherits(shp, "SpatialPointsDataFrame")) {
      shp <- sp::SpatialPointsDataFrame(
        shp,
        data = data.frame(ID = seq_len(length(shp)), stringsAsFactors = FALSE)
      )
    }
    shp <- tryCatch(
      sf::st_as_sf(shp),
      error = function(e) {
        stop(
          "writeshape: could not convert object to sf: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

  if (!is.na(crs)) {
    shp <- sf::st_set_crs(shp, sf::st_crs(crs))
  }

  sf::st_write(
    shp,
    dsn = file.path(path, paste0(fn, ".shp")),
    driver = "ESRI Shapefile",
    delete_dsn = TRUE,
    quiet = TRUE
  )
  message(msg, file, " is saved")
  invisible(shp)
}

#' Re-project coordinates betwen GCS and PCS
#' \code{ProjectCoordinate} 
#' @param  x 2-column matrix of coordinates.
#' @param  proj4string proj4string
#' @param  P2G if TRUE, a cartographic projection into lat/long, otherwise projects from lat/long into a cartographic projection.
#' @return Basic model infomation, figures and tables
#' @export
ProjectCoordinate <- function(x, proj4string, P2G = TRUE){
  # Transformed data
  x = as.matrix(x)
  y <- proj4::project(x, proj4string, inverse = P2G)
  if(P2G){
    colnames(y) = c('Lon', 'Lat')
  }else{
    colnames(y) = c('X', 'Y')
  }
  y
}

#' SpatialData to Raster
#' \code{sp2raster}
#' @param sp SpatialPolygon
#' @param mask Raster mask of mesh domain.
#' @param ngrids Number of grid along x direction.
#' @param resolution Resolution, defaul = NULL, resolution = extent / ngrids
#' @param field Index of field
#' @return Raster map
#' @export
sp2raster <- function(sp, mask = get('MASK', envir = .shud),
                      ngrids = 200, 
                      resolution = NULL, field = 1) {
  if(is.null(mask)){
    ext <- raster::extent(sp)
    xlim = ext[1:2]
    ylim = ext[3:4]
    if(resolution <= 0 || is.null(resolution)){
      dx = diff(xlim) / ngrids;
    }else{
      dx = resolution
    }
    r <- raster::raster(ext, res = dx)
  }else{
    r = mask
  }
  ## Rasterize the shapefile
  rr <- raster::rasterize(sp, r, field = field)
  return(rr)
}

#' Generate the raster mask of Mesh domain
#' \code{shud.mask}
#' @param pm \code{shud.mesh}
#' @param n Number of grid
#' @param rr Default mask in .shud environment
#' @param cellsize Resolution, defaul = NULL, resolution = extent / ngrids
#' @param proj Projection parameter
#' @return Raster map
#' @export
shud.mask  <- function (pm = readmesh(), proj = NULL,
                        rr = NULL,
                        n = 10000, cellsize = NULL){
  # mesh = readmesh(shp = TRUE); ngrids = 100; resolution = 0
  if (is.null(rr)) {
    rr <- tryCatch(get('MASK', envir = .shud), error = function(e) NULL)
  }
  if(is.null(rr)){
    spm = mesh_to_sf(pm)

    if (inherits(spm, "sf")) {
      sp0 <- sf::st_union(spm)
      bb <- sf::st_bbox(sp0)
      dx <- cellsize
      if (is.null(dx)) {
        xrange <- as.numeric(bb["xmax"] - bb["xmin"])
        yrange <- as.numeric(bb["ymax"] - bb["ymin"])
        dx <- sqrt((xrange * yrange) / n)
      }

      rr <- terra::rast(
        xmin = bb["xmin"], xmax = bb["xmax"],
        ymin = bb["ymin"], ymax = bb["ymax"],
        resolution = dx
      )
      rr[] <- 1
      rr <- terra::mask(rr, terra::vect(sp0))

      if(!is.null(proj)){
        terra::crs(rr) <- proj
      } else if (!is.na(sf::st_crs(sp0))) {
        terra::crs(rr) <- sf::st_crs(sp0)$wkt
      }
    } else {
      sp0 = rgeos::gUnaryUnion(spm)
      if(is.null(cellsize)){
        # grd <- as.data.frame(sp::spsample(spm, "regular", n = n))
        # grd <- as.data.frame(sp::spsample(spm, "regular", nsig = 2, n = n))
        grd <- sp::makegrid(sp0, n = n)
      }else{
        # grd <- as.data.frame(sp::spsample(spm, "regular", cellsize = cellsize))
        grd <- sp::makegrid(sp0, cellsize = cellsize, pretty = FALSE)
      }
      names(grd)       <- c("X", "Y")
      sp::coordinates(grd) <- c("X", "Y")
      sp::gridded(grd)     <- TRUE  # Create SpatialPixel object
      sp::fullgrid(grd)    <- TRUE  # Create SpatialGrid object
      rr = raster::raster(grd); rr[] = 1
      rr = raster::mask(rr, sp0)
      if(!is.null(proj)){
        raster::crs(rr) = proj
      }
    }
    assign('MASK', rr, envir = .shud)
  }else{
    rr = rr
  }
  rr
}

#' SpatialData to Raster
#' \code{MeshData2Raster}
#' @param x vector or matrix, length/nrow is number of cells.
#' @param rmask mask of the mesh file
#' @param stack Whether export the stack, only when the x is a matrix, i.e. (Ntime x Ncell).
#' @param proj Projejction parameter
#' @param pm shud mesh
#' @param method method for interpolation, default = 'idw'
#' @param plot Whether plot the result.
#' @return Raster map
#' @export
MeshData2Raster <- function(x = getElevation(),
                            rmask = shud.mask(proj = proj), 
                            pm = readmesh(), proj = NULL,
                            stack = FALSE, method = 'ide',
                            plot = FALSE){
  
  if(stack){
    ret <- raster::stack(apply(x, 1, FUN = MeshData2Raster) )
  }else{
    if(is.matrix(x) | is.data.frame(x)){
      x = as.numeric(x[nrow(x),])
    }
    if(any(is.na(x)) ){
      x[is.na(x)] = 0
    }
    if (any(is.infinite(x))){
      x[is.infinite(x)] = 0
    }
    xy = getCentroid(pm = pm)[, 2:3]
    
    if(grepl('idw', tolower(method))){
      val = data.frame(xy, x)
      colnames(val) = c('X', 'Y', 'Z')
      sp::coordinates(val) = c('X', 'Y')
      grd = methods::as(rmask, 'SpatialGrid')
      # if(grepl(method, 'idw')){
      # Interpolate the grid cells using a power value of 2 (idp = 2.0)
      dat <- gstat::idw(Z ~ 1, val, newdata = grd, idp = 2.0)
      r = raster::raster(dat)
    }
    if(grepl('linear', tolower(method))){
      xr = raster::rasterToPoints(rmask)
      ext = raster::extent(rmask); res = raster::res(rmask); hr = res / 2
      r0 = rmask; r0[] = 1
      xyo = raster::rasterToPoints(r0)
      xx = interp::interp(x = xy[, 1], y = xy[, 2], z = x, xo = xyo[, 1], yo = xyo[, 2] )
      r = raster::setValues(rmask, as.numeric(xx$z))
    }
    if(grepl('ide', tolower(method))){
      tps <- fields::Tps(xy, x)
      # use model to predict values at all locations
      r <- raster::interpolate(rmask, tps)
    }
    ret <- raster::mask(r, rmask)
  }
  if(plot){
    raster::plot(ret)
  }
  if(!is.null(proj)){ raster::crs(ret) <- proj }
  return(ret)
}


#' Generatue fishnet
#' \code{fishnet} Generate fishnet by the coordinates
#' @param xx  x coordinates
#' @param yy  y coordinates
#' @param crs projections parameters, defaul = epsg:4326
#' @param type option = 'polygon', 'points', 'line'
#' @return spatildata (.shp)
#' @export
#' @examples
#' library(raster)
#' ext = c(0, 8, 2, 10)
#' dx = 2; dy = 4
#' xx = seq(ext[1], ext[2], by = dx)
#' yy = seq(ext[3], ext[4], by = dy)
#' sp1 = fishnet(xx = ext[1:2], yy = ext[3:4])
#' sp2 = fishnet(xx = xx + .5 * dx, yy = yy + 0.5 * dy)
#' sp3 = fishnet(xx = xx, yy = yy, type = 'point')
#' plot(sp1, axes = TRUE, xlim = c(-1, 1) * dx + ext[1:2], ylim = c(-1, 1) * dy + ext[3:4])
#' plot(sp2, axes = TRUE, add = TRUE, border = 2)
#' plot(sp3, axes = TRUE, add = TRUE, col = 3, pch = 20)
#' grid()
fishnet <- function(xx, yy,
                    crs = sp::CRS("+init=epsg:4326"),
                    type = 'polygon'){
  nx = length(xx)
  ny = length(yy)
  if(grepl('line', tolower(type)) ){
    ymin = min(yy); ymax = max(yy)
    xmin = min(xx); xmax = max(xx)
    vline = cbind(xx, ymin, xx, ymax)
    hline =  cbind(xmin, yy, xmax, yy)
    mat = rbind(vline, hline)
    str = paste(paste('(', mat[, 1], mat[, 2], ',', mat[, 3], mat[, 4], ')'), collapse = ',')
    spl = rgeos::readWKT(paste('MULTILINESTRING(', str, ')'))
    df = as.data.frame(mat)
    colnames(df) = c('x1', 'y1', 'x2', 'y2')
    spdf = sp::SpatialLinesDataFrame(spl, data = df)
    ret = spdf
  }
  if(grepl('polygon', tolower(type)) ){
    xy = expand.grid(xx, yy)
    xm = matrix(xy[, 1], nx, ny)
    ym = matrix(xy[, 2], nx, ny)
    xloc = abind::abind(as.matrix(xm[-nx, -ny]), as.matrix(xm[-nx, -1]), as.matrix(xm[-1, -1]),
                         as.matrix(xm[-1, -ny]), as.matrix(xm[-nx, -ny]), along = 3)
    yloc = abind::abind(as.matrix(ym[-nx, -ny]), as.matrix(ym[-nx, -1]), as.matrix(ym[-1, -1]),
                         as.matrix(ym[-1, -ny]), as.matrix(ym[-nx, -ny]), along = 3)
    
    # df = as.data.frame(matrix(0, nrow = (nx - 1) * (ny - 1), 6))
    df = data.frame(as.numeric(apply(xloc, 1:2, min)),
                     as.numeric(apply(xloc, 1:2, max)),
                     as.numeric(apply(yloc, 1:2, min)),
                     as.numeric(apply(yloc, 1:2, max)))
    df = data.frame(df, rowMeans(df[, 1:2]), rowMeans(df[, 1:2 + 2]) )
    colnames(df) = c('xmin', 'xmax', 'ymin', 'ymax', 'xcenter', 'ycenter')
    str = paste('GEOMETRYCOLLECTION(',
                 paste(paste('POLYGON((',
                             paste(xm[-nx, -ny], ym[-nx, -ny], ',' ),
                             paste(xm[-nx, -1],  ym[-nx, -1], ','),
                             paste(xm[-1, -1],   ym[-1, -1], ','),
                             paste(xm[-1, -ny],  ym[-1, -ny], ','),
                             paste(xm[-nx, -ny], ym[-nx, -ny], '' ), '))' )
                 , collapse = ','),
                 ')' )
    # str = paste('MULTIPOLYGON(', paste(xt, collapse = ', '), ')')
    SRL = rgeos::readWKT(str)
    # x1 = x0@polygons[[1]]@Polygons
    # SRL = lapply(1:length(x1),  function(x, i) {Polygons(list(x[[i]]), ID = i)},  x = x1 )
    ret = sp::SpatialPolygonsDataFrame(Sr = SRL, data = df, match.ID = TRUE)
  }
  if(grepl('point', tolower(type)) ){
    xm = expand.grid(xx, yy)
    df = data.frame('X' = xm[, 1], 'Y' = xm[, 2])
    ret = sp::SpatialPointsDataFrame(coords = df, data = df, match.ID = TRUE)
  }
  raster::crs(ret) = crs
  return(ret)
}

#' Add holes into Polygons
#' \code{AddHoleToPolygon}
#' @param poly SpatialPolygons
#' @param hole Hole Polygon
#' @export
AddHoleToPolygon <- function(poly, hole){
  # https://stackoverflow.com/questions/29624895/how-to-add-a-hole-to-a-polygon-within-a-spatialpolygonsdataframe
  # invert the coordinates for Polygons to flag it as a hole
  coordsHole <-  hole@polygons[[1]]@Polygons[[1]]@coords
  newHole <- sp::Polygon(coordsHole, hole = TRUE)
  
  # punch the hole in the main poly
  listPol <- poly@polygons[[1]]@Polygons
  listPol[[length(listPol) + 1]] <- newHole
  punch <- sp::Polygons(listPol, poly@polygons[[1]]@ID)
  
  # make the polygon a SpatialPolygonsDataFrame as the entry
  new <- sp::SpatialPolygons(list(punch), proj4string = poly@proj4string)
  new <- sp::SpatialPolygonsDataFrame(new, data = as(poly, "data.frame"))
  return(new)
}
#' Cut sptialLines with threshold.
#' \code{sp.CutSptialLines}
#' @param sl SpatialLines or SpatialLineDataFrame
#' @param tol Tolerence. If the length of segment is larger than tolerance, cut the segment until the maximum segment is shorter than tolerance.
#' @export
#' @examples
#' library(rSHUD)
#' library(sp)
#' x = 1:1000 / 100
#' l1 = Lines(Line(cbind(x, sin(x)) ), ID = 'a' )
#' sl = SpatialLines( list(l1) )
#' tol1 = 5;
#' tol2 = 2
#' sl1 = sp.CutSptialLines(sl, tol1)
#' sl2 = sp.CutSptialLines(sl, tol2)
#' par(mfrow = c(1, 2))
#' plot(sl1, col = 1:length(sl1)); title(paste0('Tol=', tol1))
#' plot(sl2, col = 1:length(sl2)); title(paste0('Tol=', tol2))
#'
#' data(sh)
#' riv = sh$riv
#' x = sp.CutSptialLines(riv, tol = 5)
#' par(mfrow = c(2, 1))
#' plot(riv, col = 1:length(riv), lwd = 3);
#'
#' plot(riv, col = 'gray', lwd = 3);
#' plot(add = TRUE, x, col = 1:length(x))
sp.CutSptialLines <- function(sl, tol){
  msg = 'sp.CutSptialLines::'
  ll = rgeos::gLength(sl, byid = TRUE)
  if(all(ll < tol) ){
    ret = sl
  }else{
    nsp = length(sl)
    xsl = list(); ik = 1
    for(i in 1:nsp){
      sx = sl[i, ]
      pxy = get_coords(sx)
      np = nrow(pxy)
      dacc = cumsum( sp::LineLength(pxy, sum = FALSE))
      # dacc = getDist(pxy)
      tol = max(c(tol, min(dacc) ) )
      len = rgeos::gLength(sx)
      if(len > tol){
        nsplit = ceiling(len / tol)
      }else{
        nsplit = 1
      }
      dd = len / nsplit
      v0 = 1  # Vetex 0, Vetex 1
      message(msg, i, '/', nsp, '\t', nsplit, '\t', round(dd, 2) )
      for(k in 1:nsplit){
        if(v0 >= np){
          break
        }
        # message(msg, '\t', k, '/', nsplit)
        dk = dd * k
        v1 = order(abs(dacc - dk), decreasing = FALSE)[1] + 1
        if(v1 + 1 > np){
          v1 = np
        }
        message(msg, v0, '\t', v1)
        if(v0 == v1){
          next;
        }
        # plot(sl[i, ]); points(pxy); points(pxy[c(v0, v1), ], pch = 2, col = 2)
        xsl[[ik]] = sp::Lines(sp::Line( pxy[c(v0:v1), ]), ID = ik)
        ik = ik + 1
        # points(pxy[v0:v1,], col = k)
        v0 = v1
      }
    }
    nsl = length(xsl)
    tmp = sp::SpatialLines(xsl, proj4string = raster::crs(sl))
    ilen = rgeos::gLength(tmp, byid = TRUE)
    att = data.frame('INDEX' = 1:length(tmp), 'Length' = ilen)
    ret = sp::SpatialLinesDataFrame(tmp, data = att)
  }
  return(ret)
}

#' Extract values on Raster map along a normalized line from 0 to 1.
#' \code{extractRaster}
#' @param r Raster
#' @param xy Coordinates of the line, dim = (Npoints, 2); x and y must be in
#'   the range 0 to 1.
#' @param ext extension of value xy.
#' @param plot Whether plot result.
#' @importFrom grDevices dev.off graphics.off png rgb topo.colors
#' @importFrom graphics grid hist lines par plot points
#' @importFrom methods as
#' @importFrom stats dist rnorm time
#' @importFrom utils read.table
#' @export
#' @examples
#' library(raster)
#' r <- raster(ncol = 36, nrow = 18)
#' r[] <- 1:ncell(r)
#' extractRaster(r)
extractRaster <- function(r, xy = NULL, ext = raster::extent(r), plot = TRUE){
  if(is.null(xy)){
    ndim = dim(r)
    x = 0:ndim[2] / ndim[2]
    y = rep(0.5, length(x))
    xy = cbind(x, y)
  }
  x = ext[1] + xy[, 1] * (ext[2] - ext[1] )
  y = ext[3] + xy[, 2] * (ext[4] - ext[3] )
  if(plot){
    raster::plot(r);
    points(x, y, col = 2)
    nx = length(x)
    points(x, y)
    graphics::arrows(x[1], y[1], x[nx], y[nx], lty = 3, lwd = 1.5, col = 2)
    # lines(x,y, lwd=1.5, col=2, lty=2)
  }
  v = raster::extract(r, cbind(x, y))
  ret = cbind('x' = x, 'y' = y, 'z' = v)
  return(ret)
}

#' Simplify SpatialData.
#' @param x SpatialData
#' @return Simplified SpatialData
#' @export
SimpleSpatial <- function(x){
  # n1 = length(x@polygons)
  # nj = unlist(lapply(1:n1, function(i){ length(x@polygons[[i]]@Polygons) } ))
  # x@polygons[[1]]@Polygons[[1]]@coords
  msg = 'SimpleSpatial'
  ni = length(x@polygons)
  k = 1
  sl = list()
  for(i in 1:ni){
    nj = length(x@polygons[[1]]@Polygons)
    for(j in 1:nj){
      cd = x@polygons[[i]]@Polygons[[j]]@coords
      np = nrow(cd)
      message(msg, i, '-', j, '\t', np)
      sl[[k]] = paste('POLYGON((', paste( paste(cd[, 1], cd[, 2]), collapse = ',' ), '))')
      if(k == 1){
        str = sl[[k]]
      }else{
        str = paste(str, ',', sl[[k]] )
      }
      k = k + 1
    }
  }
  r = rgeos::readWKT(paste('GEOMETRYCOLLECTION(', str, ')'))
}

#' Find the points in distance less than tol.
#' \code{PointInDistance}
#' @param pt 2-column coordinates (x,y).
#' @param tol Tolerance
#' @return Index of points that within tol
#' @export
PointInDistance <- function(pt, tol){
  msg = 'PointInDistance'
  dm = as.matrix(stats::dist(pt, diag  = TRUE))
  dm[dm == 0] = NA
  # View(dm)
  dmin = apply(dm, 1, min, na.rm = TRUE)
  id = which(dmin < 100)
  id1 = id2 = NULL
  tmp1 = tmp = dmin[id]
  i = 2
  for(i in 1:length(id)){
    if(id[i] %in% id2){next }
    id1 = c(id1, i); id1
    tmp1[id1] = NA; tmp1
    id[which(tmp1 %in% tmp[id1])]
    id2 = unique(c(id2,  id[which(tmp1 %in% tmp[id1])]))
    id2
    tmp1 = tmp
  }
  cbind('P1' = id[id1], P2 = id2)
}

#' Conver the MULTIPOLYGONS to SINGLEPOLYGONS.
#' @param x the spatialpolygon*
#' @param id Index of the sorted (decreasing) polygons to return. default = 0;
#' @export
SinglePolygon <- function(x, id = 0){
  n1 = length(x)
  y1 = list()
  y2 = list()
  k = 1
  for(i in 1:n1){
    message('level 1:', i, '/', n1)
    x1 = x@polygons[[i]]
    
    n2 = length(x1@Polygons)
    for(j in 1:n2){
      message('level 2:', j, '/', n2)
      x2 = x1@Polygons[[j]]
      y1[[k]] = sp::Polygons(list(x2), ID = k)
      k = k + 1
    }
  }
  
  y = sp::SpatialPolygonsDataFrame(sp::SpatialPolygons(y1), data = data.frame('ID' = 2:k - 1))
  
  if(id < 1){
    return(y)
  }else{
    ia = rgeos::gArea(y, byid = TRUE)
    id = order(ia, decreasing = TRUE)[1]
    return(y[id,])
    
  }
}

#' Remove the duplicated lines, which means the FROM and TO points are identical.
#' \code{rmDuplicatedLines}
#' @param x \code{sf} data or \code{sp} ShapeLine*
#' @param ... More options in duplicated()
#' @return ShapeLine* without duplicated lines.
#' @export
rmDuplicatedLines <- function(x, ...){
  # x = spi.riv
  cd = get_coords(x)
  # dim(cd)
  ft = get_from_to_nodes(x, cd)
  id = which(duplicated(ft[, -1], MARGIN = 1, ...))
  if(length(id) > 0){
    r = x[-id,]
  }else{
    r = x
  }
  return(r)
}

#' Generate Thiesson Polygons from a point data
#' \code{voronoipolygons}
#' @param x ShapePoints* in PCS
#' @param pts Coordinates (x,y) of the points
#' @param rw The coordinates of the corners of the rectangular window enclosing the triangulation, in the order (xmin, xmax, ymin, ymax).  (default = NULL)
#' @param crs projection parameters (default = NULL)
#' @return ShapePolygon*
#' @export
#' @examples 
#' library(rgeos)
#' library(rSHUD)
#' n = 10
#' xx = rnorm(n)
#' yy = rnorm(n)
#' str = paste('MULTIPOINT(', paste(paste('(', xx, yy, ')'), collapse = ','), ')')
#' x = readWKT(str)
#' vx = voronoipolygons(x)
#' raster::plot(vx, axes = TRUE)
#' raster::plot(add = TRUE, x, col = 2)
#' #' ====END====
#' 
#' x = 1:5
#' y = 1:5
#' xy = expand.grid(x, y)
#' vx = voronoipolygons(pts = xy)
#' plot(vx, axes = TRUE)
#' points(xy)
#' #' ====END====
#' 
#' library(sf)
#' library(raster)
#' library(rgeos)
#' n = 10
#' xx = rnorm(n)
#' yy = rnorm(n)
#' str = paste('MULTIPOINT(', paste(paste('(', xx, yy, ')'), collapse = ','), ')')
#' x = readWKT(str)
#' y = readWKT(paste('MULTIPOINT(', paste(paste('(', xx + 2, yy + 2, ')'), collapse = ','), ')'))
#' e1 = extent(y)
#' e2 = extent(x)
#' rw = c(min(e1[1], e2[1]),
#'        max(e1[2], e2[2]),
#'        min(e1[3], e2[3]),
#'        max(e1[4], e2[4]) ) + c(-1, 1, -1, 1)
#' vx = voronoipolygons(x = x, rw = rw) 
#' plot(vx); plot(add = TRUE, x, col = 2); plot(add = TRUE, y, col = 3)
voronoipolygons = function(x, pts = x@coords, rw = NULL, crs = NULL) {
  z = deldir::deldir(pts[, 1], pts[, 2], rw = rw)
  w = deldir::tile.list(z)
  polys = vector(mode = 'list', length = length(w))
  for (i in seq(along = polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = sp::Polygons(list(sp::Polygon(pcrds)), ID = as.character(i))
  }
  SP = sp::SpatialPolygons(polys)
  voronoi = sp::SpatialPolygonsDataFrame(SP,
                                         data = data.frame(x = pts[, 1],
                                                           y = pts[, 2], row.names = sapply(methods::slot(SP, 'polygons'),
                                                                                            function(x) methods::slot(x, 'ID'))))
  if(!is.null(crs)){
    raster::crs(voronoi) = crs;
  }
  return(voronoi)
}



#' Generate the coverage map for forcing sites.
#' \code{ForcingCoverage}
#' @param sp.meteoSite sf POINT object or SpatialPoints*
#' @param pcs Projected Coordinate System (crs object)
#' @param gcs Geographic Coordinate System (crs object)
#' @param dem DEM raster (SpatRaster or RasterLayer)
#' @param wbd watershed boundary (sf or SpatialPolygons*)
#' @param enlarge enlarge factor for the boundary.
#' @param filenames Output filenames associated with the forcing sites.
#' @return sf object
#' @export
ForcingCoverage <- function(sp.meteoSite = NULL, filenames = paste0(sp.meteoSite$ID, '.csv'),
                            pcs, gcs = sf::st_crs(4326), 
                            dem, wbd, enlarge = 10000
                            ){
  # Convert inputs to modern sf/terra formats if needed
  if (!inherits(dem, "SpatRaster")) {
    dem <- terra::rast(dem)
  }
  if (!inherits(wbd, "sf") && !inherits(wbd, "sfc")) {
    wbd <- sf::st_as_sf(wbd)
  }
  
  if(is.null(sp.meteoSite)){
    sp.meteoSite <- sf::st_centroid(sf::st_geometry(wbd))
    sp.meteoSite <- sf::st_sf(ID = 1:length(sp.meteoSite), geometry = sp.meteoSite)
    filenames <- paste0(sp.meteoSite$ID, '.csv')
  } else if (!inherits(sp.meteoSite, "sf") && !inherits(sp.meteoSite, "sfc")) {
    sp.meteoSite <- sf::st_as_sf(sp.meteoSite)
    if (is.null(sp.meteoSite$ID)) sp.meteoSite$ID <- 1:nrow(sp.meteoSite)
  }
  
  x.pcs <- sf::st_transform(sp.meteoSite, pcs)
  x.gcs <- sf::st_transform(sp.meteoSite, gcs)
  
  ll <- sf::st_coordinates(x.gcs)
  xy <- sf::st_coordinates(x.pcs)
  
  # Extract elevation using terra
  z_ext <- terra::extract(dem, terra::vect(x.pcs))
  z <- z_ext[, 2] # terra::extract returns ID in first col
  
  att <- data.frame(1:nrow(sp.meteoSite), ll, xy, z, filenames)
  colnames(att) <- c('ID', 'Lon', 'Lat', 'X', 'Y', 'Z', 'Filename')
  
  e1 <- terra::ext(wbd)
  e2 <- terra::ext(terra::vect(x.pcs))
  rw <- c(min(e1$xmin, e2$xmin), 
         max(e1$xmax, e2$xmax),
         min(e1$ymin, e2$ymin),
         max(e1$ymax, e2$ymax)) 
         
  if(is.null(enlarge)){
    enlarge <- min(diff(rw[1:2]), diff(rw[3:4])) * 0.02
  }
  rw <- rw + c(-1, 1, -1, 1) * enlarge
  
  if(nrow(x.pcs) < 2){
    # Legacy fallback: if fishnet is needed
    sp.forc <- fishnet(xx = rw[1:2], yy = rw[3:4], crs = pcs)
    sp.forc <- sf::st_as_sf(sp.forc)
  } else {
    # Voronoi
    box <- sf::st_polygon(list(rbind(
      c(rw[1], rw[3]), c(rw[2], rw[3]), 
      c(rw[2], rw[4]), c(rw[1], rw[4]), 
      c(rw[1], rw[3])
    )))
    box_sfc <- sf::st_sfc(box, crs = pcs)
    v <- sf::st_voronoi(sf::st_union(sf::st_geometry(x.pcs)), box_sfc)
    v <- sf::st_cast(v)
    v <- sf::st_intersection(v, box_sfc)
    v_sf <- sf::st_sf(geometry = v)
    v_sf <- sf::st_join(v_sf, x.pcs)
    sf::st_crs(v_sf) <- pcs
    sp.forc <- v_sf
  }
  
  att[is.na(att)] <- -9999
  
  # For backward compatibility with scripts expecting sp
  sp.forc_sp <- sf::as_Spatial(sp.forc)
  sp.forc_sp@data <- att
  return(sp.forc_sp)
}

#' Find the subset of a extent in a grid.
#' \code{grid.subset}
#' @param ext txtent of the grid
#' @param res resolution of the grid
#' @param ext.sub the extent of subset
#' @param x x coordinates of the grids
#' @param y y coordinates of the grids
#' Subset grid by extent (DEPRECATED)
#' 
#' @description
#' \strong{This function is deprecated.} Use \code{terra::crop()} or
#' \code{terra::ext()} for raster subsetting instead.
#' 
#' This function provides grid subsetting functionality that is better handled
#' by modern spatial packages like terra, which offer more robust and efficient
#' extent-based operations.
#' 
#' @param ext Numeric vector. Full extent c(xmin, xmax, ymin, ymax)
#' @param res Numeric. Resolution (single value or c(dx, dy))
#' @param ext.sub Numeric vector. Subset extent c(xmin, xmax, ymin, ymax)
#' @param dx Numeric. X resolution (default: first element of `res`)
#' @param dy Numeric. Y resolution (default: second element of `res`)
#' @param x Numeric vector. X coordinates (auto-generated from ext and dx)
#' @param y Numeric vector. Y coordinates (auto-generated from ext and dy)
#' @return List with components: xid, yid, x, y
#' @export
#' @section Deprecated:
#' \code{grid.subset()} is deprecated. Use \code{terra::crop()} or
#' \code{terra::ext()} instead.
#' 
#' @examples
#' \dontrun{
#' # Deprecated usage
#' result <- grid.subset(ext = c(0, 100, 0, 100), res = 1, 
#'                       ext.sub = c(20, 80, 20, 80))
#' 
#' # Recommended usage with terra
#' library(terra)
#' r <- rast(xmin = 0, xmax = 100, ymin = 0, ymax = 100, res = 1)
#' r_sub <- crop(r, ext(20, 80, 20, 80))
#' }
grid.subset <- function(ext, res,
                        ext.sub,
                        dx = matrix(res, 2, 1)[1],
                        dy = matrix(res, 2, 1)[2],
                        x = seq(ext[1] + dx / 2, ext[2] - dx / 2, by = dx),
                        y = seq(ext[3] + dy / 2, ext[4] - dy / 2, by = dy)
                        ){
  .Deprecated("terra::crop", package = "terra",
              msg = paste(
                "grid.subset() is deprecated.",
                "Use terra::crop() or terra::ext() instead."
              ))
  
  xmin <- min(x - dx / 2); xmax <- max(x + dx / 2)
  ymin <- min(y - dy / 2); ymax <- max(y + dy / 2)

  if(ext.sub[1] < xmin | ext.sub[2] > xmax | ext.sub[3] < ymin | ext.sub[4] > ymax){
    warning("Extent required is larger than the bounding box of dataset")
    message(paste(ext.sub, collapse = ","))
    message(paste(c(xmin, xmax, ymin, ymax), collapse = ","))
  }
  xid <- min(which(abs(x - ext.sub[1]) <= dx / 2)):max(which(abs(x - ext.sub[2]) <= dx / 2))
  yid <- min(which(abs(y - ext.sub[3]) <= dy / 2)):max(which(abs(y - ext.sub[4]) <= dy / 2))
  nx <- length(xid); ny <- length(yid)
  x.cord <- x[xid]; y.cord <- y[yid]
  rt <- list("xid" = xid,
            "yid" = yid,
            "x" = x.cord,
            "y" = y.cord)
  return(rt)
}


#' Generate a shapefile from coordinates.
#' \code{xy2shp}
#' @param xy matrix
#' @param df attribute table
#' @param crs projection parameters
#' @param shape Shape of the result in points, lines or polygons
#' @param output Output class: \code{"sp"} (legacy) or \code{"sf"} (modern)
#' @return Spatial*DataFrame objects when \code{output = "sp"}, or an \code{sf} object when \code{output = "sf"}
#' @export
#' @examples 
#' library(raster)
#' library(sf)
#' xy = list(cbind(c(0, 2, 1), c(0, 0, 2)),  cbind(c(0, 2, 1), c(0, 0, 2)) + 2)
#' 
#' # Legacy output (default): sp
#' sp1 = xy2shp(xy = xy, shape = 'polygon')
#' raster::plot(sp1, axes = TRUE, col = 'gray')
#' 
#' sp2 = xy2shp(xy = xy, shape = 'lines')
#' raster::plot(sp2, add = TRUE, lty = 2, lwd = 3, col = 'red')
#' sp3 = xy2shp(xy = xy, shape = 'POINTS')
#' raster::plot(sp3, add = TRUE, pch = 1, cex = 2)
#' 
#' # Modern output: sf
#' sf1 = xy2shp(xy = xy[[1]], shape = 'polygon', output = 'sf')
#' plot(sf::st_geometry(sf1), border = 'blue', add = TRUE)
#' 
xy2shp <- function(xy, df = NULL, crs = NULL, shape = 'points', output = c("sp", "sf")){
  output = match.arg(output)
  shape_l = tolower(shape)

  if (!is.list(xy)) {
    xy = list(xy)
  }

  crs_sf = NA_character_
  if (!is.null(crs)) {
    crs_sf = tryCatch(sf::st_crs(crs), error = function(e) sf::st_crs(as.character(crs)))
  }

  if (grepl("point", shape_l)) {
    xy_mat = do.call(rbind, xy)
    if (is.null(df)) {
      df = data.frame("ID" = seq_len(nrow(xy_mat)))
    }
    pt_geoms = lapply(seq_len(nrow(xy_mat)), function(i) sf::st_point(as.numeric(xy_mat[i, 1:2])))
    sf_obj = sf::st_sf(df, geometry = sf::st_sfc(pt_geoms, crs = crs_sf))
  } else if (grepl("line", shape_l)) {
    if (is.null(df)) {
      df = data.frame("ID" = seq_len(length(xy)))
    }
    line_geoms = lapply(xy, function(m) sf::st_linestring(as.matrix(m)))
    sf_obj = sf::st_sf(df, geometry = sf::st_sfc(line_geoms, crs = crs_sf))
  } else if (grepl("polygon", shape_l)) {
    if (is.null(df)) {
      df = data.frame("ID" = seq_len(length(xy)))
    }
    poly_geoms = lapply(xy, function(m) {
      mm = as.matrix(m)
      if (!(mm[1, 1] == mm[nrow(mm), 1] && mm[1, 2] == mm[nrow(mm), 2])) {
        mm = rbind(mm, mm[1, ])
      }
      sf::st_polygon(list(mm))
    })
    sf_obj = sf::st_sf(df, geometry = sf::st_sfc(poly_geoms, crs = crs_sf))
  } else {
    stop("shape must be points, lines, or polygons")
  }

  if (output == "sf") {
    return(sf_obj)
  }

  spd = methods::as(sf_obj, "Spatial")
  if(!is.null(crs)){
    raster::crs(spd) = crs
  }
  return(spd)
}

