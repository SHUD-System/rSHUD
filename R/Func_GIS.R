#' Write ESRI shapefile out
#' \code{writeshape}
#' @param shp \code{sf} object. Legacy \code{sp} \code{Spatial*} objects are
#'   accepted only for deprecated compatibility and are converted with
#'   \code{sf::st_as_sf()} before writing.
#' @param crs CRS to assign before write (default: CRS taken from \code{shp} via
#'   \code{sf::st_crs()}). Use
#'   \code{NA} to skip \code{st_set_crs} (not recommended for new data).
#' @param file file path, with or without \code{.shp} extension.
#' @details
#' All outputs are written with \code{sf::st_write()}. Deprecated
#' compatibility inputs without a data slot are given a synthetic \code{ID}
#' column so conversion to \code{sf} is well-defined.
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
      sf::st_crs(sf::st_as_sf(shp))
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
    shp <- tryCatch(
      {
        if ((inherits(shp, "SpatialPolygons") &&
             !inherits(shp, "SpatialPolygonsDataFrame")) ||
            (inherits(shp, "SpatialLines") &&
             !inherits(shp, "SpatialLinesDataFrame")) ||
            (inherits(shp, "SpatialPoints") &&
             !inherits(shp, "SpatialPointsDataFrame"))) {
          geom <- sf::st_geometry(sf::st_as_sf(shp))
          sf::st_sf(
            ID = seq_len(length(shp)),
            geometry = geom
          )
        } else {
          sf::st_as_sf(shp)
        }
      },
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
  proj_str <- if (inherits(proj4string, "crs")) {
    proj4string$proj4string
  } else if (inherits(proj4string, "CRS")) {
    sf::st_crs(proj4string)$proj4string
  } else {
    as.character(proj4string)
  }
  y <- proj4::project(x, proj_str, inverse = P2G)
  if(P2G){
    colnames(y) = c('Lon', 'Lat')
  }else{
    colnames(y) = c('X', 'Y')
  }
  y
}

# sp2raster moved to R/gis_core.R as deprecated wrapper for vector_to_raster()

#' Generate the terra raster mask of a mesh domain
#' \code{shud.mask}
#' @param pm \code{shud.mesh}
#' @param n Number of grid
#' @param rr Default mask in .shud environment
#' @param cellsize Resolution, defaul = NULL, resolution = extent / ngrids
#' @param proj Projection parameter
#' @return \code{terra::SpatRaster} mask
#' @export
shud.mask  <- function (pm = read_mesh(), proj = NULL,
                        rr = NULL,
                        n = 10000, cellsize = NULL){
  if (is.null(rr)) {
    rr <- tryCatch(get('MASK', envir = .shud), error = function(e) NULL)
  }
  if (!is.null(rr) && !inherits(rr, "SpatRaster")) {
    rr <- terra::rast(rr)
  }
  if(is.null(rr)){
    spm <- mesh_to_sf(pm)
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
    assign('MASK', rr, envir = .shud)
  }
  rr
}

.mesh_data_to_raster_legacy <- function(x = getElevation(),
                                        rmask = shud.mask(proj = proj),
                                        pm = read_mesh(), proj = NULL,
                                        stack = FALSE, method = 'ide',
                                        plot = FALSE){
  if (!inherits(rmask, "SpatRaster")) {
    rmask <- terra::rast(rmask)
  }

  if(stack){
    ret <- terra::rast(lapply(seq_len(nrow(x)), function(i) {
      .mesh_data_to_raster_legacy(
        x = as.numeric(x[i, ]),
        rmask = rmask,
        pm = pm,
        proj = proj,
        stack = FALSE,
        method = method,
        plot = FALSE
      )
    }))
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
      val = data.frame(X = xy[, 1], Y = xy[, 2], Z = x)
      gs <- gstat::gstat(
        formula = Z ~ 1,
        locations = ~ X + Y,
        data = val,
        set = list(idp = 2.0)
      )
      dat <- terra::interpolate(rmask, gs)
      r <- dat[["var1.pred"]]
    }
    if(grepl('linear', tolower(method))){
      grid_xy <- terra::crds(rmask, na.rm = FALSE)
      xx <- interp::interp(
        x = xy[, 1],
        y = xy[, 2],
        z = x,
        xo = sort(unique(grid_xy[, 1])),
        yo = sort(unique(grid_xy[, 2]))
      )
      r <- rmask
      terra::values(r) <- as.numeric(xx$z)
    }
    if(grepl('ide', tolower(method))){
      tps <- fields::Tps(xy, x)
      r <- terra::interpolate(rmask, tps)
    }
    ret <- terra::mask(r, rmask)
  }
  if(plot){
    terra::plot(ret)
  }
  if(!is.null(proj)){ terra::crs(ret) <- proj }
  return(ret)
}

#' Deprecated: convert mesh data to raster
#' \code{MeshData2Raster}
#' @param ... Arguments passed to \code{\link{mesh_to_raster}}. Deprecated
#'   legacy names are mapped before forwarding: \code{x} to \code{data},
#'   \code{rmask} to \code{template}, \code{pm} to \code{mesh}, and
#'   \code{proj} to \code{crs}. The legacy \code{plot} argument is handled by
#'   plotting the returned \code{terra::SpatRaster}; it is not forwarded.
#' @return \code{terra::SpatRaster}
#' @keywords deprecated
#' @export
MeshData2Raster <- function(...){
  .Deprecated("mesh_to_raster",
              msg = "MeshData2Raster() is deprecated. Use mesh_to_raster() for modern terra workflows.")
  args <- list(...)
  has_named_arg <- function(x, name) {
    arg_names <- names(x)
    !is.null(arg_names) && name %in% arg_names
  }
  arg_value <- function(x, name) {
    x[[name, exact = TRUE]]
  }
  has_first_unnamed_arg <- function(x) {
    arg_names <- names(x)
    length(x) > 0 && (is.null(arg_names) || !nzchar(arg_names[[1]]))
  }

  if (has_named_arg(args, "x")) {
    if (!has_named_arg(args, "data") && !has_first_unnamed_arg(args)) {
      args$data <- args$x
    }
    args$x <- NULL
  }
  if (has_named_arg(args, "rmask")) {
    if (!has_named_arg(args, "template")) {
      args$template <- args$rmask
    }
    args$rmask <- NULL
  }
  if (has_named_arg(args, "pm")) {
    if (!has_named_arg(args, "mesh")) {
      args$mesh <- args$pm
    }
    args$pm <- NULL
  }
  if (has_named_arg(args, "proj")) {
    if (!has_named_arg(args, "crs")) {
      args$crs <- args$proj
    }
    args$proj <- NULL
  }

  plot_result <- isTRUE(args$plot)
  args$plot <- NULL

  if (!has_named_arg(args, "data") && !has_first_unnamed_arg(args)) {
    args$data <- getElevation()
  }
  if (!has_named_arg(args, "mesh") || is.null(arg_value(args, "mesh"))) {
    args$mesh <- read_mesh()
  }
  if (!has_named_arg(args, "template")) {
    args$template <- shud.mask(
      pm = arg_value(args, "mesh"),
      proj = arg_value(args, "crs")
    )
  }

  ret <- do.call(mesh_to_raster, args)
  if (plot_result) {
    terra::plot(ret)
  }
  ret
}


#' Generatue fishnet
#' \code{fishnet} Generate fishnet by the coordinates
#' @param xx  x coordinates
#' @param yy  y coordinates
#' @param crs projections parameters, defaul = epsg:4326
#' @param type option = 'polygon', 'points', 'line'
#' @return Legacy \code{Spatial*} object for compatibility. Prefer
#'   \code{xy2shp(..., output = "sf")} or direct \code{sf} constructors for new
#'   code.
#' @export
#' @examples
#' ext = c(0, 8, 2, 10)
#' dx = 2; dy = 4
#' xx = seq(ext[1], ext[2], by = dx)
#' yy = seq(ext[3], ext[4], by = dy)
#' sp1 = fishnet(xx = ext[1:2], yy = ext[3:4])
#' sp2 = fishnet(xx = xx + .5 * dx, yy = yy + 0.5 * dy)
#' sp3 = fishnet(xx = xx, yy = yy, type = 'point')
#' plot(sf::st_geometry(sf::st_as_sf(sp1)), axes = TRUE,
#'      xlim = c(-1, 1) * dx + ext[1:2],
#'      ylim = c(-1, 1) * dy + ext[3:4])
#' plot(sf::st_geometry(sf::st_as_sf(sp2)), axes = TRUE, add = TRUE, border = 2)
#' plot(sf::st_geometry(sf::st_as_sf(sp3)), axes = TRUE, add = TRUE, col = 3, pch = 20)
#' grid()
fishnet <- function(xx, yy,
                    crs = sf::st_crs(4326),
                    type = 'polygon'){
  crs_sf <- sf::st_crs(crs)
  nx <- length(xx)
  ny <- length(yy)
  if(grepl('line', tolower(type)) ){
    ymin <- min(yy)
    ymax <- max(yy)
    xmin <- min(xx)
    xmax <- max(xx)
    vline <- cbind(xx, ymin, xx, ymax)
    hline <- cbind(xmin, yy, xmax, yy)
    mat <- rbind(vline, hline)
    df <- as.data.frame(mat)
    colnames(df) <- c('x1', 'y1', 'x2', 'y2')
    geoms <- lapply(seq_len(nrow(mat)), function(i) {
      sf::st_linestring(matrix(mat[i, ], ncol = 2, byrow = TRUE))
    })
    ret <- sf::st_sf(df, geometry = sf::st_sfc(geoms, crs = crs_sf))
  }
  if(grepl('polygon', tolower(type)) ){
    polys <- vector("list", (nx - 1) * (ny - 1))
    df <- data.frame(
      xmin = numeric(length(polys)),
      xmax = numeric(length(polys)),
      ymin = numeric(length(polys)),
      ymax = numeric(length(polys)),
      xcenter = numeric(length(polys)),
      ycenter = numeric(length(polys))
    )
    k <- 1
    for (i in seq_len(nx - 1)) {
      for (j in seq_len(ny - 1)) {
        pxy <- rbind(
          c(xx[i], yy[j]),
          c(xx[i], yy[j + 1]),
          c(xx[i + 1], yy[j + 1]),
          c(xx[i + 1], yy[j]),
          c(xx[i], yy[j])
        )
        polys[[k]] <- sf::st_polygon(list(pxy))
        df[k, ] <- c(
          min(pxy[, 1]), max(pxy[, 1]),
          min(pxy[, 2]), max(pxy[, 2]),
          mean(range(pxy[, 1])), mean(range(pxy[, 2]))
        )
        k <- k + 1
      }
    }
    ret <- sf::st_sf(df, geometry = sf::st_sfc(polys, crs = crs_sf))
  }
  if(grepl('point', tolower(type)) ){
    xm <- expand.grid(xx, yy)
    df <- data.frame('X' = xm[, 1], 'Y' = xm[, 2])
    ret <- sf::st_as_sf(df, coords = c("X", "Y"), crs = crs_sf, remove = FALSE)
  }
  return(methods::as(ret, "Spatial"))
}

#' Add holes into polygons
#' \code{AddHoleToPolygon}
#' @param poly \code{sf} polygon object; legacy \code{SpatialPolygons} is
#'   accepted for compatibility.
#' @param hole Hole polygon as \code{sf}; legacy spatial input is accepted for
#'   compatibility.
#' @export
AddHoleToPolygon <- function(poly, hole){
  input_is_sf <- inherits(poly, "sf")
  poly_sf <- if (input_is_sf) poly else sf::st_as_sf(poly)
  hole_sf <- if (inherits(hole, "sf")) hole else sf::st_as_sf(hole)
  new <- suppressWarnings(sf::st_difference(poly_sf, sf::st_union(hole_sf)))
  if (input_is_sf) {
    return(new)
  }
  methods::as(new, "Spatial")
}
#' Deprecated: cut line geometries with a length threshold
#' \code{sp.CutSptialLines}
#' @param sl \code{sf} LINESTRING object; legacy \code{SpatialLines*} is
#'   accepted for compatibility.
#' @param tol Tolerence. If the length of segment is larger than tolerance, cut the segment until the maximum segment is shorter than tolerance.
#' @keywords deprecated
#' @export
#' @examples
#' \dontrun{
#' # Deprecated compatibility wrapper. Prefer sf-native preprocessing such as
#' # sf::st_segmentize() in new workflows.
#' line <- sf::st_sf(
#'   geometry = sf::st_sfc(sf::st_linestring(cbind(1:10, sin(1:10))))
#' )
#' cut_line <- sp.CutSptialLines(line, tol = 2)
#' }
sp.CutSptialLines <- function(sl, tol){
  msg = 'sp.CutSptialLines::'
  input_is_sp <- inherits(sl, c("Spatial", "SpatialLines", "SpatialLinesDataFrame"))
  sl_sf <- if (inherits(sl, "sf")) sl else sf::st_as_sf(sl)
  if (any(sf::st_geometry_type(sl_sf) == "MULTILINESTRING")) {
    sl_sf <- suppressWarnings(sf::st_cast(sl_sf, "LINESTRING"))
  }
  ll <- as.numeric(sf::st_length(sl_sf))
  if(all(ll < tol) ){
    ret <- sl
  }else{
    xsl <- list()
    ilen <- numeric()
    ik <- 1
    for(i in seq_len(nrow(sl_sf))){
      sx <- sl_sf[i, ]
      pxy <- sf::st_coordinates(sx)[, 1:2, drop = FALSE]
      np <- nrow(pxy)
      seglen <- sqrt(rowSums((pxy[-1, , drop = FALSE] - pxy[-np, , drop = FALSE]) ^ 2))
      dacc <- cumsum(seglen)
      tol_i <- max(c(tol, min(dacc)))
      len <- sum(seglen)
      if(len > tol_i){
        nsplit <- ceiling(len / tol_i)
      }else{
        nsplit <- 1
      }
      dd <- len / nsplit
      v0 <- 1
      message(msg, i, '/', nrow(sl_sf), '\t', nsplit, '\t', round(dd, 2) )
      for(k in seq_len(nsplit)){
        if(v0 >= np){
          break
        }
        dk <- dd * k
        v1 <- order(abs(dacc - dk), decreasing = FALSE)[1] + 1
        if(v1 + 1 > np){
          v1 <- np
        }
        message(msg, v0, '\t', v1)
        if(v0 == v1){
          next
        }
        coords_k <- pxy[v0:v1, , drop = FALSE]
        xsl[[ik]] <- sf::st_linestring(as.matrix(coords_k))
        ilen[ik] <- sum(sqrt(rowSums((coords_k[-1, , drop = FALSE] - coords_k[-nrow(coords_k), , drop = FALSE]) ^ 2)))
        ik <- ik + 1
        v0 <- v1
      }
    }
    ret_sf <- sf::st_sf(
      INDEX = seq_along(xsl),
      Length = ilen,
      geometry = sf::st_sfc(xsl, crs = sf::st_crs(sl_sf))
    )
    ret <- if (input_is_sp) methods::as(ret_sf, "Spatial") else ret_sf
  }
  return(ret)
}

#' Extract values on a terra raster along a normalized line from 0 to 1.
#' \code{extractRaster}
#' @param r \code{terra::SpatRaster}
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
#' r <- terra::rast(ncols = 36, nrows = 18)
#' terra::values(r) <- seq_len(terra::ncell(r))
#' extractRaster(r)
extractRaster <- function(r, xy = NULL, ext = NULL, plot = TRUE){
  rr <- if (inherits(r, "SpatRaster")) r else terra::rast(r)
  if (is.null(ext)) {
    ext <- terra::ext(rr)
  }
  if (inherits(ext, "SpatExtent")) {
    ext <- c(ext$xmin, ext$xmax, ext$ymin, ext$ymax)
  }
  if(is.null(xy)){
    ndim = dim(rr)
    x = 0:ndim[2] / ndim[2]
    y = rep(0.5, length(x))
    xy = cbind(x, y)
  }
  x = ext[1] + xy[, 1] * (ext[2] - ext[1] )
  y = ext[3] + xy[, 2] * (ext[4] - ext[3] )
  if(plot){
    terra::plot(rr)
    points(x, y, col = 2)
    nx = length(x)
    points(x, y)
    graphics::arrows(x[1], y[1], x[nx], y[nx], lty = 3, lwd = 1.5, col = 2)
    # lines(x,y, lwd=1.5, col=2, lty=2)
  }
  ev = terra::extract(rr, cbind(x, y))
  v = if ("ID" %in% names(ev)) ev[[2]] else ev[[1]]
  ret = cbind('x' = x, 'y' = y, 'z' = v)
  return(ret)
}

#' Simplify polygon geometry.
#' @param x \code{sf} polygon object; legacy spatial input is accepted for
#'   compatibility.
#' @return Simplified \code{sf} geometry collection
#' @export
SimpleSpatial <- function(x){
  msg = 'SimpleSpatial'
  x_sf <- if (inherits(x, "sf")) x else sf::st_as_sf(x)
  geoms <- suppressWarnings(sf::st_cast(sf::st_geometry(x_sf), "POLYGON"))
  for (i in seq_along(geoms)) {
    np <- nrow(sf::st_coordinates(geoms[[i]]))
    message(msg, i, '\t', np)
  }
  sf::st_combine(geoms)
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

#' Convert MULTIPOLYGON geometry to single POLYGON features.
#' @param x \code{sf} polygon object; legacy spatial input is accepted for
#'   compatibility.
#' @param id Index of the sorted (decreasing) polygons to return. default = 0;
#' @export
SinglePolygon <- function(x, id = 0){
  input_is_sp <- inherits(x, c("Spatial", "SpatialPolygons", "SpatialPolygonsDataFrame"))
  x_sf <- if (inherits(x, "sf")) x else sf::st_as_sf(x)
  geoms <- suppressWarnings(sf::st_cast(sf::st_geometry(x_sf), "POLYGON"))
  y_sf <- sf::st_sf(ID = seq_along(geoms), geometry = geoms)

  if(id < 1){
    return(if (input_is_sp) methods::as(y_sf, "Spatial") else y_sf)
  }

  ia <- as.numeric(sf::st_area(y_sf))
  idx <- order(ia, decreasing = TRUE)[min(id, length(ia))]
  ret <- y_sf[idx, ]
  if (input_is_sp) {
    return(methods::as(ret, "Spatial"))
  }
  ret
}

#' Remove incomplete lines and duplicated lines with identical FROM/TO nodes.
#' \code{rmDuplicatedLines}
#' @param x \code{sf} line object; legacy spatial input is accepted for
#'   compatibility.
#' @param ... More options in duplicated()
#' @return Line object without duplicated from/to nodes, preserving the input
#'   class where possible. Lines whose FROM/TO nodes cannot be determined, such
#'   as empty, one-coordinate, non-finite, or same-start/end LINESTRING
#'   geometries, are removed before duplicate FROM/TO pairs are filtered.
#' @export
rmDuplicatedLines <- function(x, ...){
  # x = spi.riv
  cd = get_coords(x)
  # dim(cd)
  ft = get_from_to_nodes(x, cd)
  valid = stats::complete.cases(ft[, c("FrNode", "ToNode"), drop = FALSE])
  keep = valid
  if (any(valid)) {
    dup = duplicated(ft[valid, c("FrNode", "ToNode"), drop = FALSE], MARGIN = 1, ...)
    keep[which(valid)[dup]] = FALSE
  }
  idx = which(keep)
  if (inherits(x, "sfc")) {
    r = x[idx]
  } else {
    r = x[idx,]
  }
  return(r)
}

#' Generate Thiessen polygons from point data
#' \code{voronoipolygons}
#' @param x \code{sf} point object in projected coordinates
#' @param pts Coordinates (x,y) of the points
#' @param rw The coordinates of the corners of the rectangular window enclosing the triangulation, in the order (xmin, xmax, ymin, ymax).  (default = NULL)
#' @param crs projection parameters (default = NULL)
#' @return Legacy \code{SpatialPolygons*} object for compatibility.
#' @export
#' @examples 
#' library(rSHUD)
#' set.seed(1)
#' n = 10
#' xx = rnorm(n)
#' yy = rnorm(n)
#' x = sf::st_as_sf(data.frame(x = xx, y = yy), coords = c("x", "y"))
#' vx = voronoipolygons(x)
#' plot(sf::st_geometry(sf::st_as_sf(vx)), axes = TRUE)
#' plot(sf::st_geometry(x), add = TRUE, col = 2)
#' #' ====END====
#' 
#' x = 1:5
#' y = 1:5
#' xy = expand.grid(x, y)
#' vx = voronoipolygons(pts = xy)
#' plot(sf::st_geometry(sf::st_as_sf(vx)), axes = TRUE)
#' points(xy)
#' #' ====END====
#' 
#' library(sf)
#' set.seed(2)
#' n = 10
#' xx = rnorm(n)
#' yy = rnorm(n)
#' x = sf::st_as_sf(data.frame(x = xx, y = yy), coords = c("x", "y"))
#' y = sf::st_as_sf(data.frame(x = xx + 2, y = yy + 2), coords = c("x", "y"))
#' e1 = sf::st_bbox(y)
#' e2 = sf::st_bbox(x)
#' rw = c(xmin = min(e1["xmin"], e2["xmin"]),
#'        xmax = max(e1["xmax"], e2["xmax"]),
#'        ymin = min(e1["ymin"], e2["ymin"]),
#'        ymax = max(e1["ymax"], e2["ymax"])) + c(xmin = -1, xmax = 1, ymin = -1, ymax = 1)
#' vx = voronoipolygons(x = x, rw = rw) 
#' plot(sf::st_geometry(sf::st_as_sf(vx)))
#' plot(sf::st_geometry(x), add = TRUE, col = 2)
#' plot(sf::st_geometry(y), add = TRUE, col = 3)
voronoipolygons = function(x = NULL, pts = NULL, rw = NULL, crs = NULL) {
  if (is.null(pts)) {
    if (is.null(x)) {
      stop("Either 'x' or 'pts' must be provided", call. = FALSE)
    }
    if (inherits(x, c("sf", "sfc"))) {
      pts <- sf::st_coordinates(sf::st_geometry(x))[, 1:2, drop = FALSE]
    } else {
      pts <- get_coords(x)
    }
  }
  z = deldir::deldir(pts[, 1], pts[, 2], rw = rw)
  w = deldir::tile.list(z)
  polys = vector(mode = 'list', length = length(w))
  for (i in seq(along = polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = sf::st_polygon(list(pcrds))
  }
  voronoi = sf::st_sf(
    data.frame(x = pts[, 1], y = pts[, 2]),
    geometry = sf::st_sfc(polys, crs = sf::st_crs(crs))
  )
  return(methods::as(voronoi, "Spatial"))
}



#' Generate the coverage map for forcing sites.
#' \code{ForcingCoverage}
#' @param sp.meteoSite \code{sf} POINT or POLYGON object; legacy
#'   \code{SpatialPoints*}/\code{SpatialPolygons*} input is accepted for
#'   compatibility. Polygon inputs are returned as one forcing zone per polygon,
#'   with representative point coordinates used for site metadata. This supports
#'   real AutoSHUD/DataPre polygon `meteoCov.shp` forcing-coverage layers as well
#'   as point station inputs.
#' @param pcs Projected Coordinate System (crs object)
#' @param gcs Geographic Coordinate System (crs object)
#' @param dem DEM as a \code{terra::SpatRaster}
#' @param wbd watershed boundary as \code{sf}; legacy \code{SpatialPolygons*}
#'   input is accepted for compatibility.
#' @param enlarge enlarge factor for the boundary.
#' @param filenames Output filenames associated with the forcing sites.
#' @return Legacy \code{SpatialPolygonsDataFrame} forcing coverage object.
#' @export
ForcingCoverage <- function(sp.meteoSite = NULL, filenames = paste0(sp.meteoSite$ID, '.csv'),
                            pcs, gcs = sf::st_crs(4326),
                            dem, wbd, enlarge = 10000
                            ){
  # Compatibility conversion for older forcing-site workflows.
  if (!inherits(dem, "SpatRaster")) {
    dem <- terra::rast(dem)
  }
  if (!inherits(wbd, "sf") && !inherits(wbd, "sfc")) {
    wbd <- sf::st_as_sf(wbd)
  }

  if(is.null(sp.meteoSite)){
    sp.meteoSite <- sf::st_centroid(sf::st_geometry(wbd))
    sp.meteoSite <- sf::st_sf(ID = 1:length(sp.meteoSite), geometry = sp.meteoSite)
    if (missing(filenames) || is.null(filenames) || length(filenames) == 0) {
      filenames <- paste0(sp.meteoSite$ID, '.csv')
    }
  } else if (inherits(sp.meteoSite, "sfc")) {
    sp.meteoSite <- sf::st_sf(geometry = sp.meteoSite)
    sp.meteoSite$ID <- seq_len(nrow(sp.meteoSite))
  } else if (!inherits(sp.meteoSite, "sf") && !inherits(sp.meteoSite, "sfc")) {
    sp.meteoSite <- sf::st_as_sf(sp.meteoSite)
  }
  if (is.null(sp.meteoSite$ID)) sp.meteoSite$ID <- seq_len(nrow(sp.meteoSite))
  if (missing(filenames) || is.null(filenames) || length(filenames) == 0) {
    filenames <- paste0(sp.meteoSite$ID, '.csv')
  }
  if (length(filenames) != nrow(sp.meteoSite)) {
    stop(
      sprintf(
        "filenames must contain exactly one value per sp.meteoSite feature; got %d for %d features",
        length(filenames), nrow(sp.meteoSite)
      ),
      call. = FALSE
    )
  }

  x.pcs <- sf::st_transform(sp.meteoSite, pcs)
  geom_types <- unique(as.character(sf::st_geometry_type(x.pcs)))
  is_point_site <- all(geom_types %in% "POINT")
  is_polygon_site <- all(geom_types %in% c("POLYGON", "MULTIPOLYGON"))
  if (!is_point_site && !is_polygon_site) {
    stop(
      paste(
        "sp.meteoSite must contain POINT sites or POLYGON/MULTIPOLYGON coverage;",
        "cast or split MULTIPOINT inputs first"
      ),
      call. = FALSE
    )
  }

  site_points.pcs <- if (is_polygon_site) {
    sf::st_point_on_surface(sf::st_geometry(x.pcs))
  } else {
    sf::st_geometry(x.pcs)
  }
  site_points.gcs <- sf::st_transform(site_points.pcs, gcs)

  ll <- sf::st_coordinates(site_points.gcs)[, 1:2, drop = FALSE]
  xy <- sf::st_coordinates(site_points.pcs)[, 1:2, drop = FALSE]

  # Extract elevation using terra
  z_ext <- terra::extract(dem, terra::vect(sf::st_as_sf(site_points.pcs)))
  z <- z_ext[, 2] # terra::extract returns ID in first col

  att <- data.frame(1:nrow(sp.meteoSite), ll, xy, z, filenames)
  colnames(att) <- c('ID', 'Lon', 'Lat', 'X', 'Y', 'Z', 'Filename')
  
  e1 <- terra::ext(wbd)
  e2 <- terra::ext(terra::vect(sf::st_as_sf(site_points.pcs)))
  rw <- c(min(e1$xmin, e2$xmin),
         max(e1$xmax, e2$xmax),
         min(e1$ymin, e2$ymin),
         max(e1$ymax, e2$ymax))

  if(is.null(enlarge)){
    enlarge <- min(diff(rw[1:2]), diff(rw[3:4])) * 0.02
  }
  rw <- rw + c(-1, 1, -1, 1) * enlarge

  if (is_polygon_site) {
    sp.forc <- x.pcs
  } else if(nrow(x.pcs) < 2){
    # Single-site fallback still uses the compatibility fishnet helper.
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
  
  # Return legacy sp object for callers that still read @data from forcing zones.
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
#' @param output Output class: \code{"sf"} (modern) or \code{"sp"} (legacy).
#' @return \code{sf} object when \code{output = "sf"}, or legacy
#'   \code{Spatial*DataFrame} when \code{output = "sp"}.
#' @export
#' @examples 
#' library(sf)
#' xy = list(cbind(c(0, 2, 1), c(0, 0, 2)),  cbind(c(0, 2, 1), c(0, 0, 2)) + 2)
#' 
#' # Modern output: sf
#' sf1 = xy2shp(xy = xy[[1]], shape = 'polygon', output = 'sf')
#' plot(sf::st_geometry(sf1), border = 'blue')
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

  methods::as(sf_obj, "Spatial")
}
