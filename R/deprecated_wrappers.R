#' Deprecated wrapper for `FromToNode()`
#'
#' Legacy compatibility wrapper for \code{\link{get_from_to_nodes}}.
#' @export
FromToNode <- function(...) {
  .Deprecated("get_from_to_nodes")
  get_from_to_nodes(...)
}

#' Deprecated wrapper for `extractCoords()`
#'
#' Legacy compatibility wrapper for \code{\link{get_coords}}.
#' @export
extractCoords <- function(...) {
  .Deprecated("get_coords")
  get_coords(...)
}

#' Deprecated wrapper for `hydrograph()`
#'
#' Legacy compatibility wrapper for \code{\link{plot_hydrograph}}.
#' @export
hydrograph <- function(...) {
  .Deprecated("plot_hydrograph")
  plot_hydrograph(...)
}

#' Deprecated wrapper for `write.riv()`
#'
#' Legacy compatibility wrapper for \code{\link{write_river}}.
#' @export
write.riv <- function(riv, file) {
  .Deprecated("write_river")
  write_river(riv, file)
}

#' Deprecated wrapper for `sp.RiverSeg()`
#'
#' Legacy compatibility wrapper for \code{\link{shud.rivseg}}.
#' @export
sp.RiverSeg <- function(...) {
  .Deprecated("shud.rivseg")
  shud.rivseg(...)
}

#' Deprecated stub for `compareMaps()`
#'
#' This legacy export is retained for transition compatibility only.
#' @export
compareMaps <- function(...) {
  .Deprecated()
  NULL
}

#' Deprecated stub for `map2d()`
#'
#' This legacy export is retained for transition compatibility only.
#' @export
map2d <- function(...) {
  .Deprecated()
  NULL
}

#' Deprecated wrapper for `plot_sp()`
#'
#' Uses modern plotting methods when possible.
#' @export
plot_sp <- function(x, field = NULL, zcol = field, ...) {
  .Deprecated("plot_polygons")

  if (is.null(field) && !is.null(zcol)) {
    field <- zcol
  }

  if (!is.null(field) && inherits(x, c("sf", "SpatVector")) && field %in% names(x)) {
    plot_polygons(x, field = field, ...)
    return(invisible(x))
  }

  if (inherits(x, "SpatVector")) {
    terra::plot(x, ...)
    return(invisible(x))
  }

  if (inherits(x, "sf")) {
    graphics::plot(sf::st_geometry(x), ...)
    return(invisible(x))
  }

  if (inherits(x, "sfc")) {
    graphics::plot(x, ...)
    return(invisible(x))
  }

  graphics::plot(x, ...)
  invisible(x)
}

#' Deprecated wrapper for `plot_tsd()`
#'
#' Uses \code{plot_timeseries()} when available, otherwise falls back to base plotting.
#' @export
plot_tsd <- function(x, ...) {
  .Deprecated("plot_timeseries")

  if (exists("plot_timeseries", mode = "function", inherits = TRUE)) {
    return(plot_timeseries(x, ...))
  }

  graphics::plot(x, ...)
  invisible(x)
}

#' Deprecated wrapper for `png.control()`
#'
#' Legacy compatibility wrapper for \code{grDevices::png()}.
#' @export
png.control <- function(...) {
  .Deprecated(msg = "png.control() is deprecated. Use grDevices::png() instead.")
  grDevices::png(...)
}

#' Deprecated wrapper for `removeholes()`
#'
#' Uses modern terra/sf polygon handling to fill interior holes.
#' @export
removeholes <- function(sp, ...) {
  .Deprecated(msg = "removeholes() is deprecated. Use terra::fillHoles() or sf workflows instead.")

  is_legacy <- inherits(sp, "Spatial")

  if (is_legacy) {
    filled <- terra::fillHoles(terra::vect(sp), ...)
    return(as(sf::st_as_sf(filled), "Spatial"))
  }

  if (inherits(sp, "SpatVector")) {
    return(terra::fillHoles(sp, ...))
  }

  if (inherits(sp, "sf")) {
    filled <- terra::fillHoles(terra::vect(sp), ...)
    return(sf::st_as_sf(filled))
  }

  warning("removeholes() only supports Spatial, sf, or SpatVector inputs.", call. = FALSE)
  NULL
}

#' Deprecated wrapper for `write.mesh()`
#'
#' Legacy compatibility wrapper for \code{\link{write_mesh}}.
#' @export
write.mesh <- function(pm, file) {
  .Deprecated("write_mesh")
  write_mesh(pm, file)
}

#' Deprecated wrapper for `write.ic()`
#'
#' Legacy compatibility wrapper for \code{\link{write_ic}}.
#' @export
write.ic <- function(x, file) {
  .Deprecated("write_ic")
  write_ic(x, file)
}

#' Deprecated wrapper for `write.df()`
#'
#' Legacy compatibility wrapper for \code{\link{write_df}}.
#' @export
write.df <- function(x, file, append = FALSE, quiet = FALSE, header = NULL) {
  .Deprecated("write_df")
  write_df(x, file, append = append, quiet = quiet, header = header)
}

#' Deprecated wrapper for `write.config()`
#'
#' Legacy compatibility wrapper for \code{\link{write_config}}.
#' @export
write.config <- function(x, file) {
  .Deprecated("write_config")
  write_config(x, file)
}

#' Deprecated wrapper for `write.forc()`
#'
#' Legacy compatibility wrapper for \code{\link{write_forc}}.
#' @export
write.forc <- function(x, file, path = "", startdate = "20000101") {
  .Deprecated("write_forc")
  write_forc(x, file, path = path, startdate = startdate)
}
