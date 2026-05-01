#' Deprecated wrapper for `FromToNode()`
#'
#' Legacy compatibility wrapper for \code{\link{get_from_to_nodes}}.
#'
#' @param ... Arguments passed to \code{\link{get_from_to_nodes}}.
#' @return The result from \code{\link{get_from_to_nodes}}.
#' @export
FromToNode <- function(...) {
  .Deprecated("get_from_to_nodes")
  get_from_to_nodes(...)
}

#' Deprecated wrapper for `extractCoords()`
#'
#' Legacy compatibility wrapper for \code{\link{get_coords}}.
#'
#' @param ... Arguments passed to \code{\link{get_coords}}.
#' @return The result from \code{\link{get_coords}}.
#' @export
extractCoords <- function(...) {
  .Deprecated("get_coords")
  get_coords(...)
}

#' Deprecated wrapper for `hydrograph()`
#'
#' Legacy compatibility wrapper for \code{\link{plot_hydrograph}}.
#'
#' @param ... Arguments passed to \code{\link{plot_hydrograph}}.
#' @return A plot object from \code{\link{plot_hydrograph}}.
#' @export
hydrograph <- function(...) {
  .Deprecated("plot_hydrograph")
  plot_hydrograph(...)
}

#' Deprecated wrapper for `write.riv()`
#'
#' Legacy compatibility wrapper for \code{\link{write_river}}.
#'
#' @param riv SHUD.RIVER object or compatible river data.
#' @param file Character. Output .riv file path.
#' @return Invisibly returns the output file path.
#' @export
write.riv <- function(riv, file) {
  .Deprecated("write_river")
  write_river(riv, file)
}

#' Deprecated wrapper for `sp.RiverSeg()`
#'
#' Legacy compatibility wrapper for \code{\link{shud.rivseg}}.
#'
#' @param ... Arguments passed to \code{\link{shud.rivseg}}.
#' @return The result from \code{\link{shud.rivseg}}.
#' @export
sp.RiverSeg <- function(...) {
  .Deprecated("shud.rivseg")
  shud.rivseg(...)
}

#' Deprecated stub for `compareMaps()`
#'
#' This legacy export is retained for transition compatibility only.
#'
#' @param ... Ignored legacy arguments.
#' @return \code{NULL}.
#' @export
compareMaps <- function(...) {
  .Deprecated()
  NULL
}

#' Deprecated stub for `map2d()`
#'
#' This legacy export is retained for transition compatibility only.
#'
#' @param ... Ignored legacy arguments.
#' @return \code{NULL}.
#' @export
map2d <- function(...) {
  .Deprecated()
  NULL
}

#' Deprecated wrapper for `plot_sp()`
#'
#' Uses modern plotting methods when possible.
#'
#' @param x Object to plot. Supported modern inputs include \code{sf},
#'   \code{sfc}, and \code{terra::SpatVector}; other inputs are forwarded to
#'   \code{\link[graphics]{plot}}.
#' @param field Character. Attribute field to use for coloring.
#' @param zcol Character. Legacy alias for \code{field}.
#' @param ... Additional arguments passed to the selected plotting function.
#' @return Invisibly returns \code{x}.
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
#' Legacy compatibility wrapper for \code{\link{plot_timeseries}}. Hydrograph
#' inputs are forwarded to \code{\link{plot_timeseries}}; other base-plottable
#' inputs use \code{\link[graphics]{plot}} to preserve legacy behavior.
#'
#' @param x Time-series object to plot.
#' @param ... Additional arguments passed to the selected plotting function.
#' @return A plot object from \code{\link{plot_timeseries}} for hydrograph
#'   inputs, otherwise \code{x} invisibly.
#' @export
plot_tsd <- function(x, ...) {
  .Deprecated("plot_timeseries")
  if (.is_hydrograph_timeseries(x)) {
    return(plot_timeseries(x, ...))
  }

  graphics::plot(x, ...)
  invisible(x)
}

.is_hydrograph_timeseries <- function(x) {
  inherits(x, c("xts", "zoo")) &&
    !is.null(ncol(x)) &&
    ncol(x) >= 2
}

#' Deprecated wrapper for `png.control()`
#'
#' Legacy compatibility wrapper for \code{grDevices::png()}.
#'
#' @param ... Arguments passed to \code{\link[grDevices]{png}}.
#' @return The result from \code{\link[grDevices]{png}}.
#' @export
png.control <- function(...) {
  .Deprecated(msg = "png.control() is deprecated. Use grDevices::png() instead.")
  grDevices::png(...)
}

#' Deprecated wrapper for `removeholes()`
#'
#' Uses modern terra/sf polygon handling to fill interior holes.
#'
#' @param sp Spatial, sf, or SpatVector polygon object.
#' @param ... Additional arguments passed to \code{\link[terra]{fillHoles}}.
#' @return Object of the same general spatial class as \code{sp}, or
#'   \code{NULL} for unsupported inputs.
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
#'
#' @param pm SHUD.MESH object containing mesh and point data.
#' @param file Character. Output .mesh file path.
#' @return Invisibly returns the output file path.
#' @export
write.mesh <- function(pm, file) {
  .Deprecated("write_mesh")
  write_mesh(pm, file)
}

#' Deprecated wrapper for `write.ic()`
#'
#' Legacy compatibility wrapper for \code{\link{write_ic}}.
#'
#' @param x List containing initial condition data.
#' @param file Character. Output .ic file path.
#' @return Invisibly returns the output file path.
#' @export
write.ic <- function(x, file) {
  .Deprecated("write_ic")
  write_ic(x, file)
}

#' Deprecated wrapper for `write.df()`
#'
#' Legacy compatibility wrapper for \code{\link{write_df}}.
#'
#' @param x Data frame or matrix to write.
#' @param file Character. Output file path.
#' @param append Logical. Whether to append to an existing file.
#' @param quiet Logical. Whether to suppress messages.
#' @param header Numeric vector. Custom header; defaults to row/column counts.
#' @return Invisibly returns the output file path.
#' @export
write.df <- function(x, file, append = FALSE, quiet = FALSE, header = NULL) {
  .Deprecated("write_df")
  write_df(x, file, append = append, quiet = quiet, header = header)
}

#' Deprecated wrapper for `write.config()`
#'
#' Legacy compatibility wrapper for \code{\link{write_config}}.
#'
#' @param x SHUD model configuration parameters or calibration data.
#' @param file Character. Output configuration file path.
#' @return Invisibly returns the output file path.
#' @export
write.config <- function(x, file) {
  .Deprecated("write_config")
  write_config(x, file)
}

#' Deprecated wrapper for `write.forc()`
#'
#' Legacy compatibility wrapper for \code{\link{write_forc}}.
#'
#' @param x Data frame of forcing-site information.
#' @param file Character. Output .forc file path.
#' @param path Character. Common path of the forcing files.
#' @param startdate Character. Start date in YYYYMMDD format.
#' @return Invisibly returns the output file path.
#' @export
write.forc <- function(x, file, path = "", startdate = "20000101") {
  .Deprecated("write_forc")
  write_forc(x, file, path = path, startdate = startdate)
}
