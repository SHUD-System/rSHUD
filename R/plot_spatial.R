#' Spatial Plotting Functions
#'
#' Functions for visualizing spatial data including mesh domains, rasters,
#' and vector data using modern spatial libraries (terra and sf).
#'
#' @name plot_spatial
#' @keywords internal
NULL

#' Plot spatial polygons with attribute coloring
#'
#' Creates a choropleth map of spatial polygon data colored by an attribute field.
#' This function uses terra or sf for spatial operations.
#'
#' @param x sf or SpatVector polygon object to plot.
#' @param field Character, name of attribute field to use for coloring.
#' @param col_fun Function to generate color palette. Default is
#'   \code{topo.colors}.
#' @param n_colors Integer, number of colors in palette. Default is 20.
#' @param main Character, plot title (optional).
#' @param ... Additional arguments passed to plotting function.
#'
#' @return Invisibly returns the input object.
#'
#' @section Migration Note:
#' This function replaces \code{plot_sp()} and uses terra/sf instead of
#' legacy spatial packages. Key differences:
#' \itemize{
#'   \item Parameter 'zcol' renamed to 'field' for consistency
#'   \item Accepts sf or SpatVector objects
#'   \item Uses terra::plot() or sf::plot() instead of legacy raster plotting
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Plot mesh polygons colored by elevation
#' mesh_sf <- mesh_to_sf()
#' plot_polygons(mesh_sf, field = "elevation")
#'
#' # Custom color palette
#' plot_polygons(mesh_sf, field = "soiltype", col_fun = heat.colors)
#' }
plot_polygons <- function(x,
                          field,
                          col_fun = topo.colors,
                          n_colors = 20,
                          main = NULL,
                          ...) {

  # Validate x parameter
  if (missing(x) || is.null(x)) {
    stop("Parameter 'x' is required", call. = FALSE)
  }

  # Check for legacy sp objects
  if (inherits(x, c("Spatial", "SpatialPolygons", "SpatialPolygonsDataFrame"))) {
    stop(
      "Legacy sp objects are not supported.\n",
      "Please convert to sf or terra format:\n",
      "  - sp -> sf: use sf::st_as_sf()\n",
      "  - sp -> terra: use terra::vect()\n",
      "See inst/MIGRATION_GUIDE.md for details.",
      call. = FALSE
    )
  }

  # Validate x is sf or SpatVector
  if (!inherits(x, c("sf", "sfc", "SpatVector"))) {
    stop(
      "Parameter 'x' must be an sf or SpatVector object, ",
      "but received ", class(x)[1],
      call. = FALSE
    )
  }

  # Validate field parameter
  if (missing(field) || is.null(field)) {
    stop("Parameter 'field' is required", call. = FALSE)
  }

  # Get field values
  if (inherits(x, c("sf", "sfc"))) {
    if (!field %in% names(x)) {
      stop(
        "Field '", field, "' not found in object. ",
        "Available fields: ", paste(names(x), collapse = ", "),
        call. = FALSE
      )
    }
    z <- x[[field]]
  } else {
    # SpatVector
    if (!field %in% names(x)) {
      stop(
        "Field '", field, "' not found in object. ",
        "Available fields: ", paste(names(x), collapse = ", "),
        call. = FALSE
      )
    }
    z <- x[[field]]
  }

  if (is.data.frame(z)) z <- z[[1]]

  # Generate colors
  colors <- col_fun(length(z))

  # Sort by field value for better visualization
  ord <- order(z)

  # Plot based on object type
  if (inherits(x, c("sf", "sfc"))) {
    # Use sf plotting
    plot(x[ord, field], col = colors, main = main, ...)
  } else {
    # Use terra plotting
    terra::plot(x[ord], col = colors, main = main, ...)
  }

  # Return invisibly
  invisible(x)
}

#' Compare multiple spatial maps side by side
#'
#' Creates a multi-panel plot comparing multiple raster or vector maps.
#' This function uses terra for plotting and supports both raster and vector
#' spatial data.
#'
#' @param maps List of SpatRaster, sf, or SpatVector objects to compare.
#'   All objects should have compatible extents and CRS.
#' @param titles Character vector of titles for each map (optional).
#'   If NULL, uses names from the list or generates default titles.
#' @param nrow Integer, number of rows in the plot layout (optional).
#'   If NULL, automatically calculated from number of maps.
#' @param ncol Integer, number of columns in the plot layout (optional).
#'   If NULL, automatically calculated from number of maps.
#' @param col Color palette function or vector of colors. Default uses
#'   \code{terrain.colors()}.
#' @param contour Logical, if TRUE adds contour lines to raster plots.
#'   Default is FALSE.
#' @param axes Logical, if TRUE shows axes. Default is FALSE.
#' @param box Logical, if TRUE shows box around plots. Default is FALSE.
#' @param mar Numeric vector of length 4 specifying plot margins.
#'   Default is c(2, 2, 2, 2).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Invisibly returns NULL.
#'
#' @section Migration Note:
#' This function replaces \code{compareMaps()} and uses terra/sf instead of
#' legacy spatial packages. Key differences:
#' \itemize{
#'   \item Parameter 'r' renamed to 'maps' for clarity
#'   \item Parameter 'mfrow' split into 'nrow' and 'ncol' for flexibility
#'   \item Accepts SpatRaster, sf, or SpatVector objects
#'   \item Uses terra::plot() instead of legacy raster plotting
#'   \item Automatically calculates layout if not specified
#'   \item Better handling of mixed raster/vector inputs
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Compare multiple rasters
#' library(terra)
#' r1 <- rast(volcano)
#' r2 <- sin(r1 / 100)
#' r3 <- cos(r1 / 100)
#' compare_maps(list(r1, r2, r3), titles = c("Original", "Sin", "Cos"))
#'
#' # With contour lines
#' compare_maps(list(r1, r2), contour = TRUE)
#'
#' # Custom layout
#' compare_maps(list(r1, r2, r3, r1 + r2), nrow = 2, ncol = 2)
#' }
compare_maps <- function(maps,
                         titles = NULL,
                         nrow = NULL,
                         ncol = NULL,
                         col = grDevices::terrain.colors,
                         contour = FALSE,
                         axes = FALSE,
                         box = FALSE,
                         mar = c(2, 2, 2, 2),
                         ...) {

  # Validate maps parameter
  if (missing(maps) || is.null(maps) || !is.list(maps)) {
    stop(
      "Parameter 'maps' must be a list of spatial objects",
      call. = FALSE
    )
  }

  n_maps <- length(maps)
  if (n_maps == 0) {
    stop("Parameter 'maps' must contain at least one map", call. = FALSE)
  }

  # Check for legacy raster/sp objects
  for (i in seq_along(maps)) {
    if (inherits(maps[[i]], c("RasterLayer", "RasterStack", "RasterBrick",
                               "Spatial", "SpatialPolygons", "SpatialLines"))) {
      stop(
        "Legacy raster/sp objects are not supported.\n",
        "Please convert to terra/sf format:\n",
        "  - raster -> terra: use terra::rast()\n",
        "  - sp -> sf: use sf::st_as_sf()\n",
        "See inst/MIGRATION_GUIDE.md for details.",
        call. = FALSE
      )
    }
  }

  # Validate all objects are spatial
  for (i in seq_along(maps)) {
    if (!inherits(maps[[i]], c("SpatRaster", "sf", "sfc", "SpatVector"))) {
      stop(
        "All elements in 'maps' must be SpatRaster, sf, or SpatVector objects. ",
        "Element ", i, " is ", class(maps[[i]])[1],
        call. = FALSE
      )
    }
  }

  # Determine layout
  if (is.null(nrow) && is.null(ncol)) {
    # Auto-calculate layout
    if (n_maps <= 2) {
      nrow <- 1
      ncol <- n_maps
    } else if (n_maps <= 4) {
      nrow <- 2
      ncol <- 2
    } else if (n_maps <= 6) {
      nrow <- 2
      ncol <- 3
    } else if (n_maps <= 9) {
      nrow <- 3
      ncol <- 3
    } else {
      nrow <- ceiling(sqrt(n_maps))
      ncol <- ceiling(n_maps / nrow)
    }
  } else if (is.null(nrow)) {
    nrow <- ceiling(n_maps / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(n_maps / nrow)
  }

  # Validate layout can fit all maps
  if (nrow * ncol < n_maps) {
    stop(
      "Layout (", nrow, " x ", ncol, ") is too small for ",
      n_maps, " maps",
      call. = FALSE
    )
  }

  # Prepare titles
  if (is.null(titles)) {
    if (!is.null(names(maps))) {
      titles <- names(maps)
    } else {
      titles <- paste("Map", seq_len(n_maps))
    }
  } else {
    if (length(titles) != n_maps) {
      warning(
        "Length of 'titles' (", length(titles), ") does not match ",
        "number of maps (", n_maps, "). Using default titles.",
        call. = FALSE
      )
      titles <- paste("Map", seq_len(n_maps))
    }
  }

  # Prepare color palette
  if (is.function(col)) {
    colors <- col(100)
  } else {
    colors <- col
  }

  # Save current par settings
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  # Set up multi-panel layout
  graphics::par(mfrow = c(nrow, ncol), mar = mar)

  # Plot each map
  for (i in seq_len(n_maps)) {
    map <- maps[[i]]

    # Plot based on object type
    if (inherits(map, "SpatRaster")) {
      # Raster plot
      suppressWarnings(terra::plot(map, col = colors, main = titles[i],
                  axes = axes, box = box, ...))

      # Add contour if requested
      if (contour) {
        suppressWarnings(terra::contour(map, add = TRUE))
      }

    } else if (inherits(map, c("sf", "sfc"))) {
      # sf plot
      suppressWarnings(plot(map, col = colors, main = titles[i],
           axes = axes, ...))

    } else if (inherits(map, "SpatVector")) {
      # SpatVector plot
      suppressWarnings(terra::plot(map, col = colors, main = titles[i],
                  axes = axes, ...))
    }

    # Add grid
    graphics::grid()
  }

  # Return invisibly
  invisible(NULL)
}
