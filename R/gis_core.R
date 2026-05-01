#' Core GIS Processing Functions
#'
#' Functions for core GIS operations including raster-vector conversions
#' and spatial interpolation using modern spatial libraries (terra and sf).
#'
#' @name gis_core
#' @keywords internal
NULL

#' Convert mesh data to raster
#'
#' Converts mesh (unstructured grid) data to a regular raster grid using
#' spatial interpolation. This function uses terra and sf for all spatial
#' operations and supports multiple interpolation methods.
#'
#' @param data Numeric vector or matrix of values to rasterize.
#'   If matrix, each row represents a time step and each column a mesh cell.
#'   Length/ncol must equal the number of mesh cells.
#' @param mesh Mesh object (optional). If NULL, uses default mesh from
#'   \code{read_mesh()}.
#' @param template SpatRaster template defining output grid (optional).
#'   If NULL, creates template from mesh extent.
#' @param method Character, interpolation method. One of:
#'   \itemize{
#'     \item "idw" - Inverse Distance Weighting (default)
#'     \item "linear" - Linear interpolation
#'     \item "nearest" - Nearest neighbor
#'   }
#' @param resolution Numeric, output raster resolution in map units (optional).
#'   Only used if template is NULL.
#' @param crs Character or CRS object, coordinate reference system (optional).
#' @param stack Logical, if TRUE and data is a matrix, returns a multi-layer
#'   SpatRaster with one layer per row (default: FALSE).
#'
#' @return SpatRaster object. If stack=TRUE and data is a matrix, returns
#'   a multi-layer raster.
#'
#' @section Migration Note:
#' This function replaces \code{MeshData2Raster()} and uses terra instead of
#' the legacy raster package. Key differences:
#' \itemize{
#'   \item Returns SpatRaster instead of RasterLayer
#'   \item Does not accept legacy spatial objects - use terra/sf only
#'   \item Method names simplified: "ide" -> "idw"
#'   \item Removed plot parameter - use terra::plot() on result
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Basic usage with vector data
#' mesh <- read_mesh("model.mesh")
#' elevation <- getElevation()
#' r <- mesh_to_raster(elevation, mesh = mesh)
#' terra::plot(r)
#'
#' # With custom resolution
#' r <- mesh_to_raster(elevation, mesh = mesh, resolution = 100)
#'
#' # Time series data (matrix)
#' timeseries_data <- matrix(rnorm(1000 * 10), nrow = 10, ncol = 1000)
#' r_stack <- mesh_to_raster(timeseries_data, mesh = mesh, stack = TRUE)
#' }
mesh_to_raster <- function(data,
                           mesh = NULL,
                           template = NULL,
                           method = c("idw", "linear", "nearest"),
                           resolution = NULL,
                           crs = NULL,
                           stack = FALSE) {

  # Match method argument
  method <- match.arg(method)

  # Validate data parameter
  if (missing(data) || is.null(data)) {
    stop("Parameter 'data' is required", call. = FALSE)
  }

  # Modern API rejects legacy raster/sp objects.
  if (inherits(data, c("Raster", "RasterLayer", "RasterStack", "RasterBrick",
                       "Spatial", "SpatialPoints", "SpatialPolygons"))) {
    stop(
      "Legacy raster/sp objects are not supported.\n",
      "Please convert to terra/sf format:\n",
      "  - raster -> terra: use terra::rast()\n",
      "  - sp -> sf: use sf::st_as_sf()\n",
      "See inst/MIGRATION_GUIDE.md for details.",
      call. = FALSE
    )
  }

  # Validate data is numeric
  if (!is.numeric(data) && !is.matrix(data) && !is.data.frame(data)) {
    stop(
      "Parameter 'data' must be a numeric vector, matrix, or data.frame, ",
      "but received ", class(data)[1],
      call. = FALSE
    )
  }

  # Handle stack processing for matrix input
  if (stack && (is.matrix(data) || is.data.frame(data))) {
    if (nrow(data) < 1) {
      stop("Matrix 'data' must have at least one row", call. = FALSE)
    }

    # Recursively process each row
    raster_list <- lapply(seq_len(nrow(data)), function(i) {
      mesh_to_raster(
        data = as.numeric(data[i, ]),
        mesh = mesh,
        template = template,
        method = method,
        resolution = resolution,
        crs = crs,
        stack = FALSE
      )
    })

    # Combine into multi-layer raster
    result <- terra::rast(raster_list)
    return(result)
  }

  # Convert matrix/data.frame to vector (use last row if matrix)
  if (is.matrix(data) || is.data.frame(data)) {
    if (nrow(data) > 1) {
      warning(
        "Matrix 'data' has multiple rows but stack=FALSE. ",
        "Using last row only. Set stack=TRUE to process all rows.",
        call. = FALSE
      )
    }
    data <- as.numeric(data[nrow(data), ])
  }

  # Handle NA and infinite values
  if (any(is.na(data))) {
    warning("NA values in 'data' will be set to 0", call. = FALSE)
    data[is.na(data)] <- 0
  }
  if (any(is.infinite(data))) {
    warning("Infinite values in 'data' will be set to 0", call. = FALSE)
    data[is.infinite(data)] <- 0
  }

  # Get mesh if not provided
  if (is.null(mesh)) {
    mesh <- read_mesh()
  }

  # Get mesh centroids
  xy <- getCentroid(pm = mesh)[, 2:3]

  # Validate data length matches mesh
  if (length(data) != nrow(xy)) {
    stop(
      "Length of 'data' (", length(data), ") does not match ",
      "number of mesh cells (", nrow(xy), ")",
      call. = FALSE
    )
  }

  # Create or validate template raster
  if (is.null(template)) {
    # Create template from mesh extent
    template <- shud.mask(pm = mesh, proj = crs)

    # Internal mask helpers now return terra; keep this guard for legacy masks.
    if (inherits(template, c("RasterLayer", "RasterStack", "RasterBrick"))) {
      template <- terra::rast(template)
    }

    # Apply resolution if specified
    if (!is.null(resolution)) {
      check_positive(resolution, "resolution")
      # Resample to new resolution
      ext <- terra::ext(template)
      new_template <- terra::rast(
        xmin = ext$xmin, xmax = ext$xmax,
        ymin = ext$ymin, ymax = ext$ymax,
        resolution = resolution
      )
      if (!is.null(crs)) {
        terra::crs(new_template) <- crs
      }
      template <- new_template
    }
  } else {
    # Validate template is a SpatRaster
    if (!inherits(template, "SpatRaster")) {
      if (inherits(template, c("RasterLayer", "RasterStack", "RasterBrick"))) {
        stop(
          "Legacy raster template not supported. ",
          "Convert to terra: template <- terra::rast(template)",
          call. = FALSE
        )
      } else {
        stop(
          "Parameter 'template' must be a SpatRaster object, ",
          "but received ", class(template)[1],
          call. = FALSE
        )
      }
    }
  }

  # Perform interpolation based on method
  if (method == "idw") {
    # Inverse Distance Weighting using gstat and terra
    if (!requireNamespace("gstat", quietly = TRUE)) {
      stop(
        "Package 'gstat' is required for IDW interpolation.\n",
        "Install it with: install.packages('gstat')",
        call. = FALSE
      )
    }

    # Create gstat model using data frame
    pts_df <- data.frame(x = xy[, 1], y = xy[, 2], z = data)
    gs <- gstat::gstat(
      formula = z ~ 1,
      locations = ~ x + y,
      data = pts_df,
      set = list(idp = 2.0)
    )

    # Perform IDW interpolation directly to raster template
    idw_result <- terra::interpolate(template, gs)

    # Extract the predicted values
    result <- idw_result[["var1.pred"]]

  } else if (method == "linear") {
    # Linear interpolation using interp package
    if (!requireNamespace("interp", quietly = TRUE)) {
      stop(
        "Package 'interp' is required for linear interpolation.\n",
        "Install it with: install.packages('interp')",
        call. = FALSE
      )
    }

    # Get grid coordinates
    grid_xy <- terra::crds(template, na.rm = FALSE)

    # Perform interpolation
    interp_result <- interp::interp(
      x = xy[, 1],
      y = xy[, 2],
      z = data,
      xo = unique(grid_xy[, 1]),
      yo = unique(grid_xy[, 2]),
      linear = TRUE,
      extrap = FALSE
    )

    # Convert to raster
    result <- template
    terra::values(result) <- as.numeric(interp_result$z)

  } else if (method == "nearest") {
    # Nearest neighbor using terra
    # Create points from mesh centroids
    pts_vect <- terra::vect(
      cbind(xy[, 1], xy[, 2]),
      type = "points",
      atts = data.frame(value = data),
      crs = terra::crs(template)
    )

    # Rasterize using nearest neighbor
    result <- terra::rasterize(
      x = pts_vect,
      y = template,
      field = "value",
      fun = "first"  # Use first (nearest) value
    )
  }

  # Apply CRS if specified
  if (!is.null(crs)) {
    terra::crs(result) <- crs
  }

  return(result)
}

#' Convert vector data to raster
#'
#' Converts vector spatial data (sf or SpatVector) to a regular raster grid.
#' This function uses terra for rasterization and supports automatic or manual
#' resolution settings and multiple aggregation functions.
#'
#' @param vector sf or SpatVector object to rasterize.
#' @param template SpatRaster template defining output grid (optional).
#'   If NULL, creates template from vector extent.
#' @param field Character or numeric, field name or index to rasterize.
#'   Default is 1 (first attribute field). Can also be a numeric vector
#'   of values matching the number of features.
#' @param resolution Numeric, output raster resolution in map units (optional).
#'   Only used if template is NULL. If NULL and template is NULL,
#'   resolution is automatically calculated from extent.
#' @param ngrids Integer, number of grid cells along x direction (optional).
#'   Only used if both template and resolution are NULL. Default is 200.
#' @param fun Character, aggregation function for multiple values per cell.
#'   Options: "first", "last", "min", "max", "mean", "sum", "count".
#'   Default is "last".
#' @param background Numeric, value for cells with no data. Default is NA.
#'
#' @return SpatRaster object with rasterized vector data.
#'
#' @section Migration Note:
#' This function replaces \code{sp2raster()} and uses terra/sf instead of
#' legacy spatial packages. Key differences:
#' \itemize{
#'   \item Returns SpatRaster instead of RasterLayer
#'   \item Accepts sf or SpatVector objects
#'   \item Does not use global MASK environment variable
#'   \item Parameter 'mask' renamed to 'template'
#'   \item More flexible field specification
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Basic usage with sf object
#' library(sf)
#' polygons <- st_read("watershed.shp")
#' r <- vector_to_raster(polygons, field = "elevation")
#' terra::plot(r)
#'
#' # With custom resolution
#' r <- vector_to_raster(polygons, field = "landcover", resolution = 100)
#'
#' # With template raster
#' template <- terra::rast(extent = terra::ext(polygons), resolution = 50)
#' r <- vector_to_raster(polygons, template = template, field = "soiltype")
#'
#' # Using aggregation function
#' points <- st_read("stations.shp")
#' r <- vector_to_raster(points, field = "rainfall", fun = "mean")
#' }
vector_to_raster <- function(vector,
                              template = NULL,
                              field = 1,
                              resolution = NULL,
                              ngrids = 200,
                              fun = "last",
                              background = NA) {

  # Validate vector parameter
  if (missing(vector) || is.null(vector)) {
    stop("Parameter 'vector' is required", call. = FALSE)
  }

  # Check for legacy sp objects and reject them
  if (inherits(vector, c("Spatial", "SpatialPoints", "SpatialPolygons",
                         "SpatialLines", "SpatialPointsDataFrame",
                         "SpatialPolygonsDataFrame", "SpatialLinesDataFrame"))) {
    stop(
      "Legacy sp objects are not supported.\n",
      "Please convert to sf or terra format:\n",
      "  - sp -> sf: use sf::st_as_sf()\n",
      "  - sp -> terra: use terra::vect()\n",
      "See inst/MIGRATION_GUIDE.md for details.",
      call. = FALSE
    )
  }

  # Validate vector is sf or SpatVector
  if (!inherits(vector, c("sf", "sfc", "SpatVector"))) {
    stop(
      "Parameter 'vector' must be an sf or SpatVector object, ",
      "but received ", class(vector)[1],
      call. = FALSE
    )
  }

  # Convert sf to SpatVector for consistent processing
  if (inherits(vector, c("sf", "sfc"))) {
    vector <- terra::vect(vector)
  }

  # Validate ngrids is positive
  if (!is.null(ngrids)) {
    check_positive(ngrids, "ngrids")
  }

  # Create or validate template raster
  if (is.null(template)) {
    # Get extent from vector
    ext <- terra::ext(vector)

    # Determine resolution
    if (is.null(resolution)) {
      # Calculate resolution from ngrids
      dx <- (ext$xmax - ext$xmin) / ngrids
    } else {
      check_positive(resolution, "resolution")
      dx <- resolution
    }

    # Create template raster
    template <- terra::rast(
      xmin = ext$xmin, xmax = ext$xmax,
      ymin = ext$ymin, ymax = ext$ymax,
      resolution = dx
    )

    # Set CRS from vector
    terra::crs(template) <- terra::crs(vector)

  } else {
    # Validate template is a SpatRaster
    if (!inherits(template, "SpatRaster")) {
      if (inherits(template, c("RasterLayer", "RasterStack", "RasterBrick"))) {
        stop(
          "Legacy raster template not supported. ",
          "Convert to terra: template <- terra::rast(template)",
          call. = FALSE
        )
      } else {
        stop(
          "Parameter 'template' must be a SpatRaster object, ",
          "but received ", class(template)[1],
          call. = FALSE
        )
      }
    }
  }

  # Validate field parameter
  if (is.numeric(field) && length(field) == 1) {
    # Field index - validate it exists
    if (field < 1 || field > ncol(vector)) {
      stop(
        "Field index ", field, " is out of range. ",
        "Vector has ", ncol(vector), " attribute fields.",
        call. = FALSE
      )
    }
    # Get field name
    field_name <- names(vector)[field]
  } else if (is.character(field) && length(field) == 1) {
    # Field name - validate it exists
    if (!field %in% names(vector)) {
      stop(
        "Field '", field, "' not found in vector attributes. ",
        "Available fields: ", paste(names(vector), collapse = ", "),
        call. = FALSE
      )
    }
    field_name <- field
  } else if (is.numeric(field) && length(field) == nrow(vector)) {
    # Vector of values - add as new field
    field_name <- "rasterize_values"
    vector[[field_name]] <- field
  } else {
    stop(
      "Parameter 'field' must be either:\n",
      "  - A single field name (character)\n",
      "  - A single field index (numeric)\n",
      "  - A numeric vector with length equal to number of features",
      call. = FALSE
    )
  }

  # Validate aggregation function
  valid_funs <- c("first", "last", "min", "max", "mean", "sum", "count")
  if (!fun %in% valid_funs) {
    stop(
      "Parameter 'fun' must be one of: ",
      paste(valid_funs, collapse = ", "),
      call. = FALSE
    )
  }

  # Perform rasterization
  result <- terra::rasterize(
    x = vector,
    y = template,
    field = field_name,
    fun = fun,
    background = background
  )

  return(result)
}

#' Deprecated: sp2raster
#'
#' @description
#' This function is deprecated. Use \code{\link{vector_to_raster}} instead.
#'
#' @param ... Arguments passed to \code{vector_to_raster}
#' @export
#' @keywords internal
sp2raster <- function(...) {
  .Deprecated("vector_to_raster")

  # Get arguments
  args <- list(...)

  # Map old parameter names to new ones
  if ("sp" %in% names(args)) {
    args$vector <- args$sp
    args$sp <- NULL
  }
  if ("mask" %in% names(args)) {
    args$template <- args$mask
    args$mask <- NULL
  }

  # Call new function
  do.call(vector_to_raster, args)
}
