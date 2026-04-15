#' Coordinate Projection Functions
#'
#' Functions for creating and working with coordinate reference systems (CRS)
#' using modern spatial libraries (sf).
#'
#' @name projection
NULL

#' Build an Albers Equal Area projection
#'
#' Creates an Albers Equal Area Conic projection based on spatial data extent.
#' This projection is suitable for areas with east-west extent and preserves
#' area measurements.
#'
#' @param spx Spatial data object (sf, SpatVector, or SpatRaster). If NULL,
#'   must provide ext parameter.
#' @param ext Numeric vector of extent: c(xmin, xmax, ymin, ymax) in
#'   longitude/latitude (EPSG:4326). If NULL, extracted from spx.
#' @return Character string with PROJ.4 CRS definition for Albers projection
#' @export
#' @examples
#' # Create Albers projection from extent
#' crs_albers <- crs.Albers(ext = c(-120, -110, 35, 45))
#'
#' # Use with spatial data
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   pts <- sf::st_as_sf(
#'     data.frame(x = c(-115, -112), y = c(40, 42)),
#'     coords = c("x", "y"),
#'     crs = 4326
#'   )
#'   crs_albers <- crs.Albers(spx = pts)
#' }
crs.Albers <- function(spx = NULL, ext = NULL) {
  # Validate inputs
  if (is.null(spx) && is.null(ext)) {
    stop("Either 'spx' or 'ext' must be provided", call. = FALSE)
  }

  # Reject legacy spatial objects
  if (!is.null(spx)) {
    if (inherits(spx, c("Spatial", "Raster", "RasterLayer", "RasterStack",
                        "RasterBrick", "SpatialPoints", "SpatialPolygons",
                        "SpatialLines", "SpatialPointsDataFrame",
                        "SpatialPolygonsDataFrame", "SpatialLinesDataFrame"))) {
      stop(
        "Legacy raster/sp objects are not supported. ",
        "Please convert to terra/sf format:\n",
        "  - For raster: use terra::rast()\n",
        "  - For sp: use sf::st_as_sf()\n",
        "See inst/MIGRATION_GUIDE.md for details.",
        call. = FALSE
      )
    }
  }

  # Extract extent from spatial object if not provided
  if (is.null(ext)) {
    # Transform to WGS84 if needed
    if (inherits(spx, c("SpatRaster", "SpatVector"))) {
      # Check if CRS is defined
      if (terra::crs(spx) == "") {
        stop("Spatial object must have a defined CRS", call. = FALSE)
      }
      # Transform to EPSG:4326 if not already
      crs_str <- terra::crs(spx)
      if (!grepl("4326", crs_str, fixed = TRUE)) {
        spx <- terra::project(spx, "EPSG:4326")
      }
      ext_obj <- terra::ext(spx)
      ext <- c(ext_obj$xmin, ext_obj$xmax, ext_obj$ymin, ext_obj$ymax)
    } else if (inherits(spx, "sf")) {
      # Check if CRS is defined
      if (is.na(sf::st_crs(spx))) {
        stop("Spatial object must have a defined CRS", call. = FALSE)
      }
      # Transform to EPSG:4326 if not already
      if (sf::st_crs(spx)$epsg != 4326) {
        spx <- sf::st_transform(spx, 4326)
      }
      bb <- sf::st_bbox(spx)
      ext <- c(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"])
    } else {
      stop(
        "'spx' must be a SpatRaster, SpatVector, or sf object",
        call. = FALSE
      )
    }
  }

  # Validate extent
  if (!is.numeric(ext) || length(ext) != 4) {
    stop("'ext' must be a numeric vector of length 4: c(xmin, xmax, ymin, ymax)",
         call. = FALSE)
  }

  # Calculate projection parameters
  x0 <- round(mean(ext[1:2]), 1)
  y0 <- round(mean(ext[3:4]), 1)
  dx <- round(diff(ext[1:2]), 2)
  dy <- round(diff(ext[3:4]), 2)

  # Calculate standard parallels
  if (dx < 1) {
    my <- 0.25
  } else {
    my <- round(dy / 4, 2)
  }
  lat1 <- y0 + my
  lat2 <- y0 - my

  # Build PROJ.4 string
  crs_str <- paste0(
    "+proj=aea",
    " +lat_1=", lat1,
    " +lat_2=", lat2,
    " +lon_0=", x0,
    " +datum=WGS84",
    " +units=m"
  )

  message("Albers Equal Area projection: ", crs_str)

  return(crs_str)
}
#' Build a Lambert Equal Area projection
#'
#' Creates a Lambert Azimuthal Equal Area projection based on spatial data
#' extent. This projection is suitable for areas with similar north-south and
#' east-west extent and preserves area measurements.
#'
#' @param spx Spatial data object (sf, SpatVector, or SpatRaster). If NULL,
#'   must provide ext parameter.
#' @param ext Numeric vector of extent: c(xmin, xmax, ymin, ymax) in
#'   longitude/latitude (EPSG:4326). If NULL, extracted from spx.
#' @return Character string with PROJ.4 CRS definition for Lambert projection
#' @export
#' @examples
#' # Create Lambert projection from extent
#' crs_lambert <- crs.Lambert(ext = c(-120, -110, 35, 45))
#'
#' # Use with spatial data
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   pts <- sf::st_as_sf(
#'     data.frame(x = c(-115, -112), y = c(40, 42)),
#'     coords = c("x", "y"),
#'     crs = 4326
#'   )
#'   crs_lambert <- crs.Lambert(spx = pts)
#' }
crs.Lambert <- function(spx = NULL, ext = NULL) {
  # Validate inputs
  if (is.null(spx) && is.null(ext)) {
    stop("Either 'spx' or 'ext' must be provided", call. = FALSE)
  }

  # Reject legacy spatial objects
  if (!is.null(spx)) {
    if (inherits(spx, c("Spatial", "Raster", "RasterLayer", "RasterStack",
                        "RasterBrick", "SpatialPoints", "SpatialPolygons",
                        "SpatialLines", "SpatialPointsDataFrame",
                        "SpatialPolygonsDataFrame", "SpatialLinesDataFrame"))) {
      stop(
        "Legacy raster/sp objects are not supported. ",
        "Please convert to terra/sf format:\n",
        "  - For raster: use terra::rast()\n",
        "  - For sp: use sf::st_as_sf()\n",
        "See inst/MIGRATION_GUIDE.md for details.",
        call. = FALSE
      )
    }
  }

  # Extract extent from spatial object if not provided
  if (is.null(ext)) {
    # Transform to WGS84 if needed
    if (inherits(spx, c("SpatRaster", "SpatVector"))) {
      # Check if CRS is defined
      if (terra::crs(spx) == "") {
        stop("Spatial object must have a defined CRS", call. = FALSE)
      }
      # Transform to EPSG:4326 if not already
      crs_str <- terra::crs(spx)
      if (!grepl("4326", crs_str, fixed = TRUE)) {
        spx <- terra::project(spx, "EPSG:4326")
      }
      ext_obj <- terra::ext(spx)
      ext <- c(ext_obj$xmin, ext_obj$xmax, ext_obj$ymin, ext_obj$ymax)
    } else if (inherits(spx, "sf")) {
      # Check if CRS is defined
      if (is.na(sf::st_crs(spx))) {
        stop("Spatial object must have a defined CRS", call. = FALSE)
      }
      # Transform to EPSG:4326 if not already
      if (sf::st_crs(spx)$epsg != 4326) {
        spx <- sf::st_transform(spx, 4326)
      }
      bb <- sf::st_bbox(spx)
      ext <- c(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"])
    } else {
      stop(
        "'spx' must be a SpatRaster, SpatVector, or sf object",
        call. = FALSE
      )
    }
  }

  # Validate extent
  if (!is.numeric(ext) || length(ext) != 4) {
    stop("'ext' must be a numeric vector of length 4: c(xmin, xmax, ymin, ymax)",
         call. = FALSE)
  }

  # Calculate projection parameters
  x0 <- round(mean(ext[1:2]), 1)
  y0 <- round(mean(ext[3:4]), 1)

  # Build PROJ.4 string
  crs_str <- paste0(
    "+proj=laea",
    " +lat_0=", y0,
    " +lon_0=", x0,
    " +datum=WGS84",
    " +units=m"
  )

  message("Lambert Azimuthal Equal Area projection: ", crs_str)

  return(crs_str)
}

#' Get UTM zone from longitude
#'
#' Calculates the UTM zone number from longitude coordinate(s).
#'
#' @param lon Numeric vector of longitude values in degrees (-180 to 180)
#' @return Integer vector of UTM zone numbers (1-60)
#' @export
#' @examples
#' # Single longitude
#' crs.long2utmZone(-75)  # Returns 18
#'
#' # Multiple longitudes
#' lons <- seq(-180, 179, 6) + 1
#' zones <- crs.long2utmZone(lons)
crs.long2utmZone <- function(lon) {
  # Validate input
  if (!is.numeric(lon)) {
    stop("'lon' must be numeric", call. = FALSE)
  }

  if (any(lon > 180 | lon < -180, na.rm = TRUE)) {
    stop(
      "Longitude values must be in range (-180, 180). ",
      "Range of provided data: [", min(lon, na.rm = TRUE), ", ",
      max(lon, na.rm = TRUE), "]",
      call. = FALSE
    )
  }

  # Calculate UTM zone
  zid <- floor((lon + 180) / 6) + 1

  return(zid)
}

#' Get UTM projection from longitude and latitude
#'
#' Creates a UTM (Universal Transverse Mercator) projection CRS string based
#' on longitude and latitude coordinates.
#'
#' @param lon Numeric, longitude in degrees (-180 to 180)
#' @param lat Numeric, latitude in degrees (-90 to 90). Default is 30 (Northern
#'   Hemisphere). Use negative values for Southern Hemisphere.
#' @return Character string with PROJ.4 CRS definition for UTM projection.
#'   If lon is a vector, returns a list of CRS strings.
#' @export
#' @examples
#' # Single location (Northern Hemisphere)
#' crs_utm <- crs.long2utm(-75, 40)
#'
#' # Southern Hemisphere
#' crs_utm_south <- crs.long2utm(150, -30)
#'
#' # Multiple longitudes
#' lons <- seq(-180, 179, 6) + 1
#' crs_list <- crs.long2utm(lons)
crs.long2utm <- function(lon, lat = 30) {
  # Validate longitude
  if (!is.numeric(lon)) {
    stop("'lon' must be numeric", call. = FALSE)
  }

  if (any(lon > 180 | lon < -180, na.rm = TRUE)) {
    stop(
      "Longitude values must be in range (-180, 180). ",
      "Range of provided data: [", min(lon, na.rm = TRUE), ", ",
      max(lon, na.rm = TRUE), "]",
      call. = FALSE
    )
  }

  # Validate latitude
  if (!is.numeric(lat)) {
    stop("'lat' must be numeric", call. = FALSE)
  }

  if (any(abs(lat) > 90, na.rm = TRUE)) {
    stop(
      "Latitude values must be in range (-90, 90). ",
      "Range of provided data: [", min(lat, na.rm = TRUE), ", ",
      max(lat, na.rm = TRUE), "]",
      call. = FALSE
    )
  }

  # Handle vector input
  if (length(lon) > 1) {
    # Recycle lat if needed
    if (length(lat) == 1) {
      lat <- rep(lat, length(lon))
    } else if (length(lat) != length(lon)) {
      stop("'lat' must be either length 1 or same length as 'lon'",
           call. = FALSE)
    }

    # Apply recursively
    result <- mapply(crs.long2utm, lon, lat, SIMPLIFY = FALSE)
    return(result)
  }

  # Calculate UTM zone
  zid <- crs.long2utmZone(lon)

  # Build PROJ.4 string
  if (lat >= 0) {
    crs_str <- paste0(
      "+proj=utm +zone=", zid,
      " +datum=WGS84 +units=m +no_defs"
    )
  } else {
    crs_str <- paste0(
      "+proj=utm +zone=", zid,
      " +south +datum=WGS84 +units=m +no_defs"
    )
  }

  return(crs_str)
}

#' Project coordinates to a new CRS
#'
#' Unified interface for projecting spatial coordinates or objects to a new
#' coordinate reference system. This function provides a consistent interface
#' for coordinate transformation across different spatial object types.
#'
#' @param x Spatial object (sf, SpatVector, SpatRaster) or numeric matrix/
#'   data.frame with coordinates (columns: x, y or lon, lat)
#' @param from_crs Source CRS. If NULL, extracted from x (for spatial objects)
#'   or assumed to be EPSG:4326 (for coordinate matrices)
#' @param to_crs Target CRS (character string, EPSG code, or PROJ.4 string)
#' @param coords Character vector of length 2 specifying coordinate column
#'   names for matrix/data.frame input. Default: c("x", "y")
#' @return Transformed spatial object or coordinate matrix/data.frame
#' @export
#' @examples
#' # Project sf object
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   pts <- sf::st_as_sf(
#'     data.frame(x = c(-75, -76), y = c(40, 41)),
#'     coords = c("x", "y"),
#'     crs = 4326
#'   )
#'   pts_utm <- project_coords(pts, to_crs = crs.long2utm(-75, 40))
#' }
#'
#' # Project coordinate matrix
#' coords <- matrix(c(-75, 40, -76, 41), ncol = 2, byrow = TRUE)
#' coords_utm <- project_coords(coords, from_crs = "EPSG:4326",
#'                               to_crs = crs.long2utm(-75, 40))
project_coords <- function(x, from_crs = NULL, to_crs,
                          coords = c("x", "y")) {
  # Validate to_crs
  if (missing(to_crs) || is.null(to_crs)) {
    stop("'to_crs' is required", call. = FALSE)
  }

  # Reject legacy spatial objects
  if (inherits(x, c("Spatial", "Raster", "RasterLayer", "RasterStack",
                    "RasterBrick", "SpatialPoints", "SpatialPolygons",
                    "SpatialLines", "SpatialPointsDataFrame",
                    "SpatialPolygonsDataFrame", "SpatialLinesDataFrame"))) {
    stop(
      "Legacy raster/sp objects are not supported. ",
      "Please convert to terra/sf format:\n",
      "  - For raster: use terra::rast()\n",
      "  - For sp: use sf::st_as_sf()\n",
      "See inst/MIGRATION_GUIDE.md for details.",
      call. = FALSE
    )
  }

  # Handle different input types
  if (inherits(x, c("SpatRaster", "SpatVector"))) {
    # Use terra::project
    if (!is.null(from_crs)) {
      warning("'from_crs' is ignored for SpatRaster/SpatVector objects; ",
              "using CRS from object", call. = FALSE)
    }
    result <- terra::project(x, to_crs)
    return(result)

  } else if (inherits(x, "sf")) {
    # Use sf::st_transform
    if (!is.null(from_crs)) {
      warning("'from_crs' is ignored for sf objects; using CRS from object",
              call. = FALSE)
    }
    result <- sf::st_transform(x, to_crs)
    return(result)

  } else if (is.matrix(x) || is.data.frame(x)) {
    # Handle coordinate matrix/data.frame
    if (is.matrix(x)) {
      if (ncol(x) != 2) {
        stop("Coordinate matrix must have exactly 2 columns (x, y)",
             call. = FALSE)
      }
      coord_mat <- x
      colnames(coord_mat) <- coords
    } else {
      # data.frame
      if (!all(coords %in% names(x))) {
        stop("Coordinate columns '", paste(coords, collapse = "', '"),
             "' not found in data.frame", call. = FALSE)
      }
      coord_mat <- as.matrix(x[, coords])
    }

    # Set default from_crs if not provided
    if (is.null(from_crs)) {
      from_crs <- "EPSG:4326"
      message("Assuming source CRS is EPSG:4326 (WGS84)")
    }

    # Create sf object for transformation
    pts_sf <- sf::st_as_sf(
      as.data.frame(coord_mat),
      coords = coords,
      crs = from_crs
    )

    # Transform
    pts_transformed <- sf::st_transform(pts_sf, to_crs)

    # Extract coordinates
    coords_transformed <- sf::st_coordinates(pts_transformed)

    # Return in same format as input
    if (is.matrix(x)) {
      return(coords_transformed)
    } else {
      # data.frame - preserve other columns
      result <- x
      result[, coords] <- coords_transformed
      return(result)
    }

  } else {
    stop(
      "'x' must be a SpatRaster, SpatVector, sf object, matrix, or data.frame",
      call. = FALSE
    )
  }
}
