#' Parameter Validation Functions
#'
#' Helper functions for validating function parameters and inputs.
#'
#' @name validators
#' @keywords internal
NULL

#' Check if a numeric value is positive
#'
#' Validates that a numeric parameter is positive (greater than zero).
#'
#' @param x Numeric value to check
#' @param name Character string, name of the parameter (for error messages)
#' @param allow_zero Logical, whether to allow zero (default: FALSE)
#' @return Invisible TRUE if valid, otherwise stops with error
#' @export
#' @examples
#' \dontrun{
#' check_positive(5, "my_param")
#' check_positive(0, "count", allow_zero = TRUE)
#' }
check_positive <- function(x, name = "value", allow_zero = FALSE) {
  # Check if numeric
  if (!is.numeric(x) || length(x) != 1) {
    stop("'", name, "' must be a single numeric value", call. = FALSE)
  }

  # Check if NA or infinite
  if (is.na(x) || is.infinite(x)) {
    stop("'", name, "' cannot be NA or infinite", call. = FALSE)
  }

  # Check if positive
  if (allow_zero) {
    if (x < 0) {
      stop("'", name, "' must be non-negative (>= 0)", call. = FALSE)
    }
  } else {
    if (x <= 0) {
      stop("'", name, "' must be positive (> 0)", call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Check if a file exists
#'
#' Validates that a file path exists and is accessible.
#'
#' @param path Character string, file path to check
#' @param name Character string, name of the parameter (for error messages)
#' @param must_be_file Logical, require path to be a file not directory
#'   (default: TRUE)
#' @return Invisible TRUE if valid, otherwise stops with error
#' @export
#' @examples
#' \dontrun{
#' check_file_exists("data/input.shp", "input_file")
#' }
check_file_exists <- function(path, name = "file", must_be_file = TRUE) {
  # Check if character
  if (!is.character(path) || length(path) != 1) {
    stop("'", name, "' must be a single character string", call. = FALSE)
  }

  # Check if empty
  if (nchar(path) == 0) {
    stop("'", name, "' cannot be an empty string", call. = FALSE)
  }

  # Check if exists
  if (!file.exists(path)) {
    stop("'", name, "' does not exist: ", path, call. = FALSE)
  }

  # Check if file (not directory)
  if (must_be_file && dir.exists(path)) {
    stop("'", name, "' must be a file, not a directory: ", path,
         call. = FALSE)
  }

  invisible(TRUE)
}

#' Check if spatial objects are compatible
#'
#' Validates that two spatial objects have compatible properties
#' (e.g., same CRS, overlapping extents).
#'
#' @param x First spatial object (SpatRaster, SpatVector, or sf)
#' @param y Second spatial object (SpatRaster, SpatVector, or sf)
#' @param check_crs Logical, check if CRS match (default: TRUE)
#' @param check_overlap Logical, check if extents overlap (default: FALSE)
#' @param name_x Character string, name of first object (for error messages)
#' @param name_y Character string, name of second object (for error messages)
#' @return Invisible TRUE if valid, otherwise stops with error
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   r1 <- terra::rast(ncol=10, nrow=10)
#'   r2 <- terra::rast(ncol=10, nrow=10)
#'   terra::crs(r1) <- terra::crs(r2) <- "EPSG:4326"
#'   check_spatial_compatible(r1, r2)
#' }
#' }
check_spatial_compatible <- function(x, y, check_crs = TRUE,
                                     check_overlap = FALSE,
                                     name_x = "x", name_y = "y") {
  # Check if both are spatial objects (terra/sf only)
  is_spatial_x <- inherits(x, c("SpatRaster", "SpatVector", "sf"))
  is_spatial_y <- inherits(y, c("SpatRaster", "SpatVector", "sf"))

  if (!is_spatial_x) {
    stop("'", name_x, "' must be a SpatRaster, SpatVector, or sf object",
         call. = FALSE)
  }
  if (!is_spatial_y) {
    stop("'", name_y, "' must be a SpatRaster, SpatVector, or sf object",
         call. = FALSE)
  }

  # Check CRS compatibility
  if (check_crs) {
    # Get CRS from both objects
    if (inherits(x, c("SpatRaster", "SpatVector"))) {
      crs_x <- terra::crs(x)
    } else {
      crs_x <- sf::st_crs(x)
    }

    if (inherits(y, c("SpatRaster", "SpatVector"))) {
      crs_y <- terra::crs(y)
    } else {
      crs_y <- sf::st_crs(y)
    }

    if (!compatible_crs(crs_x, crs_y)) {
      stop("'", name_x, "' and '", name_y, "' have incompatible CRS",
           call. = FALSE)
    }
  }

  # Check extent overlap
  if (check_overlap) {
    # Get extents
    if (inherits(x, c("SpatRaster", "SpatVector"))) {
      ext_x <- terra::ext(x)
      ext_x <- c(xmin = unname(ext_x$xmin), xmax = unname(ext_x$xmax),
                 ymin = unname(ext_x$ymin), ymax = unname(ext_x$ymax))
    } else {
      bb <- sf::st_bbox(x)
      ext_x <- c(xmin = unname(bb["xmin"]), xmax = unname(bb["xmax"]),
                 ymin = unname(bb["ymin"]), ymax = unname(bb["ymax"]))
    }

    if (inherits(y, c("SpatRaster", "SpatVector"))) {
      ext_y <- terra::ext(y)
      ext_y <- c(xmin = unname(ext_y$xmin), xmax = unname(ext_y$xmax),
                 ymin = unname(ext_y$ymin), ymax = unname(ext_y$ymax))
    } else {
      bb <- sf::st_bbox(y)
      ext_y <- c(xmin = unname(bb["xmin"]), xmax = unname(bb["xmax"]),
                 ymin = unname(bb["ymin"]), ymax = unname(bb["ymax"]))
    }

    # Check if extents overlap
    if (ext_x["xmax"] < ext_y["xmin"] || ext_x["xmin"] > ext_y["xmax"] ||
        ext_x["ymax"] < ext_y["ymin"] || ext_x["ymin"] > ext_y["ymax"]) {
      stop("'", name_x, "' and '", name_y, "' do not overlap spatially",
           call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Require projected CRS in metre units for model-building workflows
#'
#' @param x sf object to validate
#' @param name Character string, parameter name for error messages
#' @return Invisible TRUE if valid, otherwise stops with error
#' @keywords internal
#' @noRd
check_sf_projected_crs <- function(x, name = "x") {
  crs <- sf::st_crs(x)

  if (is.na(crs)) {
    stop(
      "Parameter '", name, "' must have a defined projected CRS in metres/meters ",
      "for meter-based buffer/simplify/mesh generation.",
      call. = FALSE
    )
  }

  if (isTRUE(sf::st_is_longlat(x))) {
    stop(
      "Parameter '", name, "' uses a longitude/latitude CRS, which is not ",
      "supported for meter-based buffer/simplify/mesh generation. ",
      "Transform '", name, "' to a projected CRS in metres/meters first.",
      call. = FALSE
    )
  }

  check_projected_crs_type(crs, name)
  check_crs_units_metres(crs$units_gdal, name)

  invisible(TRUE)
}

#' Require projected CRS in metre units for terra rasters
#'
#' @param x terra SpatRaster to validate
#' @param name Character string, parameter name for error messages
#' @return Invisible TRUE if valid, otherwise stops with error
#' @keywords internal
#' @noRd
check_raster_projected_crs <- function(x, name = "x") {
  crs <- terra::crs(x)

  if (is.null(crs) || is.na(crs) || identical(trimws(crs), "")) {
    stop(
      "Parameter '", name, "' must have a defined projected CRS in metres/meters ",
      "for meter-based buffer/simplify/mesh generation.",
      call. = FALSE
    )
  }

  if (isTRUE(suppressWarnings(terra::is.lonlat(x)))) {
    stop(
      "Parameter '", name, "' uses a longitude/latitude CRS, which is not ",
      "supported for meter-based buffer/simplify/mesh generation. ",
      "Transform '", name, "' to a projected CRS in metres/meters first.",
      call. = FALSE
    )
  }

  sf_crs <- tryCatch(
    sf::st_crs(crs),
    error = function(e) sf::NA_crs_
  )
  check_projected_crs_type(sf_crs, name)
  check_crs_units_metres(sf_crs$units_gdal, name)

  invisible(TRUE)
}

check_projected_crs_type <- function(crs, name) {
  wkt <- tryCatch(
    normalize_crs_wkt(crs$wkt),
    error = function(e) NA_character_
  )

  if (is.na(wkt)) {
    stop(
      "Parameter '", name, "' has a CRS that could not be parsed. ",
      "A projected CRS in metres/meters is required for meter-based ",
      "buffer/simplify/mesh generation.",
      call. = FALSE
    )
  }

  if (!grepl("^[[:space:]]*(PROJCRS|PROJCS)[[:space:]]*\\[", wkt)) {
    stop(
      "Parameter '", name, "' must have a projected CRS in metres/meters ",
      "for meter-based buffer/simplify/mesh generation.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

check_crs_units_metres <- function(units, name) {
  units <- normalize_crs_units(units)

  if (is.na(units)) {
    stop(
      "Parameter '", name, "' has CRS units that could not be determined. ",
      "A projected CRS in metres/meters is required for meter-based ",
      "buffer/simplify/mesh generation.",
      call. = FALSE
    )
  }

  if (!tolower(units) %in% c("metre", "meter")) {
    stop(
      "Parameter '", name, "' must use a projected CRS in metres/meters ",
      "for meter-based buffer/simplify/mesh generation; found CRS units '",
      units, "'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

normalize_crs_wkt <- function(wkt) {
  if (is.null(wkt) || length(wkt) == 0 || all(is.na(wkt))) {
    return(NA_character_)
  }

  wkt <- trimws(as.character(wkt[which(!is.na(wkt))[1]]))
  if (!nzchar(wkt)) {
    return(NA_character_)
  }

  wkt
}

normalize_crs_units <- function(units) {
  if (is.null(units) || length(units) == 0 || all(is.na(units))) {
    return(NA_character_)
  }

  units <- trimws(as.character(units[which(!is.na(units))[1]]))
  if (!nzchar(units)) {
    return(NA_character_)
  }

  units
}

#' Check if two CRS are compatible
#'
#' Determines if two coordinate reference systems are compatible.
#' Handles different CRS object types (character, CRS, crs, etc.).
#'
#' @param crs1 First CRS (character, CRS object, or crs object)
#' @param crs2 Second CRS (character, CRS object, or crs object)
#' @return Logical, TRUE if CRS are compatible
#' @export
#' @examples
#' \dontrun{
#' compatible_crs("EPSG:4326", "EPSG:4326")
#' compatible_crs("+proj=longlat", "+proj=longlat +datum=WGS84")
#' }
compatible_crs <- function(crs1, crs2) {
  # Handle NULL or empty CRS
  if (is.null(crs1) || is.null(crs2)) {
    return(FALSE)
  }

  # Convert to character for comparison
  crs1_str <- crs_to_string(crs1)
  crs2_str <- crs_to_string(crs2)

  # Both empty
  if (crs1_str == "" && crs2_str == "") {
    return(TRUE)
  }

  # One empty, one not
  if (crs1_str == "" || crs2_str == "") {
    return(FALSE)
  }

  # Compare strings (exact match or EPSG code match)
  if (crs1_str == crs2_str) {
    return(TRUE)
  }

  # Try to extract and compare EPSG codes
  epsg1 <- extract_epsg(crs1_str)
  epsg2 <- extract_epsg(crs2_str)

  if (!is.na(epsg1) && !is.na(epsg2)) {
    return(epsg1 == epsg2)
  }

  # If no EPSG match, consider incompatible
  return(FALSE)
}

#' Convert CRS to string representation
#'
#' @param crs CRS object (various types)
#' @return Character string
#' @keywords internal
crs_to_string <- function(crs) {
  if (is.null(crs)) {
    return("")
  }

  if (is.character(crs)) {
    return(trimws(crs))
  }

  # Handle sf crs object
  if (inherits(crs, "crs")) {
    if (!is.na(crs$input)) {
      return(trimws(crs$input))
    }
    if (!is.na(crs$wkt)) {
      return(trimws(crs$wkt))
    }
    return("")
  }

  # Try to coerce to character
  tryCatch({
    return(trimws(as.character(crs)))
  }, error = function(e) {
    return("")
  })
}

#' Extract EPSG code from CRS string
#'
#' @param crs_str Character string with CRS definition
#' @return Integer EPSG code or NA
#' @keywords internal
extract_epsg <- function(crs_str) {
  if (is.null(crs_str) || !is.character(crs_str) || nchar(crs_str) == 0) {
    return(NA_integer_)
  }

  # Try EPSG:XXXX format
  epsg_match <- regexpr("EPSG:[0-9]+", crs_str, ignore.case = TRUE)
  if (epsg_match > 0) {
    epsg_str <- regmatches(crs_str, epsg_match)
    epsg_num <- as.integer(sub("EPSG:", "", epsg_str, ignore.case = TRUE))
    return(epsg_num)
  }

  # Try +init=epsg:XXXX format
  init_match <- regexpr("\\+init=epsg:[0-9]+", crs_str, ignore.case = TRUE)
  if (init_match > 0) {
    init_str <- regmatches(crs_str, init_match)
    epsg_num <- as.integer(sub("\\+init=epsg:", "", init_str,
                                ignore.case = TRUE))
    return(epsg_num)
  }

  return(NA_integer_)
}
