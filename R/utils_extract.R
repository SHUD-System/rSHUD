#' Normalize terra::extract() output
#'
#' terra::extract() may return a vector, a single-column data frame, or an
#' ID-prefixed data frame depending on the terra version and input type. This
#' helper removes that shape variability for internal callers.
#'
#' @param values Object returned by terra::extract().
#' @param n_layers Expected number of raster layers.
#' @return Atomic vector for single-layer rasters, or a data.frame for
#'   multi-layer rasters.
#' @noRd
.unwrap_terra_extract <- function(values, n_layers = 1L) {
  if (!is.matrix(values) && !is.data.frame(values)) {
    return(values)
  }

  values_df <- as.data.frame(values)

  if (ncol(values_df) == n_layers + 1L) {
    values_df <- values_df[-1]
  } else if ("ID" %in% names(values_df) && ncol(values_df) > n_layers) {
    values_df <- values_df[setdiff(names(values_df), "ID")]
  }

  if (ncol(values_df) != n_layers) {
    stop(
      "Unexpected terra::extract() return shape: expected ", n_layers,
      " value column(s), got ", ncol(values_df),
      call. = FALSE
    )
  }

  if (n_layers == 1L) {
    return(values_df[[1]])
  }

  values_df
}

#' Extract SpatRaster values with stable output shape
#'
#' Internal wrapper around terra::extract() that normalizes return shapes
#' across terra versions and coordinate input types.
#'
#' @param raster_obj A terra SpatRaster object.
#' @param coords Coordinates or vector geometry accepted by terra::extract().
#' @param method Extraction method passed to terra::extract().
#' @param require_single_layer If TRUE, error on multilayer rasters.
#' @return Numeric vector for single-layer rasters, or a data.frame for
#'   multi-layer rasters.
#' @noRd
.extract_spatraster_values <- function(raster_obj,
                                       coords,
                                       method = "simple",
                                       require_single_layer = TRUE) {
  if (!inherits(raster_obj, "SpatRaster")) {
    stop("Parameter 'raster_obj' must be a terra SpatRaster object", call. = FALSE)
  }

  n_layers <- terra::nlyr(raster_obj)
  if (require_single_layer && n_layers != 1L) {
    stop(
      "Internal extract helper currently supports single-layer SpatRaster objects only, got ",
      n_layers, " layers",
      call. = FALSE
    )
  }

  values <- terra::extract(raster_obj, coords, method = method)
  values <- .unwrap_terra_extract(values, n_layers = n_layers)

  if (require_single_layer) {
    return(as.numeric(values))
  }

  values
}
