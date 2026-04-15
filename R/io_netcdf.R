#' NetCDF I/O Functions
#' 
#' This module provides functions for reading NetCDF files with support for
#' both legacy raster and modern terra formats. Functions maintain ncdf4
#' interface compatibility while adding terra integration.
#' 
#' @name netcdf-io
NULL

# =============================================================================
# Modern NetCDF Reading Functions
# =============================================================================

#' Read NetCDF time dimension
#' 
#' Extracts time information from NetCDF files and converts to POSIXct objects.
#' Supports various time formats and units commonly used in climate data.
#' 
#' @param ncid NetCDF connection object from ncdf4::nc_open()
#' @return POSIXct vector of time stamps
#' @family netcdf-io
#' @export
#' @examples
#' \dontrun{
#' library(ncdf4)
#' nc <- nc_open("climate_data.nc")
#' times <- read_nc_time(nc)
#' nc_close(nc)
#' }
read_nc_time <- function(ncid) {
  if (missing(ncid) || is.null(ncid)) {
    stop("Parameter 'ncid' (NetCDF connection) is required", call. = FALSE)
  }
  
  # Find time dimension
  nc_dims <- names(ncid$dim)
  time_vars <- c("time", "Time", "datetime", "Datetime", "date", "Date")
  time_var <- nc_dims[which(nc_dims %in% time_vars)[1]]
  
  if (length(time_var) == 0) {
    stop("Could not identify time variable in NetCDF file", call. = FALSE)
  }
  
  # Get time values and attributes
  times <- ncdf4::ncvar_get(ncid, time_var)
  time_att <- ncdf4::ncatt_get(ncid, time_var)
  
  # Parse time units
  time_def <- strsplit(gsub("^ *", "", time_att$units), " ")[[1]]
  time_unit <- time_def[1]
  
  # Handle timezone
  if (length(time_def) < 5) {
    tz <- "UTC"
  } else {
    tz <- time_def[5]
  }
  
  # Validate and parse start time
  time_start_parts <- strsplit(time_def[4], ":")[[1]]
  if (length(time_start_parts) != 3 || 
      any(as.numeric(time_start_parts[1]) > 24, 
          as.numeric(time_start_parts[2]) > 60, 
          as.numeric(time_start_parts[3]) > 60, 
          as.numeric(time_start_parts) < 0)) {
    warning("Invalid start time format. Using 00:00:00")
    time_def[4] <- "00:00:00"
  }
  
  time_start <- lubridate::ymd_hms(paste(time_def[3], time_def[4]), tz = tz)
  
  # Convert time units to lubridate function
  time_func <- switch(tolower(time_unit),
    "seconds" = , "second" = , "sec" = lubridate::seconds,
    "minutes" = , "minute" = , "min" = lubridate::minutes,
    "hours" = , "hour" = , "h" = lubridate::hours,
    "days" = , "day" = , "d" = lubridate::days,
    "months" = , "month" = , "m" = function(x) lubridate::period(num = x, units = "month"),
    "years" = , "year" = , "yr" = lubridate::years,
    NULL
  )
  
  if (is.null(time_func)) {
    stop("Could not understand time unit format: ", time_unit, call. = FALSE)
  }
  
  return(time_start + time_func(times))
}

#' Read NetCDF data with spatial subsetting
#' 
#' Reads NetCDF data with optional spatial and variable subsetting.
#' Supports both legacy raster and modern terra output formats.
#' 
#' @param ncid NetCDF connection object from ncdf4::nc_open()
#' @param variables Character vector of variable names to read (NULL = all)
#' @param extent Numeric vector c(xmin, xmax, ymin, ymax) for spatial subset
#' @param format Character. Output format ("modern" uses terra, "legacy" uses raster)
#' @return List containing coordinates, data array, and time information
#' @family netcdf-io
#' @export
#' @examples
#' \dontrun{
#' library(ncdf4)
#' nc <- nc_open("climate_data.nc")
#' 
#' # Read all variables
#' data <- read_nc_data(nc)
#' 
#' # Read specific variables with spatial subset
#' extent <- c(-180, -120, 30, 60)  # Western North America
#' data <- read_nc_data(nc, variables = c("temp", "precip"), extent = extent)
#' 
#' # Use modern terra format
#' data <- read_nc_data(nc, format = "modern")
#' 
#' nc_close(nc)
#' }
read_nc_data <- function(ncid, 
                        variables = NULL, 
                        extent = NULL,
                        format = c("modern", "legacy")) {
  
  format <- match.arg(format)
  msg <- "read_nc_data::"
  
  if (missing(ncid) || is.null(ncid)) {
    stop("Parameter 'ncid' (NetCDF connection) is required", call. = FALSE)
  }
  
  # Get available variables
  all_vars <- names(ncid$var)
  n_vars <- length(all_vars)
  
  # Filter out time bounds variables
  all_vars <- all_vars[!grepl("time_bnds?", all_vars)]
  
  # Determine which variables to read
  if (is.null(variables)) {
    # Read all variables except time bounds
    read_vars <- all_vars
  } else if (is.character(variables)) {
    # Validate requested variables
    missing_vars <- variables[!variables %in% all_vars]
    if (length(missing_vars) > 0) {
      stop("Variables not found in NetCDF: ", paste(missing_vars, collapse = ", "),
           call. = FALSE)
    }
    read_vars <- variables
  } else if (is.numeric(variables)) {
    # Convert indices to variable names
    if (max(variables) > n_vars || min(variables) < 1) {
      stop("Variable indices out of range: 1 to ", n_vars, call. = FALSE)
    }
    read_vars <- all_vars[variables]
    message(msg, "Reading variables: ", paste(read_vars, collapse = ", "))
  } else {
    stop("Invalid format for 'variables' parameter", call. = FALSE)
  }
  
  # Find spatial dimensions
  nc_dims <- names(ncid$dim)
  lon_var <- nc_dims[grepl("lon", tolower(nc_dims))][1]
  lat_var <- nc_dims[grepl("lat", tolower(nc_dims))][1]
  
  if (is.na(lon_var) || is.na(lat_var)) {
    stop("Could not identify longitude/latitude dimensions", call. = FALSE)
  }
  
  # Get coordinate arrays
  lon <- ncdf4::ncvar_get(ncid, varid = lon_var)
  lat <- ncdf4::ncvar_get(ncid, varid = lat_var)
  
  # Calculate grid resolution
  dx <- mean(diff(lon))
  dy <- mean(diff(lat))
  
  # Determine spatial extent
  if (is.null(extent)) {
    extent <- c(min(lon), max(lon), min(lat), max(lat))
  } else {
    # Validate extent against data bounds
    data_bounds <- c(min(lon) - dx/2, max(lon) + dx/2, 
                    min(lat) - dy/2, max(lat) + dy/2)
    
    if (extent[1] < data_bounds[1] || extent[2] > data_bounds[2] ||
        extent[3] < data_bounds[3] || extent[4] > data_bounds[4]) {
      warning("Requested extent extends beyond data bounds")
      message("Requested: ", paste(extent, collapse = ", "))
      message("Available: ", paste(data_bounds, collapse = ", "))
    }
  }
  
  # Find indices for spatial subset
  x_indices <- which(abs(lon - extent[1]) <= dx/2):which(abs(lon - extent[2]) <= dx/2)
  y_indices <- which(abs(lat - extent[3]) <= dy/2):which(abs(lat - extent[4]) <= dy/2)
  
  x_indices <- x_indices[!is.na(x_indices)]
  y_indices <- y_indices[!is.na(y_indices)]
  
  if (length(x_indices) == 0 || length(y_indices) == 0) {
    stop("No data found within specified extent", call. = FALSE)
  }
  
  nx <- length(x_indices)
  ny <- length(y_indices)
  
  x_coords <- lon[x_indices]
  y_coords <- lat[y_indices]
  
  # Prepare data array
  data_array <- array(0, 
                     dim = c(nx, ny, length(read_vars)),
                     dimnames = list(x_coords, y_coords, read_vars))
  
  # Read data for each variable
  start_indices <- c(min(x_indices), min(y_indices), 1)
  count_indices <- c(nx, ny, 1)
  
  for (i in seq_along(read_vars)) {
    var_name <- read_vars[i]
    tryCatch({
      data_array[, , i] <- ncdf4::ncvar_get(ncid, var_name, 
                                           start = start_indices, 
                                           count = count_indices)
    }, error = function(e) {
      warning("Failed to read variable '", var_name, "': ", e$message)
      data_array[, , i] <- NA
    })
  }
  
  # Get time information
  time_data <- tryCatch({
    read_nc_time(ncid)
  }, error = function(e) {
    warning("Could not read time information: ", e$message)
    NULL
  })
  
  # Create result list
  result <- list(
    x = x_coords,
    y = y_coords,
    data = data_array,
    time = time_data,
    format = format,
    extent = extent,
    variables = read_vars
  )
  
  return(result)
}

#' Convert NetCDF data to raster/terra format
#' 
#' Converts NetCDF data arrays to raster (legacy) or terra (modern) format.
#' Handles both single layers and multi-layer datasets.
#' 
#' @param nc_data List returned from read_nc_data()
#' @param x Numeric vector of x coordinates (alternative input)
#' @param y Numeric vector of y coordinates (alternative input)
#' @param data_array Array of data values (alternative input)
#' @param resolution Numeric vector c(dx, dy) for grid resolution
#' @param flip Logical. Whether to flip the array vertically
#' @param format Character. Output format ("modern" or "legacy")
#' @return SpatRaster (modern) or RasterStack/RasterLayer (legacy)
#' @family netcdf-io
#' @export
#' @examples
#' \dontrun{
#' library(ncdf4)
#' nc <- nc_open("climate_data.nc")
#' nc_data <- read_nc_data(nc)
#' 
#' # Convert to terra format (modern)
#' raster_data <- nc_to_raster(nc_data, format = "modern")
#' 
#' # Convert to raster format (legacy)
#' raster_data <- nc_to_raster(nc_data, format = "legacy")
#' 
#' nc_close(nc)
#' }
nc_to_raster <- function(nc_data = NULL,
                        x = NULL, 
                        y = NULL, 
                        data_array = NULL,
                        resolution = NULL,
                        flip = TRUE,
                        format = c("modern", "legacy")) {
  
  format <- match.arg(format)
  
  # Handle input data
  if (!is.null(nc_data)) {
    x <- nc_data$x
    y <- nc_data$y
    data_array <- nc_data$data
  } else if (is.null(x) || is.null(y) || is.null(data_array)) {
    stop("Either 'nc_data' or 'x', 'y', 'data_array' must be provided", call. = FALSE)
  }
  
  # Get dimensions
  dims <- dim(data_array)
  n_dims <- length(dims)
  
  # Calculate resolution if not provided
  if (is.null(resolution)) {
    if (length(x) == 1 || length(y) == 1) {
      stop("Resolution must be provided for single coordinate values", call. = FALSE)
    }
    dx <- abs(mean(diff(x)))
    dy <- abs(mean(diff(y)))
  } else {
    dx <- resolution[1]
    dy <- if (length(resolution) == 1) resolution[1] else resolution[2]
  }
  
  nx <- length(x)
  ny <- length(y)
  
  # Create spatial extent
  extent_vals <- c(min(x) - dx/2, max(x) + dx/2, 
                  min(y) - dy/2, max(y) + dy/2)
  
  if (format == "modern") {
    # Use terra for modern format
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("Package 'terra' is required for modern format", call. = FALSE)
    }
    
    if (n_dims > 2) {
      # Multiple layers
      raster_list <- list()
      for (i in 1:dims[3]) {
        layer_data <- data_array[, , i]
        if (flip) {
          layer_data <- layer_data[, ny:1]
        }
        
        r <- terra::rast(nrows = ny, ncols = nx, extent = extent_vals)
        terra::values(r) <- as.vector(t(layer_data))
        raster_list[[i]] <- r
      }
      result <- terra::rast(raster_list)
      if (!is.null(dimnames(data_array)[[3]])) {
        names(result) <- dimnames(data_array)[[3]]
      }
    } else {
      # Single layer
      layer_data <- data_array
      if (flip) {
        layer_data <- layer_data[, ny:1]
      }
      
      result <- terra::rast(nrows = ny, ncols = nx, extent = extent_vals)
      terra::values(result) <- as.vector(t(layer_data))
    }
    
  } else {
    # Use raster for legacy format
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("Package 'raster' is required for legacy format", call. = FALSE)
    }
    
    if (n_dims > 2) {
      # Multiple layers
      raster_list <- list()
      for (i in 1:dims[3]) {
        layer_data <- data_array[, , i]
        
        r <- raster::raster(ncols = nx, nrows = ny)
        raster::extent(r) <- extent_vals
        
        if (flip) {
          idx <- ny:1
        } else {
          idx <- 1:ny
        }
        
        raster_list[[i]] <- raster::setValues(r, t(layer_data[, idx]))
      }
      result <- raster::stack(raster_list)
      if (!is.null(dimnames(data_array)[[3]])) {
        names(result) <- dimnames(data_array)[[3]]
      }
    } else {
      # Single layer
      r <- raster::raster(ncols = nx, nrows = ny)
      raster::extent(r) <- extent_vals
      
      if (flip) {
        idx <- ny:1
      } else {
        idx <- 1:ny
      }
      
      result <- raster::setValues(r, t(data_array[, idx]))
    }
  }
  
  return(result)
}

# =============================================================================
# Utility Functions
# =============================================================================

#' Open NetCDF file with error handling
#' 
#' Wrapper around ncdf4::nc_open with improved error handling.
#' 
#' @param file Character. Path to NetCDF file
#' @return NetCDF connection object
#' @family netcdf-io
#' @export
open_nc_file <- function(file) {
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("NetCDF file does not exist: ", file, call. = FALSE)
  }
  
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required for NetCDF operations", call. = FALSE)
  }
  
  tryCatch({
    ncdf4::nc_open(file)
  }, error = function(e) {
    stop("Failed to open NetCDF file: ", e$message, call. = FALSE)
  })
}

#' Get NetCDF file information
#' 
#' Extracts basic information about NetCDF file structure.
#' 
#' @param ncid NetCDF connection object
#' @return List with file information
#' @family netcdf-io
#' @export
get_nc_info <- function(ncid) {
  if (missing(ncid) || is.null(ncid)) {
    stop("Parameter 'ncid' (NetCDF connection) is required", call. = FALSE)
  }
  
  info <- list(
    filename = ncid$filename,
    dimensions = names(ncid$dim),
    variables = names(ncid$var),
    n_dimensions = ncid$ndims,
    n_variables = ncid$nvars,
    format = ncid$format
  )
  
  # Add dimension sizes
  dim_sizes <- sapply(ncid$dim, function(d) d$len)
  names(dim_sizes) <- names(ncid$dim)
  info$dimension_sizes <- dim_sizes
  
  # Add variable information
  var_info <- lapply(ncid$var, function(v) {
    list(
      name = v$name,
      units = v$units,
      longname = v$longname,
      dimensions = names(v$dim),
      size = v$size
    )
  })
  info$variable_info <- var_info
  
  return(info)
}