#' Model Output I/O Functions
#' 
#' This module provides functions for reading SHUD model output files.
#' Functions support both legacy and modern formats with terra integration.
#' 
#' @name output-io
NULL

# =============================================================================
# Modern Output Reading Functions
# =============================================================================

#' Read SHUD model output file
#' 
#' Reads binary output files from SHUD model with support for both legacy
#' and modern formats. Returns time series data as xts objects.
#' 
#' @param keyword Character. Keyword of the output file (e.g., "eleygw", "rivqdown")
#' @param file Character. Full path to output file (optional if keyword provided)
#' @param path Character. Path to output directory (default: from environment)
#' @param version Numeric. Output file version (1.0 or 2.0+)
#' @param format Character. Output format ("modern" uses terra, "legacy" uses raster)
#' @param ascii Logical. If TRUE, convert binary to ASCII format (.csv)
#' @return xts object containing time series data
#' @family output-io
#' @export
#' @examples
#' \dontrun{
#' # Read groundwater elevation data
#' gw_data <- read_output("eleygw")
#' 
#' # Read with specific file path
#' gw_data <- read_output(file = "model.eleygw.dat")
#' 
#' # Use modern format with terra
#' gw_data <- read_output("eleygw", format = "modern")
#' }
read_output <- function(keyword = NULL, 
                       file = NULL,
                       path = get_output_path(),
                       version = get_model_version(),
                       format = c("modern", "legacy"),
                       ascii = FALSE) {
  
  format <- match.arg(format)
  msg <- "read_output::"
  
  # Determine file path
  if (is.null(file)) {
    if (is.null(keyword)) {
      stop("Either 'keyword' or 'file' must be provided", call. = FALSE)
    }
    project_name <- get_project_name()
    file <- file.path(path, paste0(project_name, ".", keyword, ".dat"))
  }
  
  # Validate file
  if (!file.exists(file)) {
    stop("Output file does not exist: ", file, call. = FALSE)
  }
  
  if (!file.access(file, mode = 4) == 0) {
    stop("Output file is not readable: ", file, call. = FALSE)
  }
  
  # Read binary data
  fid <- file(file, "rb")
  on.exit(close(fid))
  
  tmp <- readBin(fid, what = numeric(), n = 1e9, size = 8)
  
  # Parse header based on version
  if (version > 1.0) {
    # Modern version with extended header
    # header_data <- tmp[1:128]  # 1024 chars = 128 doubles (reserved for future use)
    start_time <- tmp[129]
    ncol_data <- tmp[130]
    col_indices <- tmp[130 + (1:ncol_data)]
    data_values <- tmp[-1 * (1:(130 + ncol_data))]
  } else {
    # Legacy version
    data_values <- tmp[-1 * (1:2)]
    ncol_data <- tmp[1]
    start_time <- tmp[2]
    col_indices <- 1:ncol_data
  }
  
  # Reshape data
  nrow_total <- length(data_values) / (ncol_data + 1)
  nrow_data <- round(nrow_total)
  
  if (nrow_data != nrow_total) {
    warning(msg, "File may be incomplete. Expected ", nrow_total, " rows")
  }
  
  data_matrix <- t(matrix(data_values[1:(nrow_data * (ncol_data + 1))], 
                         nrow = ncol_data + 1))
  
  # Convert to ASCII if requested
  if (ascii) {
    ascii_file <- paste0(substr(file, 1, nchar(file) - 3), "csv")
    write(paste(ncol_data, start_time), file = ascii_file, append = FALSE)
    
    temp_data <- data_matrix
    colnames(temp_data) <- c("T_MIN", paste0("X", 1:ncol_data))
    
    suppressWarnings(
      utils::write.table(temp_data, file = ascii_file, append = TRUE,
                        col.names = TRUE, row.names = FALSE, sep = "\t")
    )
    message(msg, "ASCII file written to: ", ascii_file)
  }
  
  # Create time series
  time_minutes <- data_matrix[, 1]
  time_seconds <- time_minutes * 60
  
  start_date <- as.POSIXct(as.character(start_time), format = "%Y%m%d", tz = "UTC")
  timestamps <- start_date + time_seconds
  
  # Create xts object
  if (ncol_data <= 1) {
    ts_data <- xts::as.xts(cbind(data_matrix[, -1]), order.by = timestamps)
  } else {
    ts_data <- xts::as.xts(rbind(data_matrix[, -1]), order.by = timestamps)
  }
  
  colnames(ts_data) <- paste(col_indices)
  
  # Add attributes for format tracking
  attr(ts_data, "format") <- format
  attr(ts_data, "version") <- version
  attr(ts_data, "keyword") <- keyword
  
  return(ts_data)
}

#' Read output file header information
#' 
#' Reads header information from SHUD output files for inspection.
#' 
#' @param keyword Character. Keyword of the output file
#' @param file Character. Full path to output file (optional)
#' @param path Character. Path to output directory
#' @return Character string containing header information
#' @family output-io
#' @export
#' @examples
#' \dontrun{
#' header_info <- read_output_header("eleygw")
#' }
read_output_header <- function(keyword = NULL,
                              file = NULL,
                              path = get_output_path()) {
  
  # Determine file path
  if (is.null(file)) {
    if (is.null(keyword)) {
      stop("Either 'keyword' or 'file' must be provided", call. = FALSE)
    }
    project_name <- get_project_name()
    file <- file.path(path, paste0(project_name, ".", keyword, ".dat"))
  }
  
  # Validate file
  if (!file.exists(file)) {
    stop("Output file does not exist: ", file, call. = FALSE)
  }
  
  if (!file.access(file, mode = 4) == 0) {
    stop("Output file is not readable: ", file, call. = FALSE)
  }
  
  # Read header
  fid <- file(file, "rb")
  on.exit(close(fid))
  
  header_chars <- readBin(fid, what = character(), n = 1024, size = 1)
  header_string <- paste(header_chars, collapse = "")
  
  return(header_string)
}

#' Load multiple SHUD output files
#' 
#' Loads multiple output variables and returns them as a named list.
#' Supports both modern and legacy formats.
#' 
#' @param variables Character vector. Output variable keywords to load
#' @param format Character. Output format ("modern" or "legacy")
#' @param save_rds Character. Path to save RDS file (optional)
#' @param version Numeric. Output file version
#' @return Named list of xts objects
#' @family output-io
#' @export
#' @examples
#' \dontrun{
#' # Load common variables
#' vars <- c("eleygw", "eleysurf", "rivqdown")
#' data_list <- load_output_data(vars)
#' 
#' # Use modern format
#' data_list <- load_output_data(vars, format = "modern")
#' }
load_output_data <- function(variables = get_default_variables(),
                            format = c("modern", "legacy"),
                            save_rds = NULL,
                            version = get_model_version()) {
  
  format <- match.arg(format)
  msg <- "load_output_data::"
  
  variables <- tolower(variables)
  n_vars <- length(variables)
  
  message(msg, "Loading ", n_vars, " variables with ", format, " format")
  
  # Create variable metadata
  var_types <- substr(variables, 4, 4)
  var_metadata <- data.frame(
    variable = variables,
    type = var_types,
    unit = "",
    stringsAsFactors = FALSE
  )
  
  # Assign units based on variable type
  var_metadata$unit[var_types == "y"] <- "m"      # elevation/depth
  var_metadata$unit[var_types == "v"] <- "m/d"    # velocity/flux
  var_metadata$unit[var_types == "q"] <- "m3/d"   # discharge
  
  # Load data
  result_list <- list()
  
  for (i in seq_along(variables)) {
    var_name <- variables[i]
    message(msg, i, "/", n_vars, "\t", var_name)
    
    tryCatch({
      ts_data <- read_output(keyword = var_name, format = format, version = version)
      result_list[[i]] <- ts_data
    }, error = function(e) {
      warning("Failed to load variable '", var_name, "': ", e$message)
      result_list[[i]] <- NULL
    })
  }
  
  names(result_list) <- variables
  
  # Remove NULL entries (failed loads)
  result_list <- result_list[!sapply(result_list, is.null)]
  
  # Save to RDS if requested
  if (!is.null(save_rds)) {
    saveRDS(result_list, file = save_rds)
    message(msg, "Data saved to: ", save_rds)
  }
  
  # Add metadata as attribute
  attr(result_list, "metadata") <- var_metadata
  attr(result_list, "format") <- format
  
  return(result_list)
}

# =============================================================================
# Utility Functions
# =============================================================================

#' Get default output variables
#' 
#' Returns a vector of commonly used SHUD output variables.
#' 
#' @return Character vector of variable names
#' @family output-io
#' @export
get_default_variables <- function() {
  c(
    paste0("eley", c("surf", "unsat", "gw", "snow")),
    paste0("elev", c("prcp", "infil", "rech")),
    paste0("elev", c("etp", "eta", "etev", "ettr", "etic")),
    paste0("rivq", c("down", "sub", "surf")),
    paste0("rivy", "stage")
  )
}

#' Get output path from environment
#' 
#' @return Character string with output path
#' @keywords internal
get_output_path <- function() {
  if (exists(".shud", envir = .GlobalEnv)) {
    tryCatch({
      p = get("outpath", envir = .shud)
      if(is.character(p)) as.character(p) else getwd()
    }, error = function(e) {
      getwd()  # fallback to current directory
    })
  } else {
    getwd()
  }
}

#' Get model version from environment
#' 
#' @return Numeric model version
#' @keywords internal
get_model_version <- function() {
  if (exists(".shud", envir = .GlobalEnv)) {
    tryCatch({
      v = get("version", envir = .shud)
      if(is.numeric(v)) as.numeric(v) else 2.0
    }, error = function(e) {
      2.0  # default to modern version
    })
  } else {
    2.0
  }
}

#' Get project name from environment
#' 
#' @return Character string with project name
#' @keywords internal
get_project_name <- function() {
  if (exists(".shud", envir = .GlobalEnv)) {
    tryCatch({
      p = get("PRJNAME", envir = .shud)
      if(is.character(p)) as.character(p) else "model"
    }, error = function(e) {
      "model"  # default project name
    })
  } else {
    "model"
  }
}