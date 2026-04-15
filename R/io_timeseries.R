#' Time Series I/O Functions
#' 
#' This module provides functions for reading and writing time series data
#' in SHUD format. All functions use modern snake_case naming and provide
#' comprehensive error handling.
#' 
#' @name timeseries-io
NULL

# =============================================================================
# Time Series Reading Functions
# =============================================================================

#' Read time series data file
#' 
#' Reads a time series file with header containing (nrow, ncol, timestring).
#' This function handles SHUD-specific time series format with automatic
#' time parsing and xts object creation.
#' 
#' @param file Character. Full path to time series file
#' @return List of xts objects containing time series data
#' @family timeseries-io
#' @export
#' @examples
#' \dontrun{
#' ts_data <- read_tsd("timeseries.tsd")
#' lai_data <- read_lai("model.ts.lai")
#' }
read_tsd <- function(file) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("Time-series file does not exist: ", file, call. = FALSE)
  }
  
  if (!file.access(file, mode = 4) == 0) {
    stop("Time-series file is not readable: ", file, call. = FALSE)
  }
  
  text <- readLines(file)
  r0 <- 1
  nrow <- length(text)
  xl <- list()
  
  for (i in 1:100) {
    header <- utils::read.table(text = text[r0])
    nr <- as.numeric(header[1])
    t0 <- header[3]
    tmp <- as.matrix(utils::read.table(text = text[0:nr + 1 + r0], header = TRUE))
    tsd <- xts::as.xts(tmp[, -1], 
                       order.by = as.POSIXct(paste(t0), format = "%Y%m%d") + tmp[, 1] * 86400)
    xl[[i]] <- tsd
    r0 <- r0 + nr + 2
    
    if (r0 + 1 > nrow) {
      break
    }
  }
  
  return(xl)
}

#' Read LAI time series file
#' 
#' Reads a LAI (Leaf Area Index) time series file (.ts.lai) and returns
#' properly named time series components.
#' 
#' @param file Character. Full path to LAI time series file
#' @return List of time series data with LAI and RL (Roughness Length) components
#' @family timeseries-io
#' @export
#' @examples
#' \dontrun{
#' lai_data <- read_lai("model.ts.lai")
#' plot(lai_data$LAI)
#' plot(lai_data$RL)
#' }
read_lai <- function(file = shud.filein()["ts.lai"]) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  x <- read_tsd(file)
  names(x) <- c("LAI", "RL")
  return(x)
}

# =============================================================================
# Time Series Writing Functions
# =============================================================================

#' Write time series data to file
#' 
#' Writes xts time series data to a SHUD time series file (.tsd) with proper
#' header formatting. Supports different time units and custom headers.
#' 
#' @param x xts object containing time series data
#' @param file Character. Full path to output .tsd file
#' @param dt_units Character. Time interval units ("days", "minutes", "seconds")
#' @param append Logical. Whether to append to existing file (default: FALSE)
#' @param quiet Logical. Whether to suppress messages (default: FALSE)
#' @param header Numeric vector. Custom header (default: auto-generated)
#' @family timeseries-io
#' @export
#' @examples
#' \dontrun{
#' # Create sample time series
#' dates <- seq(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
#' ts_data <- xts::xts(runif(length(dates)), order.by = dates)
#' 
#' # Write to file
#' write_tsd(ts_data, "output.tsd", dt_units = "days")
#' }
write_tsd <- function(x, file, dt_units = "days", append = FALSE, 
                      quiet = FALSE, header = NULL) {
  msg <- "write_tsd::"
  
  if (missing(x)) {
    stop("Parameter 'x' (time series data) is required", call. = FALSE)
  }
  
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!inherits(x, "xts")) {
    stop("Parameter 'x' must be an xts object", call. = FALSE)
  }
  
  # Ensure directory exists
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  
  mat <- as.matrix(rbind(x))
  nr <- nrow(x)
  nc <- ncol(x)
  
  if (!quiet) {
    message(msg, "Writing to file ", file)
  }
  
  tt <- stats::time(x)
  
  # Determine time step in seconds
  if (grepl("minute", dt_units)) {
    dt_sec <- 60  # 60 sec
  } else if (grepl("second", dt_units)) {
    dt_sec <- 1
  } else {
    dt_sec <- 3600 * 24  # days
  }
  
  time_tag <- as.numeric(difftime(tt, tt[1], units = "secs")) / dt_sec
  
  if (is.null(header)) {
    t0 <- format(time(x)[1], "%Y%m%d")
    t1 <- format(time(x)[nrow(x)], "%Y%m%d")
    header <- c(nr, nc + 1, t0, t1, dt_sec)
  }
  
  dd <- data.frame("Time_interval" = time_tag, mat)
  write(header, file = file, ncolumns = length(header), append = append, sep = "\t")
  write(colnames(dd), file = file, ncolumns = nc + 1, append = TRUE, sep = "\t")
  suppressWarnings(utils::write.table(dd, file = file, append = TRUE, sep = "\t",
                                      col.names = FALSE, row.names = FALSE))
}

#' Write xts data to file (simple format)
#' 
#' Writes xts time series data to file in simple format with TIME column.
#' This is a simpler alternative to write_tsd for basic time series output.
#' 
#' @param x xts object containing time series data
#' @param file Character. Full path to output file
#' @param append Logical. Whether to append to existing file (default: FALSE)
#' @family timeseries-io
#' @export
#' @examples
#' \dontrun{
#' # Create sample time series
#' dates <- seq(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
#' ts_data <- xts::xts(runif(length(dates)), order.by = dates)
#' 
#' # Write to simple format
#' write_xts(ts_data, "output.txt")
#' }
write_xts <- function(x, file, append = FALSE) {
  if (missing(x)) {
    stop("Parameter 'x' (time series data) is required", call. = FALSE)
  }
  
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!inherits(x, "xts")) {
    stop("Parameter 'x' must be an xts object", call. = FALSE)
  }
  
  tt <- stats::time(x)
  md <- data.frame("TIME" = tt, x)
  utils::write.table(md, file = file, append = append, sep = "\t",
                     col.names = TRUE, row.names = FALSE)
}

# =============================================================================
# Time Series Utility Functions
# =============================================================================

#' Convert time series to daily aggregation
#' 
#' Aggregates time series data to daily values using specified function.
#' 
#' @param x xts object with time series data
#' @param FUN Function to use for aggregation (default: mean)
#' @return xts object with daily aggregated data
#' @family timeseries-io
#' @export
#' @examples
#' \dontrun{
#' # Create hourly data
#' hourly_data <- xts::xts(runif(24*7), 
#'                        order.by = seq(Sys.time(), by = "hour", length.out = 24*7))
#' 
#' # Convert to daily means
#' daily_data <- ts_to_daily(hourly_data)
#' }
ts_to_daily <- function(x, FUN = mean) {
  if (!inherits(x, "xts")) {
    stop("Parameter 'x' must be an xts object", call. = FALSE)
  }
  
  xts::apply.daily(x, FUN)
}

#' Convert time series to data frame
#' 
#' Converts xts time series to data frame with time column.
#' 
#' @param x xts object with time series data
#' @param time_col Character. Name for time column (default: "Time")
#' @return data.frame with time and data columns
#' @family timeseries-io
#' @export
#' @examples
#' \dontrun{
#' dates <- seq(as.Date("2020-01-01"), as.Date("2020-01-10"), by = "day")
#' ts_data <- xts::xts(runif(length(dates)), order.by = dates)
#' df_data <- ts_to_df(ts_data)
#' }
ts_to_df <- function(x, time_col = "Time") {
  if (!inherits(x, "xts")) {
    stop("Parameter 'x' must be an xts object", call. = FALSE)
  }
  
  df <- data.frame(stats::time(x), as.data.frame(x))
  names(df)[1] <- time_col
  return(df)
}

# =============================================================================
# Deprecated Function Aliases (for backward compatibility)
# =============================================================================

#' @rdname read_tsd
#' @export
#' @section Deprecated:
#' \code{read.tsd()} is deprecated. Use \code{read_tsd()} instead.
read.tsd <- function(file) {
  .Deprecated("read_tsd", package = "rSHUD",
              msg = "read.tsd() is deprecated. Use read_tsd() instead.")
  read_tsd(file)
}

#' @rdname read_lai
#' @export
#' @section Deprecated:
#' \code{readlai()} is deprecated. Use \code{read_lai()} instead.
readlai <- function(file = shud.filein()["ts.lai"]) {
  .Deprecated("read_lai", package = "rSHUD",
              msg = "readlai() is deprecated. Use read_lai() instead.")
  read_lai(file)
}

#' @rdname write_tsd
#' @export
#' @section Deprecated:
#' \code{write.tsd()} is deprecated. Use \code{write_tsd()} instead.
write.tsd <- function(x, file, dt_units = "days", append = FALSE, 
                      quiet = FALSE, header = NULL) {
  .Deprecated("write_tsd", package = "rSHUD",
              msg = "write.tsd() is deprecated. Use write_tsd() instead.")
  write_tsd(x, file, dt_units, append, quiet, header)
}

#' @rdname write_xts
#' @export
#' @section Deprecated:
#' \code{write.xts()} is deprecated. Use \code{write_xts()} instead.
write.xts <- function(x, file, append = FALSE) {
  .Deprecated("write_xts", package = "rSHUD",
              msg = "write.xts() is deprecated. Use write_xts() instead.")
  write_xts(x, file, append)
}

#' @rdname ts_to_daily
#' @export
#' @section Deprecated:
#' \code{ts2Daily()} is deprecated. Use \code{ts_to_daily()} instead.
ts2Daily <- function(x, FUN = mean) {
  .Deprecated("ts_to_daily", package = "rSHUD",
              msg = "ts2Daily() is deprecated. Use ts_to_daily() instead.")
  ts_to_daily(x, FUN)
}

#' @rdname ts_to_df
#' @export
#' @section Deprecated:
#' \code{ts2df()} is deprecated. Use \code{ts_to_df()} instead.
ts2df <- function(x, time_col = "Time") {
  .Deprecated("ts_to_df", package = "rSHUD",
              msg = "ts2df() is deprecated. Use ts_to_df() instead.")
  ts_to_df(x, time_col)
}
