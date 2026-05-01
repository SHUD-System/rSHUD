#' Plot hydrograph with precipitation
#'
#' Creates a two-panel hydrograph plot with precipitation on top (inverted)
#' and discharge/flow on bottom. This function uses ggplot2 for improved
#' visualization quality.
#'
#' @param x xts or zoo time series matrix. The first column is precipitation
#'   (plotted on top, inverted), and remaining columns are discharge or other
#'   flow variables (plotted on bottom).
#' @param legend_position Character, position of legend for discharge plot.
#'   Options: "bottom" (default), "top", "left", "right", "none".
#' @param units Character vector of units for each variable. Length should
#'   match number of columns in x. Default is empty strings.
#' @param colors Character vector of colors for precipitation and discharge.
#'   First color for precipitation, second for discharge. Default is c(3, 4)
#'   which maps to green and blue.
#' @param heights Numeric vector of length 2 specifying relative heights of
#'   top (precipitation) and bottom (discharge) panels. Default is c(3, 7).
#' @param ylabs List or character vector of length 2 with y-axis labels for
#'   top and bottom panels (optional). If NULL, uses column names and units.
#' @param title Character, overall plot title (optional).
#'
#' @return A combined plot object (grid arrangement) that can be printed.
#'
#' @section Migration Note:
#' This function is an improved version of \code{hydrograph()} with:
#' \itemize{
#'   \item Parameter names use snake_case for consistency
#'   \item Better handling of multiple discharge variables
#'   \item Improved visual quality with ggplot2
#'   \item More flexible customization options
#'   \item Clearer parameter documentation
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' library(xts)
#' # Create sample hydrograph data
#' nday <- 365
#' dates <- as.POSIXct(as.Date('2020-01-01') + 1:nday)
#' precip <- abs(rnorm(nday, mean = 2, sd = 5))
#' discharge <- abs(rnorm(nday, mean = 10, sd = 3))
#' data <- xts(cbind(precip, discharge), order.by = dates)
#' colnames(data) <- c("Precipitation", "Discharge")
#'
#' # Basic hydrograph
#' plot_hydrograph(data)
#'
#' # With units
#' plot_hydrograph(data, units = c("mm/day", "m³/s"))
#'
#' # Multiple discharge variables
#' obs <- abs(rnorm(nday, mean = 10, sd = 3))
#' data2 <- xts(cbind(precip, discharge, obs), order.by = dates)
#' colnames(data2) <- c("Precip", "Simulated", "Observed")
#' plot_hydrograph(data2, units = c("mm", "m³/s", "m³/s"))
#' }
plot_hydrograph <- function(x,
                            legend_position = "bottom",
                            units = rep('', NCOL(x)),
                            colors = c(3, 4),
                            heights = c(3, 7),
                            ylabs = NULL,
                            title = NULL) {

  # Validate x parameter
  if (missing(x) || is.null(x)) {
    stop("Parameter 'x' is required", call. = FALSE)
  }

  # Check if x is a time series object
  if (!inherits(x, c("xts", "zoo"))) {
    stop(
      "Parameter 'x' must be an xts or zoo object, ",
      "but received ", class(x)[1],
      call. = FALSE
    )
  }

  n_columns <- NCOL(x)

  # Validate x has at least 2 columns
  if (n_columns < 2) {
    stop(
      "Parameter 'x' must have at least 2 columns ",
      "(precipitation and discharge)",
      call. = FALSE
    )
  }

  # Ensure time index is POSIXct
  zoo::index(x) <- as.POSIXct(stats::time(x))

  # Get column names
  cn <- colnames(x)
  if (is.null(cn)) {
    cn <- paste0("Var", seq_len(n_columns))
  }

  # Prepare precipitation data (first column)
  pv <- as.numeric(x[, 1])
  dfp <- data.frame(Time = stats::time(x), Precipitation = pv)

  # Prepare discharge data (remaining columns)
  dfqq <- data.frame(Time = stats::time(x), x[, -1, drop = FALSE])
  dfq <- reshape2::melt(dfqq, id = 'Time')

  # Create top plot (precipitation, inverted)
  plim_p <- range(pv, na.rm = TRUE)

  g_top <- ggplot2::ggplot()
  g_top <- g_top +
    ggplot2::coord_cartesian(ylim = plim_p) +
    ggplot2::guides(fill = "none") +
    ggplot2::geom_col(
      data = dfp,
      ggplot2::aes(x = Time, y = Precipitation),
      fill = colors[1]
    ) +
    ggplot2::scale_y_continuous(trans = "reverse") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  # Set y-label for top plot
  if (!is.null(ylabs)) {
    if (is.list(ylabs)) {
      g_top <- g_top + ggplot2::ylab(ylabs[[1]])
    } else {
      g_top <- g_top + ggplot2::ylab(ylabs[1])
    }
  } else {
    g_top <- g_top + ggplot2::labs(y = paste(cn[1], units[1]))
  }

  # Create bottom plot (discharge)
  plim_q <- range(x[, -1], na.rm = TRUE)

  g_bottom <- ggplot2::ggplot()
  g_bottom <- g_bottom +
    ggplot2::coord_cartesian(ylim = plim_q) +
    ggplot2::guides(fill = "none") +
    ggplot2::geom_line(
      data = dfq,
      ggplot2::aes(
        x = Time,
        y = value,
        linetype = variable,
        color = variable
      )
    ) +
    ggplot2::theme() +
    ggplot2::scale_colour_discrete("") +
    ggplot2::scale_linetype_manual("", values = seq_len(n_columns - 1)) +
    ggplot2::scale_shape_manual("", values = seq_len(n_columns - 1)) +
    ggplot2::labs(x = "Time")

  # Set y-label for bottom plot
  if (!is.null(ylabs)) {
    if (is.list(ylabs)) {
      g_bottom <- g_bottom + ggplot2::ylab(ylabs[[2]])
    } else {
      g_bottom <- g_bottom + ggplot2::ylab(ylabs[2])
    }
  } else {
    # Combine labels for multiple discharge variables
    discharge_labels <- paste(cn[-1], units[-1])
    g_bottom <- g_bottom + ggplot2::ylab(paste(discharge_labels, collapse = ", "))
  }

  # Handle legend position
  if (n_columns > 2) {
    g_bottom <- g_bottom +
      ggplot2::theme(
        legend.position = legend_position,
        legend.direction = 'horizontal',
        legend.title = ggplot2::element_blank()
      )
  } else {
    g_bottom <- g_bottom +
      ggplot2::theme(legend.position = 'none')
  }

  # Add overall title if provided
  if (!is.null(title)) {
    g_top <- g_top + ggplot2::ggtitle(title)
  }

  # Align plot widths
  gA <- ggplot2::ggplotGrob(g_top)
  gB <- ggplot2::ggplotGrob(g_bottom)
  maxWidth <- grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
  gA$widths[2:5] <- as.list(maxWidth)
  gB$widths[2:5] <- as.list(maxWidth)

  # Combine plots
  p <- gridExtra::grid.arrange(gA, gB, ncol = 1, heights = heights)

  return(p)
}

#' Plot generic time-series data
#'
#' Creates a ggplot line chart for ordinary \code{ts}, \code{zoo},
#' \code{xts}, matrix, or data frame time-series data. For hydrograph plots
#' with precipitation and discharge panels, use \code{\link{plot_hydrograph}}.
#'
#' @param x Time-series object or tabular data with a time column.
#' @param time_col Optional time column name or index for matrix/data frame
#'   inputs. If omitted, common names such as \code{time}, \code{date}, or
#'   \code{datetime} are used when present; otherwise the first column is used.
#' @param value_col Optional value column name(s) or index(es) for matrix/data
#'   frame inputs. Defaults to all columns except the time column.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#' @param title Character. Plot title.
#' @param legend_position Character. Position of the legend for multi-series
#'   inputs.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{geom_line}}.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @export
#' @examples
#' \dontrun{
#' library(xts)
#' dates <- as.Date("2020-01-01") + 1:10
#' x <- xts(cbind(simulated = 1:10, observed = 11:20), order.by = dates)
#' plot_timeseries(x)
#' }
plot_timeseries <- function(x,
                            time_col = NULL,
                            value_col = NULL,
                            xlab = "Time",
                            ylab = "Value",
                            title = NULL,
                            legend_position = "right",
                            ...) {
  if (missing(x) || is.null(x)) {
    stop("Parameter 'x' is required", call. = FALSE)
  }

  plot_data <- .as_timeseries_plot_data(
    x,
    time_col = time_col,
    value_col = value_col
  )

  line_args <- list(...)
  if ("main" %in% names(line_args) && is.null(title)) {
    title <- line_args$main
    line_args$main <- NULL
  }

  n_series <- length(unique(plot_data$Series))
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Time, y = Value)
  )

  if (n_series > 1L) {
    p <- p + do.call(
      ggplot2::geom_line,
      c(
        list(mapping = ggplot2::aes(color = Series, linetype = Series)),
        line_args
      )
    )
    p <- p +
      ggplot2::scale_colour_discrete("") +
      ggplot2::scale_linetype_discrete("") +
      ggplot2::theme(legend.position = legend_position)
  } else {
    p <- p + do.call(ggplot2::geom_line, line_args)
    p <- p + ggplot2::theme(legend.position = "none")
  }

  p <- p + ggplot2::labs(x = xlab, y = ylab)
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  p
}

.as_timeseries_plot_data <- function(x, time_col = NULL, value_col = NULL) {
  if (inherits(x, "ts")) {
    return(.ts_to_timeseries_plot_data(x))
  }

  if (inherits(x, c("xts", "zoo"))) {
    return(.zoo_to_timeseries_plot_data(x))
  }

  if (is.data.frame(x) || is.matrix(x)) {
    return(.tabular_to_timeseries_plot_data(
      x,
      time_col = time_col,
      value_col = value_col
    ))
  }

  stop(
    "Parameter 'x' must be a time-series object or tabular data with a time column",
    call. = FALSE
  )
}

.ts_to_timeseries_plot_data <- function(x) {
  values <- as.matrix(x)
  if (is.null(dim(values))) {
    values <- matrix(as.numeric(x), ncol = 1)
  }

  series_names <- colnames(values)
  if (is.null(series_names)) {
    series_names <- if (NCOL(values) == 1L) "Value" else paste0("Series", seq_len(NCOL(values)))
    colnames(values) <- series_names
  }

  wide <- data.frame(
    Time = as.numeric(stats::time(x)),
    values,
    check.names = FALSE
  )

  .melt_timeseries_plot_data(wide)
}

.zoo_to_timeseries_plot_data <- function(x) {
  values <- zoo::coredata(x)
  if (is.null(dim(values))) {
    values <- matrix(values, ncol = 1)
  }

  series_names <- colnames(values)
  if (is.null(series_names)) {
    series_names <- if (NCOL(values) == 1L) "Value" else paste0("Series", seq_len(NCOL(values)))
    colnames(values) <- series_names
  }

  wide <- data.frame(
    Time = zoo::index(x),
    values,
    check.names = FALSE
  )

  .melt_timeseries_plot_data(wide)
}

.tabular_to_timeseries_plot_data <- function(x,
                                             time_col = NULL,
                                             value_col = NULL,
                                             allow_first_column = TRUE) {
  df <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  if (NCOL(df) < 2L) {
    stop("Tabular time-series data must have a time column and at least one value column",
         call. = FALSE)
  }

  time_name <- .resolve_timeseries_column(
    df,
    column = time_col,
    role = "time",
    allow_first_column = allow_first_column
  )
  value_names <- .resolve_timeseries_value_columns(df, time_name, value_col)

  wide <- df[, c(time_name, value_names), drop = FALSE]
  names(wide)[1] <- "Time"

  .melt_timeseries_plot_data(wide)
}

.resolve_timeseries_column <- function(df,
                                       column = NULL,
                                       role = "time",
                                       allow_first_column = TRUE) {
  column_names <- names(df)

  if (!is.null(column)) {
    if (is.numeric(column) && length(column) == 1L &&
        column >= 1L && column <= NCOL(df)) {
      return(column_names[column])
    }

    if (is.character(column) && length(column) == 1L &&
        column %in% column_names) {
      return(column)
    }

    stop("Unknown ", role, " column: ", paste(column, collapse = ", "),
         call. = FALSE)
  }

  name_match <- match(
    c("time", "date", "datetime", "timestamp", "times", "dates"),
    tolower(column_names),
    nomatch = 0L
  )
  name_match <- name_match[name_match > 0L]
  if (length(name_match) > 0L) {
    return(column_names[name_match[1]])
  }

  time_like <- vapply(
    df,
    function(col) inherits(col, c("Date", "POSIXct", "POSIXlt")),
    logical(1)
  )
  if (any(time_like)) {
    return(column_names[which(time_like)[1]])
  }

  if (allow_first_column) {
    return(column_names[1])
  }

  NULL
}

.resolve_timeseries_value_columns <- function(df, time_name, value_col = NULL) {
  if (is.null(value_col)) {
    value_names <- setdiff(names(df), time_name)
  } else if (is.numeric(value_col)) {
    if (any(value_col < 1L | value_col > NCOL(df))) {
      stop("Unknown value column index", call. = FALSE)
    }
    value_names <- names(df)[value_col]
  } else if (is.character(value_col)) {
    missing_values <- setdiff(value_col, names(df))
    if (length(missing_values) > 0L) {
      stop("Unknown value column: ", paste(missing_values, collapse = ", "),
           call. = FALSE)
    }
    value_names <- value_col
  } else {
    stop("Parameter 'value_col' must be a column name or index", call. = FALSE)
  }

  value_names <- setdiff(value_names, time_name)
  if (length(value_names) == 0L) {
    stop("At least one value column is required", call. = FALSE)
  }

  value_names
}

.melt_timeseries_plot_data <- function(wide) {
  value_names <- setdiff(names(wide), "Time")
  numeric_values <- vapply(
    wide[value_names],
    function(col) is.numeric(col) || is.integer(col) || is.logical(col),
    logical(1)
  )

  if (!all(numeric_values)) {
    stop("Time-series value columns must be numeric", call. = FALSE)
  }

  plot_data <- reshape2::melt(
    wide,
    id.vars = "Time",
    variable.name = "Series",
    value.name = "Value"
  )
  plot_data$Series <- factor(plot_data$Series, levels = value_names)
  plot_data
}

.is_generic_timeseries_input <- function(x) {
  if (inherits(x, c("ts", "xts", "zoo"))) {
    return(TRUE)
  }

  if (!(is.data.frame(x) || is.matrix(x)) || NCOL(x) < 2L) {
    return(FALSE)
  }

  df <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  allow_first_column <- is.matrix(x) && is.null(colnames(x))
  !is.null(.resolve_timeseries_column(
    df,
    allow_first_column = allow_first_column
  ))
}
