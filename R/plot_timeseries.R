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

#' Plot time-series hydrograph
#'
#' Compatibility alias for \code{\link{plot_hydrograph}} retained for code that
#' used the transitional \code{plot_timeseries()} name.
#'
#' @param x xts or zoo time-series matrix.
#' @param ... Additional arguments passed to \code{\link{plot_hydrograph}}.
#' @return A combined plot object from \code{\link{plot_hydrograph}}.
#' @export
#' @examples
#' \dontrun{
#' library(xts)
#' dates <- as.POSIXct(as.Date("2020-01-01") + 1:10)
#' x <- xts(cbind(precip = 1:10, discharge = 11:20), order.by = dates)
#' plot_timeseries(x)
#' }
plot_timeseries <- function(x, ...) {
  plot_hydrograph(x, ...)
}
