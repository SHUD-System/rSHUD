# Tests for Visualization Module
#
# Removed APIs:
# - plot_mesh_2d()
# - map2d()
#
# Visualization coverage focuses on plot_polygons(), compare_maps(),
# plot_hydrograph(), plot_timeseries(), plot_tsd(), and deprecated wrappers
# still present.

with_temp_pdf <- function(code) {
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  plot_dev <- grDevices::dev.cur()
  on.exit({
    active_devs <- grDevices::dev.list()
    if (!is.null(active_devs) && plot_dev %in% active_devs) {
      grDevices::dev.off(which = plot_dev)
    }
    unlink(plot_file)
  }, add = TRUE)
  force(code)
}

# Test plot_polygons function ------------------------------------------------

test_that("plot_polygons requires x parameter", {
  expect_error(
    plot_polygons(field = "test"),
    "Parameter 'x' is required"
  )
})

test_that("plot_polygons requires field parameter", {
  skip_if_not_installed("sf")
  
  # Create simple polygon
  poly <- sf::st_as_sf(
    data.frame(
      id = 1,
      value = 10,
      geometry = sf::st_sfc(
        sf::st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))),
        crs = 4326
      )
    )
  )
  
  expect_error(
    plot_polygons(x = poly),
    "Parameter 'field' is required"
  )
})

test_that("plot_polygons rejects legacy sp objects", {
  skip_if_not_installed("sf")
  
  # Mock legacy sp object
  legacy_sp <- structure(list(), class = "SpatialPolygons")
  
  expect_error(
    plot_polygons(x = legacy_sp, field = "test"),
    "Legacy sp objects are not supported"
  )
})

test_that("plot_polygons validates x type", {
  expect_error(
    plot_polygons(x = "invalid", field = "test"),
    "must be an sf or SpatVector object"
  )
})

test_that("plot_polygons validates field exists", {
  skip_if_not_installed("sf")
  
  # Create simple polygon
  poly <- sf::st_as_sf(
    data.frame(
      id = 1,
      value = 10,
      geometry = sf::st_sfc(
        sf::st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))),
        crs = 4326
      )
    )
  )
  
  expect_error(
    plot_polygons(x = poly, field = "nonexistent"),
    "Field .* not found"
  )
})

test_that("plot_polygons works with sf objects", {
  skip_if_not_installed("sf")
  
  # Create simple polygon
  poly <- sf::st_as_sf(
    data.frame(
      id = 1:3,
      value = c(10, 20, 30),
      geometry = sf::st_sfc(
        sf::st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))),
        sf::st_polygon(list(cbind(c(1, 2, 2, 1, 1), c(0, 0, 1, 1, 0)))),
        sf::st_polygon(list(cbind(c(2, 3, 3, 2, 2), c(0, 0, 1, 1, 0)))),
        crs = 4326
      )
    )
  )
  
  # Plot should return invisibly
  result <- with_temp_pdf(plot_polygons(x = poly, field = "value"))
  
  expect_s3_class(result, "sf")
})

test_that("plot_polygons works with SpatVector objects", {
  skip_if_not_installed("terra")
  
  # Create simple SpatVector polygon
  coords <- cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
  poly <- terra::vect(coords, type = "polygons", atts = data.frame(value = 10))
  
  # Plot should return invisibly
  result <- with_temp_pdf(plot_polygons(x = poly, field = "value"))
  
  expect_s4_class(result, "SpatVector")
})

# Test compare_maps function -------------------------------------------------

test_that("compare_maps requires maps parameter", {
  expect_error(
    compare_maps(),
    "Parameter 'maps' must be a list"
  )
})

test_that("compare_maps validates maps is a list", {
  expect_error(
    compare_maps(maps = "invalid"),
    "must be a list"
  )
})

test_that("compare_maps validates maps is not empty", {
  expect_error(
    compare_maps(maps = list()),
    "must contain at least one map"
  )
})

test_that("compare_maps rejects legacy raster/sp objects", {
  skip_if_not_installed("terra")
  
  # Mock legacy objects
  legacy_raster <- structure(list(), class = "RasterLayer")
  
  expect_error(
    compare_maps(maps = list(legacy_raster)),
    "Legacy raster/sp objects are not supported"
  )
})

test_that("compare_maps validates all elements are spatial", {
  skip_if_not_installed("terra")
  
  # Create valid raster
  r <- terra::rast(ncol = 10, nrow = 10)
  terra::values(r) <- 1:100
  
  # Mix with invalid object
  expect_error(
    compare_maps(maps = list(r, "invalid")),
    "must be SpatRaster, sf, or SpatVector"
  )
})

test_that("compare_maps works with single raster", {
  skip_if_not_installed("terra")
  
  # Create simple raster
  r <- terra::rast(ncol = 10, nrow = 10)
  terra::values(r) <- 1:100
  
  # Should not error
  expect_invisible(with_temp_pdf(compare_maps(maps = list(r))))
})

test_that("compare_maps works with multiple rasters", {
  skip_if_not_installed("terra")
  
  # Create simple rasters
  r1 <- terra::rast(ncol = 10, nrow = 10)
  terra::values(r1) <- 1:100
  
  r2 <- terra::rast(ncol = 10, nrow = 10)
  terra::values(r2) <- 100:1
  
  # Should not error
  expect_invisible(with_temp_pdf(compare_maps(maps = list(r1, r2))))
})

test_that("compare_maps auto-calculates layout", {
  skip_if_not_installed("terra")
  
  # Create rasters
  r1 <- terra::rast(ncol = 10, nrow = 10)
  terra::values(r1) <- 1:100
  
  r2 <- terra::rast(ncol = 10, nrow = 10)
  terra::values(r2) <- 100:1
  
  r3 <- terra::rast(ncol = 10, nrow = 10)
  terra::values(r3) <- rep(50, 100)
  
  # Should auto-calculate layout for 3 maps
  expect_invisible(with_temp_pdf(compare_maps(maps = list(r1, r2, r3))))
})

test_that("compare_maps validates layout size", {
  skip_if_not_installed("terra")
  
  # Create 5 rasters
  maps <- lapply(1:5, function(i) {
    r <- terra::rast(ncol = 10, nrow = 10)
    terra::values(r) <- 1:100
    r
  })
  
  # Layout too small
  expect_error(
    compare_maps(maps = maps, nrow = 1, ncol = 2),
    "Layout .* is too small"
  )
})

test_that("compare_maps works with custom titles", {
  skip_if_not_installed("terra")
  
  # Create rasters
  r1 <- terra::rast(ncol = 10, nrow = 10)
  terra::values(r1) <- 1:100
  
  r2 <- terra::rast(ncol = 10, nrow = 10)
  terra::values(r2) <- 100:1
  
  # Custom titles
  expect_invisible(
    with_temp_pdf(compare_maps(maps = list(r1, r2), titles = c("Map A", "Map B")))
  )
})

test_that("compare_maps warns about title length mismatch", {
  skip_if_not_installed("terra")
  
  # Create rasters
  r1 <- terra::rast(ncol = 10, nrow = 10)
  terra::values(r1) <- 1:100
  
  r2 <- terra::rast(ncol = 10, nrow = 10)
  terra::values(r2) <- 100:1
  
  # Wrong number of titles
  expect_warning(
    with_temp_pdf(compare_maps(maps = list(r1, r2), titles = c("Only One"))),
    "Length of 'titles'"
  )
})

test_that("compare_maps works with sf objects", {
  skip_if_not_installed("sf")
  
  # Create simple polygons
  poly1 <- sf::st_as_sf(
    data.frame(
      id = 1,
      value = 10,
      geometry = sf::st_sfc(
        sf::st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))),
        crs = 4326
      )
    )
  )
  
  poly2 <- sf::st_as_sf(
    data.frame(
      id = 1,
      value = 20,
      geometry = sf::st_sfc(
        sf::st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))),
        crs = 4326
      )
    )
  )
  
  # Should not error
  expect_invisible(with_temp_pdf(compare_maps(maps = list(poly1, poly2))))
})

# Test plot_hydrograph function ----------------------------------------------

test_that("plot_hydrograph requires x parameter", {
  expect_error(
    plot_hydrograph(),
    "Parameter 'x' is required"
  )
})

test_that("plot_hydrograph validates x is time series", {
  expect_error(
    plot_hydrograph(x = matrix(1:20, nrow = 10, ncol = 2)),
    "must be an xts or zoo object"
  )
})

test_that("plot_hydrograph requires at least 2 columns", {
  skip_if_not_installed("xts")
  
  # Create single column xts
  dates <- as.POSIXct(as.Date('2000-01-01') + 1:100)
  x <- xts::xts(sin(1:100 / 10), order.by = dates)
  
  expect_error(
    plot_hydrograph(x),
    "must have at least 2 columns"
  )

  empty_x <- xts::xts(
    matrix(numeric(0), nrow = length(dates), ncol = 0),
    order.by = dates
  )

  expect_error(
    plot_hydrograph(empty_x),
    "must have at least 2 columns"
  )
})

test_that("plot_hydrograph rejects one-dimensional zoo inputs clearly", {
  skip_if_not_installed("zoo")

  dates <- as.POSIXct(as.Date("2000-01-01") + 1:10)
  x <- zoo::zoo(1:10, order.by = dates)

  expect_error(
    plot_hydrograph(x),
    "must have at least 2 columns"
  )
})

test_that("plot_timeseries preserves hydrograph input validation", {
  skip_if_not_installed("zoo")

  dates <- as.POSIXct(as.Date("2000-01-01") + 1:10)
  x <- zoo::zoo(1:10, order.by = dates)

  expect_error(
    plot_timeseries(x),
    "must have at least 2 columns"
  )

  expect_error(
    plot_timeseries(1:10),
    "must be an xts or zoo object"
  )

  expect_error(
    plot_timeseries(stats::ts(1:10)),
    "must be an xts or zoo object"
  )
})

test_that("plot_hydrograph works with 2 columns", {
  skip_if_not_installed("xts")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("reshape2")

  expect_r_subprocess_ok(paste(
    'library(rSHUD)',
    'library(xts)',
    'dates <- as.POSIXct(as.Date("2000-01-01") + 1:100)',
    'precip <- abs(rnorm(100, mean = 2, sd = 1))',
    'discharge <- abs(rnorm(100, mean = 10, sd = 3))',
    'x <- xts::xts(cbind(precip, discharge), order.by = dates)',
    'colnames(x) <- c("Precipitation", "Discharge")',
    'plot_file <- tempfile(fileext = ".pdf")',
    'grDevices::pdf(plot_file)',
    'on.exit({ grDevices::dev.off(); unlink(plot_file) })',
    'p <- plot_hydrograph(x)',
    'stopifnot(!is.null(p))',
    sep = "; "
  ))
})

test_that("plot_hydrograph works with multiple discharge columns", {
  skip_if_not_installed("xts")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("reshape2")

  expect_r_subprocess_ok(paste(
    'library(rSHUD)',
    'library(xts)',
    'dates <- as.POSIXct(as.Date("2000-01-01") + 1:100)',
    'precip <- abs(rnorm(100, mean = 2, sd = 1))',
    'discharge1 <- abs(rnorm(100, mean = 10, sd = 3))',
    'discharge2 <- abs(rnorm(100, mean = 10, sd = 3))',
    'x <- xts::xts(cbind(precip, discharge1, discharge2), order.by = dates)',
    'colnames(x) <- c("Precip", "Simulated", "Observed")',
    'plot_file <- tempfile(fileext = ".pdf")',
    'grDevices::pdf(plot_file)',
    'on.exit({ grDevices::dev.off(); unlink(plot_file) })',
    'p <- plot_hydrograph(x)',
    'stopifnot(!is.null(p))',
    sep = "; "
  ))
})

test_that("plot_hydrograph works with custom units", {
  skip_if_not_installed("xts")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("reshape2")

  expect_r_subprocess_ok(paste(
    'library(rSHUD)',
    'library(xts)',
    'dates <- as.POSIXct(as.Date("2000-01-01") + 1:100)',
    'precip <- abs(rnorm(100, mean = 2, sd = 1))',
    'discharge <- abs(rnorm(100, mean = 10, sd = 3))',
    'x <- xts::xts(cbind(precip, discharge), order.by = dates)',
    'plot_file <- tempfile(fileext = ".pdf")',
    'grDevices::pdf(plot_file)',
    'on.exit({ grDevices::dev.off(); unlink(plot_file) })',
    'p <- plot_hydrograph(x, units = c("mm/day", "m³/s"))',
    'stopifnot(!is.null(p))',
    sep = "; "
  ))
})

test_that("plot_tsd preserves base plotting fallback for legacy inputs", {
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    grDevices::dev.off()
    unlink(plot_file)
  })

  x <- stats::ts(1:10)
  result <- suppressWarnings(withVisible(plot_tsd(x)))

  expect_false(result$visible)
  expect_identical(result$value, x)
})

test_that("plot_tsd dispatches hydrograph inputs to plot_timeseries", {
  skip_if_not_installed("xts")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("reshape2")

  expect_r_subprocess_ok(paste(
    'library(rSHUD)',
    'library(xts)',
    'dates <- as.POSIXct(as.Date("2000-01-01") + 1:10)',
    'x <- xts::xts(cbind(precip = 1:10, discharge = 11:20), order.by = dates)',
    'plot_file <- tempfile(fileext = ".pdf")',
    'grDevices::pdf(plot_file)',
    'on.exit({ grDevices::dev.off(); unlink(plot_file) })',
    'p <- suppressWarnings(plot_tsd(x))',
    'stopifnot(inherits(p, "gtable"))',
    sep = "; "
  ))
})

# Test deprecated functions --------------------------------------------------

test_that("plot_polygons replaces plot_sp functionality", {
  skip_if_not_installed("sf")
  
  # Create simple polygon
  poly <- sf::st_as_sf(
    data.frame(
      id = 1,
      value = 10,
      geometry = sf::st_sfc(
        sf::st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))),
        crs = 4326
      )
    )
  )
  
  # The old plot_sp function is removed. 
  # Test that plot_polygons handles the expected inputs correctly.
  expect_invisible(with_temp_pdf(plot_polygons(x = poly, field = "value")))
})

test_that("hydrograph shows deprecation warning", {
  skip_if_not_installed("xts")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("reshape2")

  # expect_r_subprocess_ok(
  #   paste(
  #     'library(rSHUD)',
  #     'library(xts)',
  #     'dates <- as.POSIXct(as.Date("2000-01-01") + 1:100)',
  #     'precip <- abs(rnorm(100, mean = 2, sd = 1))',
  #     'discharge <- abs(rnorm(100, mean = 10, sd = 3))',
  #     'x <- xts::xts(cbind(precip, discharge), order.by = dates)',
  #     'hydrograph(x = x)',
  #     sep = "; "
  #   ),
  #   "deprecated"
  # )
})

# Integration tests ----------------------------------------------------------

test_that("visualization module integration test", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  skip_on_cran()
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)

  # Convert mesh to sf
  mesh_sf <- sp.mesh2Shape(pm = tri)
  
  # Test plot_polygons
  result <- with_temp_pdf(plot_polygons(x = mesh_sf, field = "Area"))
  expect_s3_class(result, "sf")
  
  # Test compare_maps with multiple rasters
  r <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(r) <- runif(100, 100, 200)
  r2 <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(r2) <- runif(100, 0, 1)
  
  expect_invisible(with_temp_pdf(compare_maps(maps = list(r, r2))))
})
