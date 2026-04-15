# Test Helper Functions and Data
#
# This file contains helper functions and test data setup for the test suite.

#' Create a simple test raster (terra)
#'
#' @return SpatRaster object
create_test_raster_terra <- function() {
  if (!requireNamespace("terra", quietly = TRUE)) {
    testthat::skip("terra not available")
  }
  r <- terra::rast(ncol = 10, nrow = 10,
                   xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  terra::values(r) <- 1:100
  terra::crs(r) <- "EPSG:4326"
  return(r)
}

#' Create a simple test polygon (sf)
#'
#' @return sf object
create_test_polygon_sf <- function() {
  if (!requireNamespace("sf", quietly = TRUE)) {
    testthat::skip("sf not available")
  }
  coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0),
                   ncol = 2, byrow = TRUE)
  poly <- sf::st_polygon(list(coords))
  sfc <- sf::st_sfc(poly, crs = 4326)
  sf_obj <- sf::st_sf(id = 1, geometry = sfc)
  return(sf_obj)
}

#' Create a simple test polygon (terra)
#'
#' @return SpatVector object
create_test_polygon_terra <- function() {
  if (!requireNamespace("terra", quietly = TRUE)) {
    testthat::skip("terra not available")
  }
  coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0),
                   ncol = 2, byrow = TRUE)
  poly <- terra::vect(coords, type = "polygons", crs = "EPSG:4326")
  return(poly)
}

#' Create temporary test file
#'
#' @param content Character vector of file content
#' @param ext File extension (default: "txt")
#' @return Character string with file path
create_temp_file <- function(content = "test", ext = "txt") {
  tmp_file <- tempfile(fileext = paste0(".", ext))
  writeLines(content, tmp_file)
  return(tmp_file)
}

#' Clean up temporary files
#'
#' @param files Character vector of file paths to remove
cleanup_temp_files <- function(files) {
  for (f in files) {
    if (file.exists(f)) {
      unlink(f)
    }
  }
}

#' Stable replacement for testthat::skip_if_not_installed
#'
#' In this project, namespace loading can crash some local installed-package
#' checks. Probe package availability in a child R process so the main test
#' process can safely skip instead of aborting.
namespace_probe_cache <- new.env(parent = emptyenv())

can_load_namespace <- function(pkg) {
  if (pkg %in% loadedNamespaces()) {
    return(TRUE)
  }

  if (exists(pkg, envir = namespace_probe_cache, inherits = FALSE)) {
    return(get(pkg, envir = namespace_probe_cache, inherits = FALSE))
  }

  expr <- sprintf(
    'quit(status = if (requireNamespace("%s", quietly = TRUE)) 0 else 1)',
    pkg
  )
  r_bin <- file.path(R.home("bin"), "R")
  stdout_file <- tempfile()
  stderr_file <- tempfile()
  on.exit(unlink(c(stdout_file, stderr_file)), add = TRUE)
  status <- tryCatch(
    suppressWarnings(system2(
      r_bin,
      c("--vanilla", "-q", "-e", shQuote(expr)),
      stdout = stdout_file,
      stderr = stderr_file
    )),
    error = function(e) 1L
  )
  ok <- isTRUE(as.integer(status) == 0L)
  assign(pkg, ok, envir = namespace_probe_cache)
  ok
}

skip_if_not_installed <- function(pkg, ...) {
  if (!can_load_namespace(pkg)) {
    testthat::skip(paste0("Package '", pkg, "' not available"))
  }
}

run_r_subprocess <- function(expr) {
  r_bin <- file.path(R.home("bin"), "R")
  stdout_file <- tempfile()
  stderr_file <- tempfile()
  on.exit(unlink(c(stdout_file, stderr_file)), add = TRUE)

  status <- tryCatch(
    suppressWarnings(system2(
      r_bin,
      c("--vanilla", "-q", "-e", shQuote(expr)),
      stdout = stdout_file,
      stderr = stderr_file
    )),
    error = function(e) 1L
  )

  output <- paste(
    c(
      readLines(stdout_file, warn = FALSE),
      readLines(stderr_file, warn = FALSE)
    ),
    collapse = "\n"
  )

  list(
    status = as.integer(status),
    output = output
  )
}

expect_r_subprocess_ok <- function(expr, pattern = NULL) {
  result <- run_r_subprocess(expr)
  testthat::expect_equal(result$status, 0L, info = result$output)
  if (!is.null(pattern)) {
    testthat::expect_match(result$output, pattern)
  }
  invisible(result)
}

#' Skip tests when s2 runtime is unavailable
#'
#' Some sf workflows load s2 dynamically. On machines with a broken s2/abseil
#' runtime, these tests should be skipped rather than fail unrelated package
#' checks.
skip_if_s2_unavailable <- function() {
  if (!can_load_namespace("s2")) {
    testthat::skip("s2 runtime not available")
  }
}
