# Tests for Validation Functions

test_that("check_positive validates positive numbers", {
  expect_true(check_positive(5, "test"))
  expect_true(check_positive(0.1, "test"))
  expect_true(check_positive(1000, "test"))
})

test_that("check_positive rejects non-positive numbers", {
  expect_error(check_positive(0, "test"), "positive")
  expect_error(check_positive(-1, "test"), "positive")
  expect_error(check_positive(-0.5, "test"), "positive")
})

test_that("check_positive allow_zero parameter works", {
  expect_true(check_positive(0, "test", allow_zero = TRUE))
  expect_true(check_positive(5, "test", allow_zero = TRUE))
  expect_error(check_positive(-1, "test", allow_zero = TRUE), "non-negative")
})

test_that("check_positive validates input type", {
  expect_error(check_positive("5", "test"), "numeric")
  expect_error(check_positive(c(1, 2), "test"), "single")
  expect_error(check_positive(NA, "test"), "numeric value")
  expect_error(check_positive(Inf, "test"), "infinite")
})

test_that("check_file_exists validates existing files", {
  tmp_file <- create_temp_file()
  expect_true(check_file_exists(tmp_file, "test"))
  cleanup_temp_files(tmp_file)
})

test_that("check_file_exists rejects non-existent files", {
  expect_error(check_file_exists("/nonexistent/file.txt", "test"),
               "does not exist")
})

test_that("check_file_exists rejects directories when must_be_file=TRUE", {
  tmp_dir <- tempdir()
  expect_error(check_file_exists(tmp_dir, "test", must_be_file = TRUE),
               "must be a file")
})

test_that("check_file_exists allows directories when must_be_file=FALSE", {
  tmp_dir <- tempdir()
  expect_true(check_file_exists(tmp_dir, "test", must_be_file = FALSE))
})

test_that("check_file_exists validates input type", {
  expect_error(check_file_exists(123, "test"), "character")
  expect_error(check_file_exists(c("a", "b"), "test"), "single")
  expect_error(check_file_exists("", "test"), "empty")
})

test_that("check_spatial_compatible works with compatible objects", {
  skip_if_not_installed("terra")

  r1 <- create_test_raster_terra()
  r2 <- create_test_raster_terra()

  expect_true(check_spatial_compatible(r1, r2))
})

test_that("check_spatial_compatible detects CRS mismatch", {
  skip_if_not_installed("terra")

  r1 <- create_test_raster_terra()
  r2 <- create_test_raster_terra()
  terra::crs(r2) <- "EPSG:3857"

  expect_error(check_spatial_compatible(r1, r2), "incompatible CRS")
})

test_that("check_spatial_compatible can skip CRS check", {
  skip_if_not_installed("terra")

  r1 <- create_test_raster_terra()
  r2 <- create_test_raster_terra()
  terra::crs(r2) <- "EPSG:3857"

  expect_true(check_spatial_compatible(r1, r2, check_crs = FALSE))
})

test_that("check_spatial_compatible detects non-overlapping extents", {
  skip_if_not_installed("terra")

  r1 <- terra::rast(ncol = 10, nrow = 10,
                    xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  r2 <- terra::rast(ncol = 10, nrow = 10,
                    xmin = 20, xmax = 30, ymin = 20, ymax = 30)
  terra::crs(r1) <- terra::crs(r2) <- "EPSG:4326"

  expect_error(check_spatial_compatible(r1, r2, check_overlap = TRUE),
               "do not overlap")
})

test_that("projected CRS validators accept projected sf and raster inputs", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0),
                   ncol = 2, byrow = TRUE)
  domain <- sf::st_sf(
    geometry = sf::st_sfc(sf::st_polygon(list(coords)), crs = 3857)
  )
  dem <- terra::rast(ncol = 10, nrow = 10,
                     xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                     crs = "EPSG:3857")

  expect_true(check_sf_projected_crs(domain, "domain"))
  expect_true(check_raster_projected_crs(dem, "dem"))
})

test_that("projected CRS validators reject missing CRS", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0),
                   ncol = 2, byrow = TRUE)
  domain <- sf::st_sf(
    geometry = sf::st_sfc(sf::st_polygon(list(coords)))
  )
  dem <- terra::rast(ncol = 10, nrow = 10,
                     xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  terra::crs(dem) <- ""

  expect_error(check_sf_projected_crs(domain, "domain"),
               "domain.*projected CRS in metres/meters")
  expect_error(check_raster_projected_crs(dem, "dem"),
               "dem.*projected CRS in metres/meters")
})

test_that("projected CRS validators reject longitude latitude CRS", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  coords <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0),
                   ncol = 2, byrow = TRUE)
  domain <- sf::st_sf(
    geometry = sf::st_sfc(sf::st_polygon(list(coords)), crs = 4326)
  )
  dem <- terra::rast(ncol = 10, nrow = 10,
                     xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                     crs = "EPSG:4326")

  expect_error(check_sf_projected_crs(domain, "domain"),
               "domain.*longitude/latitude CRS.*Transform 'domain'")
  expect_error(check_raster_projected_crs(dem, "dem"),
               "dem.*longitude/latitude CRS.*Transform 'dem'")
})

test_that("projected CRS validators reject projected feet units", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0),
                   ncol = 2, byrow = TRUE)
  domain <- sf::st_sf(
    geometry = sf::st_sfc(sf::st_polygon(list(coords)), crs = 2272)
  )
  dem <- terra::rast(ncol = 10, nrow = 10,
                     xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                     crs = "EPSG:2272")

  expect_error(check_sf_projected_crs(domain, "domain"),
               "domain.*metres/meters.*US survey foot")
  expect_error(check_raster_projected_crs(dem, "dem"),
               "dem.*metres/meters.*US survey foot")
})

test_that("projected CRS unit check rejects unknown units", {
  expect_error(check_crs_units_metres(NA_character_, "domain"),
               "domain.*CRS units.*could not be determined")
  expect_error(check_crs_units_metres("", "dem"),
               "dem.*CRS units.*could not be determined")
})

test_that("compatible_crs detects matching CRS", {
  expect_true(compatible_crs("EPSG:4326", "EPSG:4326"))
  expect_true(compatible_crs("+proj=longlat", "+proj=longlat"))
})

test_that("compatible_crs detects different CRS", {
  expect_false(compatible_crs("EPSG:4326", "EPSG:3857"))
  expect_false(compatible_crs("+proj=longlat", "+proj=utm"))
})

test_that("compatible_crs handles NULL inputs", {
  expect_false(compatible_crs(NULL, "EPSG:4326"))
  expect_false(compatible_crs("EPSG:4326", NULL))
  expect_false(compatible_crs(NULL, NULL))
})

test_that("compatible_crs handles empty strings", {
  expect_true(compatible_crs("", ""))
  expect_false(compatible_crs("", "EPSG:4326"))
  expect_false(compatible_crs("EPSG:4326", ""))
})

test_that("compatible_crs extracts and compares EPSG codes", {
  expect_true(compatible_crs("EPSG:4326", "+init=epsg:4326"))
  expect_false(compatible_crs("EPSG:4326", "+init=epsg:3857"))
})
