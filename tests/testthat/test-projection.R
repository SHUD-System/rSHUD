# Tests for Coordinate Projection Functions

# Test crs.Albers function ----------------------------------------------------

test_that("crs.Albers requires either spx or ext", {
  expect_error(
    crs.Albers(),
    "Either 'spx' or 'ext' must be provided"
  )
})

test_that("crs.Albers rejects legacy spatial objects", {
  # Mock legacy objects
  legacy_sp <- structure(list(), class = "SpatialPoints")
  legacy_raster <- structure(list(), class = "RasterLayer")

  expect_error(
    crs.Albers(spx = legacy_sp),
    "Legacy raster/sp objects are not supported"
  )
  expect_error(
    crs.Albers(spx = legacy_raster),
    "Legacy raster/sp objects are not supported"
  )
})

test_that("crs.Albers validates extent format", {
  expect_error(
    crs.Albers(ext = c(1, 2, 3)),
    "must be a numeric vector of length 4"
  )
  expect_error(
    crs.Albers(ext = "invalid"),
    "must be a numeric vector of length 4"
  )
})

test_that("crs.Albers works with extent vector", {
  crs_str <- crs.Albers(ext = c(-120, -110, 35, 45))

  expect_type(crs_str, "character")
  expect_match(crs_str, "\\+proj=aea")
  expect_match(crs_str, "\\+lat_1=")
  expect_match(crs_str, "\\+lat_2=")
  expect_match(crs_str, "\\+lon_0=")
  expect_match(crs_str, "\\+datum=WGS84")
  expect_match(crs_str, "\\+units=m")
})

test_that("crs.Albers works with sf object", {
  skip_if_not_installed("sf")

  pts <- sf::st_as_sf(
    data.frame(x = c(-115, -112), y = c(40, 42)),
    coords = c("x", "y"),
    crs = 4326
  )

  crs_str <- crs.Albers(spx = pts)

  expect_type(crs_str, "character")
  expect_match(crs_str, "\\+proj=aea")
})

test_that("crs.Albers works with SpatVector", {
  skip_if_not_installed("terra")

  pts <- terra::vect(
    cbind(x = c(-115, -112), y = c(40, 42)),
    type = "points",
    crs = "EPSG:4326"
  )

  crs_str <- crs.Albers(spx = pts)

  expect_type(crs_str, "character")
  expect_match(crs_str, "\\+proj=aea")
})

test_that("crs.Albers requires CRS on spatial objects", {
  skip_if_not_installed("terra")

  pts <- terra::vect(
    cbind(x = c(-115, -112), y = c(40, 42)),
    type = "points"
  )

  expect_error(
    crs.Albers(spx = pts),
    "must have a defined CRS"
  )
})

# Test crs.Lambert function ---------------------------------------------------

test_that("crs.Lambert requires either spx or ext", {
  expect_error(
    crs.Lambert(),
    "Either 'spx' or 'ext' must be provided"
  )
})

test_that("crs.Lambert rejects legacy spatial objects", {
  # Mock legacy objects
  legacy_sp <- structure(list(), class = "SpatialPolygons")

  expect_error(
    crs.Lambert(spx = legacy_sp),
    "Legacy raster/sp objects are not supported"
  )
})

test_that("crs.Lambert works with extent vector", {
  crs_str <- crs.Lambert(ext = c(-120, -110, 35, 45))

  expect_type(crs_str, "character")
  expect_match(crs_str, "\\+proj=laea")
  expect_match(crs_str, "\\+lat_0=")
  expect_match(crs_str, "\\+lon_0=")
  expect_match(crs_str, "\\+datum=WGS84")
  expect_match(crs_str, "\\+units=m")
})

test_that("crs.Lambert works with sf object", {
  skip_if_not_installed("sf")

  pts <- sf::st_as_sf(
    data.frame(x = c(-115, -112), y = c(40, 42)),
    coords = c("x", "y"),
    crs = 4326
  )

  crs_str <- crs.Lambert(spx = pts)

  expect_type(crs_str, "character")
  expect_match(crs_str, "\\+proj=laea")
})

# Test crs.long2utmZone function ----------------------------------------------

test_that("crs.long2utmZone validates input type", {
  expect_error(
    crs.long2utmZone("invalid"),
    "'lon' must be numeric"
  )
})

test_that("crs.long2utmZone validates longitude range", {
  expect_error(
    crs.long2utmZone(200),
    "must be in range \\(-180, 180\\)"
  )
  expect_error(
    crs.long2utmZone(-200),
    "must be in range \\(-180, 180\\)"
  )
})

test_that("crs.long2utmZone calculates correct zones", {
  # Test known values
  expect_equal(crs.long2utmZone(-75), 18)  # Eastern US
  expect_equal(crs.long2utmZone(0), 31)    # Prime meridian
  expect_equal(crs.long2utmZone(150), 56)  # Eastern Australia

  # Test vector input
  lons <- c(-75, 0, 150)
  zones <- crs.long2utmZone(lons)
  expect_equal(zones, c(18, 31, 56))
})

test_that("crs.long2utmZone handles edge cases", {
  expect_equal(crs.long2utmZone(-180), 1)
  expect_equal(crs.long2utmZone(179), 60)
})

# Test crs.long2utm function --------------------------------------------------

test_that("crs.long2utm validates longitude", {
  expect_error(
    crs.long2utm("invalid"),
    "'lon' must be numeric"
  )
  expect_error(
    crs.long2utm(200),
    "must be in range \\(-180, 180\\)"
  )
})

test_that("crs.long2utm validates latitude", {
  expect_error(
    crs.long2utm(-75, lat = "invalid"),
    "'lat' must be numeric"
  )
  expect_error(
    crs.long2utm(-75, lat = 100),
    "must be in range \\(-90, 90\\)"
  )
})

test_that("crs.long2utm creates correct CRS for Northern Hemisphere", {
  crs_str <- crs.long2utm(-75, 40)

  expect_type(crs_str, "character")
  expect_match(crs_str, "\\+proj=utm")
  expect_match(crs_str, "\\+zone=18")
  expect_match(crs_str, "\\+datum=WGS84")
  expect_false(grepl("\\+south", crs_str))
})

test_that("crs.long2utm creates correct CRS for Southern Hemisphere", {
  crs_str <- crs.long2utm(150, -30)

  expect_type(crs_str, "character")
  expect_match(crs_str, "\\+proj=utm")
  expect_match(crs_str, "\\+zone=56")
  expect_match(crs_str, "\\+south")
})

test_that("crs.long2utm handles vector input", {
  lons <- c(-75, 0, 150)
  lats <- c(40, 50, -30)

  result <- crs.long2utm(lons, lats)

  expect_type(result, "list")
  expect_length(result, 3)
  expect_match(result[[1]], "\\+zone=18")
  expect_match(result[[2]], "\\+zone=31")
  expect_match(result[[3]], "\\+zone=56")
  expect_match(result[[3]], "\\+south")
})

test_that("crs.long2utm recycles single latitude", {
  lons <- c(-75, -76, -77)
  result <- crs.long2utm(lons, lat = 40)

  expect_type(result, "list")
  expect_length(result, 3)
})

test_that("crs.long2utm validates lat/lon length mismatch", {
  expect_error(
    crs.long2utm(c(-75, -76), c(40, 41, 42)),
    "must be either length 1 or same length as 'lon'"
  )
})

# Test project_coords function ------------------------------------------------

test_that("project_coords requires to_crs", {
  skip_if_not_installed("sf")

  pts <- sf::st_as_sf(
    data.frame(x = c(-75, -76), y = c(40, 41)),
    coords = c("x", "y"),
    crs = 4326
  )

  expect_error(
    project_coords(pts),
    "'to_crs' is required"
  )
})

test_that("project_coords rejects legacy spatial objects", {
  legacy_sp <- structure(list(), class = "SpatialPoints")

  expect_error(
    project_coords(legacy_sp, to_crs = "EPSG:3857"),
    "Legacy raster/sp objects are not supported"
  )
})

test_that("project_coords works with sf objects", {
  skip_if_not_installed("sf")

  pts <- sf::st_as_sf(
    data.frame(x = c(-75, -76), y = c(40, 41)),
    coords = c("x", "y"),
    crs = 4326
  )

  pts_utm <- project_coords(pts, to_crs = crs.long2utm(-75, 40))

  expect_s3_class(pts_utm, "sf")
  expect_equal(nrow(pts_utm), 2)
})

test_that("project_coords works with SpatVector", {
  skip_if_not_installed("terra")

  pts <- terra::vect(
    cbind(x = c(-75, -76), y = c(40, 41)),
    type = "points",
    crs = "EPSG:4326"
  )

  pts_utm <- project_coords(pts, to_crs = crs.long2utm(-75, 40))

  expect_s4_class(pts_utm, "SpatVector")
  expect_equal(nrow(pts_utm), 2)
})

test_that("project_coords works with SpatRaster", {
  skip_if_not_installed("terra")

  r <- terra::rast(
    xmin = -76, xmax = -74, ymin = 39, ymax = 41,
    resolution = 0.1, crs = "EPSG:4326"
  )
  terra::values(r) <- 1:terra::ncell(r)

  r_utm <- project_coords(r, to_crs = crs.long2utm(-75, 40))

  expect_s4_class(r_utm, "SpatRaster")
})

test_that("project_coords works with coordinate matrix", {
  skip_if_not_installed("sf")

  coords <- matrix(c(-75, 40, -76, 41), ncol = 2, byrow = TRUE)

  coords_utm <- project_coords(coords, from_crs = "EPSG:4326",
                                to_crs = crs.long2utm(-75, 40))

  expect_true(is.matrix(coords_utm))
  expect_equal(nrow(coords_utm), 2)
  expect_equal(ncol(coords_utm), 2)
})

test_that("project_coords validates matrix dimensions", {
  coords <- matrix(1:9, ncol = 3)

  expect_error(
    project_coords(coords, to_crs = "EPSG:3857"),
    "must have exactly 2 columns"
  )
})

test_that("project_coords works with data.frame", {
  skip_if_not_installed("sf")

  df <- data.frame(
    x = c(-75, -76),
    y = c(40, 41),
    value = c(10, 20)
  )

  df_utm <- project_coords(df, from_crs = "EPSG:4326",
                           to_crs = crs.long2utm(-75, 40),
                           coords = c("x", "y"))

  expect_s3_class(df_utm, "data.frame")
  expect_equal(nrow(df_utm), 2)
  expect_true("value" %in% names(df_utm))
  expect_equal(df_utm$value, c(10, 20))
})

test_that("project_coords validates coordinate columns in data.frame", {
  df <- data.frame(a = 1:5, b = 1:5)

  expect_error(
    project_coords(df, to_crs = "EPSG:3857", coords = c("x", "y")),
    "Coordinate columns .* not found"
  )
})

test_that("project_coords assumes EPSG:4326 for matrix without from_crs", {
  skip_if_not_installed("sf")

  coords <- matrix(c(-75, 40, -76, 41), ncol = 2, byrow = TRUE)

  expect_message(
    coords_utm <- project_coords(coords, to_crs = crs.long2utm(-75, 40)),
    "Assuming source CRS is EPSG:4326"
  )

  expect_true(is.matrix(coords_utm))
})

test_that("project_coords warns when from_crs ignored for spatial objects", {
  skip_if_not_installed("sf")

  pts <- sf::st_as_sf(
    data.frame(x = c(-75, -76), y = c(40, 41)),
    coords = c("x", "y"),
    crs = 4326
  )

  expect_warning(
    project_coords(pts, from_crs = "EPSG:3857", to_crs = crs.long2utm(-75, 40)),
    "'from_crs' is ignored for sf objects"
  )
})

test_that("project_coords validates input type", {
  expect_error(
    project_coords("invalid", to_crs = "EPSG:3857"),
    "must be a SpatRaster, SpatVector, sf object, matrix, or data.frame"
  )
})
