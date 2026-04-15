# Integration Tests for Main Interface Functions
# Tests for auto_build_model() and quick_model()

test_that("auto_build_model rejects legacy raster objects", {
  skip_if_not_installed("raster")
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # Create legacy raster object
  r_old <- raster::raster(ncol = 10, nrow = 10)
  raster::values(r_old) <- 1:100
  
  # Create valid sf domain
  coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0), ncol = 2, byrow = TRUE)
  poly <- sf::st_polygon(list(coords))
  domain <- sf::st_sf(geometry = sf::st_sfc(poly))
  
  # Should reject legacy raster with clear error
  expect_error(
    auto_build_model(
      project_name = "test",
      domain = domain,
      dem = r_old
    ),
    "Legacy raster objects are not supported"
  )
  
  expect_error(
    auto_build_model(
      project_name = "test",
      domain = domain,
      dem = r_old
    ),
    "terra::rast"
  )
})

test_that("auto_build_model rejects legacy sp objects", {
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # Create legacy sp object
  coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0), ncol = 2, byrow = TRUE)
  poly <- sp::Polygon(coords)
  polys <- sp::Polygons(list(poly), ID = "1")
  sp_obj <- sp::SpatialPolygons(list(polys))
  
  # Create valid terra DEM
  dem <- terra::rast(ncol = 10, nrow = 10)
  terra::values(dem) <- 1:100
  
  # Should reject legacy sp with clear error
  expect_error(
    auto_build_model(
      project_name = "test",
      domain = sp_obj,
      dem = dem
    ),
    "Legacy sp objects are not supported"
  )
  
  expect_error(
    auto_build_model(
      project_name = "test",
      domain = sp_obj,
      dem = dem
    ),
    "sf::st_read"
  )
})

test_that("auto_build_model accepts terra/sf objects", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  skip_if_s2_unavailable()
  
  # Create small test data
  coords <- matrix(c(0, 0, 100, 0, 100, 100, 0, 100, 0, 0), ncol = 2, byrow = TRUE)
  poly <- sf::st_polygon(list(coords))
  domain <- sf::st_sf(id = 1, geometry = sf::st_sfc(poly, crs = 4326))
  
  # Create DEM
  dem <- terra::rast(ncol = 10, nrow = 10, 
                     xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                     crs = "EPSG:4326")
  terra::values(dem) <- seq(100, 200, length.out = 100)
  
  # Create temporary output directory
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  # Should work with terra/sf
  expect_no_error({
    result <- auto_build_model(
      project_name = "test_model",
      domain = domain,
      dem = dem,
      output_dir = temp_dir,
      years = 2000:2001,
      verbose = FALSE
    )
  })
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("auto_build_model validates required parameters", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # Create valid objects
  coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0), ncol = 2, byrow = TRUE)
  poly <- sf::st_polygon(list(coords))
  domain <- sf::st_sf(geometry = sf::st_sfc(poly))
  
  dem <- terra::rast(ncol = 10, nrow = 10)
  terra::values(dem) <- 1:100
  
  # Missing project_name
  expect_error(
    auto_build_model(domain = domain, dem = dem),
    "project_name|missing, with no default"
  )
  
  # Missing domain
  expect_error(
    auto_build_model(project_name = "test", dem = dem),
    "domain"
  )
  
  # Missing dem
  expect_error(
    auto_build_model(project_name = "test", domain = domain),
    "dem"
  )
})

test_that("quick_model works with minimal inputs", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  skip_if_s2_unavailable()
  
  # Create small test data
  coords <- matrix(c(0, 0, 100, 0, 100, 100, 0, 100, 0, 0), ncol = 2, byrow = TRUE)
  poly <- sf::st_polygon(list(coords))
  domain <- sf::st_sf(id = 1, geometry = sf::st_sfc(poly, crs = 4326))
  
  # Create DEM
  dem <- terra::rast(ncol = 10, nrow = 10,
                     xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                     crs = "EPSG:4326")
  terra::values(dem) <- seq(100, 200, length.out = 100)
  
  # Create temporary output directory
  temp_dir <- tempfile()
  dir.create(temp_dir)
  result <- NULL
  
  # Should work with minimal inputs
  expect_no_error({
    result <- quick_model(
      project_name = "quick_test",
      domain = domain,
      dem = dem,
      output_dir = temp_dir,
      years = 2000:2001,
      verbose = FALSE
    )
  })
  if (is.null(result)) {
    return(invisible())
  }
  
  # Check result structure
  expect_type(result, "list")
  expect_true("mesh" %in% names(result))
  expect_true("files" %in% names(result))
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("quick_model validates spatial object types", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # Create valid objects
  coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0), ncol = 2, byrow = TRUE)
  poly <- sf::st_polygon(list(coords))
  domain <- sf::st_sf(geometry = sf::st_sfc(poly))
  
  dem <- terra::rast(ncol = 10, nrow = 10)
  terra::values(dem) <- 1:100
  
  # Invalid domain type
  expect_error(
    quick_model(
      project_name = "test",
      domain = "not_sf",
      dem = dem
    ),
    "domain.*must be an sf object"
  )
  
  # Invalid dem type
  expect_error(
    quick_model(
      project_name = "test",
      domain = domain,
      dem = "not_raster"
    ),
    "dem.*must be a terra SpatRaster"
  )
})

test_that("autoBuildModel deprecated function shows migration guidance", {
  expect_error(suppressWarnings(autoBuildModel()), "deprecated|MIGRATION|auto_build_model")
})
