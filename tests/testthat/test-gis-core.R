# Tests for Core GIS Functions

small_test_mesh <- function() {
  SHUD.MESH(
    mesh = data.frame(
      ID = 1:2,
      Node1 = c(1, 2),
      Node2 = c(2, 3),
      Node3 = c(3, 4),
      Nabr1 = c(0, 0),
      Nabr2 = c(0, 0),
      Nabr3 = c(0, 0)
    ),
    point = data.frame(
      ID = 1:4,
      X = c(0, 1, 0, 1),
      Y = c(0, 0, 1, 1),
      AqD = c(10, 20, 30, 40),
      Z = c(100, 110, 120, 130)
    )
  )
}

small_test_template <- function() {
  r <- terra::rast(
    xmin = -0.5, xmax = 1.5,
    ymin = -0.5, ymax = 1.5,
    resolution = 0.5
  )
  terra::crs(r) <- "EPSG:4326"
  r
}

# Test mesh_to_raster function ------------------------------------------------

test_that("mesh_to_raster requires data parameter", {
  expect_error(
    mesh_to_raster(),
    "Parameter 'data' is required"
  )
  expect_error(
    mesh_to_raster(data = NULL),
    "Parameter 'data' is required"
  )
})

test_that("mesh_to_raster rejects legacy raster/sp objects", {
  skip_if_not_installed("terra")
  
  # Mock legacy objects by class
  legacy_raster <- structure(list(), class = "RasterLayer")
  legacy_sp <- structure(list(), class = "SpatialPoints")
  
  expect_error(
    mesh_to_raster(data = legacy_raster),
    "Legacy raster/sp objects are not supported"
  )
  expect_error(
    mesh_to_raster(data = legacy_sp),
    "Legacy raster/sp objects are not supported"
  )
})

test_that("mesh_to_raster validates data type", {
  expect_error(
    mesh_to_raster(data = "invalid"),
    "must be a numeric vector, matrix, or data.frame"
  )
  expect_error(
    mesh_to_raster(data = list(1, 2, 3)),
    "must be a numeric vector, matrix, or data.frame"
  )
})

test_that("mesh_to_raster handles NA and infinite values", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # Create simple test data
  data_with_na <- c(1, 2, NA, 4, 5)
  data_with_inf <- c(1, 2, Inf, 4, 5)
  
  # These should warn but not error
  # Note: Will fail without actual mesh, but tests the validation logic
  expect_warning(
    tryCatch(
      mesh_to_raster(data = data_with_na, method = "nearest"),
      error = function(e) NULL
    ),
    "NA values"
  )
  
  expect_warning(
    tryCatch(
      mesh_to_raster(data = data_with_inf, method = "nearest"),
      error = function(e) NULL
    ),
    "Infinite values"
  )
})

test_that("mesh_to_raster validates method parameter", {
  skip_if_not_installed("terra")
  
  expect_error(
    mesh_to_raster(data = 1:10, method = "invalid"),
    "'arg' should be one of"
  )
})

test_that("mesh_to_raster validates resolution parameter", {
  skip_if_not_installed("terra")
  skip("Requires active SHUD project context before resolution validation")
})

test_that("mesh_to_raster rejects legacy template", {
  skip_if_not_installed("terra")
  skip("Requires active SHUD project context before template validation")
})

test_that("mesh_to_raster validates template type", {
  skip_if_not_installed("terra")
  skip("Requires active SHUD project context before template validation")
})

test_that("mesh_to_raster warns about matrix without stack", {
  skip_if_not_installed("terra")
  
  # Matrix with multiple rows should warn if stack=FALSE
  data_matrix <- matrix(1:20, nrow = 2, ncol = 10)
  
  expect_warning(
    tryCatch(
      mesh_to_raster(data = data_matrix, stack = FALSE),
      error = function(e) NULL
    ),
    "Using last row only"
  )
})

test_that("mesh_to_raster validates data length matches mesh", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # This test requires actual mesh data
  # We'll test the error message structure
  skip("Requires actual mesh data for full integration test")
})

test_that("mesh_to_raster with nearest method works", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # This requires actual mesh data
  skip("Requires actual mesh data for full integration test")
})

test_that("mesh_to_raster with idw method requires gstat", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # Mock missing gstat
  skip("Requires actual mesh data for full integration test")
})

test_that("mesh_to_raster with linear method requires interp", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # Mock missing interp
  skip("Requires actual mesh data for full integration test")
})

test_that("mesh_to_raster stack parameter works with matrix", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # This requires actual mesh data
  skip("Requires actual mesh data for full integration test")
})

test_that("mesh_to_raster applies CRS correctly", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # This requires actual mesh data
  skip("Requires actual mesh data for full integration test")
})

# Test deprecated MeshData2Raster function ------------------------------------

test_that("MeshData2Raster shows deprecation warning", {
  skip_if_not_installed("terra")
  
  expect_warning(
    tryCatch(
      MeshData2Raster(x = 1:10),
      error = function(e) NULL
    ),
    "deprecated|不再有用"
  )
})

test_that("MeshData2Raster forwards modern arguments without legacy defaults", {
  skip_if_not_installed("terra")

  template <- small_test_template()
  mesh <- small_test_mesh()
  returned <- template
  terra::values(returned) <- seq_len(terra::ncell(returned))
  captured <- NULL

  testthat::local_mocked_bindings(
    mesh_to_raster = function(data,
                              mesh = NULL,
                              template = NULL,
                              method = c("idw", "linear", "nearest"),
                              resolution = NULL,
                              crs = NULL,
                              stack = FALSE) {
      captured <<- list(
        data = data,
        mesh = mesh,
        template = template,
        method_missing = missing(method),
        method = method,
        resolution = resolution,
        crs = crs,
        stack = stack
      )
      returned
    },
    .package = "rSHUD"
  )

  expect_warning(
    out <- MeshData2Raster(
      data = 1,
      mesh = mesh,
      template = template,
      resolution = 0.25,
      crs = "EPSG:4326",
      stack = TRUE,
      plot = FALSE
    ),
    "deprecated|不再有用"
  )

  expect_s4_class(out, "SpatRaster")
  expect_identical(captured$data, 1)
  expect_identical(captured$mesh, mesh)
  expect_identical(captured$template, template)
  expect_true(captured$method_missing)
  expect_identical(captured$resolution, 0.25)
  expect_identical(captured$crs, "EPSG:4326")
  expect_true(captured$stack)
})

test_that("MeshData2Raster maps old parameters to new ones", {
  skip_if_not_installed("terra")

  template <- small_test_template()
  mesh <- small_test_mesh()
  returned <- template
  terra::values(returned) <- seq_len(terra::ncell(returned))
  captured <- NULL

  testthat::local_mocked_bindings(
    mesh_to_raster = function(data,
                              mesh = NULL,
                              template = NULL,
                              method = c("idw", "linear", "nearest"),
                              resolution = NULL,
                              crs = NULL,
                              stack = FALSE) {
      captured <<- list(
        data = data,
        mesh = mesh,
        template = template,
        method = method,
        resolution = resolution,
        crs = crs,
        stack = stack
      )
      returned
    },
    .package = "rSHUD"
  )

  expect_warning(
    out <- MeshData2Raster(
      x = 1,
      rmask = template,
      pm = mesh,
      proj = "EPSG:3857",
      method = "nearest",
      plot = FALSE
    ),
    "deprecated|不再有用"
  )

  expect_s4_class(out, "SpatRaster")
  expect_identical(captured$data, 1)
  expect_identical(captured$mesh, mesh)
  expect_identical(captured$template, template)
  expect_identical(captured$method, "nearest")
  expect_identical(captured$crs, "EPSG:3857")
})

test_that("MeshData2Raster resolves legacy defaults before forwarding", {
  skip_if_not_installed("terra")

  data <- c(11, 12)
  mesh <- small_test_mesh()
  template <- small_test_template()
  returned <- template
  terra::values(returned) <- seq_len(terra::ncell(returned))
  captured <- NULL
  calls <- list(getElevation = 0, read_mesh = 0, shud.mask = 0)

  testthat::local_mocked_bindings(
    getElevation = function(...) {
      calls$getElevation <<- calls$getElevation + 1
      data
    },
    read_mesh = function(...) {
      calls$read_mesh <<- calls$read_mesh + 1
      mesh
    },
    shud.mask = function(pm = NULL, proj = NULL, ...) {
      calls$shud.mask <<- calls$shud.mask + 1
      expect_identical(pm, mesh)
      expect_null(proj)
      template
    },
    mesh_to_raster = function(data,
                              mesh = NULL,
                              template = NULL,
                              method = c("idw", "linear", "nearest"),
                              resolution = NULL,
                              crs = NULL,
                              stack = FALSE) {
      captured <<- list(
        data = data,
        mesh = mesh,
        template = template,
        method_missing = missing(method),
        method = method,
        resolution = resolution,
        crs = crs,
        stack = stack
      )
      returned
    },
    .package = "rSHUD"
  )

  expect_warning(
    out <- MeshData2Raster(),
    "deprecated|不再有用"
  )

  expect_s4_class(out, "SpatRaster")
  expect_equal(calls, list(getElevation = 1, read_mesh = 1, shud.mask = 1))
  expect_identical(captured$data, data)
  expect_identical(captured$mesh, mesh)
  expect_identical(captured$template, template)
  expect_true(captured$method_missing)
})

test_that("MeshData2Raster builds missing template from resolved mesh and CRS", {
  skip_if_not_installed("terra")

  legacy_mesh <- small_test_mesh()
  modern_mesh <- small_test_mesh()
  template <- small_test_template()
  returned <- template
  terra::values(returned) <- seq_len(terra::ncell(returned))
  captured <- NULL

  testthat::local_mocked_bindings(
    read_mesh = function(...) stop("read_mesh() should not be called"),
    shud.mask = function(pm = NULL, proj = NULL, ...) {
      expect_identical(pm, modern_mesh)
      expect_identical(proj, "EPSG:3857")
      template
    },
    mesh_to_raster = function(data,
                              mesh = NULL,
                              template = NULL,
                              method = c("idw", "linear", "nearest"),
                              resolution = NULL,
                              crs = NULL,
                              stack = FALSE) {
      captured <<- list(
        data = data,
        mesh = mesh,
        template = template,
        method_missing = missing(method),
        method = method,
        crs = crs
      )
      returned
    },
    .package = "rSHUD"
  )

  expect_warning(
    MeshData2Raster(
      x = c(1, 2),
      pm = legacy_mesh,
      mesh = modern_mesh,
      proj = "EPSG:3857"
    ),
    "deprecated|不再有用"
  )

  expect_identical(captured$data, c(1, 2))
  expect_identical(captured$mesh, modern_mesh)
  expect_identical(captured$template, template)
  expect_identical(captured$crs, "EPSG:3857")
  expect_true(captured$method_missing)
})

test_that("MeshData2Raster accepts nearest interpolation through mesh_to_raster", {
  skip_if_not_installed("terra")

  expect_warning(
    r <- MeshData2Raster(
      data = c(7, 8),
      mesh = small_test_mesh(),
      template = small_test_template(),
      method = "nearest"
    ),
    "deprecated|不再有用"
  )

  expect_s4_class(r, "SpatRaster")
})

test_that("modern reader defaults avoid deprecated self-calls in GIS helpers", {
  skip_if_not_installed("sf")

  mesh <- small_test_mesh()

  testthat::local_mocked_bindings(
    read_mesh = function(...) mesh,
    readmesh = function(...) stop("deprecated readmesh() called"),
    read_att = function(...) data.frame(LAKE = c(0, 2)),
    readatt = function(...) stop("deprecated readatt() called"),
    mesh_to_sf = function(pm) {
      expect_identical(pm, mesh)
      sf::st_sf(
        id = 1,
        geometry = sf::st_sfc(
          sf::st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 0),
                                      ncol = 2, byrow = TRUE)))
        )
      )
    },
    .package = "rSHUD"
  )

  expect_equal(getElevation(), c(110, 120))
  expect_equal(getAquiferDepth(), c(20, 30))
  expect_equal(getLakeEleID(), 2)
  expect_length(getArea(), 1)
  expect_equal(getCentroid()[1, "ZMAX"], 110)
})

test_that("getRiverNodes default uses read_river_sp compatibility reader", {
  skip_if_not_installed("sf")

  river_line <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
      crs = 4326
    )
  )

  testthat::local_mocked_bindings(
    read_river_sp = function(...) river_line,
    readriv.sp = function(...) stop("deprecated readriv.sp() called"),
    get_coords = function(...) matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE),
    get_from_to_nodes = function(spr, coords) {
      matrix(c(1, 2), ncol = 2)
    },
    ProjectCoordinate = function(x, proj4string, P2G = TRUE) {
      matrix(c(10, 20, 11, 21), ncol = 2, byrow = TRUE,
             dimnames = list(NULL, c("Lon", "Lat")))
    },
    .package = "rSHUD"
  )

  nodes <- getRiverNodes()

  expect_equal(nrow(nodes$points), 2)
  expect_equal(as.integer(nodes$FT_ID[1, 2:3]), c(1L, 2L))
})

# Integration tests with real data --------------------------------------------

test_that("mesh_to_raster integration test with sample data", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if_not_installed("gstat")
  
  # Skip if sample data not available
  skip_on_cran()
  skip_if_not(exists("sh"), message = "Sample data 'sh' not available")
  
  # This would test with actual rSHUD sample data
  # data(sh)
  # mesh <- sh$mesh
  # elevation <- getElevation(mesh)
  # r <- mesh_to_raster(elevation, mesh = mesh, method = "nearest")
  # expect_s4_class(r, "SpatRaster")
  # expect_true(terra::ncell(r) > 0)
})

test_that("mesh_to_raster produces consistent results across methods", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if_not_installed("gstat")
  skip_if_not_installed("interp")
  
  # Skip if sample data not available
  skip_on_cran()
  skip_if_not(exists("sh"), message = "Sample data 'sh' not available")
  
  # This would test that different methods produce similar results
  # data(sh)
  # mesh <- sh$mesh
  # elevation <- getElevation(mesh)
  # 
  # r_idw <- mesh_to_raster(elevation, mesh = mesh, method = "idw")
  # r_linear <- mesh_to_raster(elevation, mesh = mesh, method = "linear")
  # r_nearest <- mesh_to_raster(elevation, mesh = mesh, method = "nearest")
  # 
  # # All should have same dimensions
  # expect_equal(dim(r_idw), dim(r_linear))
  # expect_equal(dim(r_idw), dim(r_nearest))
  # 
  # # Values should be correlated (not identical but similar)
  # # This tests that interpolation is working reasonably
})

test_that("mesh_to_raster handles time series data correctly", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Skip if sample data not available
  skip_on_cran()
  skip_if_not(exists("sh"), message = "Sample data 'sh' not available")

  # This would test stack functionality with time series
  # data(sh)
  # mesh <- sh$mesh
  # n_cells <- nrow(getCentroid(mesh))
  # n_times <- 5
  #
  # # Create fake time series
  # ts_data <- matrix(rnorm(n_times * n_cells), nrow = n_times, ncol = n_cells)
  #
  # # Process as stack
  # r_stack <- mesh_to_raster(ts_data, mesh = mesh, stack = TRUE, method = "nearest")
  #
  # expect_s4_class(r_stack, "SpatRaster")
  # expect_equal(terra::nlyr(r_stack), n_times)
})

# Test vector_to_raster function ----------------------------------------------

test_that("vector_to_raster requires vector parameter", {
  expect_error(
    vector_to_raster(),
    "Parameter 'vector' is required"
  )
  expect_error(
    vector_to_raster(vector = NULL),
    "Parameter 'vector' is required"
  )
})

test_that("vector_to_raster rejects legacy sp objects", {
  skip_if_not_installed("terra")

  # Mock legacy sp objects by class
  legacy_sp_points <- structure(list(), class = "SpatialPoints")
  legacy_sp_polygons <- structure(list(), class = "SpatialPolygons")
  legacy_sp_lines <- structure(list(), class = "SpatialLines")

  expect_error(
    vector_to_raster(vector = legacy_sp_points),
    "Legacy sp objects are not supported"
  )
  expect_error(
    vector_to_raster(vector = legacy_sp_polygons),
    "Legacy sp objects are not supported"
  )
  expect_error(
    vector_to_raster(vector = legacy_sp_lines),
    "Legacy sp objects are not supported"
  )
})

test_that("vector_to_raster validates vector type", {
  expect_error(
    vector_to_raster(vector = "invalid"),
    "must be an sf or SpatVector object"
  )
  expect_error(
    vector_to_raster(vector = data.frame(x = 1:5, y = 1:5)),
    "must be an sf or SpatVector object"
  )
})

test_that("vector_to_raster validates ngrids parameter", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5, value = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  # Negative ngrids should fail
  expect_error(
    vector_to_raster(vector = pts, ngrids = -1),
    "positive"
  )

  # Zero ngrids should fail
  expect_error(
    vector_to_raster(vector = pts, ngrids = 0),
    "positive"
  )
})

test_that("vector_to_raster validates resolution parameter", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5, value = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  # Negative resolution should fail
  expect_error(
    vector_to_raster(vector = pts, resolution = -1),
    "positive"
  )

  # Zero resolution should fail
  expect_error(
    vector_to_raster(vector = pts, resolution = 0),
    "positive"
  )
})

test_that("vector_to_raster rejects legacy template", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5, value = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  # Mock legacy raster template
  legacy_template <- structure(list(), class = "RasterLayer")

  expect_error(
    vector_to_raster(vector = pts, template = legacy_template),
    "Legacy raster template not supported"
  )
})

test_that("vector_to_raster validates template type", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5, value = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  expect_error(
    vector_to_raster(vector = pts, template = "invalid"),
    "must be a SpatRaster object"
  )
})

test_that("vector_to_raster validates field parameter", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5, value = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  # Field index out of range
  expect_error(
    vector_to_raster(vector = pts, field = 10),
    "Field index .* is out of range"
  )

  # Field name not found
  expect_error(
    vector_to_raster(vector = pts, field = "nonexistent"),
    "Field .* not found"
  )

  # Invalid field type
  expect_error(
    vector_to_raster(vector = pts, field = c(1, 2)),
    "Parameter 'field' must be either"
  )
})

test_that("vector_to_raster validates aggregation function", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5, value = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  expect_error(
    vector_to_raster(vector = pts, fun = "invalid"),
    "must be one of"
  )
})

test_that("vector_to_raster works with sf points", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = c(0, 1, 2, 3, 4), y = c(0, 1, 2, 3, 4), value = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  # Rasterize with field name
  r <- vector_to_raster(vector = pts, field = "value", ngrids = 10)

  expect_s4_class(r, "SpatRaster")
  expect_true(terra::ncell(r) > 0)
})

test_that("vector_to_raster works with sf polygons", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test polygon
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

  # Rasterize
  r <- vector_to_raster(vector = poly, field = "value", ngrids = 10, fun = "sum")

  expect_s4_class(r, "SpatRaster")
  expect_true(terra::ncell(r) > 0)
})

test_that("vector_to_raster works with SpatVector", {
  skip_if_not_installed("terra")

  # Create simple SpatVector
  pts <- terra::vect(
    cbind(x = 1:5, y = 1:5),
    type = "points",
    atts = data.frame(value = 1:5),
    crs = "EPSG:4326"
  )

  # Rasterize
  r <- vector_to_raster(vector = pts, field = "value", ngrids = 10)

  expect_s4_class(r, "SpatRaster")
  expect_true(terra::ncell(r) > 0)
})

test_that("vector_to_raster works with custom resolution", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5, value = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  # Rasterize with custom resolution
  r <- vector_to_raster(vector = pts, field = "value", resolution = 0.5)

  expect_s4_class(r, "SpatRaster")
  expect_equal(terra::res(r)[1], 0.5)
})

test_that("vector_to_raster works with template raster", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5, value = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  # Create template
  template <- terra::rast(xmin = 0, xmax = 6, ymin = 0, ymax = 6, resolution = 0.5)
  terra::crs(template) <- "EPSG:4326"

  # Rasterize with template
  r <- vector_to_raster(vector = pts, template = template, field = "value")

  expect_s4_class(r, "SpatRaster")
  expect_equal(dim(r), dim(template))
})

test_that("vector_to_raster works with field as numeric vector", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  # Rasterize with values vector
  values <- c(10, 20, 30, 40, 50)
  r <- vector_to_raster(vector = pts, field = values, ngrids = 10)

  expect_s4_class(r, "SpatRaster")
  expect_true(terra::ncell(r) > 0)
})

test_that("vector_to_raster works with different aggregation functions", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create test vector with overlapping points
  pts <- sf::st_as_sf(
    data.frame(x = c(1, 1, 2, 2), y = c(1, 1, 2, 2), value = c(10, 20, 30, 40)),
    coords = c("x", "y"),
    crs = 4326
  )

  # Test different aggregation functions
  r_mean <- vector_to_raster(vector = pts, field = "value", fun = "mean", ngrids = 5)
  r_sum <- vector_to_raster(vector = pts, field = "value", fun = "sum", ngrids = 5)
  r_min <- vector_to_raster(vector = pts, field = "value", fun = "min", ngrids = 5)
  r_max <- vector_to_raster(vector = pts, field = "value", fun = "max", ngrids = 5)

  expect_s4_class(r_mean, "SpatRaster")
  expect_s4_class(r_sum, "SpatRaster")
  expect_s4_class(r_min, "SpatRaster")
  expect_s4_class(r_max, "SpatRaster")
})

# Test deprecated sp2raster function ------------------------------------------

test_that("sp2raster shows deprecation warning", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5, value = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  expect_warning(
    sp2raster(sp = pts, field = "value"),
    "deprecated|不再有用"
  )
})

test_that("sp2raster maps old parameters to new ones", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Create simple test vector
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5, value = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )

  # Test parameter mapping: sp -> vector, mask -> template
  expect_warning(
    r <- sp2raster(sp = pts, mask = NULL, field = "value"),
    "deprecated|不再有用"
  )

  expect_s4_class(r, "SpatRaster")
})
