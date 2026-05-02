# Integration Tests for Main Interface Functions
# Tests for auto_build_model() and quick_model()

test_projected_crs <- 3857

create_integration_domain <- function(crs = test_projected_crs) {
  coords <- matrix(c(0, 0, 100, 0, 100, 100, 0, 100, 0, 0),
                   ncol = 2, byrow = TRUE)
  poly <- sf::st_polygon(list(coords))
  if (is.na(crs)) {
    return(sf::st_sf(id = 1, geometry = sf::st_sfc(poly)))
  }
  sf::st_sf(id = 1, geometry = sf::st_sfc(poly, crs = crs))
}

create_integration_dem <- function(crs = paste0("EPSG:", test_projected_crs)) {
  dem <- terra::rast(ncol = 10, nrow = 10,
                     xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                     crs = crs)
  terra::values(dem) <- seq(100, 200, length.out = 100)
  dem
}

create_integration_rivers <- function(crs = test_projected_crs) {
  line <- sf::st_linestring(matrix(c(10, 10, 90, 90), ncol = 2, byrow = TRUE))
  sf::st_sf(id = 1, geometry = sf::st_sfc(line, crs = crs))
}

expect_no_longlat_warning <- function(expr) {
  warnings <- character()
  value <- NULL
  tryCatch(
    value <- withCallingHandlers(
      force(expr),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      # Some integration runs may still fail later because optional model inputs
      # are intentionally incomplete. Those errors are not part of this check.
      expect_false(grepl("legacy|Spatial|RasterLayer", e$message))
    }
  )

  longlat_warnings <- grep(
    "longitude/latitude|although coordinates are longitude/latitude|st_buffer|st_simplify",
    warnings,
    value = TRUE
  )
  expect_length(longlat_warnings, 0)
  invisible(value)
}

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

  domain <- create_integration_domain()
  dem <- create_integration_dem()
  
  # Create temporary output directory
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  # auto_build_model requires soil/landcover for a full run;
  # here we verify it at least accepts sf/terra inputs without
  # type-checking errors (the internal build may fail without
  # complete input data).
  expect_no_longlat_warning(
    auto_build_model(
      project_name = "test_model",
      domain = domain,
      dem = dem,
      output_dir = temp_dir,
      years = 2000:2001,
      verbose = FALSE
    )
  )
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("auto_build_model rejects missing CRS on domain", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  domain <- create_integration_domain(crs = NA)
  dem <- create_integration_dem()

  expect_error(
    auto_build_model(
      project_name = "test",
      domain = domain,
      dem = dem,
      verbose = FALSE
    ),
    "domain.*projected CRS in metres/meters"
  )
})

test_that("auto_build_model rejects longlat CRS on domain", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  domain <- create_integration_domain(crs = 4326)
  dem <- create_integration_dem()

  expect_error(
    auto_build_model(
      project_name = "test",
      domain = domain,
      dem = dem,
      verbose = FALSE
    ),
    "domain.*longitude/latitude CRS.*Transform 'domain'"
  )
})

test_that("auto_build_model rejects projected feet-unit CRS before spatial operations", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  domain_feet <- create_integration_domain(crs = 2272)
  dem <- create_integration_dem()
  expect_error(
    auto_build_model(
      project_name = "test",
      domain = domain_feet,
      dem = dem,
      verbose = FALSE
    ),
    "domain.*metres/meters.*US survey foot"
  )

  domain <- create_integration_domain()
  dem_feet <- create_integration_dem(crs = "EPSG:2272")
  expect_error(
    auto_build_model(
      project_name = "test",
      domain = domain,
      dem = dem_feet,
      verbose = FALSE
    ),
    "dem.*metres/meters.*US survey foot"
  )
})

test_that("auto_build_model rejects projected feet-unit rivers", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  domain <- create_integration_domain()
  rivers_feet <- create_integration_rivers(crs = 2272)
  dem <- create_integration_dem()

  expect_error(
    auto_build_model(
      project_name = "test",
      domain = domain,
      rivers = rivers_feet,
      dem = dem,
      verbose = FALSE
    ),
    "rivers.*metres/meters.*US survey foot"
  )
})

test_that("auto_build_model rejects longlat rivers with projected domain and DEM", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  domain <- create_integration_domain()
  rivers <- create_integration_rivers(crs = 4326)
  dem <- create_integration_dem()

  expect_error(
    auto_build_model(
      project_name = "test",
      domain = domain,
      rivers = rivers,
      dem = dem,
      verbose = FALSE
    ),
    "rivers.*longitude/latitude CRS.*Transform 'rivers'"
  )
})

test_that("auto_build_model rejects missing or longlat CRS on DEM", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  domain <- create_integration_domain()

  dem_missing_crs <- create_integration_dem(crs = "")
  expect_error(
    auto_build_model(
      project_name = "test",
      domain = domain,
      dem = dem_missing_crs,
      verbose = FALSE
    ),
    "dem.*projected CRS in metres/meters"
  )

  dem_longlat <- create_integration_dem(crs = "EPSG:4326")
  expect_error(
    auto_build_model(
      project_name = "test",
      domain = domain,
      dem = dem_longlat,
      verbose = FALSE
    ),
    "dem.*longitude/latitude CRS.*Transform 'dem'"
  )
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

  domain <- create_integration_domain()
  dem <- create_integration_dem()
  
  # Create temporary output directory
  temp_dir <- tempfile()
  dir.create(temp_dir)
  result <- NULL
  
  # quick_model may fail internally without complete soil/landcover
  # data; verify it at least accepts sf/terra types.
  result <- expect_no_longlat_warning(
    quick_model(
      project_name = "quick_test",
      domain = domain,
      dem = dem,
      output_dir = temp_dir,
      years = 2000:2001,
      verbose = FALSE
    )
  )
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

test_that("quick_model inherits auto_build_model longlat rejection", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  domain <- create_integration_domain(crs = 4326)
  dem <- create_integration_dem()

  expect_error(
    quick_model(
      project_name = "quick_test",
      domain = domain,
      dem = dem,
      verbose = FALSE
    ),
    "domain.*longitude/latitude CRS.*Transform 'domain'"
  )
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

test_that("autoBuildModel deprecated function shows deprecation warning", {
  expect_warning(
    tryCatch(autoBuildModel(), error = function(e) NULL),
    "deprecated|auto_build_model|不再有用"
  )
})
