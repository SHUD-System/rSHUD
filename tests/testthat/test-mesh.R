# Tests for Mesh Generation Module

if (can_load_namespace("sf")) {
  sf::sf_use_s2(FALSE)
}

# Test shud.triangle function -------------------------------------------------

test_that("shud.triangle requires wb parameter", {
  expect_error(
    shud.triangle(),
    "Parameter 'wb'"
  )
})

test_that("shud.triangle accepts sf objects", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  
  # Create simple boundary
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  
  # Generate mesh
  mesh <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  expect_true(is.list(mesh))
  expect_true(all(c("P", "T") %in% names(mesh)))
  expect_true(nrow(mesh$P) > 0)
  expect_true(nrow(mesh$T) > 0)
})

test_that("shud.triangle validates geometry type", {
  skip_if_not_installed("sf")
  
  # Create POINT geometry (invalid)
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5),
    coords = c("x", "y")
  )
  
  expect_error(
    shud.triangle(wb = pts),
    "must have POLYGON or MULTIPOLYGON geometry"
  )
})

test_that("shud.triangle validates q parameter", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  
  # Create simple boundary
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  
  # Negative angle should fail
  expect_error(
    shud.triangle(wb = boundary, q = -1),
    "must be positive"
  )
  
  # Angle > 35 should warn and adjust
  expect_warning(
    shud.triangle(wb = boundary, q = 40),
    "may cause triangulation issues"
  )
})

test_that("shud.triangle handles holes", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  
  # Create boundary with hole
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  
  hole <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(3,3, 7,3, 7,7, 3,7, 3,3),
                                           ncol=2, byrow=TRUE))))
  )
  
  # Generate mesh with hole
  mesh <- shud.triangle(wb = boundary, hole = hole, q = 30)
  
  expect_true(is.list(mesh))
  expect_true(nrow(mesh$T) > 0)
})

test_that("shud.triangle works with rivers", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  
  # Create boundary
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  
  # Create river
  river <- sf::st_as_sf(
    sf::st_sfc(sf::st_linestring(matrix(c(2,2, 8,8), ncol=2, byrow=TRUE)))
  )
  
  # Generate mesh with river
  mesh <- shud.triangle(wb = boundary, riv = river, q = 30, a = 5)
  
  expect_true(is.list(mesh))
  expect_true(nrow(mesh$T) > 0)
})

test_that("shud.triangle works with additional constraint points", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  
  # Create boundary
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  
  # Additional points
  pts <- matrix(c(5, 5, 3, 7, 7, 3), ncol=2, byrow=TRUE)
  
  # Generate mesh with constraint points
  mesh <- shud.triangle(wb = boundary, pts = pts, q = 30, a = 5)
  
  expect_true(is.list(mesh))
  expect_true(nrow(mesh$T) > 0)
})

# Test sp.mesh2Shape function -------------------------------------------------

test_that("sp.mesh2Shape handles triangulation objects", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  
  # Create simple mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  # Convert to sf
  mesh_sf <- sp.mesh2Shape(pm = tri)
  
  expect_true(inherits(mesh_sf, "sf"))
  expect_true("Area" %in% colnames(mesh_sf))
})

test_that("sp.mesh2Shape validates input type", {
  expect_error(
    sp.mesh2Shape(pm = "invalid"),
    "must be either a SHUD.MESH object"
  )
})

test_that("sp.mesh2Shape applies CRS correctly", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  
  # Create simple mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  # Convert with CRS
  mesh_sf <- sp.mesh2Shape(pm = tri, crs = 4326)
  
  expect_equal(sf::st_crs(mesh_sf)$epsg, 4326)
})

test_that("sp.mesh2Shape handles SHUD.MESH objects", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create simple mesh and DEM
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  # Create DEM
  dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(dem) <- runif(100, 100, 200)
  
  # Create SHUD.MESH object
  mesh_obj <- shud.mesh(tri, dem = dem, AqDepth = 10)
  
  # Convert to sf
  mesh_sf <- sp.mesh2Shape(pm = mesh_obj)
  
  expect_true(inherits(mesh_sf, "sf"))
  expect_true("Area" %in% colnames(mesh_sf))
  expect_true("AqDepth" %in% colnames(mesh_sf))
  expect_true("Zsurf" %in% colnames(mesh_sf))
})

test_that("sp.mesh2Shape with custom attributes", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  
  # Create simple mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  # Custom attributes
  attrs <- data.frame(
    ID = seq_len(nrow(tri$T)),
    value = runif(nrow(tri$T))
  )
  
  # Convert with attributes
  mesh_sf <- sp.mesh2Shape(pm = tri, dbf = attrs)
  
  expect_true(inherits(mesh_sf, "sf"))
  expect_true("value" %in% colnames(mesh_sf))
})

# Test Tri2Centroid function --------------------------------------------------

test_that("Tri2Centroid validates input", {
  expect_error(
    Tri2Centroid("invalid"),
    "must have T .* and P .* components"
  )
})

test_that("Tri2Centroid computes centroids correctly", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  
  # Create simple mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  # Calculate centroids
  centroids <- Tri2Centroid(tri)
  
  expect_true(is.matrix(centroids))
  expect_equal(ncol(centroids), 2)
  expect_equal(nrow(centroids), nrow(tri$T))
})

test_that("Tri2Centroid centroids are within bounds", {
  skip_if_not_installed("sf")
  skip_if_not_installed("RTriangle")
  
  # Create simple mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  # Calculate centroids
  centroids <- Tri2Centroid(tri)
  
  # Centroids should be within boundary
  expect_true(all(centroids[, 1] >= 0 & centroids[, 1] <= 10))
  expect_true(all(centroids[, 2] >= 0 & centroids[, 2] <= 10))
})

# Test shud.mesh function -----------------------------------------------------

test_that("shud.mesh requires tri and dem parameters", {
  expect_error(
    shud.mesh(),
    "Parameter 'tri'|missing|tri|缺少"
  )
})

test_that("shud.mesh creates SHUD.MESH object", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh and DEM
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(dem) <- runif(100, 100, 200)
  
  # Create SHUD.MESH
  mesh_obj <- shud.mesh(tri, dem = dem, AqDepth = 10)
  
  expect_s4_class(mesh_obj, "Untructure Domain")
  expect_true(is.data.frame(mesh_obj@mesh))
  expect_true(is.data.frame(mesh_obj@point))
})

test_that("shud.mesh validates triangulation structure", {
  skip_if_not_installed("terra")
  
  # Invalid triangulation
  dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(dem) <- runif(100, 100, 200)
  
  expect_error(
    shud.mesh(tri = list(T = matrix(1:9, 3, 3)), dem = dem),
    "must be a triangulation object"
  )
})

test_that("shud.mesh mesh data frame has correct structure", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(dem) <- runif(100, 100, 200)
  
  mesh_obj <- shud.mesh(tri, dem = dem, AqDepth = 10)
  
  # Check mesh structure
  expect_true(all(c("ID", "Node1", "Node2", "Node3", "Nabr1", "Nabr2", "Nabr3", "Zmax") 
                  %in% colnames(mesh_obj@mesh)))
  expect_equal(nrow(mesh_obj@mesh), nrow(tri$T))
})

test_that("shud.mesh point data frame has correct structure", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(dem) <- runif(100, 100, 200)
  
  mesh_obj <- shud.mesh(tri, dem = dem, AqDepth = 10)
  
  # Check point structure
  expect_true(all(c("ID", "X", "Y", "AqDepth", "Elevation") 
                  %in% colnames(mesh_obj@point)))
  expect_equal(nrow(mesh_obj@point), nrow(tri$P))
})

test_that("shud.mesh handles spatially variable aquifer depth", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(dem) <- runif(100, 100, 200)
  
  # Variable aquifer depth
  aq_raster <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(aq_raster) <- runif(100, 5, 20)
  
  mesh_obj <- shud.mesh(tri, dem = dem, r.aq = aq_raster)
  
  # Check aquifer depth varies
  expect_true(length(unique(mesh_obj@point$AqDepth)) > 1)
})

# Test shud.att function ------------------------------------------------------

test_that("shud.att creates attribute data frame", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  # Create attributes
  att <- shud.att(tri)
  
  expect_true(is.data.frame(att))
  expect_true(all(c("INDEX", "SOIL", "GEOL", "LC", "FORC", "MF", "BC", "SS", "LAKE") 
                  %in% colnames(att)))
  expect_equal(nrow(att), nrow(tri$T))
})

test_that("shud.att extracts soil attributes", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  # Create soil raster
  soil <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(soil) <- sample(1:5, 100, replace=TRUE)
  
  # Extract attributes
  att <- shud.att(tri, r.soil = soil)
  
  expect_true(all(att$SOIL %in% 1:5))
})

test_that("shud.att handles multiple attribute rasters", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  # Create attribute rasters
  soil <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(soil) <- sample(1:5, 100, replace=TRUE)
  
  geol <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(geol) <- sample(1:3, 100, replace=TRUE)
  
  lc <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(lc) <- sample(1:10, 100, replace=TRUE)
  
  # Extract attributes
  att <- shud.att(tri, r.soil = soil, r.geol = geol, r.lc = lc)
  
  expect_true(all(att$SOIL %in% 1:5))
  expect_true(all(att$GEOL %in% 1:3))
  expect_true(all(att$LC %in% 1:10))
})

test_that("shud.att handles lake polygons", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  # Create lake polygon
  lake <- sf::st_as_sf(
    data.frame(
      ID = 1,
      geometry = sf::st_sfc(
        sf::st_polygon(list(matrix(c(3,3, 7,3, 7,7, 3,7, 3,3), ncol=2, byrow=TRUE)))
      )
    )
  )
  sf::st_crs(lake) <- NA
  
  # Extract attributes
  att <- shud.att(tri, sp.lake = lake)
  
  # Some cells should be marked as lake
  expect_true(any(att$LAKE > 0))
})

# Test S4 class compatibility -------------------------------------------------

test_that("SHUD.MESH S4 class is properly defined", {
  expect_true(methods::isClass("Untructure Domain"))
})

test_that("SHUD.MESH objects can be created", {
  mesh_df <- data.frame(
    ID = 1:3,
    Node1 = c(1, 2, 3),
    Node2 = c(2, 3, 4),
    Node3 = c(3, 4, 5),
    Nabr1 = c(0, 1, 2),
    Nabr2 = c(2, 3, 0),
    Nabr3 = c(3, 0, 1),
    Zmax = c(100, 110, 120)
  )
  
  point_df <- data.frame(
    ID = 1:5,
    X = c(0, 1, 2, 3, 4),
    Y = c(0, 1, 2, 3, 4),
    AqDepth = rep(10, 5),
    Elevation = c(100, 105, 110, 115, 120)
  )
  
  mesh_obj <- SHUD.MESH(mesh = mesh_df, point = point_df)
  
  expect_s4_class(mesh_obj, "Untructure Domain")
  expect_equal(mesh_obj@mesh, mesh_df)
  expect_equal(mesh_obj@point, point_df)
})

test_that("SHUD.MESH slots can be accessed", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(dem) <- runif(100, 100, 200)
  
  mesh_obj <- shud.mesh(tri, dem = dem, AqDepth = 10)
  
  # Access slots
  expect_true(is.data.frame(mesh_obj@mesh))
  expect_true(is.data.frame(mesh_obj@point))
  expect_true(nrow(mesh_obj@mesh) > 0)
  expect_true(nrow(mesh_obj@point) > 0)
})

# Test getter functions -------------------------------------------------------

test_that("getElevation extracts elevation from SHUD.MESH", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(dem) <- runif(100, 100, 200)
  
  mesh_obj <- shud.mesh(tri, dem = dem, AqDepth = 10)
  
  # Get elevation
  elev <- getElevation(mesh_obj)
  
  expect_true(is.numeric(elev))
  expect_equal(length(elev), nrow(mesh_obj@mesh))
  expect_true(all(elev >= 100 & elev <= 200))
})

test_that("getAquiferDepth extracts aquifer depth from SHUD.MESH", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(dem) <- runif(100, 100, 200)
  
  mesh_obj <- shud.mesh(tri, dem = dem, AqDepth = 15)
  
  # Get aquifer depth
  aq_depth <- getAquiferDepth(mesh_obj)
  
  expect_true(is.numeric(aq_depth))
  expect_equal(length(aq_depth), nrow(mesh_obj@mesh))
  expect_true(all(abs(aq_depth - 15) < 0.01))
})

test_that("getCentroid computes mesh centroids", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(dem) <- runif(100, 100, 200)
  
  mesh_obj <- shud.mesh(tri, dem = dem, AqDepth = 10)
  
  # Get centroids
  centroids <- getCentroid(mesh_obj)
  
  expect_true(is.matrix(centroids) || is.data.frame(centroids))
  expect_equal(nrow(centroids), nrow(mesh_obj@mesh))
  expect_true(all(c("X", "Y") %in% colnames(centroids)))
})

test_that("getVertex extracts mesh vertices", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("RTriangle")
  skip_if_not_installed("abind")
  
  # Create mesh
  boundary <- sf::st_as_sf(
    sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
                                           ncol=2, byrow=TRUE))))
  )
  tri <- shud.triangle(wb = boundary, q = 30, a = 5)
  
  dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
  terra::values(dem) <- runif(100, 100, 200)
  
  mesh_obj <- shud.mesh(tri, dem = dem, AqDepth = 10)
  
  # Get vertices
  vertices <- getVertex(mesh_obj)
  
  expect_true(is.array(vertices))
  expect_equal(dim(vertices)[1], nrow(mesh_obj@mesh))
  expect_equal(dim(vertices)[2], 3)  # 3 vertices per triangle
  expect_equal(dim(vertices)[3], 4)  # X, Y, AqD, ZMAX
})
