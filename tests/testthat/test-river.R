# Tests for River Processing Module

if (can_load_namespace("sf")) {
  sf::sf_use_s2(FALSE)
}

# Helper function to create test river network
create_test_river_network <- function() {
  skip_if_not_installed("sf")
  
  # Create a simple dendritic river network
  # Main stem with two tributaries
  coords_list <- list(
    # Main stem (outlet segment)
    matrix(c(5, 5, 5, 0), ncol = 2, byrow = TRUE),
    # Main stem (upstream segment)
    matrix(c(5, 10, 5, 5), ncol = 2, byrow = TRUE),
    # Left tributary
    matrix(c(0, 10, 5, 10), ncol = 2, byrow = TRUE),
    # Right tributary
    matrix(c(10, 10, 5, 10), ncol = 2, byrow = TRUE)
  )
  
  geoms <- lapply(coords_list, sf::st_linestring)
  rivers <- sf::st_sf(
    id = 1:4,
    geometry = sf::st_sfc(geoms, crs = 4326)
  )
  
  rivers
}

# Test calc_river_order function ---------------------------------------------

test_that("calc_river_order requires sf input", {
  expect_error(
    calc_river_order("invalid"),
    "must be an sf object"
  )
  expect_error(
    calc_river_order(data.frame(x = 1:5, y = 1:5)),
    "must be an sf object"
  )
})

test_that("calc_river_order rejects legacy sp objects", {
  skip_if_not_installed("sf")
  
  # Mock legacy sp object
  legacy_sp <- structure(list(), class = "SpatialLines")
  
  expect_error(
    calc_river_order(legacy_sp),
    "must be an sf object"
  )
})

test_that("calc_river_order validates geometry type", {
  skip_if_not_installed("sf")
  
  # Create POINT geometry (invalid)
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )
  
  expect_error(
    calc_river_order(pts),
    "must contain LINESTRING"
  )
})

test_that("calc_river_order calculates stream order correctly", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Calculate order
  order <- calc_river_order(rivers)
  
  expect_true(is.numeric(order))
  expect_equal(length(order), nrow(rivers))
  expect_true(all(order > 0))
  
  # Tributaries should be order 1
  expect_equal(order[3], 1)
  expect_equal(order[4], 1)
  
  # Main stem should be higher order
  expect_true(order[1] >= order[3])
  expect_true(order[2] >= order[3])
})

test_that("calc_river_order handles MULTILINESTRING", {
  skip_if_not_installed("sf")
  
  # Create MULTILINESTRING
  coords1 <- matrix(c(0, 0, 5, 5), ncol = 2, byrow = TRUE)
  coords2 <- matrix(c(5, 5, 10, 10), ncol = 2, byrow = TRUE)
  
  multi_line <- sf::st_multilinestring(list(coords1, coords2))
  rivers <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(multi_line, crs = 4326)
  )
  
  # Should convert and process
  order <- calc_river_order(rivers)
  
  expect_true(is.numeric(order))
  expect_equal(length(order), 2)  # Should be split into 2 linestrings
})

# Test calc_river_downstream function ----------------------------------------

test_that("calc_river_downstream requires sf input", {
  expect_error(
    calc_river_downstream("invalid"),
    "must be an sf object"
  )
})

test_that("calc_river_downstream validates geometry type", {
  skip_if_not_installed("sf")
  
  # Create POINT geometry (invalid)
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )
  
  expect_error(
    calc_river_downstream(pts),
    "must contain LINESTRING"
  )
})

test_that("calc_river_downstream identifies outlets correctly", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Calculate downstream segments
  downstream <- calc_river_downstream(rivers)
  
  expect_true(is.numeric(downstream))
  expect_equal(length(downstream), nrow(rivers))
  
  # First segment should be outlet (no downstream)
  expect_equal(downstream[1], -1)
  
  # Other segments should have downstream
  expect_true(all(downstream[2:4] > 0 | downstream[2:4] == -1))
})

test_that("calc_river_downstream with provided coords", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Get coordinates
  coords <- get_coords(rivers)
  
  # Calculate downstream with coords
  downstream <- calc_river_downstream(rivers, coords = coords)
  
  expect_true(is.numeric(downstream))
  expect_equal(length(downstream), nrow(rivers))
})

test_that("calc_river_downstream handles MULTILINESTRING", {
  skip_if_not_installed("sf")
  
  # Create MULTILINESTRING
  coords1 <- matrix(c(0, 0, 5, 5), ncol = 2, byrow = TRUE)
  coords2 <- matrix(c(5, 5, 10, 10), ncol = 2, byrow = TRUE)
  
  multi_line <- sf::st_multilinestring(list(coords1, coords2))
  rivers <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(multi_line, crs = 4326)
  )
  
  # Should convert and process
  downstream <- calc_river_downstream(rivers)
  
  expect_true(is.numeric(downstream))
  expect_equal(length(downstream), 2)
})

# Test calc_river_path function ----------------------------------------------

test_that("calc_river_path requires sf input", {
  expect_error(
    calc_river_path("invalid"),
    "must be an sf object"
  )
})

test_that("calc_river_path validates geometry type", {
  skip_if_not_installed("sf")
  
  # Create POINT geometry (invalid)
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )
  
  expect_error(
    calc_river_path(pts),
    "must contain LINESTRING"
  )
})


test_that("calc_river_path creates dissolved paths", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Calculate paths
  paths <- calc_river_path(rivers)
  
  expect_true(is.list(paths))
  expect_true(all(c("seg_ids", "point_ids", "paths") %in% names(paths)))
  
  # Check seg_ids
  expect_true(is.list(paths$seg_ids))
  expect_true(length(paths$seg_ids) > 0)
  
  # Check point_ids
  expect_true(is.list(paths$point_ids))
  expect_equal(length(paths$point_ids), length(paths$seg_ids))
  
  # Check paths sf object
  expect_s3_class(paths$paths, "sf")
  expect_true("path_id" %in% colnames(paths$paths))
  expect_true("n_segments" %in% colnames(paths$paths))
})

test_that("calc_river_path with provided downstream", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Calculate downstream first
  downstream <- calc_river_downstream(rivers)
  
  # Calculate paths with downstream
  paths <- calc_river_path(rivers, downstream = downstream)
  
  expect_true(is.list(paths))
  expect_s3_class(paths$paths, "sf")
})



# Test calc_river_properties function ----------------------------------------

test_that("calc_river_properties requires sf input", {
  expect_error(
    calc_river_properties("invalid"),
    "must be an sf object"
  )
})

test_that("calc_river_properties validates geometry type", {
  skip_if_not_installed("sf")
  
  # Create POINT geometry (invalid)
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )
  
  expect_error(
    calc_river_properties(pts),
    "must contain LINESTRING"
  )
})

test_that("calc_river_properties validates properties parameter", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Invalid property
  expect_error(
    calc_river_properties(rivers, properties = "invalid"),
    "Invalid properties"
  )
})

test_that("calc_river_properties requires DEM for slope", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Slope requires DEM
  expect_error(
    calc_river_properties(rivers, properties = "slope"),
    "dem.*required for slope"
  )
})

test_that("calc_river_properties validates DEM type", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Mock legacy raster
  legacy_raster <- structure(list(), class = "RasterLayer")
  
  expect_error(
    calc_river_properties(rivers, dem = legacy_raster, properties = "slope"),
    "must be a SpatRaster"
  )
})

test_that("calc_river_properties calculates order", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Calculate order only
  result <- calc_river_properties(rivers, properties = "order")
  
  expect_s3_class(result, "sf")
  expect_true("order" %in% colnames(result))
  expect_true(is.numeric(result$order))
  expect_equal(nrow(result), nrow(rivers))
})

test_that("calc_river_properties calculates downstream", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Calculate downstream only
  result <- calc_river_properties(rivers, properties = "downstream")
  
  expect_s3_class(result, "sf")
  expect_true("downstream" %in% colnames(result))
  expect_true(is.numeric(result$downstream))
})

test_that("calc_river_properties calculates length", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Calculate length only
  result <- calc_river_properties(rivers, properties = "length")
  
  expect_s3_class(result, "sf")
  expect_true("length" %in% colnames(result))
  expect_true(is.numeric(result$length))
  expect_true(all(result$length > 0))
})

test_that("calc_river_properties calculates slope", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  # Create DEM
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Calculate slope
  result <- calc_river_properties(rivers, dem = dem, properties = "slope")
  
  expect_s3_class(result, "sf")
  expect_true("slope" %in% colnames(result))
  expect_true(is.numeric(result$slope))
  expect_true(all(result$slope > 0))  # All slopes should be positive
})

test_that("calc_river_properties calculates all properties", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  # Create DEM
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Calculate all properties
  result <- calc_river_properties(rivers, dem = dem, properties = "all")
  
  expect_s3_class(result, "sf")
  expect_true(all(c("order", "downstream", "length", "slope") %in% colnames(result)))
})

test_that("calc_river_properties calculates multiple properties", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Calculate order and downstream
  result <- calc_river_properties(rivers, properties = c("order", "downstream"))
  
  expect_s3_class(result, "sf")
  expect_true("order" %in% colnames(result))
  expect_true("downstream" %in% colnames(result))
  expect_false("length" %in% colnames(result))
  expect_false("slope" %in% colnames(result))
})

test_that("calc_river_properties handles MULTILINESTRING", {
  skip_if_not_installed("sf")
  
  # Create MULTILINESTRING
  coords1 <- matrix(c(0, 0, 5, 5), ncol = 2, byrow = TRUE)
  coords2 <- matrix(c(5, 5, 10, 10), ncol = 2, byrow = TRUE)
  
  multi_line <- sf::st_multilinestring(list(coords1, coords2))
  rivers <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(multi_line, crs = 4326)
  )
  
  # Should convert and process
  result <- calc_river_properties(rivers, properties = "length")
  
  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 2)  # Should be split
})

# Test generate_river_types function -----------------------------------------

test_that("generate_river_types validates n parameter", {
  expect_error(
    generate_river_types(-1),
    "positive integer"
  )
  expect_error(
    generate_river_types(0),
    "positive integer"
  )
  expect_error(
    generate_river_types("invalid"),
    "positive integer"
  )
})

test_that("generate_river_types creates correct structure", {
  river_types <- generate_river_types(5)
  
  expect_true(is.data.frame(river_types))
  expect_equal(nrow(river_types), 5)
  
  # Check required columns
  required_cols <- c("Index", "Depth", "BankSlope", "Width", "Sinuosity",
                     "Manning", "Cwr", "KsatH", "BedThick")
  expect_true(all(required_cols %in% colnames(river_types)))
})

test_that("generate_river_types with custom widths", {
  widths <- c(2, 5, 10, 20, 40)
  river_types <- generate_river_types(5, width = widths)
  
  expect_equal(river_types$Width, widths)
})

test_that("generate_river_types with custom depths", {
  depths <- c(1, 2, 3, 4, 5)
  river_types <- generate_river_types(5, depth = depths)
  
  expect_equal(river_types$Depth, depths)
})

test_that("generate_river_types with custom manning", {
  manning <- c(0.03, 0.04, 0.05, 0.06, 0.07)
  river_types <- generate_river_types(5, manning = manning)
  
  expect_equal(river_types$Manning, manning)
})

test_that("generate_river_types with scalar parameters", {
  river_types <- generate_river_types(5, width = 10, depth = 5, manning = 0.05)
  
  expect_true(all(river_types$Width == 10))
  expect_true(all(river_types$Depth == 5))
  expect_true(all(river_types$Manning == 0.05))
})

test_that("generate_river_types validates parameter lengths", {
  expect_error(
    generate_river_types(5, width = c(1, 2)),
    "must have length 1 or n"
  )
  expect_error(
    generate_river_types(5, depth = c(1, 2, 3)),
    "must have length 1 or n"
  )
})

# Test calc_river_width_from_area function -----------------------------------

test_that("calc_river_width_from_area validates area parameter", {
  expect_error(
    calc_river_width_from_area(-1),
    "positive"
  )
  expect_error(
    calc_river_width_from_area(0),
    "positive"
  )
})

test_that("calc_river_width_from_area validates n_types parameter", {
  expect_error(
    calc_river_width_from_area(1000, n_types = -1),
    "positive integer"
  )
  expect_error(
    calc_river_width_from_area(1000, n_types = 0),
    "positive integer"
  )
})

test_that("calc_river_width_from_area returns correct length", {
  widths <- calc_river_width_from_area(1000, n_types = 5)
  
  expect_true(is.numeric(widths))
  expect_equal(length(widths), 5)
  expect_true(all(widths > 0))
})

test_that("calc_river_width_from_area scales with area", {
  widths_small <- calc_river_width_from_area(100, n_types = 5)
  widths_large <- calc_river_width_from_area(10000, n_types = 5)
  
  # Larger area should give larger widths
  expect_true(all(widths_large > widths_small))
})

test_that("calc_river_width_from_area increases with order", {
  widths <- calc_river_width_from_area(1000, n_types = 5)
  
  # Widths should increase with order (index)
  expect_true(all(diff(widths) > 0))
})

# Test build_river_network function ------------------------------------------

test_that("build_river_network requires sf input", {
  skip_if_not_installed("terra")
  
  dem <- terra::rast(ncol = 10, nrow = 10)
  
  expect_error(
    build_river_network("invalid", dem),
    "must be an sf object"
  )
})

test_that("build_river_network requires SpatRaster DEM", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Mock legacy raster
  legacy_raster <- structure(list(), class = "RasterLayer")
  
  expect_error(
    build_river_network(rivers, legacy_raster),
    "must be a SpatRaster"
  )
})

test_that("build_river_network validates geometry type", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  # Create POINT geometry (invalid)
  pts <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5),
    coords = c("x", "y"),
    crs = 4326
  )
  
  dem <- terra::rast(ncol = 10, nrow = 10)
  
  expect_error(
    build_river_network(pts, dem),
    "must contain LINESTRING"
  )
})

test_that("build_river_network creates complete network", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  # Create DEM
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Build network
  network <- build_river_network(rivers, dem)
  
  expect_true(is.list(network))
  expect_true(all(c("network", "river_types", "points") %in% names(network)))
  
  # Check network
  expect_s3_class(network$network, "sf")
  expect_true(all(c("Index", "Down", "Type", "Slope", "Length", "BC") %in% 
                  colnames(network$network)))
  
  # Check river_types
  expect_true(is.data.frame(network$river_types))
  expect_true(nrow(network$river_types) > 0)
  
  # Check points
  expect_true(is.data.frame(network$points))
  expect_equal(nrow(network$points), nrow(rivers))
  expect_true(all(c("From.x", "From.y", "From.z", "To.x", "To.y", "To.z") %in% 
                  colnames(network$points)))
})

test_that("build_river_network with area parameter", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  # Create DEM
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Build network with area
  network <- build_river_network(rivers, dem, area = 1000)
  
  expect_true(is.list(network))
  expect_true(is.data.frame(network$river_types))
  
  # Widths should be calculated from area
  expect_true(all(network$river_types$Width > 0))
})

test_that("build_river_network with provided river_order", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  # Create DEM
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Provide river order
  river_order <- c(2, 2, 1, 1)
  
  # Build network
  network <- build_river_network(rivers, dem, river_order = river_order)
  
  expect_equal(network$network$Type, river_order)
})

test_that("build_river_network with provided downstream", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  # Create DEM
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Provide downstream
  downstream <- c(-1, 1, 2, 2)
  
  # Build network
  network <- build_river_network(rivers, dem, downstream = downstream)
  
  expect_equal(network$network$Down, downstream)
})

test_that("build_river_network validates river_order length", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Wrong length
  river_order <- c(1, 2)
  
  expect_error(
    build_river_network(rivers, dem, river_order = river_order),
    "must have length equal to number of river segments"
  )
})

test_that("build_river_network validates downstream length", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Wrong length
  downstream <- c(-1, 1)
  
  expect_error(
    build_river_network(rivers, dem, downstream = downstream),
    "must have length equal to number of river segments"
  )
})

test_that("build_river_network handles MULTILINESTRING", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  # Create MULTILINESTRING
  coords1 <- matrix(c(0, 0, 5, 5), ncol = 2, byrow = TRUE)
  coords2 <- matrix(c(5, 5, 10, 10), ncol = 2, byrow = TRUE)
  
  multi_line <- sf::st_multilinestring(list(coords1, coords2))
  rivers <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(multi_line, crs = 4326)
  )
  
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Should convert and process
  network <- build_river_network(rivers, dem)
  
  expect_true(is.list(network))
  expect_equal(nrow(network$network), 2)  # Should be split
})

# Test get_river_outlets function --------------------------------------------

test_that("get_river_outlets from sf object", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Add downstream column
  rivers$Down <- c(-1, 1, 2, 2)
  
  outlets <- get_river_outlets(rivers)
  
  expect_true(is.numeric(outlets))
  expect_equal(outlets, 1)
})

test_that("get_river_outlets from numeric vector", {
  downstream <- c(-1, 1, 2, 2, -1)
  
  outlets <- get_river_outlets(downstream)
  
  expect_true(is.numeric(outlets))
  expect_equal(outlets, c(1, 5))
})

test_that("get_river_outlets requires Down or downstream column", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  expect_error(
    get_river_outlets(rivers),
    "must have 'Down' or 'downstream' column"
  )
})

test_that("get_river_outlets validates input type", {
  expect_error(
    get_river_outlets("invalid"),
    "must be an sf object or numeric vector"
  )
})

test_that("get_river_outlets with downstream column name", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  # Add downstream column (lowercase)
  rivers$downstream <- c(-1, 1, 2, 2)
  
  outlets <- get_river_outlets(rivers)
  
  expect_equal(outlets, 1)
})

# Test helper functions ------------------------------------------------------

test_that("get_coords extracts unique coordinates", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  
  coords <- get_coords(rivers)
  
  expect_true(is.matrix(coords))
  expect_equal(ncol(coords), 2)
  expect_true(nrow(coords) > 0)
  
  # Should have unique coordinates
  expect_equal(nrow(coords), nrow(unique(coords)))
})

test_that("get_from_to_nodes identifies nodes correctly", {
  skip_if_not_installed("sf")
  
  rivers <- create_test_river_network()
  coords <- get_coords(rivers)
  
  ft <- get_from_to_nodes(rivers, coords)
  
  expect_true(is.matrix(ft))
  expect_equal(ncol(ft), 3)
  expect_equal(nrow(ft), nrow(rivers))
  expect_true(all(c("ID", "FrNode", "ToNode") %in% colnames(ft)))
  
  # Node indices should be valid
  expect_true(all(ft[, "FrNode"] > 0))
  expect_true(all(ft[, "ToNode"] > 0))
  expect_true(all(ft[, "FrNode"] <= nrow(coords)))
  expect_true(all(ft[, "ToNode"] <= nrow(coords)))
})

test_that("sp.RiverOrder from/to points keep exact upstream-downstream match", {
  skip_if_not_installed("sp")
  skip_if_not_installed("raster")
  skip_if_not_installed("rgeos")
  
  line_a <- sp::Line(cbind(1:5, rep(0, 5)))
  line_b <- sp::Line(cbind(5:9, rep(0, 5)))
  reaches <- sp::SpatialLines(
    list(
      sp::Lines(list(line_a), ID = "1"),
      sp::Lines(list(line_b), ID = "2")
    ),
    proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
  )
  
  order <- sp.RiverOrder(reaches)
  coords_original <- get_coords(reaches)
  ft <- get_from_to_nodes(reaches, coords_original)
  
  from_xy <- coords_original[ft[, "FrNode"], , drop = FALSE]
  to_xy <- coords_original[ft[, "ToNode"], , drop = FALSE]
  
  expect_equal(length(order), 2)
  expect_equal(unname(from_xy), rbind(c(1, 0), c(5, 0)))
  expect_equal(unname(to_xy), rbind(c(5, 0), c(9, 0)))
  expect_identical(unname(to_xy[1, ]), unname(from_xy[2, ]))
})

test_that("get_from_to_nodes returns NA nodes for invalid line geometries", {
  skip_if_not_installed("sf")

  rivers <- sf::st_sf(
    id = 1:5,
    geometry = sf::st_sfc(
      sf::st_linestring(),
      sf::st_linestring(matrix(c(2, 2), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(3, 3, 3, 3), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(0, 0, Inf, 1), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
      crs = 4326
    )
  )

  coords <- get_coords(rivers)
  ft <- get_from_to_nodes(rivers, coords)

  expect_equal(nrow(ft), 5)
  expect_true(all(is.na(ft[1:4, c("FrNode", "ToNode")])))
  expect_false(anyNA(ft[5, c("FrNode", "ToNode")]))
})

test_that("rmDuplicatedLines removes duplicate reaches and invalid lines", {
  skip_if_not_installed("sf")

  rivers <- sf::st_sf(
    id = 1:6,
    geometry = sf::st_sfc(
      sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)),
      sf::st_linestring(),
      sf::st_linestring(matrix(c(2, 2), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(3, 3, 3, 3), ncol = 2, byrow = TRUE)),
      sf::st_linestring(matrix(c(4, 4, 5, 5), ncol = 2, byrow = TRUE)),
      crs = 4326
    )
  )

  out <- rmDuplicatedLines(rivers)

  expect_s3_class(out, "sf")
  expect_equal(nrow(out), 2)
  expect_equal(out$id, c(1, 6))
  ft <- get_from_to_nodes(out)
  expect_false(anyNA(ft[, c("FrNode", "ToNode")]))
})

# Test S4 class conversion functions -----------------------------------------

test_that("as_shud_river validates input", {
  expect_error(
    as_shud_river("invalid"),
    "must be a list"
  )
  
  expect_error(
    as_shud_river(list(network = NULL)),
    "Missing required components"
  )
})

test_that("as_shud_river creates SHUD.RIVER object", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  # Create DEM
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Build network
  network_list <- build_river_network(rivers, dem)
  
  # Convert to SHUD.RIVER
  shud_river <- as_shud_river(network_list)
  
  expect_s4_class(shud_river, "SHUD River")
  expect_true(is.data.frame(shud_river@river))
  expect_true(is.data.frame(shud_river@rivertype))
  expect_true(is.data.frame(shud_river@point))
})

test_that("shud_river_to_sf validates input", {
  expect_error(
    shud_river_to_sf("invalid"),
    "must be a SHUD.RIVER object"
  )
})

test_that("shud_river_to_sf extracts sf from modern format", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  # Create DEM
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Build network and convert
  network_list <- build_river_network(rivers, dem)
  shud_river <- as_shud_river(network_list)
  
  # Extract sf
  rivers_sf <- shud_river_to_sf(shud_river)
  
  expect_s3_class(rivers_sf, "sf")
  expect_true(all(c("Index", "Down", "Type", "Slope", "Length", "BC") %in% 
                  colnames(rivers_sf)))
})

test_that("is_modern_river identifies modern format", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  # Create DEM
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Build network and convert
  network_list <- build_river_network(rivers, dem)
  shud_river <- as_shud_river(network_list)
  
  # Check if modern
  expect_true(is_modern_river(shud_river))
})

test_that("is_modern_river validates input", {
  expect_error(
    is_modern_river("invalid"),
    "must be a SHUD.RIVER object"
  )
})

# Integration tests with sample data -----------------------------------------

test_that("river module integration test with sample data", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_on_cran()
  
  # This would test with actual rSHUD sample data if available
  # skip_if_not(exists("sh"), message = "Sample data 'sh' not available")
  
  # For now, use synthetic data
  rivers <- create_test_river_network()
  
  # Create DEM
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Full workflow
  # 1. Calculate properties
  rivers_props <- calc_river_properties(rivers, dem, properties = "all")
  
  expect_s3_class(rivers_props, "sf")
  expect_true(all(c("order", "downstream", "length", "slope") %in% 
                  colnames(rivers_props)))
  
  # 2. Build complete network
  network <- build_river_network(rivers, dem, area = 1000)
  
  expect_true(is.list(network))
  expect_s3_class(network$network, "sf")
  
  # 3. Convert to SHUD.RIVER
  shud_river <- as_shud_river(network)
  
  expect_s4_class(shud_river, "SHUD River")
  
  # 4. Convert back to sf
  rivers_sf <- shud_river_to_sf(shud_river)
  
  expect_s3_class(rivers_sf, "sf")
  
  # 5. Calculate paths
  paths <- calc_river_path(rivers)
  
  expect_true(is.list(paths))
  expect_s3_class(paths$paths, "sf")
  
  # 6. Get outlets
  outlets <- get_river_outlets(network$network)
  
  expect_true(is.numeric(outlets))
  expect_true(length(outlets) > 0)
})

test_that("river module handles complex network topology", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  # Create more complex network with multiple branches
  coords_list <- list(
    # Main stem outlet
    matrix(c(5, 3, 5, 0), ncol = 2, byrow = TRUE),
    # Main stem middle
    matrix(c(5, 6, 5, 3), ncol = 2, byrow = TRUE),
    # Main stem upper
    matrix(c(5, 10, 5, 6), ncol = 2, byrow = TRUE),
    # Left tributary lower
    matrix(c(0, 3, 5, 3), ncol = 2, byrow = TRUE),
    # Left tributary upper
    matrix(c(0, 6, 5, 6), ncol = 2, byrow = TRUE),
    # Right tributary lower
    matrix(c(10, 3, 5, 3), ncol = 2, byrow = TRUE),
    # Right tributary upper
    matrix(c(10, 6, 5, 6), ncol = 2, byrow = TRUE),
    # Headwater 1
    matrix(c(0, 10, 5, 10), ncol = 2, byrow = TRUE),
    # Headwater 2
    matrix(c(10, 10, 5, 10), ncol = 2, byrow = TRUE)
  )
  
  geoms <- lapply(coords_list, sf::st_linestring)
  rivers <- sf::st_sf(
    id = 1:9,
    geometry = sf::st_sfc(geoms, crs = 4326)
  )
  
  # Create DEM
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Calculate order
  order <- calc_river_order(rivers)
  
  expect_equal(length(order), 9)
  expect_true(max(order) >= 2)  # Should have at least order 2
  
  # Calculate downstream
  downstream <- calc_river_downstream(rivers)
  
  expect_equal(length(downstream), 9)
  expect_equal(sum(downstream < 0), 1)  # Should have exactly 1 outlet
  
  # Build network
  network <- build_river_network(rivers, dem)
  
  expect_equal(nrow(network$network), 9)
  expect_true(max(network$network$Type) >= 2)
})

test_that("river module preserves CRS through workflow", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  
  rivers <- create_test_river_network()
  
  # Create DEM with specific CRS
  dem <- terra::rast(ncol = 20, nrow = 20, xmin = -1, xmax = 11, ymin = -1, ymax = 11)
  terra::values(dem) <- seq(100, 200, length.out = 400)
  terra::crs(dem) <- "EPSG:4326"
  
  # Calculate properties
  rivers_props <- calc_river_properties(rivers, dem, properties = "all")
  
  expect_equal(sf::st_crs(rivers_props), sf::st_crs(rivers))
  
  # Build network
  network <- build_river_network(rivers, dem)
  
  expect_equal(sf::st_crs(network$network), sf::st_crs(rivers))
  
  # Calculate paths
  paths <- calc_river_path(rivers)
  
  expect_equal(sf::st_crs(paths$paths), sf::st_crs(rivers))
})
