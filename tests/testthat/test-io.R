# I/O Module Tests
# 
# This file tests all I/O functions including:
# - Modern snake_case functions
# - Backward compatibility with legacy functions
# - Deprecation warnings
# - Format conversion
# - Error handling and validation

# Load required packages for testing
library(testthat)

# =============================================================================
# SHUD File I/O Tests
# =============================================================================

test_that("new read functions have proper parameter validation", {
  # Test missing file parameter
  expect_error(
    read_mesh(NULL),
    "Parameter 'file' is required"
  )
  
  expect_error(
    read_river(NA),
    "Parameter 'file' is required"
  )
  
  expect_error(
    read_att(""),
    "Parameter 'file' is required"
  )
  
  # Test non-existent file
  expect_error(
    read_mesh("nonexistent.mesh"),
    "Mesh file does not exist"
  )
  
  expect_error(
    read_river("nonexistent.riv"),
    "River file does not exist"
  )
  
  expect_error(
    read_att("nonexistent.att"),
    "Attributes file does not exist"
  )
})

test_that("read_df function validates parameters correctly", {
  # Test missing parameters
  expect_error(
    read_df(),
    "Either 'file' or 'text' parameter is required"
  )
  
  # Test non-existent file
  expect_error(
    read_df("nonexistent.txt"),
    "File does not exist"
  )
})

test_that("read_tsd function validates parameters correctly", {
  # Test missing file parameter
  expect_error(
    read_tsd(),
    "Parameter 'file' is required"
  )
  
  expect_error(
    read_tsd(NULL),
    "Parameter 'file' is required"
  )
  
  # Test non-existent file
  expect_error(
    read_tsd("nonexistent.tsd"),
    "Time-series file does not exist"
  )
})

test_that("all new read functions exist and are exported", {
  # Test that all new functions exist
  expect_true(exists("read_mesh"))
  expect_true(exists("read_river"))
  expect_true(exists("read_att"))
  expect_true(exists("read_rivseg"))
  expect_true(exists("read_para"))
  expect_true(exists("read_calib"))
  expect_true(exists("read_config"))
  expect_true(exists("read_ic"))
  expect_true(exists("read_soil"))
  expect_true(exists("read_geol"))
  expect_true(exists("read_lc"))
  expect_true(exists("read_forc_fn"))
  expect_true(exists("read_forc_csv"))
  expect_true(exists("read_df"))
  expect_true(exists("read_tsd"))
  expect_true(exists("read_lai"))
  expect_true(exists("read_river_sp"))
})

# =============================================================================
# Write Function Tests
# =============================================================================

test_that("new write functions have proper parameter validation", {
  temp_dir <- tempdir()
  
  # Test missing data parameter
  expect_error(
    write_mesh(file = file.path(temp_dir, "test.mesh")),
    "Parameter 'pm' \\(SHUD\\.MESH object\\) is required"
  )
  
  expect_error(
    write_river(file = file.path(temp_dir, "test.riv")),
    "Parameter 'pr' \\(SHUD\\.RIVER object\\) is required"
  )
  
  expect_error(
    write_ic(file = file.path(temp_dir, "test.ic")),
    "Parameter 'x' \\(initial condition data\\) is required"
  )
  
  expect_error(
    write_config(file = file.path(temp_dir, "test.para")),
    "Parameter 'x' \\(configuration data\\) is required"
  )
  
  expect_error(
    write_forc(file = file.path(temp_dir, "test.forc")),
    "Parameter 'x' \\(forcing sites data\\) is required"
  )
  
  expect_error(
    write_df(file = file.path(temp_dir, "test.txt")),
    "Parameter 'x' \\(data to write\\) is required"
  )
  
  # Test missing file parameter
  test_mesh <- new("Untructure Domain", 
                   mesh = data.frame(node1 = 1:3, node2 = 2:4, node3 = 3:5),
                   point = data.frame(x = c(0, 1, 2), y = c(0, 1, 2)))
  
  expect_error(
    write_mesh(test_mesh),
    "Parameter 'file' is required"
  )
  
  expect_error(
    write_mesh(test_mesh, NULL),
    "Parameter 'file' is required"
  )
  
  expect_error(
    write_mesh(test_mesh, ""),
    "Parameter 'file' is required"
  )
  
  # Test wrong data type
  expect_error(
    write_mesh("not_a_mesh", file.path(temp_dir, "test.mesh")),
    "Parameter 'pm' must be a SHUD\\.MESH object"
  )
  
  expect_error(
    write_river("not_a_river", file.path(temp_dir, "test.riv")),
    "Parameter 'pr' must be a SHUD\\.RIVER object"
  )
  
  expect_error(
    write_ic("not_a_list", file.path(temp_dir, "test.ic")),
    "Parameter 'x' must be a list with at least 2 elements"
  )
  
  expect_error(
    write_forc("not_a_dataframe", file.path(temp_dir, "test.forc")),
    "Parameter 'x' must be a data frame"
  )
})

test_that("write_tsd and write_xts validate parameters correctly", {
  temp_dir <- tempdir()
  
  # Test missing data parameter
  expect_error(
    write_tsd(file = file.path(temp_dir, "test.tsd")),
    "Parameter 'x' \\(time series data\\) is required"
  )
  
  expect_error(
    write_xts(file = file.path(temp_dir, "test.txt")),
    "Parameter 'x' \\(time series data\\) is required"
  )
  
  # Test missing file parameter
  expect_error(
    write_tsd(data.frame(a = 1:3)),
    "Parameter 'file' is required"
  )
  
  expect_error(
    write_xts(data.frame(a = 1:3)),
    "Parameter 'file' is required"
  )
  
  # Test wrong data type
  expect_error(
    write_tsd("not_xts", file.path(temp_dir, "test.tsd")),
    "Parameter 'x' must be an xts object"
  )
  
  expect_error(
    write_xts("not_xts", file.path(temp_dir, "test.txt")),
    "Parameter 'x' must be an xts object"
  )
})

test_that("all new write functions exist and are exported", {
  # Test that all new write functions exist
  expect_true(exists("write_mesh"))
  expect_true(exists("write_river"))
  expect_true(exists("write_ic"))
  expect_true(exists("write_config"))
  expect_true(exists("write_forc"))
  expect_true(exists("write_df"))
  expect_true(exists("write_tsd"))
  expect_true(exists("write_xts"))
})

# =============================================================================
# Time Series I/O Tests
# =============================================================================

test_that("time series functions validate parameters correctly", {
  temp_dir <- tempdir()
  
  # Test read_tsd parameter validation
  expect_error(
    read_tsd(),
    "Parameter 'file' is required"
  )
  
  expect_error(
    read_tsd("nonexistent.tsd"),
    "Time-series file does not exist"
  )
  
  # Test write_tsd parameter validation
  expect_error(
    write_tsd(file = file.path(temp_dir, "test.tsd")),
    "Parameter 'x' \\(time series data\\) is required"
  )
  
  expect_error(
    write_tsd("not_xts", file.path(temp_dir, "test.tsd")),
    "Parameter 'x' must be an xts object"
  )
})

test_that("time series utility functions work correctly", {
  skip_if_not_installed("xts")
  
  # Create test time series
  dates <- seq(as.Date("2020-01-01"), as.Date("2020-01-10"), by = "day")
  ts_data <- xts::xts(runif(length(dates)), order.by = dates)
  
  # Test ts_to_df
  df_result <- ts_to_df(ts_data)
  expect_s3_class(df_result, "data.frame")
  expect_equal(nrow(df_result), length(dates))
  expect_true("Time" %in% names(df_result))
  
  # Test ts_to_daily (should be identity for daily data)
  daily_result <- ts_to_daily(ts_data)
  expect_s3_class(daily_result, "xts")
  expect_equal(nrow(daily_result), nrow(ts_data))
})

test_that("time series functions exist and are exported", {
  expect_true(exists("read_tsd"))
  expect_true(exists("write_tsd"))
  expect_true(exists("read_lai"))
  expect_true(exists("write_xts"))
  expect_true(exists("ts_to_daily"))
  expect_true(exists("ts_to_df"))
})

# =============================================================================
# Output I/O Tests
# =============================================================================

test_that("output reading functions validate parameters correctly", {
  # Test read_output parameter validation
  expect_error(
    read_output(),
    "Either 'keyword' or 'file' must be provided"
  )
  
  expect_error(
    read_output(file = "nonexistent.dat"),
    "Output file does not exist"
  )
  
  # Test read_output_header parameter validation
  expect_error(
    read_output_header(),
    "Either 'keyword' or 'file' must be provided"
  )
})

test_that("output utility functions work correctly", {
  # Test get_default_variables
  default_vars <- get_default_variables()
  expect_type(default_vars, "character")
  expect_true(length(default_vars) > 0)
  expect_true(any(grepl("^eley", default_vars)))
  expect_true(any(grepl("^rivq", default_vars)))
})

test_that("output functions exist and are exported", {
  expect_true(exists("read_output"))
  expect_true(exists("read_output_header"))
  expect_true(exists("load_output_data"))
  expect_true(exists("get_default_variables"))
})

# =============================================================================
# NetCDF I/O Tests
# =============================================================================

test_that("netcdf functions validate parameters correctly", {
  # Test read_nc_time parameter validation
  expect_error(
    read_nc_time(),
    "Parameter 'ncid' \\(NetCDF connection\\) is required"
  )
  
  # Test read_nc_data parameter validation
  expect_error(
    read_nc_data(),
    "Parameter 'ncid' \\(NetCDF connection\\) is required"
  )
  
  # Test open_nc_file parameter validation
  expect_error(
    open_nc_file(),
    "Parameter 'file' is required"
  )
  
  expect_error(
    open_nc_file("nonexistent.nc"),
    "NetCDF file does not exist"
  )
})

test_that("nc_to_raster validates parameters correctly", {
  # Test parameter validation
  expect_error(
    nc_to_raster(),
    "Either 'nc_data' or 'x', 'y', 'data_array' must be provided"
  )
  
  expect_error(
    nc_to_raster(x = 1:5),
    "Either 'nc_data' or 'x', 'y', 'data_array' must be provided"
  )
})

test_that("netcdf functions exist and are exported", {
  expect_true(exists("read_nc_time"))
  expect_true(exists("read_nc_data"))
  expect_true(exists("nc_to_raster"))
  expect_true(exists("open_nc_file"))
  expect_true(exists("get_nc_info"))
})

# =============================================================================
# Backward Compatibility Tests
# =============================================================================

test_that("legacy function names still exist and show deprecation warnings", {
  # Test that legacy functions exist (they should be defined in the original files)
  legacy_functions <- c(
    "readmesh", "readriv", "readatt", "readpara", "readcalib", 
    "readconfig", "readic", "readsoil", "readgeol", "readlc",
    "readforc.fn", "readforc.csv", "readlai", "readout", "readout.header"
  )
  
  for (func_name in legacy_functions) {
    expect_true(exists(func_name), 
                info = paste("Legacy function", func_name, "should exist"))
  }
})

test_that("deprecated GIS functions show warnings", {
  skip_if_not_installed("terra")
  
  # Test deprecated sp2raster function
  expect_warning(
    tryCatch(
      sp2raster(),
      error = function(e) NULL
    ),
    "sp2raster.*deprecated.*vector_to_raster|sp2raster.*不再有用|vector_to_raster"
  )
  
  # Test deprecated MeshData2Raster function  
  expect_warning(
    tryCatch(
      MeshData2Raster(),
      error = function(e) NULL
    ),
    "MeshData2Raster.*deprecated.*mesh_to_raster|MeshData2Raster.*不再有用|mesh_to_raster"
  )
})

test_that("modern functions handle legacy input gracefully", {
  # Test that modern functions provide helpful error messages for legacy inputs
  temp_dir <- tempdir()
  
  # Test with invalid legacy-style input
  expect_error(
    read_mesh(NULL),
    "Parameter 'file' is required"
  )
  
  expect_error(
    write_mesh("not_a_mesh_object", file.path(temp_dir, "test.mesh")),
    "Parameter 'pm' must be a SHUD\\.MESH object"
  )
})

# =============================================================================
# Format Conversion Tests
# =============================================================================

test_that("format conversion works correctly", {
  skip_if_not_installed("xts")
  
  # Create test data
  dates <- seq(as.Date("2020-01-01"), as.Date("2020-01-03"), by = "day")
  ts_data <- xts::xts(runif(length(dates)), order.by = dates)
  
  # Test conversion to data frame
  df_data <- ts_to_df(ts_data, time_col = "Date")
  expect_s3_class(df_data, "data.frame")
  expect_equal(nrow(df_data), length(dates))
  expect_true("Date" %in% names(df_data))
  
  # Test daily aggregation (should be identity for daily data)
  daily_data <- ts_to_daily(ts_data, FUN = mean)
  expect_s3_class(daily_data, "xts")
  expect_equal(nrow(daily_data), nrow(ts_data))
})

test_that("time series format conversion handles different units", {
  skip_if_not_installed("xts")
  
  temp_dir <- tempdir()
  
  # Create hourly test data
  hourly_dates <- seq(as.POSIXct("2020-01-01 00:00:00"), 
                     as.POSIXct("2020-01-01 23:00:00"), by = "hour")
  hourly_data <- xts::xts(runif(length(hourly_dates)), order.by = hourly_dates)
  
  # Test writing with different time units
  test_file <- file.path(temp_dir, "test_hourly.tsd")
  
  expect_silent(write_tsd(hourly_data, test_file, dt_units = "minutes", quiet = TRUE))
  expect_true(file.exists(test_file))
  
  # Test conversion to daily
  daily_converted <- ts_to_daily(hourly_data, FUN = mean)
  expect_s3_class(daily_converted, "xts")
  expect_true(nrow(daily_converted) <= nrow(hourly_data))
})

test_that("NetCDF format conversion works", {
  skip_if_not_installed("terra")
  
  # Create mock NetCDF-like data
  x_coords <- seq(-180, 180, by = 1)
  y_coords <- seq(-90, 90, by = 1)
  data_array <- array(runif(length(x_coords) * length(y_coords)), 
                     dim = c(length(x_coords), length(y_coords)))
  
  # Test conversion to raster format
  raster_result <- nc_to_raster(x = x_coords, y = y_coords, 
                               data_array = data_array, format = "modern")
  
  expect_s4_class(raster_result, "SpatRaster")
  expect_equal(terra::ncell(raster_result), length(x_coords) * length(y_coords))
})

# =============================================================================
# Integration Tests
# =============================================================================

test_that("I/O modules work together", {
  skip_if_not_installed("xts")
  
  temp_dir <- tempdir()
  
  # Create test time series data
  dates <- seq(as.Date("2020-01-01"), as.Date("2020-01-05"), by = "day")
  ts_data <- xts::xts(matrix(runif(length(dates) * 3), ncol = 3), 
                     order.by = dates)
  colnames(ts_data) <- c("var1", "var2", "var3")
  
  # Test write and read cycle
  test_file <- file.path(temp_dir, "test_integration.tsd")
  
  # Write data
  expect_silent(write_tsd(ts_data, test_file, quiet = TRUE))
  expect_true(file.exists(test_file))
  
  # Read data back
  read_data <- read_tsd(test_file)
  expect_type(read_data, "list")
  expect_true(length(read_data) >= 1)
  expect_s3_class(read_data[[1]], "xts")
})

test_that("complete I/O workflow with different formats", {
  skip_if_not_installed("xts")
  
  temp_dir <- tempdir()
  
  # Create comprehensive test data
  dates <- seq(as.Date("2020-01-01"), as.Date("2020-01-10"), by = "day")
  
  # Multi-variable time series
  ts_data <- xts::xts(
    data.frame(
      temp = runif(length(dates), 10, 30),
      precip = runif(length(dates), 0, 10),
      flow = runif(length(dates), 1, 100)
    ),
    order.by = dates
  )
  
  # Test TSD format
  tsd_file <- file.path(temp_dir, "test_complete.tsd")
  expect_silent(write_tsd(ts_data, tsd_file, quiet = TRUE))
  expect_true(file.exists(tsd_file))
  
  # Test XTS format
  xts_file <- file.path(temp_dir, "test_complete.txt")
  expect_silent(write_xts(ts_data, xts_file))
  expect_true(file.exists(xts_file))
  
  # Read back and verify
  tsd_read <- read_tsd(tsd_file)
  expect_type(tsd_read, "list")
  expect_true(length(tsd_read) >= 1)
  
  # Test format conversion
  df_converted <- ts_to_df(ts_data)
  expect_s3_class(df_converted, "data.frame")
  expect_equal(nrow(df_converted), length(dates))
  expect_equal(ncol(df_converted), ncol(ts_data) + 1)  # +1 for time column
})

# =============================================================================
# Error Handling and Edge Cases
# =============================================================================

test_that("I/O functions handle edge cases gracefully", {
  temp_dir <- tempdir()
  
  # Test empty file paths
  expect_error(read_tsd(""), "Time-series file does not exist")
  expect_error(write_tsd(NULL, ""), "Parameter 'file' is required")
  
  # Test non-existent directories for write operations
  deep_path <- file.path(temp_dir, "deep", "nested", "path", "test.tsd")
  
  skip_if_not_installed("xts")
  dates <- seq(as.Date("2020-01-01"), as.Date("2020-01-03"), by = "day")
  ts_data <- xts::xts(runif(length(dates)), order.by = dates)
  
  # Should create directories automatically
  expect_silent(write_tsd(ts_data, deep_path, quiet = TRUE))
  expect_true(file.exists(deep_path))
})

test_that("NetCDF functions handle missing dependencies gracefully", {
  # Test parameter validation without requiring ncdf4
  expect_error(
    read_nc_time(),
    "Parameter 'ncid' \\(NetCDF connection\\) is required"
  )
  
  expect_error(
    open_nc_file(""),
    "Parameter 'file' is required"
  )
  
  expect_error(
    get_nc_info(),
    "Parameter 'ncid' \\(NetCDF connection\\) is required"
  )
})

test_that("output functions handle missing environment gracefully", {
  # Test functions that depend on .shud environment
  expect_type(get_default_variables(), "character")
  expect_true(length(get_default_variables()) > 0)
  
  # These are internal helpers and should still be callable without .shud.
  expect_type(rSHUD:::get_output_path(), "character")
  expect_type(rSHUD:::get_model_version(), "double")
  expect_type(rSHUD:::get_project_name(), "character")
})