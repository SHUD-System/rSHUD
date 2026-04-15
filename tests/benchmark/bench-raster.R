# Performance Benchmarks for rSHUD Modernization
#
# This script compares the performance of new terra/sf implementations
# against legacy raster/sp implementations to verify >= 20% performance
# improvement as required by Requirement 12.
#
# Run this script with: source("tests/benchmark/bench-raster.R")

# Check for required packages
if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop(
    "Package 'microbenchmark' is required for benchmarks.\n",
    "Install it with: install.packages('microbenchmark')",
    call. = FALSE
  )
}

if (!requireNamespace("terra", quietly = TRUE)) {
  stop(
    "Package 'terra' is required for benchmarks.\n",
    "Install it with: install.packages('terra')",
    call. = FALSE
  )
}

if (!requireNamespace("sf", quietly = TRUE)) {
  stop(
    "Package 'sf' is required for benchmarks.\n",
    "Install it with: install.packages('sf')",
    call. = FALSE
  )
}

# Load required packages
library(microbenchmark)
library(terra)
library(sf)

# Try to load legacy packages for comparison (optional)
has_legacy <- requireNamespace("raster", quietly = TRUE) &&
              requireNamespace("sp", quietly = TRUE) &&
              requireNamespace("rgeos", quietly = TRUE)

if (!has_legacy) {
  message("Legacy packages (raster/sp) not available.")
  message("Benchmarks will only test modern implementations.")
}

# Helper function to format benchmark results
format_benchmark <- function(mb_result) {
  summary_df <- summary(mb_result)
  summary_df$mean_ms <- summary_df$mean / 1e6  # Convert to milliseconds
  summary_df$median_ms <- summary_df$median / 1e6
  return(summary_df)
}

# Helper function to calculate speedup
calculate_speedup <- function(legacy_time, modern_time) {
  speedup <- (legacy_time - modern_time) / legacy_time * 100
  return(speedup)
}

cat("\n=== rSHUD Performance Benchmarks ===\n")
cat("Testing terra/sf vs raster/sp implementations\n")
cat("Target: >= 20% performance improvement\n\n")

# Store all results
all_results <- list()

# ==============================================================================
# Benchmark 1: Raster Creation and Basic Operations
# ==============================================================================

cat("Benchmark 1: Raster Creation (1000x1000 grid)\n")
cat(strrep("-", 60), "\n")

n_rows <- 1000
n_cols <- 1000
values <- matrix(rnorm(n_rows * n_cols), nrow = n_rows, ncol = n_cols)

if (has_legacy) {
  bench1 <- microbenchmark(
    legacy_raster = {
      r <- raster::raster(nrows = n_rows, ncols = n_cols,
                         xmn = 0, xmx = 100, ymn = 0, ymx = 100)
      raster::values(r) <- values
    },
    modern_terra = {
      r <- terra::rast(nrows = n_rows, ncols = n_cols,
                      xmin = 0, xmax = 100, ymin = 0, ymax = 100)
      terra::values(r) <- values
    },
    times = 50
  )
  
  result1 <- format_benchmark(bench1)
  print(result1[, c("expr", "mean_ms", "median_ms")])
  
  speedup1 <- calculate_speedup(
    result1$mean_ms[result1$expr == "legacy_raster"],
    result1$mean_ms[result1$expr == "modern_terra"]
  )
  
  cat(sprintf("\nSpeedup: %.1f%%", speedup1))
  if (speedup1 >= 20) {
    cat(" ✓ PASS (>= 20%)\n")
  } else {
    cat(" ✗ FAIL (< 20%)\n")
  }
  
  all_results$raster_creation <- list(
    speedup = speedup1,
    pass = speedup1 >= 20,
    details = result1
  )
} else {
  # Just test modern implementation
  bench1 <- microbenchmark(
    modern_terra = {
      r <- terra::rast(nrows = n_rows, ncols = n_cols,
                      xmin = 0, xmax = 100, ymin = 0, ymax = 100)
      terra::values(r) <- values
    },
    times = 50
  )
  result1 <- format_benchmark(bench1)
  print(result1[, c("expr", "mean_ms", "median_ms")])
  cat("\nLegacy comparison not available\n")
}

# ==============================================================================
# Benchmark 2: Raster Cropping
# ==============================================================================

cat("\n\nBenchmark 2: Raster Cropping\n")
cat(strrep("-", 60), "\n")

# Create test raster
r_terra <- terra::rast(nrows = 1000, ncols = 1000,
                       xmin = 0, xmax = 100, ymin = 0, ymax = 100)
terra::values(r_terra) <- rnorm(1000 * 1000)

# Create crop extent
crop_ext <- c(25, 75, 25, 75)

if (has_legacy) {
  r_raster <- raster::raster(r_terra)
  crop_ext_raster <- raster::extent(crop_ext)
  crop_ext_terra <- terra::ext(crop_ext)
  
  bench2 <- microbenchmark(
    legacy_raster = {
      r_crop <- raster::crop(r_raster, crop_ext_raster)
    },
    modern_terra = {
      r_crop <- terra::crop(r_terra, crop_ext_terra)
    },
    times = 100
  )
  
  result2 <- format_benchmark(bench2)
  print(result2[, c("expr", "mean_ms", "median_ms")])
  
  speedup2 <- calculate_speedup(
    result2$mean_ms[result2$expr == "legacy_raster"],
    result2$mean_ms[result2$expr == "modern_terra"]
  )
  
  cat(sprintf("\nSpeedup: %.1f%%", speedup2))
  if (speedup2 >= 20) {
    cat(" ✓ PASS (>= 20%)\n")
  } else {
    cat(" ✗ FAIL (< 20%)\n")
  }
  
  all_results$raster_crop <- list(
    speedup = speedup2,
    pass = speedup2 >= 20,
    details = result2
  )
} else {
  crop_ext_terra <- terra::ext(crop_ext)
  bench2 <- microbenchmark(
    modern_terra = {
      r_crop <- terra::crop(r_terra, crop_ext_terra)
    },
    times = 100
  )
  result2 <- format_benchmark(bench2)
  print(result2[, c("expr", "mean_ms", "median_ms")])
  cat("\nLegacy comparison not available\n")
}

# ==============================================================================
# Benchmark 3: Raster Resampling
# ==============================================================================

cat("\n\nBenchmark 3: Raster Resampling (bilinear)\n")
cat(strrep("-", 60), "\n")

# Create source and target rasters
r_source_terra <- terra::rast(nrows = 500, ncols = 500,
                               xmin = 0, xmax = 100, ymin = 0, ymax = 100)
terra::values(r_source_terra) <- rnorm(500 * 500)

r_target_terra <- terra::rast(nrows = 200, ncols = 200,
                               xmin = 0, xmax = 100, ymin = 0, ymax = 100)

if (has_legacy) {
  r_source_raster <- raster::raster(r_source_terra)
  r_target_raster <- raster::raster(r_target_terra)
  
  bench3 <- microbenchmark(
    legacy_raster = {
      r_resamp <- raster::resample(r_source_raster, r_target_raster,
                                   method = "bilinear")
    },
    modern_terra = {
      r_resamp <- terra::resample(r_source_terra, r_target_terra,
                                  method = "bilinear")
    },
    times = 30
  )
  
  result3 <- format_benchmark(bench3)
  print(result3[, c("expr", "mean_ms", "median_ms")])
  
  speedup3 <- calculate_speedup(
    result3$mean_ms[result3$expr == "legacy_raster"],
    result3$mean_ms[result3$expr == "modern_terra"]
  )
  
  cat(sprintf("\nSpeedup: %.1f%%", speedup3))
  if (speedup3 >= 20) {
    cat(" ✓ PASS (>= 20%)\n")
  } else {
    cat(" ✗ FAIL (< 20%)\n")
  }
  
  all_results$raster_resample <- list(
    speedup = speedup3,
    pass = speedup3 >= 20,
    details = result3
  )
} else {
  bench3 <- microbenchmark(
    modern_terra = {
      r_resamp <- terra::resample(r_source_terra, r_target_terra,
                                  method = "bilinear")
    },
    times = 30
  )
  result3 <- format_benchmark(bench3)
  print(result3[, c("expr", "mean_ms", "median_ms")])
  cat("\nLegacy comparison not available\n")
}

# ==============================================================================
# Benchmark 4: Vector Operations - Buffer
# ==============================================================================

cat("\n\nBenchmark 4: Vector Buffer Operation (1000 points)\n")
cat(strrep("-", 60), "\n")

# Create test points
n_pts <- 1000
coords <- data.frame(
  x = runif(n_pts, 0, 100),
  y = runif(n_pts, 0, 100),
  id = 1:n_pts
)

pts_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 4326)

if (has_legacy) {
  pts_sp <- as(pts_sf, "Spatial")
  
  bench4 <- microbenchmark(
    legacy_sp = {
      buf <- rgeos::gBuffer(pts_sp, width = 1, byid = TRUE)
    },
    modern_sf = {
      buf <- sf::st_buffer(pts_sf, dist = 1)
    },
    times = 50
  )
  
  result4 <- format_benchmark(bench4)
  print(result4[, c("expr", "mean_ms", "median_ms")])
  
  speedup4 <- calculate_speedup(
    result4$mean_ms[result4$expr == "legacy_sp"],
    result4$mean_ms[result4$expr == "modern_sf"]
  )
  
  cat(sprintf("\nSpeedup: %.1f%%", speedup4))
  if (speedup4 >= 20) {
    cat(" ✓ PASS (>= 20%)\n")
  } else {
    cat(" ✗ FAIL (< 20%)\n")
  }
  
  all_results$vector_buffer <- list(
    speedup = speedup4,
    pass = speedup4 >= 20,
    details = result4
  )
} else {
  bench4 <- microbenchmark(
    modern_sf = {
      buf <- sf::st_buffer(pts_sf, dist = 1)
    },
    times = 50
  )
  result4 <- format_benchmark(bench4)
  print(result4[, c("expr", "mean_ms", "median_ms")])
  cat("\nLegacy comparison not available\n")
}

# ==============================================================================
# Benchmark 5: Vector Operations - Intersection
# ==============================================================================

cat("\n\nBenchmark 5: Vector Intersection (100 polygons)\n")
cat(strrep("-", 60), "\n")

# Create test polygons
n_poly <- 100
polys_list <- lapply(1:n_poly, function(i) {
  x <- runif(1, 0, 90)
  y <- runif(1, 0, 90)
  sf::st_polygon(list(cbind(
    c(x, x + 10, x + 10, x, x),
    c(y, y, y + 10, y + 10, y)
  )))
})

polys_sf <- sf::st_sf(
  id = 1:n_poly,
  geometry = sf::st_sfc(polys_list, crs = 4326)
)

# Create clip polygon
clip_poly_sf <- sf::st_sf(
  geometry = sf::st_sfc(
    sf::st_polygon(list(cbind(c(25, 75, 75, 25, 25), c(25, 25, 75, 75, 25)))),
    crs = 4326
  )
)

if (has_legacy) {
  polys_sp <- as(polys_sf, "Spatial")
  clip_poly_sp <- as(clip_poly_sf, "Spatial")
  
  bench5 <- microbenchmark(
    legacy_sp = {
      result <- rgeos::gIntersection(polys_sp, clip_poly_sp, byid = TRUE)
    },
    modern_sf = {
      result <- sf::st_intersection(polys_sf, clip_poly_sf)
    },
    times = 30
  )
  
  result5 <- format_benchmark(bench5)
  print(result5[, c("expr", "mean_ms", "median_ms")])
  
  speedup5 <- calculate_speedup(
    result5$mean_ms[result5$expr == "legacy_sp"],
    result5$mean_ms[result5$expr == "modern_sf"]
  )
  
  cat(sprintf("\nSpeedup: %.1f%%", speedup5))
  if (speedup5 >= 20) {
    cat(" ✓ PASS (>= 20%)\n")
  } else {
    cat(" ✗ FAIL (< 20%)\n")
  }
  
  all_results$vector_intersection <- list(
    speedup = speedup5,
    pass = speedup5 >= 20,
    details = result5
  )
} else {
  bench5 <- microbenchmark(
    modern_sf = {
      result <- sf::st_intersection(polys_sf, clip_poly_sf)
    },
    times = 30
  )
  result5 <- format_benchmark(bench5)
  print(result5[, c("expr", "mean_ms", "median_ms")])
  cat("\nLegacy comparison not available\n")
}

# ==============================================================================
# Benchmark 6: Rasterization (Vector to Raster)
# ==============================================================================

cat("\n\nBenchmark 6: Rasterization (1000 points to raster)\n")
cat(strrep("-", 60), "\n")

# Create test points with values
n_pts <- 1000
coords_val <- data.frame(
  x = runif(n_pts, 0, 100),
  y = runif(n_pts, 0, 100),
  value = rnorm(n_pts)
)

pts_sf <- sf::st_as_sf(coords_val, coords = c("x", "y"), crs = 4326)
pts_vect <- terra::vect(pts_sf)

# Create template raster
template_terra <- terra::rast(nrows = 100, ncols = 100,
                               xmin = 0, xmax = 100, ymin = 0, ymax = 100)
terra::crs(template_terra) <- "EPSG:4326"

if (has_legacy) {
  pts_sp <- as(pts_sf, "Spatial")
  template_raster <- raster::raster(template_terra)
  
  bench6 <- microbenchmark(
    legacy_raster = {
      r <- raster::rasterize(pts_sp, template_raster, field = "value")
    },
    modern_terra = {
      r <- terra::rasterize(pts_vect, template_terra, field = "value")
    },
    times = 50
  )
  
  result6 <- format_benchmark(bench6)
  print(result6[, c("expr", "mean_ms", "median_ms")])
  
  speedup6 <- calculate_speedup(
    result6$mean_ms[result6$expr == "legacy_raster"],
    result6$mean_ms[result6$expr == "modern_terra"]
  )
  
  cat(sprintf("\nSpeedup: %.1f%%", speedup6))
  if (speedup6 >= 20) {
    cat(" ✓ PASS (>= 20%)\n")
  } else {
    cat(" ✗ FAIL (< 20%)\n")
  }
  
  all_results$rasterization <- list(
    speedup = speedup6,
    pass = speedup6 >= 20,
    details = result6
  )
} else {
  bench6 <- microbenchmark(
    modern_terra = {
      r <- terra::rasterize(pts_vect, template_terra, field = "value")
    },
    times = 50
  )
  result6 <- format_benchmark(bench6)
  print(result6[, c("expr", "mean_ms", "median_ms")])
  cat("\nLegacy comparison not available\n")
}

# ==============================================================================
# Summary Report
# ==============================================================================

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("PERFORMANCE BENCHMARK SUMMARY\n")
cat(strrep("=", 70), "\n\n")

if (has_legacy && length(all_results) > 0) {
  # Calculate overall statistics
  speedups <- sapply(all_results, function(x) x$speedup)
  passes <- sapply(all_results, function(x) x$pass)
  
  cat("Individual Test Results:\n")
  cat(strrep("-", 70), "\n")
  for (test_name in names(all_results)) {
    result <- all_results[[test_name]]
    status <- if (result$pass) "✓ PASS" else "✗ FAIL"
    cat(sprintf("%-30s: %6.1f%% speedup  %s\n",
                test_name, result$speedup, status))
  }
  
  cat("\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("Overall Statistics:\n"))
  cat(sprintf("  Mean speedup:     %.1f%%\n", mean(speedups)))
  cat(sprintf("  Median speedup:   %.1f%%\n", median(speedups)))
  cat(sprintf("  Min speedup:      %.1f%%\n", min(speedups)))
  cat(sprintf("  Max speedup:      %.1f%%\n", max(speedups)))
  cat(sprintf("  Tests passed:     %d / %d (%.0f%%)\n",
              sum(passes), length(passes), sum(passes) / length(passes) * 100))
  
  cat("\n")
  if (all(passes)) {
    cat("✓ ALL TESTS PASSED: Performance improvement >= 20% achieved!\n")
  } else if (mean(speedups) >= 20) {
    cat("⚠ PARTIAL PASS: Mean speedup >= 20% but some tests below threshold\n")
  } else {
    cat("✗ TESTS FAILED: Performance improvement < 20%\n")
  }
  
  # Save results to file
  results_file <- "tests/benchmark/benchmark_results.rds"
  saveRDS(all_results, results_file)
  cat(sprintf("\nDetailed results saved to: %s\n", results_file))
  
} else {
  cat("Legacy packages not available - comparison benchmarks not run.\n")
  cat("Modern implementations tested successfully.\n")
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("Benchmark complete!\n")
cat(strrep("=", 70), "\n\n")
