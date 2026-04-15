# rSHUD Performance Benchmarks

This directory contains performance benchmarks for the rSHUD modernization project. The benchmarks compare the performance of new terra/sf implementations against legacy raster/sp implementations to verify the >= 20% performance improvement requirement (Requirement 12).

## Files

- `bench-raster.R` - Main benchmark script comparing raster and vector operations
- `check_setup.R` - Setup verification script to check package availability
- `benchmark_results.rds` - Saved benchmark results (generated after running)

## Quick Start

```r
# 1. Check if your environment is ready
source("tests/benchmark/check_setup.R")

# 2. Install missing packages if needed
install.packages("microbenchmark")  # Required for benchmarks

# 3. Run benchmarks
source("tests/benchmark/bench-raster.R")
```

---

## Requirements

### Required Packages
- `microbenchmark` - For accurate performance measurements
- `terra` - Modern raster operations
- `sf` - Modern vector operations

### Optional Packages (for comparison)
- `raster` - Legacy raster operations
- `sp` - Legacy vector operations  
- `rgeos` - Legacy geometry operations

**Note**: If legacy packages are not installed, the benchmarks will only test modern implementations without comparison.

---

## Performance Targets

### Overall Target

**Minimum Performance Improvement: >= 20%**

All core operations should demonstrate at least 20% performance improvement when using terra/sf compared to legacy raster/sp implementations.

### Operation-Specific Targets

Based on the design document and terra/sf benchmarks from the literature:

| Operation Type | Target Improvement | Test Data Scale | Notes |
|---------------|-------------------|-----------------|-------|
| Raster Read/Write | 2-5x (100-400%) | 1000x1000 | terra uses more efficient I/O |
| Raster Cropping | 2-5x (100-400%) | 1000x1000 | terra's C++ backend is faster |
| Raster Resampling | 1.5-3x (50-200%) | 500x500 → 200x200 | Bilinear interpolation |
| Vector Buffer | 2-5x (100-400%) | 1000 points | sf's buffer is highly optimized |
| Vector Intersection | 2-5x (100-400%) | 100 polygons | sf uses spatial indexing |
| Rasterization | 2-5x (100-400%) | 1000 points | terra's rasterize is faster |
| Overall Workflow | 20-50% | Complete model | End-to-end improvement |

---

## Benchmark Tests

The benchmark suite includes the following tests:

1. **Raster Creation** (1000x1000 grid)
   - Tests basic raster object creation and value assignment
   - Compares `raster::raster()` vs `terra::rast()`

2. **Raster Cropping**
   - Tests spatial subsetting of raster data
   - Compares `raster::crop()` vs `terra::crop()`

3. **Raster Resampling** (bilinear interpolation)
   - Tests resampling to different resolution
   - Compares `raster::resample()` vs `terra::resample()`

4. **Vector Buffer** (1000 points)
   - Tests buffer operation on point geometries
   - Compares `rgeos::gBuffer()` vs `sf::st_buffer()`

5. **Vector Intersection** (100 polygons)
   - Tests spatial intersection operations
   - Compares `rgeos::gIntersection()` vs `sf::st_intersection()`

6. **Rasterization** (1000 points to raster)
   - Tests vector-to-raster conversion
   - Compares `raster::rasterize()` vs `terra::rasterize()`

---

## Running Benchmarks

### From R Console

```r
# Run all benchmarks
source("tests/benchmark/bench-raster.R")
```

### From Command Line

```bash
# Check setup first
Rscript tests/benchmark/check_setup.R

# Run benchmarks
Rscript tests/benchmark/bench-raster.R
```

---

## Interpreting Results

### Output Format

For each benchmark, the script reports:
- Mean execution time (milliseconds)
- Median execution time (milliseconds)
- Speedup percentage (improvement over legacy)
- Pass/Fail status (✓ PASS if >= 20% improvement)

### Summary Statistics

The final summary includes:
- Mean speedup across all tests
- Median speedup
- Min/Max speedup
- Number of tests passed
- Overall pass/fail status

### Success Criteria

#### Individual Tests
- Each benchmark test must show >= 20% improvement
- Tests that fail individual threshold should be investigated
- Acceptable reasons for < 20%:
  - Very small datasets where overhead dominates
  - Operations where legacy implementation was already optimal
  - System-specific factors (disk I/O, memory)

#### Overall Performance
- Mean speedup across all tests >= 20%
- Median speedup >= 20%
- At least 80% of tests pass individual threshold

#### Real-World Performance
- Complete model building workflow >= 20% faster
- Large dataset operations (>1GB) show >= 30% improvement
- Memory usage reduced or comparable

### Example Output

```
=== rSHUD Performance Benchmarks ===
Testing terra/sf vs raster/sp implementations
Target: >= 20% performance improvement

Benchmark 1: Raster Creation (1000x1000 grid)
------------------------------------------------------------
           expr  mean_ms median_ms
1 legacy_raster   125.3     123.1
2  modern_terra    45.2      44.8

Speedup: 63.9% ✓ PASS (>= 20%)

...

======================================================================
PERFORMANCE BENCHMARK SUMMARY
======================================================================

Individual Test Results:
----------------------------------------------------------------------
raster_creation               :   63.9% speedup  ✓ PASS
raster_crop                   :   78.2% speedup  ✓ PASS
raster_resample               :   45.1% speedup  ✓ PASS
vector_buffer                 :   52.3% speedup  ✓ PASS
vector_intersection           :   38.7% speedup  ✓ PASS
rasterization                 :   71.5% speedup  ✓ PASS

----------------------------------------------------------------------
Overall Statistics:
  Mean speedup:     58.3%
  Median speedup:   57.8%
  Min speedup:      38.7%
  Max speedup:      78.2%
  Tests passed:     6 / 6 (100%)

✓ ALL TESTS PASSED: Performance improvement >= 20% achieved!

Detailed results saved to: tests/benchmark/benchmark_results.rds
```

---

## Measurement Methodology

### Timing
- Use `microbenchmark` package for accurate measurements
- Run each test 30-100 times for statistical reliability
- Report mean and median times
- Calculate speedup as: `(legacy_time - modern_time) / legacy_time * 100`

### Test Data
- Use synthetic data with realistic sizes
- Test with multiple data scales (small, medium, large)
- Include edge cases (very small, very large)

### Environment
- Run on clean R session
- Disable parallel processing for consistency
- Document system specifications
- Run multiple times to account for variability

---

## Known Performance Characteristics

### terra Advantages
- **Memory efficiency**: Uses pointers instead of copying data
- **C++ backend**: Core operations in compiled code
- **Lazy evaluation**: Defers computation until needed
- **Chunked processing**: Handles large files efficiently

### sf Advantages
- **GEOS integration**: Efficient geometry operations
- **Spatial indexing**: Fast spatial queries
- **Tidy data**: Works well with dplyr/tidyverse
- **Standards compliance**: Follows OGC standards

### Potential Bottlenecks
- **Interpolation**: IDW and linear methods may be slower than nearest neighbor
- **Small datasets**: Overhead may dominate for very small data
- **Disk I/O**: File operations may be limited by disk speed
- **Memory**: Very large datasets may require chunked processing

---

## Expected Results

Based on terra/sf documentation and community benchmarks:

- **Best case**: 100-400% improvement (2-5x faster)
- **Typical case**: 50-100% improvement (1.5-2x faster)
- **Worst case**: 20-30% improvement (minimum acceptable)

The overall workflow should show 20-50% improvement as specified in the design document.

---

## Troubleshooting

### Legacy Packages Not Available

If you see "Legacy packages not available", install them with:

```r
install.packages(c("raster", "sp", "rgeos"))
```

**Note**: These are only needed for comparison benchmarks, not for rSHUD functionality.

### Memory Issues

For large datasets, you may need to increase R's memory limit:

```r
# Increase memory limit (Windows)
memory.limit(size = 8000)

# Or run with more memory from command line
R --max-mem-size=8G
```

### Slow Performance

Benchmarks can take several minutes to complete. To run faster (with less accuracy):

1. Reduce the `times` parameter in each `microbenchmark()` call
2. Reduce the size of test datasets (e.g., 500x500 instead of 1000x1000)

---

## Further Analysis

To analyze saved results:

```r
# Load results
results <- readRDS("tests/benchmark/benchmark_results.rds")

# Extract speedups
speedups <- sapply(results, function(x) x$speedup)
print(speedups)

# Plot results
barplot(speedups, 
        main = "Performance Improvements",
        ylab = "Speedup (%)",
        las = 2,
        col = ifelse(speedups >= 20, "green", "red"))
abline(h = 20, col = "blue", lty = 2)
```

---

## Validation Process

1. **Run Benchmarks**: Execute `bench-raster.R`
2. **Review Results**: Check individual and overall speedups
3. **Investigate Failures**: Analyze any tests < 20%
4. **Optimize**: Improve implementations if needed
5. **Document**: Record results and any issues
6. **Verify**: Test with real rSHUD workflows

---

## Reporting

After running benchmarks, document:

1. System specifications (CPU, RAM, OS, R version)
2. Package versions (terra, sf, raster, sp)
3. Individual test results (speedup %)
4. Overall statistics (mean, median, min, max)
5. Pass/fail status
6. Any anomalies or issues
7. Recommendations for optimization

Save results to `benchmark_results.rds` for future reference and comparison.

---

## References

- terra package: https://rspatial.org/terra/
- sf package: https://r-spatial.github.io/sf/
- Performance comparisons: https://geocompr.robinlovelace.net/
- rSHUD design document: `.kiro/specs/rshud-modernization/design.md`

---

## Notes

- Benchmarks use `microbenchmark` package for accurate timing
- Each test runs multiple iterations (30-100 times) for statistical reliability
- Results may vary based on system specifications and load
- Legacy packages are optional - benchmarks will run without them but won't show comparisons
- Detailed results are saved to `benchmark_results.rds` for further analysis
