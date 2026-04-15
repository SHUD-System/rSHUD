# Migration Guide: raster/sp to terra/sf

## Overview

Starting from rSHUD version 2.1.0, the package has fully migrated to modern spatial data packages:
- **terra** (replaces raster)
- **sf** (replaces sp)

The old `raster` and `sp` packages are **no longer supported**.

## Why This Change?

1. **Performance**: terra is 10-100x faster than raster
2. **Memory**: Better memory management for large datasets
3. **Maintenance**: raster/sp are no longer actively maintained
4. **Modern API**: Cleaner, more consistent interface

## Quick Migration

### For Raster Data

**Old code (raster package):**
```r
library(raster)
dem <- raster("dem.tif")
slope <- terrain(dem, "slope")
```

**New code (terra package):**
```r
library(terra)
dem <- rast("dem.tif")
slope <- terrain(dem, "slope")
```

### For Vector Data

**Old code (sp package):**
```r
library(sp)
library(rgdal)
watershed <- readOGR("watershed.shp")
area <- gArea(watershed)
```

**New code (sf package):**
```r
library(sf)
watershed <- st_read("watershed.shp")
area <- st_area(watershed)
```

## Common Function Replacements

### Raster Operations

| raster package | terra package |
|----------------|---------------|
| `raster()` | `rast()` |
| `stack()` | `rast()` |
| `brick()` | `rast()` |
| `extent()` | `ext()` |
| `crs()` | `crs()` (same name) |
| `values()` | `values()` (same name) |
| `crop()` | `crop()` (same name) |
| `mask()` | `mask()` (same name) |
| `projectRaster()` | `project()` |
| `writeRaster()` | `writeRaster()` (same name) |

### Vector Operations

| sp/rgeos | sf |
|----------|-----|
| `readOGR()` | `st_read()` |
| `writeOGR()` | `st_write()` |
| `spTransform()` | `st_transform()` |
| `gBuffer()` | `st_buffer()` |
| `gIntersection()` | `st_intersection()` |
| `gUnion()` | `st_union()` |
| `gArea()` | `st_area()` |
| `gLength()` | `st_length()` |
| `gSimplify()` | `st_simplify()` |
| `coordinates()` | `st_coordinates()` |
| `bbox()` | `st_bbox()` |

## Using rSHUD Functions

All rSHUD functions now expect terra/sf objects. **Use terra and sf functions directly** - rSHUD does not provide wrapper functions:

```r
library(rSHUD)
library(terra)
library(sf)

# Load raster data with terra
dem <- rast("dem.tif")

# Load vector data with sf
watershed <- st_read("watershed.shp")

# Use terra/sf functions for spatial operations
dem_crop <- crop(dem, watershed)
slope <- terrain(dem_crop, "slope")

# Use rSHUD functions for hydrological modeling
mesh <- autoBuildModel(
  dem = dem_crop,
  boundary = watershed,
  # ... other parameters
)
```

**Important**: rSHUD does not wrap terra/sf functions. Learn and use terra/sf directly:
- terra documentation: https://rspatial.github.io/terra/
- sf documentation: https://r-spatial.github.io/sf/

## Error Messages

If you try to use old format objects, you'll see clear error messages:

```r
# Using old raster object
library(raster)
dem_old <- raster("dem.tif")

# This will fail with helpful message:
result <- some_rshud_function(dem_old)
# Error: Input must be a SpatRaster object.
# Old raster/sp formats are no longer supported.
# Please use terra::rast() to load your data.
```

## Installation

Make sure you have the required packages:

```r
install.packages(c("terra", "sf"))
```

For Ubuntu/Debian users, you may need system dependencies:

```bash
sudo apt install gdal-bin libgdal-dev
```

## Need Help?

- terra documentation: https://rspatial.github.io/terra/
- sf documentation: https://r-spatial.github.io/sf/
- Report issues: https://github.com/SHUD-System/rSHUD/issues

## Benefits You'll See

1. **Faster processing**: Especially for large rasters
2. **Lower memory usage**: Handle bigger datasets
3. **Better error messages**: Clearer debugging
4. **Future-proof**: Active development and support
5. **Consistent API**: Easier to learn and use

## Example: Complete Workflow Migration

**Old workflow:**
```r
library(raster)
library(sp)
library(rgeos)

# Load data
dem <- raster("dem.tif")
watershed <- readOGR("watershed.shp")

# Process
dem_crop <- crop(dem, watershed)
dem_mask <- mask(dem_crop, watershed)
slope <- terrain(dem_mask, "slope")

# Buffer watershed
watershed_buf <- gBuffer(watershed, width = 1000)
area <- gArea(watershed_buf)

# Save
writeRaster(slope, "slope.tif")
writeOGR(watershed_buf, ".", "watershed_buffer", driver = "ESRI Shapefile")
```

**New workflow:**
```r
library(terra)
library(sf)

# Load data
dem <- rast("dem.tif")
watershed <- st_read("watershed.shp")

# Process
dem_crop <- crop(dem, watershed)
dem_mask <- mask(dem_crop, watershed)
slope <- terrain(dem_mask, "slope")

# Buffer watershed
watershed_buf <- st_buffer(watershed, dist = 1000)
area <- st_area(watershed_buf)

# Save
writeRaster(slope, "slope.tif")
st_write(watershed_buf, "watershed_buffer.shp")
```

Notice how the new code is:
- Cleaner (fewer packages)
- More consistent (st_* prefix for sf functions)
- Faster (especially for large files)
