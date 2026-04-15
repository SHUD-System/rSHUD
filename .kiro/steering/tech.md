# Technology Stack

## Language & Build System

- **Language**: R (≥ 3.5.0)
- **Build System**: Standard R package structure with Rcpp integration
- **Package Manager**: devtools for development, standard R package installation

## Core Dependencies

### Spatial/GIS Libraries
- `raster` (≥ 2.1.0) - Raster data manipulation
- `sp` - Spatial data classes and methods
- `rgeos` - Geometry operations (simplification, buffering, unions)
- `proj4` - Coordinate projection transformations
- `gstat` - Geostatistical modeling
- `terra` - Modern raster/vector processing (suggested)

### Mesh Generation
- `RTriangle` - Constrained Delaunay triangulation (custom fork: shulele/RTriangle/pkg)
- `geometry` - Computational geometry
- `deldir` - Delaunay triangulation and Voronoi diagrams

### Time Series & Analysis
- `xts` - Extensible time series
- `zoo` - Time series infrastructure
- `hydroGOF` - Hydrological goodness-of-fit metrics
- `lubridate` - Date/time manipulation

### Data Processing
- `reshape2` - Data reshaping
- `abind` - Array binding
- `doParallel` - Parallel processing support

### Visualization
- `ggplot2` - Graphics
- `gridExtra`, `grid` - Grid graphics
- `fields` - Spatial data visualization

### I/O & Utilities
- `ncdf4` - NetCDF file handling
- `Rcpp` - C++ integration for performance-critical functions

## System Dependencies (Ubuntu/Debian)

```bash
sudo apt -y install gdal-bin libgdal-dev gcc g++ gfortran
sudo apt -y install r-cran-systemfonts r-cran-textshaping
```

## Common Commands

### Installation

```r
# Install from GitHub
devtools::install_github("SHUD-System/rSHUD")

# Install custom RTriangle dependency
devtools::install_github("shulele/RTriangle", subdir="pkg")
```

### Building Package (Development)

```r
# Build and check package
devtools::check()

# Build documentation
devtools::document()

# Run tests
devtools::test()

# Load package for development
devtools::load_all()
```

### Running Examples

```r
# Load package
library(rSHUD)

# List available demos
list.files(system.file("demo", package = "rSHUD"))

# Run specific demo
demo(demo_autoBuild, package = "rSHUD")
demo(demo_sac, package = "rSHUD")
```

### Building from Source

```bash
# Standard R package build
R CMD build .
R CMD check rSHUD_*.tar.gz
R CMD INSTALL rSHUD_*.tar.gz
```

## C++ Integration

The package uses Rcpp for performance-critical operations:
- `src/polygonArea.cpp` - Polygon area calculations
- `src/triTopology.cpp` - Triangle topology operations
- Compiled objects: `.o` files and shared library `.so`

## Documentation System

- **roxygen2** (7.2.1) - Documentation generation from inline comments
- **Markdown support** enabled in roxygen
- Man pages in `man/` directory (auto-generated, do not edit manually)

## File Formats

SHUD uses custom text-based formats with tab-separated values:
- `.mesh` - Mesh domain (triangles and points)
- `.riv` - River network
- `.att` - Attributes (soil, geology, land cover, forcing)
- `.ic` - Initial conditions
- `.para` - Model parameters
- `.calib` - Calibration parameters
- `.forc` - Forcing data file list
- `.tsd` - Time series data (with header: nrow, ncol, start_date, end_date, dt)
