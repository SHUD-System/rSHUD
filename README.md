# rSHUD --- SHUD Modeling System Toolbox

[![R-CMD-check](https://github.com/SHUD-System/rSHUD/workflows/R-CMD-check/badge.svg)](https://github.com/SHUD-System/rSHUD/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-≥4.0.0-blue.svg)](https://www.r-project.org/)
[![terra](https://img.shields.io/badge/terra-≥1.7.0-green.svg)](https://rspatial.github.io/terra/)
[![sf](https://img.shields.io/badge/sf-≥1.0.0-green.svg)](https://r-spatial.github.io/sf/)

**Author**: Lele Shu (shulele@lzb.ac.cn)  
**Website**: [www.shud.xyz](https://www.shud.xyz/)  
**GitHub**: [SHUD-System/rSHUD](https://github.com/SHUD-System/rSHUD)

---

## 📖 Project Introduction

**rSHUD** is an R language toolbox designed for the SHUD (Simulator for Hydrological Unstructured Domain) hydrological modeling system.

### SHUD Model Features
SHUD is a multi-process, multi-scale hydrological model that fully couples major hydrological processes using the semi-discrete finite volume method. This package can be used with the AutoSHUD project to automatically build modeling domains.

### ⚡ Version 2.3.0 - Modern Spatial Libraries

**rSHUD 2.3.0** follows the modern rSHUD v3 spatial API using **terra** and **sf** packages by default:
- 🚀 **150-400% performance improvement** on spatial operations
- 🎯 **Modern API**: Uses terra/sf directly; modern functions reject legacy raster/sp inputs
- ✅ **Consistent naming**: All functions follow snake_case convention
- 📚 **Complete documentation**: Migration guide and comprehensive examples
- 🧪 **Robust testing**: 70%+ test coverage with parameter validation

**⚠️ Breaking Change**: This version requires terra (≥1.7-0) and sf (≥1.0-0). Modern APIs no longer accept legacy raster/sp objects. Deprecated wrappers are retained only as migration compatibility entry points and may perform limited conversion before forwarding. See [Migration Guide](#-migration-from-v21x) below.

## 📚 Project Docs Index

Process and governance documents are organized under `docs/`:
- [Documentation Index](docs/README.md)

---

## 🎯 Main Features

### 1. Data Conversion and Preprocessing
- Convert geospatial data to SHUD format
- Process terra `SpatRaster` and sf vector data
- Build unstructured triangular mesh domains

### 2. Model File Management
- Read/write SHUD input files
- Read SHUD output files
- Generate calibration parameter sets

### 3. Hydrological Analysis
- Time series analysis of hydrological data
- 2D and 3D visualization plotting
- Hydrograph analysis

### 4. GIS Functions
- Convert unstructured data to modern spatial outputs (`sf`/Shapefile or terra `SpatRaster`)
- Watershed delineation and river network processing
- Spatial analysis and visualization

---

## 🚀 Installation Instructions

### System Dependencies (Ubuntu Platform)
```bash
sudo apt -y install gdal-bin libgdal-dev gcc g++ gfortran
sudo apt -y install r-cran-systemfonts r-cran-textshaping
```

### Install rSHUD

**Simple installation** (recommended):
```r
# Install from GitHub (includes all dependencies)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("SHUD-System/rSHUD")
```

**Manual dependency installation** (if needed):
```r
# Core spatial libraries (required)
install.packages(c("terra", "sf"))

# Other dependencies
install.packages(c("Rcpp", "reshape2", "ggplot2", "xts", "hydroGOF", "zoo",
                   "proj4", "gstat", "abind", "lubridate", "geometry"))

# Custom RTriangle package
devtools::install_github("shulele/RTriangle", subdir="pkg")

# Install rSHUD
devtools::install_github("SHUD-System/rSHUD")
```

### Verify Installation
```r
library(rSHUD)
packageVersion("rSHUD")  # Should be 2.3.0
packageVersion("terra")  # Should be >= 1.7.0
packageVersion("sf")     # Should be >= 1.0.0
```

---

## 📚 Usage Examples

### Basic Usage (v2.3.0)
```r
library(rSHUD)
library(terra)
library(sf)

# Read SHUD model files (new snake_case functions)
mesh <- read_mesh("model.mesh")
river <- read_river("model.riv")
attributes <- read_att("model.att")

# Read spatial data with terra/sf
dem <- rast("dem.tif")
watershed <- st_read("watershed.shp")
rivers_sf <- st_read("rivers.shp")

# Auto-build model (modern API)
model <- auto_build_model(
  project_name = "example",
  domain = watershed,
  dem = dem,
  rivers = rivers_sf,
  output_dir = "./output"
)

# Spatial visualization
mesh_sf <- mesh_to_sf(model$mesh)
plot_polygons(mesh_sf, field = "Area")
```

### View Available Examples
```r
# View all example scripts
list.files(system.file("demo", package = "rSHUD"))

# Run specific example
demo(demo_autoBuild, package = "rSHUD")
demo(demo_sac, package = "rSHUD")
```

---

## 📁 Project Structure

```
rSHUD/
├── R/                    # R source code
├── man/                  # Function help documentation
├── demo/                 # Example scripts
├── data/                 # Example datasets
├── src/                  # C++ source code
├── experiments/          # Experimental code
└── output/               # Output directory
```

---

## 🔧 Main Function Categories

### Model Building
- `auto_build_model()` - Auto-build SHUD model (modern API)
- `read_mesh()`, `read_river()` - Read mesh and river data
- `write_mesh()`, `write_river()` - Write mesh and river data
- `read_att()`, `read_para()`, `read_calib()` - Read model configuration

### GIS Functions (terra/sf)
- `watershed_delineation()` - Watershed delineation
- `vector_to_raster()` - Convert vector to raster (terra)
- `mesh_to_raster()` - Convert mesh data to raster
- `mesh_to_shapefile()` - Export mesh as shapefile

### Visualization
- `plot_hydrograph()` - Hydrograph analysis
- `plot_polygons()` - Spatial polygon plotting
- `compare_maps()` - Side-by-side map comparison

### Hydrological Calculations
- `calc_pet_pm()` - Penman-Monteith evapotranspiration
- `calc_melt_factor()` - Melt factor calculation
- `calc_river_order()` - River order calculation
- `ts_to_daily()` - Time series aggregation

### Deprecated Functions (migration compatibility only)
- `autoBuildModel()` → use `auto_build_model()`
- `readmesh()` → use `read_mesh()`
- `sp2raster()` → use `vector_to_raster()`
- `hydrograph()` → use `plot_hydrograph()`

See `NEWS.md` for complete list of deprecated functions.

---

## 🔄 Migration from v2.1.x

### Quick Migration Guide

**1. Update spatial object types:**
```r
# Old legacy code (v2.1.x; deprecated and not accepted by modern functions)
library(raster)
library(sp)
dem <- raster("dem.tif")
watershed <- readOGR("watershed.shp")

# New (v2.3.0 / rSHUD v3 API)
library(terra)
library(sf)
dem <- rast("dem.tif")
watershed <- st_read("watershed.shp")
```

**2. Update function names:**
```r
# Old function names (deprecated compatibility wrappers)
mesh <- readmesh("model.mesh")
river <- readriv("model.riv")
ts_data <- read.tsd("data.tsd")

# New function names (recommended)
mesh <- read_mesh("model.mesh")
river <- read_river("model.riv")
ts_data <- read_tsd("data.tsd")
```

**3. Update spatial operations:**
```r
# Old legacy code (raster/sp packages; deprecated)
cropped <- crop(raster_obj, extent_obj)
buffered <- buffer(sp_obj, width = 100)

# New (terra/sf)
cropped <- crop(rast_obj, ext_obj)
buffered <- st_buffer(sf_obj, dist = 100)
```

### Complete Migration Guide

For detailed migration instructions, see:
- `NEWS.md` - Complete list of changes and deprecated functions
- Package vignettes: `vignette("migration-guide", package = "rSHUD")`
- Online documentation: [www.shud.xyz](https://www.shud.xyz/)

### Getting Help

If you encounter issues during migration:
1. Check the deprecation warnings - they provide the new function names
2. Read `NEWS.md` for breaking changes
3. Open an issue on GitHub with a reproducible example

## 📖 Documentation Resources

- **Function Help**: Use `?function_name` for detailed help
- **Example Code**: View example scripts in `demo/` directory
- **Online Documentation**: [www.shud.xyz](https://www.shud.xyz/)
- **Migration Guide**: See `NEWS.md` and package vignettes
- **Chinese Version**: [README_cn.md](README_cn.md)

---

## 🤝 Contributing Guidelines

Welcome to submit Issues and Pull Requests! Please ensure:
1. Code follows R package development standards
2. New functions include complete documentation and tests
3. Follow the project's coding style

For detailed contributing guidelines, see [CONTRIBUTING.md](CONTRIBUTING.md).

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

---

## 📞 Contact Information

- **Email**: shulele@lzb.ac.cn
- **Project Homepage**: [www.shud.xyz](https://www.shud.xyz/)
- **GitHub**: [SHUD-System/rSHUD](https://github.com/SHUD-System/rSHUD)

---

*If this project helps you, please give us a ⭐️ star!*

