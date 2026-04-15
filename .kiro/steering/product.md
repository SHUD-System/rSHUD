# Product Overview

rSHUD is an R package toolbox for the SHUD (Simulator for Hydrologic Unstructured Domains) hydrological modeling system.

## Purpose

SHUD is a multi-process, multi-scale hydrological model that fully couples major hydrological processes using the semi-discrete finite volume method on unstructured triangular meshes. rSHUD provides comprehensive tools to:

- Convert geospatial data (raster/vector) into SHUD model format
- Build unstructured triangular mesh domains for watershed modeling
- Process river networks and perform watershed delineation
- Read and write SHUD model input/output files
- Perform hydrological analysis and water balance calculations
- Visualize spatial and temporal hydrological data (2D/3D)

## Key Capabilities

- **Automated Model Building**: `autoBuildModel()` function automates the entire workflow from raw geospatial data to a complete SHUD model
- **GIS Processing**: Watershed delineation, river network simplification, mesh generation using constrained Delaunay triangulation
- **Hydrological Calculations**: PET (Penman-Monteith), melt factors, LAI/roughness length time series
- **Data I/O**: Comprehensive read/write functions for all SHUD file formats (.mesh, .riv, .att, .ic, .para, .calib, etc.)
- **Visualization**: Spatial plotting of mesh domains, river networks, and time series analysis

## Target Users

Hydrologists, environmental scientists, and researchers working with distributed hydrological modeling, particularly those using the SHUD modeling system or working with unstructured mesh-based models.
