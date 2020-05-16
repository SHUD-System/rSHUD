# rSHUD --- A toolbox for SHUD modeling system

- Lele Shu (lele.shu@gmail.com)
- www.shud.xyz

## SHUD
The Simulator for Hydrological Unstructured Domain (SHUD) is a multiprocess, multi-scale hydrologic model where the major hydrological processes are fully coupled using the semi-discrete finite volume method. 

This package can be used with the AutoSHUD project, that can build modeling domain automatically.

## Purpose of the package:
1. Convert the geospatial data into SHUD format. The tool kit is able to process the raster and vector data, then building the unstructured triangular mesh domain for PIHM.
2. Write/read the SHUD input files.
3. Read the SHUD output files.
4. Generate the calibration parameter set.
5. Time-Series analysis on hydrologic data
6. Two-dimensional and 3-Dimentional plot.
7. GIS analysis. Convert the unstructure data into spatial data (Shapefile or Raster)


## Installation
```
install.packages("devtools")
devtools::install_github("SHUD-System/rSHUD")
```

## Note:
Current rSHUD requires different version of RTriangle package. you must install that via github(Dec 2019):
```
install.packages("devtools")
devtools::install_github("davidcsterratt/RTriangle", subdir="pkg")
```

