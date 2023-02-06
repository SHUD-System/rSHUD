# rSHUD --- A toolbox for SHUD modeling system

- Lele Shu (shulele@lzb.ac.cn)
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


*If you work on Ubuntu platform, you must install software below before installing rSHUD and packages.*

```
sudo apt -y install gdal-bin libgdal-dev  gcc g++ gfortran
sudo apt -y install r-cran-systemfonts r-cran-textshaping
```

We suggest user install the dependant packages before install rSHUD.

```
libs = c('reshape2','xts','hydroGOF','zoo','RTriangle','proj4','gstat','abind','lubridate','interp','geometry',
         'testthat', 'rmarkdown', 'ncdf4', 'blogdown', 'doParallel', 'knitr', 'rmarkdown', 'deldir',
         'ncdf4', 'devtools')
nx = length(libs)
for(i in 1:nx){
  message(i, '/', nx, '\t',  libs[i])
  if(require(libs[i], character.only = TRUE)){
  }else{
    install.packages(libs[i], dependencies = TRUE, INSTALL_opts = '--no-lock')
  }
}

if(!require(devtools)){
install.packages("devtools", dependencies = TRUE, INSTALL_opts = '--no-lock')
}
devtools::install_github("shulele/RTriangle", subdir="pkg")

```
Then you may install the rSHUD directly from github.

```
if(!require(devtools)){ install.packages("devtools") }
devtools::install_github("SHUD-System/rSHUD")
```




