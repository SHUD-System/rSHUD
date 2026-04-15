# =============================================================================
# Demo: Auto-build SHUD Model
# =============================================================================
# This demo shows how to automatically build a SHUD model using the rSHUD package
# 
# Author: rSHUD Team
# Date: 2023
# =============================================================================

# Clear workspace and set random seed for reproducibility
rm(list=ls())
set.seed(123)

# Load required libraries
clib=c('raster', 'terra', 'sf', 'fields', 'xts')
x=lapply(clib, library, character.only=T)

# Load rSHUD package
library(rSHUD)

# Load example data (Sacramento watershed)
data(sac)
indata = sac

# Extract forcing data information
sp.forc = indata[['forc']]
forc.fns = paste0(sp.forc@data[, 'NLDAS_ID'], '.csv')
message('Forcing files: ', paste(forc.fns, collapse=', '))

# Auto-build SHUD model
# Note: Using tempdir() for demonstration, use actual path in production
outdir = file.path(tempdir(), 'demo_autobuild')
message('Output directory: ', outdir)

pm = auto_build_model(
  project_name = "sac",
  domain = sf::st_as_sf(indata[['wbd']]),
  dem = terra::rast(indata[['dem']]),
  rivers = sf::st_as_sf(indata[['riv']]),
  soil = terra::rast(indata[['rsoil']]),
  geology = terra::rast(indata[['rsoil']]),
  landcover = terra::rast(indata[['rlc']]),
  forcing_sites = sf::st_as_sf(sp.forc),
  forcing_files = sp.forc@data,
  output_dir = outdir,
  verbose = TRUE
)

# Convert mesh to shapefile for visualization
spm = mesh_to_sf(pm$mesh)

# Visualize the results
plot(indata$dem, main='DEM with Watershed and Mesh')
plot(sf::st_as_sf(indata$wbd), add=T, lwd=2, border=2)
plot(spm, add=T, lwd=.2)
plot(sf::st_as_sf(indata$riv), add=T, col=4)

message('Demo completed successfully!')
message('Model files saved to: ', outdir)
