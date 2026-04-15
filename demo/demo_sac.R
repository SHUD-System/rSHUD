# =============================================================================
# Demo: SHUD Model Construction for Sacramento Watershed
# =============================================================================
# This demo shows how to manually construct a SHUD model for the Sacramento watershed
# including mesh generation, river processing, and file writing
# 
# Author: rSHUD Team
# Date: 2023
# =============================================================================

# Clear workspace and set random seed for reproducibility
rm(list=ls())
set.seed(123)

# Load required libraries
clib=c('rgeos', 'raster', 'sp', 'sf')
x=lapply(clib, library, character.only=T)

# Load rSHUD package
library(rSHUD)

# Set up project paths
prjname = 'sac'
PRJNAME = prjname

# Note: Using tempdir() for demonstration, use actual path in production
inpath <- file.path('../demo/sac', 'input', prjname)
outpath <- file.path('../demo/sac', 'output', paste0(prjname, '.out'))
fin <- shud.filein(prjname, inpath = inpath, outpath = outpath)

# Create output directories
pngout = file.path(inpath, 'fig')
gisout = file.path(inpath, 'gis')
dir.create(inpath, showWarnings = F, recursive = T)
dir.create(pngout, showWarnings = F, recursive = T)
dir.create(gisout, showWarnings = F, recursive = T)

message('Input directory: ', inpath)
message('Output directory: ', outpath)

# Model parameters
a.max = 1e6 * .2;    # Maximum triangle area
q.min = 33;           # Minimum angle
tol.riv = 200;        # River simplification tolerance
tol.wb = 200;         # Watershed boundary simplification tolerance
tol.len = 500;        # Length simplification tolerance
AqDepth = 10;         # Aquifer depth
years = 2000:2010;    # Simulation years
ny = length(years)
nday = 365*ny + round(ny/4)

message('Model period: ', min(years), ' to ', max(years))
message('Total days: ', nday)

# Load Sacramento watershed data
data(sac)
wbd = sac[['wbd']]
riv = sac[['riv']]
dem = sac[['dem']]
rsoil = sac[['rsoil']]
rgeol = sac[['rsoil']]
rlc = sac[['rlc']]
asoil = sac[['asoil']]
ageol = sac[['asoil']]
alc = sac[['alc']]
sp.forc = sac[['forc']]

# Visualize input data
if(interactive()) {
  plot(riv, main='River Network')
  plot(sp.forc, main='Forcing Sites')
}

# Create buffer around watershed
wbbuf = rgeos::gBuffer(wbd, width = 2000)
dem = raster::crop(dem, wbbuf)

# Save initial data visualization
png(file = file.path(pngout, 'data_0.png'), height=11, width=11, res=100, unit='in')
plot(dem, main='DEM with Watershed and Rivers')
plot(wbd, add=T, border=2, lwd=2)
plot(riv, add=T, lwd=2, col=4)
dev.off()

# Simplify river network (skipping simplification to preserve topology)
riv.s2 = sf::st_as_sf(riv)

if(interactive()) {
  plot(riv, main='River Network')
}

# Simplify watershed boundary
wb.dis = sf::st_union(sf::st_as_sf(wbd))
wb.s1 = sf::st_simplify(wb.dis, dTolerance=tol.wb)
wb.s2 = sf::st_cast(wb.s1, "POLYGON")
wb.s2 = sf::st_sf(geometry = wb.s2)

# Save simplified data visualization
png(file = file.path(pngout, 'data_1.png'), height=11, width=11, res=100, unit='in')
plot(dem, main='DEM with Simplified Watershed and Rivers')
plot(wb.s2, add=T, border=2, lwd=2)
plot(riv.s2, add=T, lwd=2, col=4)
dev.off()

# Set final simplified geometries
wb.simp = wb.s2
riv.simp = riv.s2

# Generate triangular mesh
tri = shud.triangle(wb=wb.simp, q=q.min, a=a.max)
message('Generated ', nrow(tri$T), ' triangles')

if(interactive()) {
  plot(tri, asp=1, main='Triangular Mesh')
}

# Generate SHUD mesh
pm = shud.mesh(tri, dem=dem, AqDepth = AqDepth)
spm = mesh_to_sf(pm, crs = sf::st_crs(riv))
writeshape(spm, file=file.path(gisout, 'domain.shp'))

message('Mesh generated with ', nrow(pm@mesh), ' cells')

# Generate mesh attributes
pa = shud.att(tri, r.soil = rsoil, r.geol = rgeol, r.lc = rlc, r.forc = sp.forc)

# Generate river network
pr = shud.river(riv, dem)
message('River network generated with ', nrow(pr@river), ' reaches')

# Generate IC
pic = shud.ic(nrow(pm@mesh), nrow(pr@river))

# Generate topological relation between river and mesh
spm = mesh_to_sf(pm)
sf::st_crs(spm) = sf::st_crs(riv)
prs = shud.rivseg(spm, sf::st_as_sf(riv))

# Generate parameter configurations
cfg.para = shud.para(nday=nday)
cfg.calib = shud.calib()

# Write SHUD input files
write_mesh(pm, file = fin['md.mesh'])
write_river(pr, file = fin['md.riv'])
write_ic(pic, file = fin['md.ic'])
write_df(pa, file=fin['md.att'])
write_df(prs, file=fin['md.rivseg'])
write_config(cfg.para, fin['md.para'])
write_config(cfg.calib, fin['md.calib'])

# Generate forcing data
forc.fns = paste0(sp.forc@data[, 'NLDAS_ID'], '.csv')
message('Forcing files: ', paste(forc.fns, collapse=', '))

# Note: Using tempdir() for demonstration, use actual path in production
write_forc(sp.forc@data, path=file.path('../demo/sac', 'forcing'), file=fin['md.forc'])

# Save river shapefile
spr = sf::st_as_sf(riv)
writeshape(spr, file=file.path(gisout, 'river.shp'))

message('Demo completed successfully!')
message('All files written to: ', inpath)
