rm(list=ls())

# === 1. load library ============
clib=c('rgdal', 'rgeos', 'raster', 'sp', 'fields', 'xts')
x=lapply(clib, library, character.only=T)
library(rSHUD)

# === 2. create directories  ============
dir.prj = '~/Documents/Ex_waerma'
dir.forc = file.path(dir.prj, 'forc')
dir.fig = file.path(dir.prj, 'figure')
dir.create(dir.forc, showWarnings = FALSE, recursive = TRUE)
dir.create(dir.fig, showWarnings = FALSE, recursive = TRUE)

# === 3. setup the project ============
prjname = 'waerma'
model.in <- file.path(dir.prj, 'input', prjname)
model.out <- file.path(dir.prj, 'output', paste0(prjname, '.out'))
fin=shud.filein(prjname, inpath = model.in, outpath = model.out )
if (dir.exists(model.in)){
  unlink(model.in, recursive = T, force = T)
}
dir.create(model.in, showWarnings = F, recursive = T)

# === 4. load and reproject data ============
data(waerma)
wbd=waerma[['wbd']]
meteosite = waerma[['meteosite']] # This is in GCS

crs.pcs = crs.Albers(wbd)
dem = projectRaster(waerma[['dem']], crs=crs.pcs)
wbd= spTransform(waerma[['wbd']], CRSobj = crs.pcs)
riv= spTransform(waerma[['riv']], CRSobj = crs.pcs)

# sl=mask(terra::terrain(dem, opt='slope', unit='tangent'), wbd)
# plot(sl)

r0.soil = waerma[['soil']]
att.soil = waerma[['att']]$soil
rcl.soil=cbind(att.soil[, 1], 1:nrow(att.soil))
r.soil = projectRaster(reclassify(r0.soil, rcl.soil), crs = crs.pcs, method="ngb")

r0.geol = waerma[['geol']]
att.geol = waerma[['att']]$geol
rcl.geol=cbind(att.geol[, 1], 1:nrow(att.geol))
r.geol = projectRaster(reclassify(r0.geol, rcl.geol), crs = crs.pcs, method="ngb")

r0.lc = waerma[['lc']]
att.lc = waerma[['att']]$lc
rcl.lc=cbind(att.lc[, 1], 1:nrow(att.lc))
r.lc = projectRaster(reclassify(r0.lc, rcl.lc), crs = crs.pcs, method="ngb")

tsd.forc = waerma$tsd$forc
tsd.lai = waerma$tsd$lai[[1]]

# === 5. some threshold for model deployment ============
AREA = 9853260 # KNOWN Area
a.max = 150*150;
q.min = 33;
tol.riv = 50
tol.wb = 50
aqd = 6
NX = AREA / a.max
years = seq(as.numeric(format(min(time(tsd.forc[[1]])), '%Y')), 
            as.numeric(format(max(time(tsd.forc[[1]])), '%Y')))
ndays = days_in_year(years)

# === 6. domain decomposition ============
# simplify the river network.
riv.simp = rgeos::gSimplify(riv, tol=tol.riv, topologyPreserve = T)
# riv.simp = sp.CutSptialLines(sl=riv.simp, tol=1000)

# desolve and simplify the watershed boundary
wb.dis = rgeos::gUnionCascaded(wbd)
wb.simp = rgeos::gSimplify(wb.dis, tol=tol.wb, topologyPreserve = T)

# !! Triangulation
tri = shud.triangle(wb=wb.simp,q=q.min, a=a.max, S=NX)
# generate  .sp.mesh 
pm=shud.mesh(tri,dem=dem, AqDepth = aqd)
sp.mesh=sp.mesh2Shape(pm=pm)
ncell = nrow(pm@mesh)
print(ncell)

# generate  .riv
pr=shud.river(riv.simp, dem = dem)
pr@rivertype$Width=c(3, 4, 5)
pr@rivertype$Depth=c(2, 3, 3)
pr@rivertype$BankSlope=c(3, 3, 3)

# === 7. TSD DATA ============
fns.meteo  =  paste0(meteosite$FILENAME, '.csv')
range(time(tsd.forc[[1]]))
for(i in 1:length(fns.meteo)){
  write.tsd(tsd.forc[[i]], file = file.path(dir.forc, fns.meteo[i]))
}
tsd.mf = MeltFactor(years = seq(as.numeric(format(min(time(tsd.forc[[1]])), '%Y')), 
                                as.numeric(format(max(time(tsd.forc[[1]])), '%Y'))))

# Coverage of meteorological sites.
sp.forc = ForcingCoverage(sp.meteoSite = meteosite,  
                                 filenames= fns.meteo,
                                 pcs=crs.pcs, gcs=crs(meteosite), 
                                 dem=dem, wbd=wbd)
write.forc(sp.forc@data, 
           path = './forc',
           # path = normalizePath(dir.forc),
           startdate = format(min(time(tsd.forc[[1]])), '%Y%m%d'), 
           file=fin['md.forc'])


# === 8. attributes ============
# generate  .sp.att
pa=shud.att(tri, r.soil = r.soil, r.geol = r.geol, r.lc=r.lc, r.forc = sp.forc)
head(pa)


# === 9. toplogical relation between river and triangle  ============

# Cut the rivers with triangles
spm = sp.mesh2Shape(pm)
crs(spm) =crs(riv)
spr=riv
sp.seg = sp.RiverSeg(spm, spr)
# Generate the River segments table
prs = shud.rivseg(sp.seg)

# Generate initial condition
pic = shud.ic(nrow(pm@mesh), nrow(pr@river), AqD = aqd)


go.plot <- function(){
  z=getElevation(pm = pm)
  loc = getCentroid(pm=pm)
  idx.ord = order(z)
  col=colorspace::diverge_hcl(length(z))
  plot(sp.mesh[idx.ord, ], axes=TRUE, col=col, lwd=.5); plot(spr , col='blue', add=T, lwd=3); 
  image.plot( legend.only=TRUE, zlim= range(z), col=col, horizontal = F,
              legend.args = list('text'='Elevation (m)', side=3, line=.05, font=2, adj= .2),
              smallplot= c(.79,.86, 0.20,0.4))
}
ia = getArea(pm=pm)
png(filename = file.path(dir.fig, paste0(prjname, '_mesh.png')), height = 9, width = 6, res = 400, units = 'in')
par(mfrow=c(2,1), mar=c(3, 3.5, 1.5,1) )
go.plot(); 
mtext(side=3, text = '(a)')
mtext(side=3, text = paste0('Ncell = ', ncell), line=-1)
hist(ia, xlab='', nclass=20, main='', ylab='')
mtext(side=3, text = '(b)')
mtext(side=2, text = 'Frequency', line=2)
mtext(side=1, text = expression(paste("Area (", m^2, ")")), line=2)
box(); 
# grid()
dev.off()


# === 10. configuration files  ============
# model configuration, parameter
cfg.para = shud.para(nday=ndays)

# calibration file
cfg.calib = shud.calib()

para.lc = lc.GLC()
para.soil = PTF.soil(att.soil[, -1])  # only 4 columns (Silt, clay, OM, bulk density) as input
para.geol = PTF.geol(att.geol[, -1])

# === 11. write  input files. ============
write.mesh(pm, file = fin['md.mesh'])
write.riv(pr, file = fin['md.riv'])
write.ic(pic, file = fin['md.ic'])

write.df(pa, file=fin['md.att'])
write.df(prs, file=fin['md.rivseg'])
write.config(cfg.para, fin['md.para'])
write.config(cfg.calib, fin['md.calib'])

write.tsd(tsd.lai, fin['md.lai'] )
write.tsd(tsd.mf, fin['md.mf'] )

write.df(para.lc, file=fin['md.lc'])
write.df(para.soil, file=fin['md.soil'])
write.df(para.geol, file=fin['md.geol'])
writeshape(riv.simp, file=file.path(dirname(fin['md.att']), 'riv'))
print(nrow(pm@mesh))

