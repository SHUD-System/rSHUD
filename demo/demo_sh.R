rm(list=ls())
clib=c('rgdal', 'rgeos', 'raster', 'sp', 'fields')
x=lapply(clib, library, character.only=T)

library(rSHUD)
prjname = 'sh'
model.in <- file.path('../demo/input', prjname)
model.out <- file.path('../demo/output', paste0(prjname, '.out'))
fin=shud.filein(prjname, inpath = model.in, outpath = model.out )

dir.create(model.in, showWarnings = F, recursive = T)

load("./data/sh.rda")
wbd=sh[['wbd']]
riv=sh[['riv']]
dem=sh[['dem']]
tsd.forc=sh[['forc']]

# sl = terrain(dem, v=slope, unit='tangent')
# cellStats(sl, quantile)

a.max = 200;
q.min = 33;
tol.riv = 5
tol.wb = 5
aqd=3
NX = 800
years = seq(as.numeric(format(min(time(tsd.forc)), '%Y')), 
            as.numeric(format(max(time(tsd.forc)), '%Y')))
ndays = days_in_year(years)


riv.simp = rgeos::gSimplify(riv, tol=tol.riv, topologyPreserve = T)
riv.simp = sp.CutSptialLines(sl=riv.simp, tol=20)

wb.dis = rgeos::gUnionCascaded(wbd)
wb.simp = rgeos::gSimplify(wb.dis, tol=tol.wb, topologyPreserve = T)

# shp.riv =raster::crop(riv.simp, wb.simp)
# shp.wb = raster::intersect( wb.simp, riv.simp)

tri = shud.triangle(wb=wb.simp,q=q.min, a=a.max, S=NX)
# generate  .sp.mesh 
pm=shud.mesh(tri,dem=dem, AqDepth = aqd)
sp.mesh=sp.mesh2Shape(pm=pm)
ncell = nrow(pm@mesh)
print(ncell)
# generate  .sp.att

pa=shud.att(tri)

# generate  .riv
pr=shud.river(riv.simp, dem = dem)

# Cut the rivers with triangles
spm = sp.mesh2Shape(pm)
crs(spm) =crs(riv)
spr=riv
sp.seg = sp.RiverSeg(spm, spr)
# Generate the River segments table
prs = shud.rivseg(sp.seg)

# Generate initial condition
pic = shud.ic(nrow(pm@mesh), nrow(pr@river), AqD = aqd, p1 = 0.2, p2=0.2)

go.plot <- function(){
  z=getElevation(pm = pm)
  loc = getCentroid(pm=pm)
  idx.ord = order(z)
  col=colorspace::diverge_hcl(length(z))
  plot(sp.mesh[idx.ord, ], axes=TRUE, col=col, lwd=.5); plot(spr , col='blue', add=T, lwd=3); 
  # image.plot( legend.only=TRUE, zlim= range(z), col=col, horizontal = T,legend.lab="Elevation (m)", 
  #             smallplot= c(.6,.9, 0.22,0.26))
  image.plot( legend.only=TRUE, zlim= range(z), col=col, horizontal = F,
              legend.args = list('text'='Elevation (m)', side=3, line=.05, font=2, adj= .2),
              smallplot= c(.79,.86, 0.20,0.4))
}
ia = getArea(pm=pm)
png(filename = '~/sh_mesh.png', height = 9, width = 6, res = 400, units = 'in')
par(mfrow=c(2,1), mar=c(3, 3.5, 1.5,1) )
go.plot(); 
mtext(side=3, text = '(a)')
hist(ia, xlab='', nclass=20, main='', ylab='')
mtext(side=3, text = '(b)')
mtext(side=2, text = 'Frequency', line=2)
mtext(side=1, text = expression(paste("Area (", m^2, ")")), line=2)
box(); 
# grid()
dev.off()


sp.c = SpatialPointsDataFrame(gCentroid(wb.simp, byid = TRUE), 
                              data=data.frame('ID' = 'forcing'), match.ID = FALSE)
sp.forc = ForcingCoverage(sp.meteoSite = sp.c,  pcs=crs(wb.simp), dem=dem, wbd=wbd)
write.forc(sp.forc@data, 
           path = file.path('./input', prjname),
           startdate = format(min(time(tsd.forc)), '%Y%m%d'), 
           file=fin['md.forc'])
write.tsd(tsd.forc, file = file.path(fin['inpath'], 'forcing.csv'))

# model configuration, parameter
cfg.para = shud.para(nday=ndays)
# calibration
cfg.calib = shud.calib()

para.lc = lc.NLCD(lc=42) # 42 is the forest in NLCD classes
para.soil = PTF.soil()
para.geol = PTF.geol()

tsd.lai =  LaiRf.NLCD(lc=42, years=years)
write.tsd(tsd.lai$LAI, file = fin['md.lai'])

tsd.mf = MeltFactor(years=years)
write.tsd(tsd.mf, file = fin['md.mf'])
# write  input files.
write.mesh(pm, file = fin['md.mesh'])
write.riv(pr, file=fin['md.riv'])
write.ic(pic, file=fin['md.ic'])

write.df(pa, file=fin['md.att'])
write.df(prs, file=fin['md.rivseg'])
write.config(cfg.para, fin['md.para'])
write.config(cfg.calib, fin['md.calib'])

write.df(para.lc, file=fin['md.lc'])
write.df(para.soil, file=fin['md.soil'])
write.df(para.geol, file=fin['md.geol'])
writeshape(riv.simp, file=file.path(dirname(fin['md.att']), 'riv'))
print(nrow(pm@mesh))

