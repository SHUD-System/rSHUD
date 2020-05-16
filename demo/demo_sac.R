rm(list=ls())
clib=c('rgdal', 'rgeos', 'raster', 'sp')
x=lapply(clib, library, character.only=T)
# fns=list.files('R/', glob2rx('*.R'), full.names = T)
# x=lapply(fns, source)
# fns=list.files('src/', glob2rx('*.cpp'), full.names = T)
# x=lapply(fns, Rcpp::sourceCpp)

library(rSHUD)
prjname = 'sac'
PRJNAME=prjname
inapth = file.path(prjname)

inpath <- file.path('../demo/input', prjname)
outpath <- file.path('../demo/outpath', paste0(prjname, '.out'))
fin <- shud.filein(prjname, inpath = inpath, outpath = outpath)
x=list.files(inpath, pattern = glob2rx(paste0(prjname, '.*.*')), full.names = T)
file.remove(x)

pngout = file.path(inpath, 'fig')
gisout = file.path(inpath, 'gis')
dir.create(inpath, showWarnings = F, recursive = T)
dir.create(pngout, showWarnings = F, recursive = T)
dir.create(gisout, showWarnings = F, recursive = T)

a.max = 1e6 * .2;
q.min = 33;
tol.riv = 200
tol.wb = 200
tol.len = 500
AqDepth = 10
years = 2000:2010
ny=length(years)
nday = 365*ny +round(ny/4)

data(sac)
wbd=sac[['wbd']]
riv=sac[['riv']]; plot(riv)
dem=sac[['dem']]
rsoil=sac[['rsoil']]
rgeol=sac[['rsoil']]
rlc=sac[['rlc']]
asoil=sac[['asoil']]
ageol=sac[['asoil']]
alc =sac[['alc']]
sp.forc =sac[['forc']]
plot(sp.forc)

wbbuf = rgeos::gBuffer(wbd, width = 2000)
dem = raster::crop(dem, wbbuf)

png(file = file.path(pngout, 'data_0.png'), height=11, width=11, res=100, unit='in')
plot(dem); plot(wbd, add=T, border=2, lwd=2); plot(riv, add=T, lwd=2, col=4)
dev.off()

riv.s1 = rgeos::gSimplify(riv, tol=tol.riv, topologyPreserve = T)
riv.s2 = sp.simplifyLen(riv, tol.len)
plot(riv.s1); plot(riv.s2, add=T, col=3)

wb.dis = rgeos::gUnionCascaded(wbd)
wb.s1 = rgeos::gSimplify(wb.dis, tol=tol.wb, topologyPreserve = T)
wb.s2 = sp.simplifyLen(wb.s1, tol.len)

png(file = file.path(pngout, 'data_1.png'), height=11, width=11, res=100, unit='in')
plot(dem); plot(wb.s2, add=T, border=2, lwd=2); 
plot(riv.s2, add=T, lwd=2, col=4)
dev.off()


# shp.riv =raster::crop(riv.simp, wb.simp)
# shp.wb = raster::intersect( wb.simp, riv.simp)
wb.simp = wb.s2
riv.simp = riv.s2

tri = shud.triangle(wb=wb.simp,q=q.min, a=a.max)

plot(tri, asp=1)

# generate shud .mesh 
pm=shud.mesh(tri,dem=dem, AqDepth = AqDepth)
spm = sp.mesh2Shape(pm, crs = crs(riv))
writeshape(spm, crs(wbd), file=file.path(gisout, 'domain'))

# generate shud .att
pa=shud.att(tri, r.soil = rsoil, r.geol = rgeol, r.lc = rlc, r.forc = sp.forc )
forc.fns = paste0(sp.forc@data[, 'NLDAS_ID'], '.csv')
forc.fns
write.forc(forc.fns, path='/Users/leleshu/Dropbox/PIHM/Projects/SAC/forcing/csv2000-2017', file=fin['md.forc'])

# generate shud .riv
# undebug(shud.river)
pr=shud.river(riv.simp, dem)
# Correct river slope to avoid negative slope
# pr = correctRiverSlope(pr)

# shudriver to Shapefile
# spr = sp.riv2shp(pr)
spr = riv
writeshape(spr, crs(wbd), file=file.path(gisout, 'river'))

# Cut the rivers with triangles
# sp.seg = sp.RiverSeg(pm, pr)
sp.seg=sp.RiverSeg(spm, spr)
writeshape(sp.seg, crs(wbd), file=file.path(gisout, 'seg'))

# Generate the River segments table
prs = shud.rivseg(sp.seg)

# Generate initial condition
pic = shud.ic(nrow(pm@mesh), nrow(pr@river))

# Generate shapefile of river
# spr = sp.riv2shp(pr); 

# Generate shapefile of mesh domain
sp.dm = sp.mesh2Shape(pm)
png(file = file.path(pngout, 'data_2.png'), height=11, width=11, res=100, unit='in')
zz = sp.dm@data[,'Zsurf']
ord=order(zz)
col=terrain.colors(length(sp.dm))
plot(sp.dm[ord, ], col = col)
plot(spr, add=T, lwd=3)
dev.off()

# model configuration, parameter
cfg.para = shud.para(nday = nday)
# calibration
cfg.calib = shud.calib()
#soil/geol/landcover
lc = unlist(alc)
para.lc = PTF.lc(lc)
para.soil = PTF.soil(asoil)
para.geol = PTF.geol(asoil)

# 43-mixed forest in NLCD classification
# 23-developed, medium           
# 81-crop land
# 11-water
lr=fun.lairl(lc, years=years)
png(file = file.path(pngout, 'data_lairl.png'), height=11, width=11, res=100, unit='in')
par(mfrow=c(2,1))
col=1:length(lc)
plot(lr$LAI, col=col, main='LAI'); legend('top', paste0(lc), col=col, lwd=1)
plot(lr$RL, col=col, main='Roughness Length'); legend('top', paste0(lc), col=col, lwd=1)
dev.off()
write.tsd(lr$LAI, file = fin['md.lai'])
write.tsd(lr$RL, file = fin['md.rl'])

#MeltFactor
mf = MeltFactor(years = years)
write.tsd(mf, file=fin['md.mf'])

# write shud input files.
write.mesh(pm, file = fin['md.mesh'])
write.riv(pr, file=fin['md.riv'])
write.ic(pic, file=fin['md.ic'])

write.df(pa, file=fin['md.att'])
write.df(prs, file=fin['md.rivseg'])
write.df(para.lc, file=fin['md.lc'])
write.df(para.soil, file=fin['md.soil'])
write.df(para.geol, file=fin['md.geol'])

write.config(cfg.para, fin['md.para'])
write.config(cfg.calib, fin['md.calib'])
# x=xts::xts(sin(1:100/10)*10+100+rnorm(100)*10, order.by=as.Date('2000-01-01')+1:100-1); colnames(x)='Q'
# write.tsd(x, file=file.path(inpath, paste0(prjname, '.tsd.obs')), backup = FALSE)
print(nrow(pm@mesh))
