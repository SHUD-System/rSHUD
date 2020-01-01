rm(list=ls())
clib=c('rgdal', 'rgeos', 'raster', 'sp')
x=lapply(clib, library, character.only=T)

library(PIHMgisR)
#test_check("PIHMgisR")
prjname = 'snp'
PRJNAME=prjname
inapth = file.path(prjname)

pihmout <- file.path('../demo', prjname)
fin <- PIHM.filein(prjname, indir = pihmout)
x=list.files(pihmout, pattern = glob2rx(paste0(prjname, '.*.*')), full.names = T)
file.remove(x)

pngout = file.path(pihmout, 'fig')
gisout = file.path(pihmout, 'gis')
dir.create(pihmout, showWarnings = F, recursive = T)
dir.create(pngout, showWarnings = F, recursive = T)
dir.create(gisout, showWarnings = F, recursive = T)

a.max = 1e6 * .2;
q.min = 33;
tol.riv = 200
tol.wb = 200
tol.lake = 200
tol.len = 500
AqDepth = 30
years = 2000:2010
ny=length(years)
nday = 365*ny +round(ny/4)
 # source('../PIHMgisR.code/Data_Sunapee.R')

data("sunapee")
indata = sunapee
wbd=indata[['wbd']]
riv=indata[['riv']]; plot(riv)
dem=indata[['dem']]
rsoil=indata[['rsoil']]
rgeol=indata[['rsoil']]
rlc=indata[['rlc']]
asoil=indata[['asoil']]
ageol=indata[['asoil']]
alc =indata[['alc']]
lake = indata[['lake']]
plot(dem)
plot(lake, add=T)

wbbuf = rgeos::gBuffer(wbd, width = 2000)
dem = raster::crop(dem, wbbuf)


png(file = file.path(pngout, 'data_0.png'), height=11, width=11, res=100, unit='in')
plot(dem); plot(wbd, add=T, border=2, lwd=2); plot(riv, add=T, lwd=2, col=4)
dev.off()

riv.s1 = rgeos::gSimplify(riv, tol=tol.riv, topologyPreserve = T)
riv.s2 = sp.simplifyLen(riv, tol.len)
plot(riv.s1); plot(riv.s2, add=T, col=3)
# lake.simp = rgeos::gSimplify(lake, tol=tol.lake, topologyPreserve = T)
# lk.sim = rgeos::gSimplify(lake, tol=tol.lake, topologyPreserve = T)

wblk = AddHoleToPolygon(wbd, lake)
riv.s2=raster::crop(riv.s2, wblk)
plot(riv.s2)

wb.dis = rgeos::gUnionCascaded(wblk)
wb.s1 = rgeos::gSimplify(wb.dis, tol=tol.wb, topologyPreserve = T)
wb.s2 = sp.simplifyLen(wb.s1, tol.len)

png(file = file.path(pngout, 'data_1.png'), height=11, width=11, res=100, unit='in')
plot(dem); 
plot(wb.s2, add=T, border=2, lwd=2); 
plot(riv.s2, add=T, lwd=2, col=4)
# plot(lk.sim, add=T, lwd=2, border=2)
dev.off()


# shp.riv =raster::crop(riv.simp, wb.simp)
# shp.wb = raster::intersect( wb.simp, riv.simp)
wb.simp = wb.s2
riv.simp = riv.s2

stop()
tri = m.DomainDecomposition(wb=wb.simp,q=q.min, a=a.max)
plot(tri, asp=1)

# generate PIHM .mesh 
pm=pihmMesh(tri,dem=dem, AqDepth = AqDepth)
sm = sp.mesh2Shape(pm)
writeshape(sm, crs(wbd), file=file.path(gisout, 'domain'))

# generate PIHM .att
# debug(pihmAtt)
pa=pihmAtt(tri, r.soil = rsoil, r.geol = rgeol, r.lc = rlc)
forc.fns = paste0(sp.forc@data[, 'NLDAS_ID'], '.csv')
forc.fns
writeforc(forc.fns, path='/Users/leleshu/Dropbox/PIHM/Projects/SAC/forcing/csv2000-2017', file=fin['md.forc'])

# generate PIHM .riv
pr=pihmRiver(riv.simp, dem)
# Correct river slope to avoid negative slope
pr = correctRiverSlope(pr)
# PIHMriver to Shapefile
sriv = sp.riv2shp(pr)
writeshape(sriv, crs(wbd), file=file.path(gisout, 'river'))

# Cut the rivers with triangles
sp.seg = sp.RiverSeg(pm, pr)
writeshape(sp.seg, crs(wbd), file=file.path(gisout, 'seg'))

# Generate the River segments table
prs = pihmRiverSeg(sp.seg)

# Generate initial condition
pic = pihm.init(nrow(pm@mesh), nrow(pr@river))

# Generate shapefile of river
spp.riv = sp.riv2shp(pr); 

# Generate shapefile of mesh domain
sp.dm = sp.mesh2Shape(pm)
png(file = file.path(pngout, 'data_2.png'), height=11, width=11, res=100, unit='in')
zz = sp.dm@data[,'Zsurf']
ord=order(zz)
col=terrain.colors(length(sp.dm))
plot(sp.dm[ord, ], col = col)
plot(spp.riv, col=spp.riv@data[,5] + 1 , add=T, lwd=3)
dev.off()

# model configuration, parameter
cfg.para = pihmpara(nday = nday)
# calibration
cfg.calib = pihmcalib()
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

# write PIHM input files.
writemesh(pm, file = fin['md.mesh'])
writeriv(pr, file=fin['md.riv'])
writeinit(pic, file=fin['md.ic'])

write.df(pa, file=fin['md.att'])
write.df(prs, file=fin['md.rivseg'])
write.df(para.lc, file=fin['md.lc'])
write.df(para.soil, file=fin['md.soil'])
write.df(para.geol, file=fin['md.geol'])

write.config(cfg.para, fin['md.para'])
write.config(cfg.calib, fin['md.calib'])
print(nrow(pm@mesh))
