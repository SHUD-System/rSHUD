rm(list=ls())
clib=c('rgdal', 'rgeos', 'raster', 'sp')
x=lapply(clib, library, character.only=T)
# fns=list.files('R/', glob2rx('*.R'), full.names = T)
# x=lapply(fns, source)
# fns=list.files('src/', glob2rx('*.cpp'), full.names = T)
# x=lapply(fns, Rcpp::sourceCpp)

library(PIHMgisR)
#test_check("PIHMgisR")
prjname = 'sh'

pihmout <- file.path('../demo', prjname)
fin <- PIHM.filein(prjname, indir = pihmout)
if (dir.exists(pihmout)){
  unlink(pihmout, recursive = T, force = T)
}
dir.create(pihmout, showWarnings = F, recursive = T)

a.max = 200;
q.min = 33;
tol.riv=50
tol.wb=50
aqd=3


data(sh)
wbd=sh[['wbd']]
riv=sh[['riv']]
dem=sh[['dem']]


riv.simp = rgeos::gSimplify(riv, tol=10, topologyPreserve = T)

wb.dis = rgeos::gUnionCascaded(wbd)
length(wb.dis)
wb.simp = rgeos::gSimplify(wb.dis, tol=10, topologyPreserve = T)

# shp.riv =raster::crop(riv.simp, wb.simp)
# shp.wb = raster::intersect( wb.simp, riv.simp)

tri = m.DomainDecomposition(wb=wb.simp,q=q.min, a=a.max)
# generate PIHM .mesh 
pm=pihmMesh(tri,dem=dem, AqDepth = aqd)

# generate PIHM .att
# debug(pihmAtt)
pa=pihmAtt(tri)

# generate PIHM .riv
pr=pihmRiver(riv.simp, dem = dem)
stop()
# Cut the rivers with triangles
spm = sp.mesh2Shape(pm)
crs(spm) =crs(riv)
spr=riv
sp.seg = sp.RiverSeg(spm, spr)
# Generate the River segments table
prs = pihmRiverSeg(sp.seg)

# Generate initial condition
pic = pihm.init(nrow(pm@mesh), nrow(pr@river), AqD = aqd)

plot(tri, asp=1, type='n'); plot(spr , add=T, lwd=3)

# model configuration, parameter
cfg.para = pihmpara(nday=365)
# calibration
cfg.calib = pihmcalib()

para.lc = PTF.lc(43)
para.soil = PTF.soil()
para.geol = PTF.geol()

# write PIHM input files.
writemesh(pm, file = fin['md.mesh'])
writeriv(pr, file=fin['md.riv'])
writeinit(pic, file=fin['md.ic'])

write.df(pa, file=fin['md.att'])
write.df(prs, file=fin['md.rivseg'])
write.config(cfg.para, fin['md.para'])
write.config(cfg.calib, fin['md.calib'])

write.df(para.lc, file=fin['md.lc'])
write.df(para.soil, file=fin['md.soil'])
write.df(para.geol, file=fin['md.geol'])
writeshape(riv.simp, file=file.path(dirname(fin['md.att']), 'riv'))
print(nrow(pm@mesh))
