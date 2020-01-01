
rm(list=ls())
clib=c('rgdal', 'rgeos', 'raster', 'sp')
x=lapply(clib, library, character.only=T)

library(PIHMgisR)
data(sac)
indata = sac
sp.forc =indata[['forc']]
forc.fns = paste0(sp.forc@data[, 'NLDAS_ID'], '.csv')
forc.fns
pm=autoPIHMgis(sac, forcfiles = forc.fns, outdir='../demo')


spm=sp.mesh2Shape(pm)
plot(indata$dem)
plot(indata$wbd, add=T, lwd=2, border=2)
plot(spm, add=T, lwd=.2)
plot(indata$riv, add=T, col=4)
