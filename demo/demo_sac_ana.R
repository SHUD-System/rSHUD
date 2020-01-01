
rm(list=ls())
clib=c('rgdal', 'rgeos', 'raster', 'sp')
x=lapply(clib, library, character.only=T)

library(PIHMgisR)
#test_check("PIHMgisR")
prjname = 'sac'
PRJNAME=prjname
inapth = file.path(prjname)
pdir = getwd()
cdir <- file.path('../demo')
setwd(cdir)
PIHM(prjname = prjname)

spr=sp.riv2shp()
varname=c(paste0('eleq',c( 'sub', 'surf') ), 
          paste0('eley',c( 'surf','unsat', 'gw') ), 
          paste0('elev',c('prcp','etp','infil', 'rech') ),
          paste0('elev',paste0('et',0:2) ),
          paste0('rivq',c('down', 'sub', 'surf') ),
          paste0('rivy','stage') )
# undebug(BasicPlot)
x=BasicPlot(varname = varname, imap = TRUE, sp.riv = spr)

setwd(pdir)