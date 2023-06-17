rm(list=ls())

# === pre1. load library ============
clib=c('rgdal', 'rgeos', 'raster', 'sp', 'fields', 'xts', 'ggplot2')
x=lapply(clib, library, character.only=T)
library(rSHUD)

# === pre2. create directories  ============
dir.prj = '../demo'
dir.forc = file.path(dir.prj, 'forc')
dir.fig = file.path(dir.prj, 'figure')
dir.create(dir.forc, showWarnings = FALSE, recursive = TRUE)
dir.create(dir.fig, showWarnings = FALSE, recursive = TRUE)

# === pre3. setup the project ============
prjname = 'sh'
model.in <- file.path(dir.prj, 'input', prjname)
model.out <- file.path(dir.prj, 'output', paste0(prjname, '.out'))

fin=shud.filein(prjname, inpath = model.in, outpath = model.out )
shud.env(prjname, inpath = model.in, outpath = model.out )
dir.create(model.in, showWarnings = F, recursive = T)

ia=getArea()
ncell=length(ia)
spm=sp.mesh2Shape()

gplotfun <- function(r, leg.lab='value'){
  map.p <- rasterToPoints(r)
  #Make the points a dataframe for ggplot
  df <- data.frame(map.p)
  #Make appropriate column headings
  colnames(df) <- c('X', 'Y', 'MAP')
  #Now make the map
  g= ggplot(data=df, aes(y=Y, x=X)) +
    geom_raster(aes(fill=MAP)) +
    # geom_point(data=sites, aes(x=x, y=y), color=”white”, size=3, shape=4) +
    theme_bw() +
    coord_equal() +
    scale_fill_continuous(leg.lab) +
    theme(
      # axis.title.x = element_text(size=16),
      # axis.title.y = element_text(size=16, angle=90),
      # axis.text.x = element_text(size=14),
      # axis.text.y = element_text(size=14),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = 'right',
      legend.key = element_blank()
    )
  return(g)
}

gl=list()
# === 1. plot Q (discharge) data ============
oid = getOutlets()
qdown = readout('rivqdown')
prcp = readout('elevprcp')
xt = 1:(365*1)+365*2
q=qdown[xt, oid]
pq = cbind(q, rowMeans(prcp[xt,]))[,2:1]
# gl[[1]] = autoplot(q)+xlab('')+ylab('Discharge (m^3/day)')+theme_bw()
gl[[1]] = hydrograph(pq, ylabs = c('Preciptation (mm)', 'Discharge (cmd)'))
gl[[1]] 

# === 2. plot Water Balance ============
xl=loaddata(varname=c('rivqdown', 'eleveta', 'elevetp', 'elevprcp', 'eleygw'))
wb=wb.all(xl=xl, plot=F)[-(1:12)]*1000
gl[[2]] = hydrograph(wb, ylabs = c('Storage (mm)', 'Flux (mm/mon)'))
gl[[2]]

# === 3. plot groundwater data ============
eleygw = readout('eleygw')[xt, ]
ts.gw=apply.daily(eleygw, sum)/ncell
# plot(ts.gw)  
gw.mean = apply(eleygw, 2, mean)
r.gw = MeshData2Raster(gw.mean)
gl[[3]] =gplotfun(r.gw, leg.lab='Depth (m)')
gl[[3]]

# === 4. plot ETa data ============
eleveta = readout('eleveta')[xt, ]
ts.eta=apply.monthly(eleveta, sum)
# plot(ts.eta)  
eta.mean = apply(eleveta, 2, mean)
r.eta = MeshData2Raster(eta.mean)*1000 # mm/day
# plot(spm, axes=TRUE)
# plot(add=T, r.eta)
gl[[4]]=gplotfun(r.eta, leg.lab='Rate (mm/day) ')
gl[[4]]

# === Saving the plots ============
gg=gridExtra::arrangeGrob(grobs=gl, nrow=2, ncol=2)
ggsave(plot = gg, filename = file.path(dir.fig, 'waerma_res.png'), width = 7, height=7, dpi=400, units = 'in')

