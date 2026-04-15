rm(list=ls())

# ===# === pre1. load library ============
clib=c('rgeos', 'raster', 'sp', 'sf', 'fields', 'xts', 'ggplot2')
x=lapply(clib, library, character.only=T)
library(rSHUD)

# === pre2. create directories  ============
dir.prj = '../demo/waerma'
dir.forc = file.path(dir.prj, 'forc')
dir.fig = file.path(dir.prj, 'figure')
dir.create(dir.forc, showWarnings = FALSE, recursive = TRUE)
dir.create(dir.fig, showWarnings = FALSE, recursive = TRUE)

# === pre3. setup the project ============
prjname = 'waerma'
model.in <- file.path(dir.prj, 'input', prjname)
model.out <- file.path(dir.prj, 'output', paste0(prjname, '.out'))

fin=shud.filein(prjname, inpath = model.in, outpath = model.out )
shud.env(prjname, inpath = model.in, outpath = model.out )
dir.create(model.in, showWarnings = F, recursive = T)

ia=getArea()
ncell=length(ia)

gplotfun <- function(r, leg.lab='value'){
  if (inherits(r, "SpatRaster")) {
    df <- as.data.frame(r, xy = TRUE)
  } else {
    map.p <- rasterToPoints(r)
    df <- data.frame(map.p)
  }
  colnames(df) <- c('X', 'Y', 'MAP')
  g= ggplot(data=df, aes(y=Y, x=X)) +
    geom_raster(aes(fill=MAP)) +
    # geom_point(data=sites, aes(x=x, y=y), color="white", size=3, shape=4) +
    theme_bw() + coord_equal() +
    # scale_fill_continuous(leg.lab) +
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
xt = 1:(365*2)+365*1
# xt=1:nrow(qdown)
q=qdown[xt, oid]
pq = cbind(q, rowMeans(prcp[xt,]))[,2:1]
# gl[[1]] = autoplot(q)+xlab('')+ylab('Discharge (m^3/day)')+theme_bw()
gl[[1]] = plot_hydrograph(pq, ylabs = c('Preciptation (mm)', 'Discharge (cmd)'))
plot(gl[[1]] )

# === 2. plot Water Balance ============
xl=loaddata(varname=c('rivqdown', 'eleveta', 'elevetp', 'elevprcp', 'eleygw'))
wb=wb.all(xl=xl, plot=F)[(1:24)+12, ]*1000
gl[[2]] = plot_hydrograph(wb, ylabs = c('Storage (mm)', 'Flux (mm/mon)'), legend.position='top')
gl[[2]]

# === 3. plot groundwater data ============
eleygw = readout('eleygw')[xt, ]
ts.gw=apply.daily(eleygw, sum)/ncell
# plot(ts.gw)  
gw.mean = apply(eleygw, 2, mean)
aqd =getAquiferDepth()
r.gw = mesh_to_raster(gw.mean)
d.gw = aqd - r.gw
d.gw[d.gw<0]=0
gl[[3]] =gplotfun(d.gw, leg.lab='Depth (m)')+
  scale_fill_gradient(low = "darkblue", high = "yellow")
gl[[3]]

# === 4. plot ETa data ============
eleveta = readout('eleveta')[xt, ]
ts.eta=apply.monthly(eleveta, sum)
# plot(ts.eta)  
eta.mean = apply(eleveta, 2, mean)*365
r.eta = mesh_to_raster(eta.mean)*1000 # mm/day
# plot(spm, axes=TRUE)
# plot(add=T, r.eta)
gl[[4]]=gplotfun(r.eta, leg.lab='Rate (mm/day) ')+
  scale_fill_gradient(low = "white", high = "blue")
gl[[4]]

# === Saving the plots ============
gg=gridExtra::arrangeGrob(grobs=gl, nrow=2, ncol=2)
ggsave(plot = gg, filename = file.path(dir.fig, 'waerma_res.png'), 
       width = 9, height=9, dpi=400, units = 'in')

for(i in 1:4){
  ggsave(plot = gl[[i]], filename = file.path(dir.fig, paste0('waerma_res_', i, '.png')),
         width = 4.5, height=5, dpi=400, units = 'in')
}
shud.calib()
