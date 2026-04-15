# install.packages('png')
library(SHUDtoolbox)
library(raster)
library(png) 
library(rgeos)
library(rasterVis)

rm(list=ls())

makeLogo <- function(fn, dx, nbuf, ratio, ncell, ncolor, str='SHUD', wd, ht,
                     fout){
  graphics.off()
  go.png <- function(fn, 
                     str='SHUD', ht, wd, family='mono'){
    # str='SHUD'; ht=55;family='mono'
    nc = nchar(str)
    png(fn, width=wd, height=ht)
    par(mar=c(0,0,0,0))
    plot(0,0, type='n', axes=F, xlim=c(0, wd), ylim=c(0,ht))
    text(wd/2, ht/2, str, cex=wd/35, font=2, family=family)
    # grid()
    dev.off()
  }; 
  # go.png('shud.png', ht=50, ratio=ratio, family='mono')
  
  zplane <- function(r){
    m=as.matrix(r)
    nr=nrow(m); nc=ncol(m)
    rb=raster( t(matrix(rep((nc:1/nc *2), nr), nrow=nc)), xmx=nc*dx, ymx=nr*dx, xmn=0, ymn=0);
    rb
  }
  go.raster <- function(fn, dx, ratio){
    # dx=100
    img2 <- readPNG(fn) 
    m=apply(img2[, , 1], 1:2, mean)
    nr=nrow(m); nc=ncol(m)
    range(m)
    m[m<.9]=-1; 
    m[m>=.9]=-2;
    rm=raster(m+2, xmn=0, xmx=nc*dx,  ymn=0, ymx=nr*dx ); rs=rm
    # x=t(matrix(rep(1:nc, nr), nc, nr))
    # y=matrix(rep(1:nr, nc), nr, nc)
    # m0=sqrt(x^2+y^2)
    # r0=raster(m0, xmn=0, xmx=nc*dx,  ymn=0, ymx=nr*dx/ratio )
    # rs=rm+r0; rs
    # plot(rm+r0)
    rs
  }
  w.fun <- function(n,simple=FALSE){
    if(simple) return(matrix(1, n,n))
    m=round(n/2) +1
    x=expand.grid(1:n - m, 1:n -m)
    z=sqrt(x[,1]^2 + x[,2]^2)
    zm=max(z)
    mat=matrix((zm-z)/zm, n,n)
    # mat[lower.tri(mat)] <- 0
    mat
  }
  message('make png')
  go.png(fn = fn, wd=wd, ht=ht)
  message('make raster')
  re = go.raster(fn = fn, dx = dx, ratio=ratio)
  # w=w.fun(nbuf, simple=T)
  w=matrix(nbuf*nbuf, nbuf, nbuf)
  message('make DEM')
  dem=focal(re, w=w, fun=mean, na.rm=TRUE)
  
  ext=extent(dem)
  
  bnd=fishnet(ext=ext+c(1,-1,1,-1)*6*dx, n=1);
  AA=diff(ext[1:2]) * diff(ext[3:4])
  
  message('Build Mesh')
  x=shud.triangle(wb=bnd, a=AA/ncell, q=33); #plot(x, type='n')
  pm=shud.mesh(x, dem=dem)
  spm=mesh_to_sf(pm=pm)
  z=spm@data$Zmax; 
  col=(colorspace::diverge_hcl(ncolor, h=c(246,100), c=96, alpha =c(.85, 1)) )
  # col=(colorspace::diverge_hcl(ncolor, h=c(246,40), c=96, alpha =c(.85, 1)) )
  v <- seq(round(min(z)), round(max(z)), length.out = 100) 
  cols=col[findInterval(z, v)]
  # plot(spm, col=col, lwd=0.2)
  # stop()
  # dev.off()
  go.plot <- function(){
    par(mar=rep(0,4),bg=NA); 
    plot(spm, col=cols, lwd=0.4) 
    # plot(dem, legend.only=TRUE, col=col, horizontal=TRUE,
    #      smallplot=c(0.1, 0.9, 0.1, 0.20),
    #      axis.args=list(col.axis='blue', lwd = 0,
    #                     font.axis=1, cex.axis=0.1,tck = 0, line=-.85),
    #      legend.args=list(text='') )
  }; 
  go.plot()
  # fout='/Users/leleshu/Dropbox/SHUD/github/website/static/img/shud.png'
  png(fout, width=wd, height = ht)
  go.plot()
  dev.off()
  
}

# fn= file.path(tmpDir(), 'shud.png')
fn= file.path('shud.png')
dx=100
nbuf = 13
ncell=2550
ncolor = 100

wd=800
ratio = 3/4
makeLogo(fn=fn, dx=dx, nbuf=nbuf, ratio=ratio, ncell=ncell,
         ncolor=ncolor, str='SHUD', 
         fout='/Users/leleshu/Dropbox/SHUD/github/shud.xyz/static/img/shud.png',
         wd=wd, ht=ratio*wd)

wd=400
ratio = 2/5
makeLogo(fn=fn, dx=dx, nbuf=nbuf, ratio=ratio, ncell=ncell,
         ncolor=ncolor, str='SHUD', 
         fout='/Users/leleshu/Dropbox/SHUD/github/shud.xyz/static/img/logo.png',
         wd=wd, ht=ratio*wd)


# plot(dem)
stop()
library(rgl)
tri3d <-function(pm, zscale){
  arr=getVertex(pm)
  x = t(arr[, , 1])
  y = t(arr[, , 2])
  z = t(arr[, , 4])*zscale
  nz=length(z)
  # col=terrain.colors(length(z))
  col=(colorspace::diverge_hcl(nz, h=c(246,40), c=96))
  cm = matrix(col[order(z)], nrow=3)
  rgl::rgl.clear();
  # rgl::par3d(windowRect = c(10, 10, 800, 800))
  Sys.sleep(0.1)
  rgl::rgl.pop("lights")
  rgl::light3d(theta=40, phi=45, specular="gray")
  rgl::triangles3d(x,y,z,color=col,box=TRUE)
}
tri3d(pm, zscale=20)


go.shape <- function(msk, fn=file.path(tmpDir(), 'shud.png')){
  img2 <- readPNG(fn) 
  m=apply(img2[, , 1], 1:2, mean)
  # image(m)
  nr=nrow(m)
  nc=ncol(m)
  m[m>=.9]=2;
  m[m<.9]=1;
  # r=raster(m, xmx=nc*dx, ymx=nr*dx, xmn=0, ymn=0);
  r=msk
  r[] = m
  spr=rasterToPolygons(r, dissolve = TRUE)
  wb = spr[1, ]
  writeshape(wb, file=file.path(tmpDir(), 'shud') )
  # save(wb, file = 'data/shud.rda')
  wb
}