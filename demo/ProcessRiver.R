# clib=c('rgdal', 'rgeos', 'raster', 'sp')
# x=lapply(clib, library, character.only=T)
# # fns=list.files('R/', glob2rx('*.R'), full.names = T)
# # x=lapply(fns, source)
# # fns=list.files('src/', glob2rx('*.cpp'), full.names = T)
# # x=lapply(fns, Rcpp::sourceCpp)
# 
# library(PIHMgisR)
# #test_check("PIHMgisR")
# prjname = 'sh'
# 
# pihmout <- file.path('../demo', prjname)
# fin <- PIHM.filein(prjname, indir = pihmout)
# if (dir.exists(pihmout)){
#   unlink(pihmout, recursive = T, force = T)
# }
# dir.create(pihmout, showWarnings = F, recursive = T)
# 
# a.max = 1000;
# q.min = 33;
# tol.riv=50 
# tol.wb=50
# # riv = readOGR('Experiments/riv.shp')
# # saveRDS(riv, file = 'data/River.rda')
# riv = readRDS(file = 'data/River.rda')
# # riv = rgeos::gSimplify(riv, tol = 200)
# # ord=PIHMgisR::sp.RiverOrder(riv)
# par(mfrow=c(1,1))
# 
# 
id2Line <- function(id, key){
  lp=key
  ct = count(as.numeric(id))
  cn = as.numeric( names(ct))
  # ctf = count(id[,1]) #count of frNode
  # ctt = count(id[,2])
  lfr = which(id[,1] == key)
  pto = which( cn == id[lfr,2] )
  if( ct[pto] == 2 ) {
    # 1 - header/outlet, >2 - joint point
    # message(key, '\t', pto)
    ret = id2Line(id, key = pto)
    return (c(lp, ret) )
  }else{
    return(c(key, cn[pto]) )
  }
}
library(PIHMgisR);
# rm(list=ls())
library(raster); library(sp); library(rgeos)
riv = readRDS(file = 'Experiments/River.rda')
# riv -- orig
# riv1 = simplified by length
# riv2 = gsimplified
# riv3 = Dissolved.
plotpp <- function(x, cols=c(2,3), pch=2,  add=T, ...){
  plot(x, col=cols[1], add = add, ...);
  points(extractCoords(x), col=cols[2], pch=pch)
}
pc <- function(x){
  plot(x, col=1:length(x))
}
par(mfrow=c(1,1))
# plotpp(riv, col=c(1,1), add=F, pch=1)
riv1 = SimplifybyLen(riv, split_length=500, plot.results = F)
plotpp(riv1, col=c(2,2), pch=2)

riv2 = gSimplify(riv1, tol=1500)
plotpp(riv2)
riv2 = riv1

riv3 = sp.DissolveLines(riv2)
pc(riv3)

cdu = extractCoords(sp, unique = T)
ft=FromToNode(sp)
ct= count(ft)
pnode = as.numeric(names(ct))[which(ct != 2)]

sp=riv3
sl=list()
k=0
nnode = length(pnode)
for(i in 1:nnode ){
  key = pnode[i]
  points(cdu[key, 1], cdu[key, 2], pch=20)
  if( length(which(ft[,1] == key)) > 0){
    k = k +1
    # message(k, '\t', i, '/', nnode )
    ids = id2Line(ft, key)
    # cdf[ids] = k
    points(cdu[ids, ], col=k, pch=k)
    # text(cdu[ids[1], 1], cdu[ids[1], 2], paste(k), col=2 )
    ln = sp::Line( cdu[ids,] )
    sl[[k]] = sp::Lines(ln, ID=paste(k) )
  }else{
    # print(key)
  }
}
sll = sp::SpatialLines(sl)
# ord= sp.RiverOrder(sll)
df = data.frame('INDEX'=1:length(sl))
sld = sp::SpatialLinesDataFrame(sll, data = df)
sld
sld
sp.RiverPath <- function(sp, idown = sp.RiverDown(sp)){
  coord = extractCoords(sp)
  ft = FromToNode(sp, coord = coord)
  nsp = length(sp)
  jid = as.numeric(names(count(ft, n = c(2:100 )) ) ) #id of joint points
  keys = c(which(idown<0), which(ft[,2] %in% jid) )
  updown = cbind(1:nsp, idown)
  goUp <- function(updown, id0){
    uid = which(idown == id0) # id of my upstream
    if(length(uid) == 1){ #up stream exist
      ret = c(goUp(updown, uid), id0)
    }else{
      ret = id0
    }
    ret
  }
  StreamPath = lapply(keys, function(x) goUp(updown, x))
  #========Point ID==========
  nstr = length(StreamPath)
  pl = list()
  for(i in 1:nstr){
    sid = StreamPath[[i]]
    pid = c(ft[sid, 1], ft[sid[length(sid)], 2])
    pl[[i]] = pid
  }
  #=========Spatial Lines===============
  ll = list()
  for(i in 1:nstr){
    sid = StreamPath[[i]]
    pid = c(ft[sid, 1], ft[sid[length(sid)], 2])
    ll[[i]] = sp::Lines( sp::Line(coord[pid, ] ), ID = i)
  }
  spx = sp::SpatialLines(ll, proj4string = raster::crs(sp))
  ret <- list(SegIDs = StreamPath,
              PointIDs= pl,
              sp = spx)
}

ord = sp.RiverOrder(sld)
rp = sp.RiverPath(sld)

# plot(sld, col=ord)
# uord = sort( unique(ord) )
# nord = length(uord)
# for(i in 1:nord){
#   iord = uord[i]
#   
# }

# plot(riv3, col=1:length(riv3), add=T);
# points(extractCoords(riv2), col=3, pch=3)
stop()
ord = sp.RiverOrder(riv2)
plot(riv2, col=ord)

# plot(riv3, col = ord)
#
uord = sort(unique(ord), decreasing = TRUE)
# nu = length(uord)
stop()
# for(i in 1:nu){
i=2
ro = uord[i]
sl = riv3[which(ord == ro),]
plot(sl, col=1:length(sl))
m.sl = rgeos::gLineMerge(sl)
m.sl
if(i == 1){
  sout = m.sl
}else{
  sout = gUnion(sout, m.sl)
}
# }
# plot(sout, col = 1:length(sout))
# m.sl = rgeos::gLineMerge(riv3, id=ord)
# plot(m.sl, col=1:length(m.sl))
#
# sp.RiverOrder(m.sl)