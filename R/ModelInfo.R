
#' Summary of the basic information of model input data.
#' \code{ModelInfo} 
#' @param  path Where to export the basic model infomation.
#' @param  crs.pcs The crs4string that defines the (x,y) in the model.
#' @return Basic model infomation, figures and tables
#' @export
ModelInfo <- function(path = shud.filein()['outpath'],
                      spr = readriv.sp(),
                      crs.pcs=raster::crs(spr) ){
  # path=pp$outpath
  outdir = file.path(path, 'ModelInfo')
  figdir = file.path(outdir, 'Figure')
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
  files=shud.filein()
  # colnames(files)='shud_files'
  
  pm=readmesh(); pr=readriv(); pseg=readchannel()
  ia=getArea(pm = pm); 
  
  #======Basic Info ============
  iRivLen = rgeos::gLength(spr, byid=TRUE)
  iSegLen = pseg$Length
  RL = sum(iRivLen)
  RivEle = sort(unique(pseg$iEle) )
  AA  =sum(ia); 
  z=getElevation()
  da = cbind('Ele_Area'=summary(ia), 
             'Riv_Len'=summary(iRivLen), 
             'Ele_Elv'=summary(z), 
             'Seg_Len'=summary(iSegLen)
  )
  x=data.frame('Stats'=rownames(da), da)
  write.table(x, file = file.path(outdir, 'Basic_Sp_stat.txt'), 
              col.names = TRUE, row.names = FALSE)
  x=c('Ncell' = nrow(pm@mesh),
      'Nriv' = nrow(pr@river),
      'NrivType' = nrow(pr@rivertype),
      'NrivSeg' = nrow(pseg),
      'NwetEle' = length(RivEle),
      'Z_max_m' = max(z),
      'Z_min_m' = min(z),
      'Area_km2' = AA * 1E-6,
      'Area_mean_km2' = mean(ia) * 1E-6,
      # 'Area_min'= min(ia),
      # 'Area_max' = max(ia)
      'RivLength_km' = RL * 1E-3,
      'RivLen_min_m' = min(iRivLen),
      'RivSlope_mean' = mean(pr@river$Slope)
  )
  write.table(x, file = file.path(outdir, 'Basic_Info.txt'), col.names = FALSE)
  ret = x
  # ============================
  ma=MeshAtt()
  ma.st=apply(ma, 2, summary)
  ma.st
  write.table(t(ma.st), file = file.path(outdir, 'Basic_mesh_states.csv'), 
              col.names = TRUE, row.names = TRUE, append = FALSE)
  
  ################
  rl=getRiverNodes(spr=spr)
  ################
  centroid = getCentroid()
  triangle = pm@mesh[, 1:7]
  vertex = pm@point[, -4]
  
  if(!is.null(crs.pcs)) {
    ll = ProjectCoordinate(centroid[, 2:3], crs.pcs, P2G=TRUE)
    centroid = cbind(centroid, ll)
    
    ll = ProjectCoordinate(vertex[, 2:3], crs.pcs, P2G=TRUE)
    vertex = cbind(vertex, ll)
  }else{
  }
  xl = list('files' = files,
             'vertex' = vertex,
             'centroid' = centroid,
             'triangle' = triangle,
             'area' = ia,
            'riv.pts' = rl$points,
            'riv.ft' = rl$FT_ID
             )
  fns = c('Files.txt',
          'sp_Ele_vertex.csv',
          'sp_Ele_centroid.csv',
          'sp_Ele_triangle.csv',
          'sp_Ele_area.csv',
          'sp_Riv_points.csv',
          'sp_Riv_FrTo.csv'
          )
  
  for(i in 1:length(xl)){
    fn=file.path(outdir, fns[i])
    write.table(file = fn, xl[[i]], quote=F, row.names = FALSE, col.names = TRUE)
  }
  md = data.frame( triangle,
                  Area=ia,
                  Centroid = centroid
                  )
  
  #=======Area Histgram ===============
  md = data.frame(ia); colnames(md)='Area'
  nb=20
  p <- ggplot2::ggplot(md, ggplot2::aes(x=Area)) + 
    ggplot2::labs( x = "Area (m2)" ) +
    ggplot2::geom_histogram(bins=nb) +
    ggplot2::geom_freqpoly(bins=nb, col=2, lty=2)
  ggplot2::ggsave(file.path(figdir, paste0('Sp_AreaHist','.png')), p)
  
  myhist <- function(x, icols, pp){
    par(mfrow = pp)
    cn=colnames(x)
    for(i in 1:length(icols)){
      hist(x[, icols[i]], xlab=cn[ icols[i] ], main='')
    }
  }
  #======Soil============
  s=readsoil()
  png.control(fn = paste0('Para_soil','.png'), path=figdir, ht=6, wd=8, res=144) 
  myhist(s, icols=c(2, 3, 6,7,9), pp=c(2,3))
  dev.off()
  
  #======Geol============
  g=readgeol()
  png.control(fn = paste0('Para_geol','.png'), path=figdir, ht=6, wd=8, res=144) 
  myhist(g, icols=c(2,3,4,4,7,8), pp=c(2,3))
  dev.off()
  return(ret)
}
#' Summary of the attributions on Mesh
#' \code{MeshAtt} 
#' @param  pm SHUD.MESH, .sp.mesh
#' @param  att .sp.att of input
#' @param  soil .para.soil of input
#' @param  geol .para.geol of input
#' @param  lc .para.lc of input
#' @return Attributes of each element.
#' @export
MeshAtt<- function(pm=readmesh(), att=readatt(), 
                    soil=readsoil(), geol=readgeol(), lc=readlc() ){
  # pm=readmesh();att=readatt();
  # soil=readsoil(); geol=readgeol();
  # lc=readlc()
  y=data.frame('MESH'=pm@mesh,
             'ATT'=att,
             'SOIL'=soil[att$SOIL,],
             'GEOL'=geol[att$GEOL,],
             'LC'=lc[att$LC,]
             )
  return(y)
}
