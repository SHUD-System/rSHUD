#' 
#' #' Smooth river slope with IDW algorithm
#' #' \code{getOutlets} 
#' #' @param pr SHUD.RIVER class
#' #' @return updated SHUD.RIVER
#' #' @export
#' smoothRiver <- function(pr){
#'   idw <- function(ll, v, p=1){
#'     w = (1/ll)^p
#'     rt = sum(v * w) / sum(w)
#'     rt
#'   }
#'   nflag = 1
#'   pfr = pr@river[,'FrNode']
#'   pto = pr@river[,'ToNode']
#'   idown = pr@river[,'Down']
#'   len = RiverLength(pr)
#'   oid = getOutlets(pr)
#'   
#'   np = length(pz)
#'   i=1
#'   
#'   pz = pr@point[,'Zmax']
#'   i = 89
#'   for(i in 1:np){
#'     s1 = which(pfr == i)
#'     s2 = which(pto == i)
#'     p1 = pto[s1 ]
#'     p2 = pfr[s2 ]
#'     
#'     s = c(s1,s2)
#'     p = c(p1,p2)
#'     if(length(s) > 1){
#'       ll = len[s]
#'       zz = pz[p]
#'       nz = idw(ll, zz)
#'       message(nz - pz[i])
#'       pz[i] = nz
#'     }
#'   }
#'   
#'   pr@point[,'Zmax'] = pz
#'   sl2 = RiverSlope(pr)
#'   sl1 = RiverSlope(ppr)
#'   # matplot(type='l', sl2)
#'   
#'   # ppz = ppr@point[,'Zmax'] 
#'   
#'   dz = pz - ppz
#'   dz[order(dz)]
#'   order(dz) 
#'   pr
#' }



