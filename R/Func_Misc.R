#' Find the centroid of the triangles
#' \code{toCentroid}
#' @param tri Triangle, N x 3 dimension
#' @param point Coordinates.
#' @return Coordinates of triangles centroids
#' @export
toCentroid <- function(tri, point){
  x=cbind(point[tri[,1],1],
          point[tri[,2],1],
          point[tri[,3],1])
  y=cbind(point[tri[,1],2],
          point[tri[,2],2],
          point[tri[,3],2])
  xc=rowMeans(x)
  yc=rowMeans(y)
  ret <- cbind('x'=xc,'y'=yc)
}

#' Eudilic Distance
#' \code{Eudist}
#' @param p1 coordinates of point 1
#' @param p2 coordinates of point 2
#' @return Distance between \code{p1} and \code{p2}
#' @export
#' @examples
#' p1 = c(1,1)
#' p2 = c(2,2)
#' Eudist(p1,p2)
Eudist <- function(p1,p2){
  d = sqrt( sum( (p1-p2) ^2) )
  d
}

#' extract Coordinates of SpatialLines or  SpatialPolygons
#' \code{extractCoords}
#' @param x SpatialLines or SpatialPolygons
#' @param unique Whether return unique coordinates or duplicated coordinates
#' @param aslist Whether return as list of coordinates
#' @return coordinates of points, m x 2.
#' @export
extractCoords<-function(x, unique=TRUE, aslist = FALSE){
  spl <- methods::as(x, "SpatialLines")
  spp <- methods::as(spl, "SpatialPoints")
  pts = sp::coordinates(spp)
  if (unique){
    ret = unique(pts)
  }else{
    ret = pts
  }
  return(ret)
}

#' Convert the \code{xy} (X,Y) into Index in \code{coords}
#' \code{xy2ID}
#' @param xy coordinate (x,y), subset of \code{coords}
#' @param coord coordinate set
#' @return Index in coodinate set
#' @export
#' @examples
#' coord = matrix(rnorm(10), 5,2)
#' xy = coord[c(5,1,3),]
#' xy2ID(xy, coord)
xy2ID <- function(xy, coord){
  if(is.matrix(xy) | is.data.frame(xy)){
  }else{
    xy = matrix(xy, ncol=2)
  }
  ng=nrow(xy)
  id=rep(0,ng)
  for(i in 1:ng){
    dd = which(rowMatch(xy[i,], coord))
    # print(dd)
    if(length(dd) > 0){
      id[i]= dd
    }
  }
  id
}

#' Count the frequency of numbers
#' \code{count}
#' @param x target value
#' @param n specific frequency
#' @return count of specific frequency
#' @export
#' @examples
#' x=c(round(rnorm(100, 2)+1), 0,0)
#' count(x)
#' count(x, sum(x==0))
count <- function(x, n=NULL){
  if(is.null(n)){
    rt = table(x)
  }else{
    ct = table(x)
    rt = ct[which(ct %in% n)]
  }
  rt
}

#' Find the outliers in a dataset
#' \code{which_outliers}
#' @param x dataset
#' @param na.rm Whether ignore the na value.
#' @param probs Range of non-outlier range
#' @param ... More options in quatile();
#' @return Index of the outliers
#' @export
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' x <- c(-10, x, 10)
#' id <- which_outliers(x, probs=c(0.05,0.95))
#' y = x
#' y[id]=NA
#' par(mfrow = c(1, 2))
#' boxplot(x, ylim=range(x))
#' boxplot(y, ylim=range(x))
#' par(mfrow=c(1,1))
which_outliers <- function(x, na.rm = TRUE, probs=c(.5, .95),...) {
  qnt <- stats::quantile(x, probs=probs, na.rm = na.rm, ...)
  H <- 1.5 * stats::IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  id=which(is.na(y))
  id
}

#' Convert SpatialLines or SpatialPolygons a Planar Straight Line Graph object
#' \code{sp2PSLG}
#' @param sp SpatialLines or SpatialPolygons
#' @return pslg class
#' @importFrom methods as is
#' @export
#' @examples
#' library(rgeos)
#' sl = readWKT("MULTIPOLYGON(((1 1,5 1,5 5,1 5,1 1),(2 2,2 3,3 3,3 2,2 2)),((6 3,9 2,9 4,6 3)))")
#' x = sp2PSLG(sl)
sp2PSLG<-function(sp){
  if(methods::is(sp, 'sf')){
    sp=methods::as(sp, 'Spatial')
  }
  sl = methods::as(sp, 'SpatialLines')
  nsl = length(sl)
  P = extractCoords(sl, unique = TRUE)
  S=NULL
  for(i in 1:nsl){
    cc = extractCoords(sl[i,], unique = FALSE)
    e = xy2ID(cc,P)
    S = rbind(S, cbind(e[-length(e)], e[-1]))
  }
  ret <- RTriangle::pslg(P=P,S=unique(S) )
}

#' Find the common coordinates in two groups of coordinates.
#' @param x,y Matrix of coodinates.
#' @return Matrix (m*2). Columns are index of x rows and y rows.
#' @export
CommonPoints <- function(x, y){
  x = rbind(x)
  y = rbind(y)
  k=0
  nx = nrow(x)
  ny = nrow(y)
  mx = ncol(x)
  my = ncol(y)
  ret = matrix(0, nrow=nx, ncol=2)
  if(mx!=my){
    return(ret);
  }
  for(ix in 1:nx){
    d1 = x[ix, 1] - y[, 1] ;
    d2 = x[ix, 2] - y[, 2] ;
    rowy = which(d1 ==0 & d2 ==0)
    if(length(rowy) > 0){
      k = k+1
      ret[k, 2] = rowy
      ret[k, 1] = ix
    }
    # ret[ix] = any( (x[ix, 1] - y[, 1] == 0) & (x[ix, 2] - y[, 2] == 0) )
  }
  if(k>0){
    ret = rbind(ret[1:k, ])
  }else{
    ret = NULL
  }
  return(ret)
}
