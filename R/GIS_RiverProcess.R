#' Deprecated: determine downstream index for legacy river lines
#' \code{sp.RiverDown}
#' @param sp Legacy \code{SpatialLines*} or \code{sf} LINESTRING object.
#' @param coord coordinates of \code{sp}.
#' @return Index of downstream for each segements
#' @keywords deprecated
#' @export
sp.RiverDown <- function(sp, coord = get_coords(sp)){
  msg = 'sp.RiverDown::'
  ft = rbind(get_from_to_nodes(sp, coords = coord)[, 2:3])
  nsp = if (inherits(sp, "sf")) nrow(sp) else length(sp)
  idown = rep(-3, nsp)
  for(i in 1:nsp){
    pto = ft[i,2]
    id = which(ft[,1] == pto)
    nid = length(id)
    if(nid == 1){
      idown[i] = id
    }else if(length(id) > 1){
      message(msg, i, '/', nsp, '\t', nid, ' Downstream found')
      message('Downstream IDs: ', paste(id, collapse=', '))
      idown[i] = id[1]
      # xid = c(id, i)
      # plot(sp[xid, ]); points(coord[ft[i, ], ], col=1:2);
      # points(coord[ft[id[1],], ], col=1:2, pch=2)
      # points(coord[ft[id[2],], ], col=1:2, pch=3)
    }else{
      # VOID
    }
  }
  idown
}
# sp.slt = spr
# rivdown = sp.RiverDown(sp.slt)
#' Deprecated: determine river path from downstream relationships
#' \code{sp.RiverPath}
#' @param sp Legacy \code{SpatialLines*} or \code{sf} LINESTRING object.
#' @param idown downstream index of \code{sp}.
#' @param coord coordinates of \code{sp}.
#' @param tol.simplify tolerance for simplification of river lines.
#' @return List of River Path(SegIDs, PointIDs, sp)
#' @keywords deprecated
#' @export
#' @examples
#' \dontrun{
#' # Deprecated compatibility example. Prefer calc_river_path() with sf input.
#' data(sh)
#' riv = sh$riv
#' x = sp.CutSptialLines(riv, tol = 200)
#' xx = sp.RiverPath(x)
#' names(xx)
#' }
sp.RiverPath <- function(sp, 
                         coord = get_coords(sp), 
                         idown = sp.RiverDown(sp, coord = coord),
                         tol.simplify = 0){
  msg = 'sp.RiverPath::'
  sp_sf <- if (inherits(sp, "sf")) sp else sf::st_as_sf(sp)
  if(tol.simplify > 0){
    sp_sf <- sf::st_simplify(sp_sf, dTolerance = tol.simplify, preserveTopology = TRUE)
  }
  pt.list <- lapply(sf::st_geometry(sp_sf), function(line) {
    sf::st_coordinates(line)[, 1:2, drop = FALSE]
  })
  id.list <- lapply(pt.list, function(x) { xy2ID(xy = x, coord = coord)})

  ft = get_from_to_nodes(sp_sf, coords = coord)[, 2:3]
  nsp = nrow(sp_sf)
  p.frq = data.frame(table(unlist(id.list)) )
  p0 = ft[!(ft[, 1] %in% ft[,2]), 1]
  tmp =  p.frq[p.frq[,2]>2, 1] # more than twice.
  p1 = as.numeric(levels(tmp))[tmp]
  p2 = ft[!(ft[, 2] %in% ft[,1]), 2] # TO that not in FROM
  
  # jid = which(table(as.matrix(ft)) > 2) #id of joint points
  riv.keyid = c(which(idown<0), which(ft[,2] %in% p1) )
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
  message(msg, 'to Upstream')
  StreamPath = lapply(riv.keyid, function(x) goUp(updown, x))
  #========Point ID==========
  nstr = length(StreamPath)
  pl = list()
  message(msg, 'searching from/to nodes')
  for(i in 1:nstr){
    sid = StreamPath[[i]]
    pid = unique(unlist(lapply(1:length(sid), function(x) id.list[[sid[x]]])))
    # print(sid)
    # print(pid)
    pl[[i]] = pid
  }
  #=========Spatial Lines===============
  ll = vector("list", nstr)
  message(msg, 'build new spatialLines')
  for(i in 1:nstr){
    sid = StreamPath[[i]]
    pid = unique(unlist(lapply(1:length(sid), function(x) id.list[[sid[x]]])))
    ll[[i]] = sf::st_linestring(as.matrix(coord[pid, ]))
  }
  spx_sf <- sf::st_sf(
    ID = seq_len(length(ll)),
    geometry = sf::st_sfc(ll, crs = sf::st_crs(sp_sf))
  )
  ret <- list(SegIDs = StreamPath,
              PointIDs= pl,
              sp = methods::as(spx_sf, "Spatial"))
}


#' Deprecated: calculate river order for legacy river lines
#' \code{sp.RiverOrder}
#' @param sp Legacy \code{SpatialLines*} or \code{sf} LINESTRING object.
#' @param coord coordinates of \code{sp}.
#' @return Stream order vector.
#' @keywords deprecated
#' @export
#' @examples 
#' \dontrun{
#' if (requireNamespace("sp", quietly = TRUE)) {
#'   library(sp)
#'   data(sac)
#'   riv=sac$riv
#'   ord = sp.RiverOrder(riv)
#'   print(ord)
#'   idx = sort(unique(ord))
#'   plot(riv, col=ord)
#'   legend('topleft', legend=idx, col=idx, lwd=1, lty=1)
#' }
#' }
sp.RiverOrder <- function(sp, coord = get_coords(sp)){
  msg='sp.RiverOrder::'
  nsp <- if (inherits(sp, "sf")) nrow(sp) else length(sp)
  get1st <- function(x){
    fr = x[,2]
    to = x[,3]
    tb= table(x[,2:3])
    pid = as.numeric(names(tb))
    ncount = as.numeric(tb)
    p.out = to[which(! to %in% fr)]
    p.jnt = pid[which(ncount> 2)]
    p.key = c(p.out, p.jnt)
    # print(p.key)
    tid=fr[ !fr %in% to] #target from-id
    # plot(riv.simp)
    # points(coord[ ])
    # points(coord[p.key,], col=2)
    nd = length(tid)
    ret = NULL
    # message(msg, 'Number of Streams = ', nd)
    for(i in 1:nd){
      cr = tid[i]
      for(j in 1:100000){
        # message(msg, 'Stream ',i,'\tSeg ', j, '\tFrom Node ', cr)
        rid = which(fr %in% cr)
        ret=c(ret, rid)
        sid= to[rid] #ID of to-point;
        # if(length(sid) < 1 | length(sid)>1){
        #   print(sid)
        # }
        if( any(sid  %in% p.key) ){ #if sid IS IN key-points(outlets or joint point).
          break;
        }else{
          cr = sid;
        }
      }
    }
    ret;
  }

  ft0 = get_from_to_nodes(sp, coords = coord)
  ft = rbind(unique(ft0[, 2:3]))
  if(nsp != nrow(ft)){
    message(msg, 'ERROR: duplicated river reches extis.')
    stop('STOP WITH ERROR')
  }
  x = cbind(seq_len(nsp), ft)
  y =x
  x.ord =x[,1]*0
  for(i in 1:100000){
    ids = y[,1]
    message(msg, 'Order = ', i)
    id = get1st(y)
    x.ord[ ids[id] ] = i
    y=rbind(y[-1 * id ,])
    ny = length(y)
    # message(msg, 'Left Segments =', ny)
    if(ny<=0){ break }
  }
  x.ord
}
#' Return node indices for line geometry
#' \code{NodeIDList}
#' @param sp \code{sf} LINESTRING object; legacy \code{SpatialLines*} is accepted
#'   for compatibility.
#' @param coord Coordinate of vertex in \code{sp}
#' @return list of node indices for each line feature.
#' @export
NodeIDList <- function(sp, 
            coord = get_coords(sp) ){
  force(coord)
  sp_sf <- if (inherits(sp, "sf")) sp else sf::st_as_sf(sp)
  pt.list <- lapply(sf::st_geometry(sp_sf), function(line) {
    sf::st_coordinates(line)[, 1:2, drop = FALSE]
  })
  id.list <- lapply(pt.list, function(x) { xy2ID(x, coord = coord)})
}


#' parameters for river types
#' \code{RiverType}
#' @param n number of types
#' @param width Width of rivers
#' @return data.frame of parameters
#' @export
RiverType <- function(n, width = 2 * (1:n) ){
  cn = c('Index', 'Depth', 'BankSlope',
         'Width', 'Sinuosity', 'Manning',
         'Cwr', 'KsatH', 'BedThick')
  nc = length(cn)
  rtype = cbind(1:n,
                5.0 + 0.5 * 1:n, #Depth
                0, #bankslope
                width, #Width
                1.0, #Sinuosity
                0.04, # manning's n, s/m^(1/3)
                # 4.63e-07, #manning's n, day/m^(1/3)
                0.6, #CWR
                0.1, #KsatH
                .1 # Bed Sediment Thickness.
  )
  colnames(rtype) = cn
  rtype
}
#' Cut rivers by mesh and build river--mesh segment layer
#'
#' Intersects reach lines with mesh polygons, then attaches the tabular fields
#' expected for SHUD \code{.rivseg} output: \code{Index} (1\ldots{}n),
#' intersection attributes (e.g. \code{iRiv}, \code{iEle}), and \code{Length}
#' from \code{sf::st_length()}.
#'
#' @param sf_mesh Mesh polygons as \code{sf} (e.g. from \code{\link{mesh_to_sf}}).
#' @param sf_riv River reaches as \code{sf} line geometries.
#' @return \code{sf} line object with columns \code{Index}, original attributes,
#'   and \code{Length} (geometry retained for mapping and \code{st_write}).
#' @export
shud.rivseg <- function(sf_mesh, sf_riv) {
  if (!inherits(sf_mesh, "sf")) {
    stop(
      "'sf_mesh' must be an sf object (e.g. mesh_to_sf(); use sf::st_as_sf() to convert sp)",
      call. = FALSE
    )
  }
  if (!inherits(sf_riv, "sf")) {
    stop(
      "'sf_riv' must be an sf object",
      call. = FALSE
    )
  }
  mesh_sf <- sf_mesh
  riv_sf  <- sf_riv

  n_mesh <- nrow(mesh_sf)
  n_riv  <- nrow(riv_sf)
  if (!"iEle" %in% names(mesh_sf)) mesh_sf$iEle <- seq_len(n_mesh)
  if (!"iRiv" %in% names(riv_sf))  riv_sf$iRiv  <- seq_len(n_riv)

  crs_mesh <- sf::st_crs(mesh_sf)
  crs_riv  <- sf::st_crs(riv_sf)
  if (is.na(crs_riv) && !is.na(crs_mesh)) {
    riv_sf <- sf::st_set_crs(riv_sf, crs_mesh)
  } else if (!is.na(crs_mesh) && !is.na(crs_riv) && crs_mesh != crs_riv) {
    riv_sf <- sf::st_transform(riv_sf, crs_mesh)
  }

  # Only repair features that are actually invalid (avoid full-scan overhead)
  inv_mesh <- which(!sf::st_is_valid(mesh_sf))
  if (length(inv_mesh) > 0L)
    mesh_sf[inv_mesh, ] <- suppressWarnings(sf::st_make_valid(mesh_sf[inv_mesh, ]))
  inv_riv <- which(!sf::st_is_valid(riv_sf))
  if (length(inv_riv) > 0L)
    riv_sf[inv_riv, ] <- suppressWarnings(sf::st_make_valid(riv_sf[inv_riv, ]))

  # Fast path: geos direct GEOS bindings (sparse pre-filter + C-level intersection)
  # Fallback: sf::st_intersection (always available)
  if (requireNamespace("geos", quietly = TRUE)) {
    seg <- .rivseg_geos(riv_sf, mesh_sf)
  } else {
    seg <- tryCatch(
      sf::st_intersection(riv_sf, mesh_sf),
      error = function(e) {
        stop(
          "shud.rivseg: st_intersection failed: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

  seg <- seg[!sf::st_is_empty(seg), , drop = FALSE]
  if (nrow(seg) == 0L) {
    warning("shud.rivseg: no features after intersection", call. = FALSE)
    return(seg)
  }

  gt  <- sf::st_geometry_type(seg)
  seg <- seg[gt %in% c(
    "LINESTRING", "MULTILINESTRING", "GEOMETRYCOLLECTION"
  ), , drop = FALSE]
  if (nrow(seg) == 0L) {
    warning("shud.rivseg: no line features after intersection", call. = FALSE)
    return(seg)
  }
  if (any(sf::st_geometry_type(seg) == "GEOMETRYCOLLECTION"))
    seg <- sf::st_collection_extract(seg, type = "LINESTRING", warn = FALSE)

  seg <- seg[!sf::st_is_empty(seg) & sf::st_dimension(seg) == 1L, , drop = FALSE]
  if (nrow(seg) == 0L) {
    warning("shud.rivseg: no line segments retained", call. = FALSE)
    return(seg)
  }

  n   <- nrow(seg)
  len <- as.numeric(sf::st_length(seg))
  d   <- sf::st_drop_geometry(seg)
  d$Index  <- NULL
  d$Length <- NULL
  tab <- cbind(
    data.frame(Index = seq_len(n), stringsAsFactors = FALSE),
    d,
    Length = len,
    stringsAsFactors = FALSE
  )
  sf::st_sf(tab, geometry = sf::st_geometry(seg), crs = sf::st_crs(seg))
}

# Internal helper: geos fast path for shud.rivseg
# Uses sf::st_intersects() for STRtree-based candidate filtering, then
# geos::geos_intersection() for C-level geometry computation without R-object overhead.
.rivseg_geos <- function(riv_sf, mesh_sf) {
  riv_attrs  <- sf::st_drop_geometry(riv_sf)
  mesh_attrs <- sf::st_drop_geometry(mesh_sf)

  # Sparse candidate pairs via STRtree (bbox pre-filter)
  candidates <- sf::st_intersects(riv_sf, mesh_sf)
  riv_idx    <- rep(seq_along(candidates), lengths(candidates))
  mesh_idx   <- unlist(candidates)

  empty_result <- function() {
    sf::st_sf(
      cbind(riv_attrs[integer(0), , drop = FALSE],
            mesh_attrs[integer(0), , drop = FALSE]),
      geometry = sf::st_sfc(crs = sf::st_crs(riv_sf))
    )
  }

  if (length(riv_idx) == 0L) return(empty_result())

  # Convert candidate geometries to geos once (avoid repeated conversion)
  riv_g  <- geos::as_geos_geometry(sf::st_geometry(riv_sf))
  mesh_g <- geos::as_geos_geometry(sf::st_geometry(mesh_sf))

  # Vectorised exact intersection at C level
  isect <- geos::geos_intersection(riv_g[riv_idx], mesh_g[mesh_idx])
  keep  <- !geos::geos_is_empty(isect)

  if (!any(keep)) return(empty_result())

  riv_idx  <- riv_idx[keep]
  mesh_idx <- mesh_idx[keep]
  isect    <- isect[keep]

  geom_sfc <- sf::st_as_sfc(isect)
  if (is.na(sf::st_crs(geom_sfc))) sf::st_crs(geom_sfc) <- sf::st_crs(riv_sf)

  sf::st_sf(
    cbind(riv_attrs[riv_idx,  , drop = FALSE],
          mesh_attrs[mesh_idx, , drop = FALSE]),
    geometry = geom_sfc
  )
}

#' Find shared points in line geometries
#' \code{SharedPoints}
#' @param sp \code{sf} LINESTRING object; legacy \code{SpatialLines} is accepted
#'   for compatibility.
#' @param plot Whether to plot the results for debugging
#' @return Topologic relationship between line components.
#' @export
SharedPoints <- function(sp, plot=FALSE){
  msg='SharedPoints::'
#  sp=riv
  if (inherits(sp, "sf") || inherits(sp, "sfc")) {
    pl = lapply(sf::st_geometry(sp), function(line) sf::st_coordinates(line)[, 1:2, drop = FALSE])
    nsp = nrow(sp)
  } else {
    sp_sf <- sf::st_as_sf(sp)
    pl = lapply(sf::st_geometry(sp_sf), function(line) sf::st_coordinates(line)[, 1:2, drop = FALSE])
    nsp = length(sp)
  }
  ret = matrix(0, nrow = nsp, ncol=nsp)
  for(i in 1:(nsp-1) ){
    for(j in (i+1):nsp ){
      if(plot) {
        plot(rbind(pl[[i]], pl[[j]]), asp=1); 
        points(pl[[i]], col=2)
        points(pl[[j]], col=3)
      }
      cp = CommonPoints(pl[[i]], pl[[j]]) 
      if( is.null(cp) ){
        
      }else{
        message(msg, i, '/', j)
        # print(ret)
        # print(cp)
        # readline("go")
        ret[i,j] = cp[,1]
        ret[j,i] = cp[,2]
      }
    }
  }
  ret
}
