#' Generate the mesh data from the triangulation
#' \code{shud.mesh}
#' @param tri Triangles defination
#' @param dem Elevation. Projection of the DEM raster must be same as watershed boundary.
#' @param AqDepth Aquifer depth, numeric.
#' @param r.aq  Aquifer Thickness. Raster object
#' @return Triangle mesh of model domain.
#' @noRd
#' @examples
#' data(sac)
#' wbd = sac[['wbd']]
#' dem = sac[['dem']]
#' a.max = 1e6 * 1;
#' q.min = 33;
#' tol.riv = 200
#' tol.wb = 200
#' tol.len = 500
#' AqDepth = 10
#'
#' wbbuf = sf::st_buffer(sf::st_as_sf(wbd), dist = 2000)
#' dem = terra::crop(dem, terra::vect(wbbuf))
#'
#' wb.dis = sf::st_union(sf::st_as_sf(wbd))
#' wb.simp = sf::st_simplify(wb.dis, dTolerance = tol.wb, preserveTopology = TRUE)
#'
#'
#' tri = shud.triangle(wb = wb.simp, q = q.min, a = a.max)
#' plot(tri, asp = 1)
#'
#' # generate SHUD .mesh
#' pm = shud.mesh(tri, dem = dem, AqDepth = AqDepth)
#' sm = mesh_to_sf(pm)
#' plot(sm)
shud.mesh_legacy <- function(tri, dem, AqDepth = 10, r.aq = dem * 0 + AqDepth){
  shud.mesh(tri = tri, dem = dem, AqDepth = AqDepth, r.aq = r.aq)
}

#' Centroids of the triangulation
#' \code{Tri2Centroid}
#' @param tri Triangles defination
#' @return Centroids of the triangles, m x 2;
#' @noRd
Tri2Centroid_legacy <- function(tri){
  Tri2Centroid(tri)
}

#' extract Coordinates of SpatialLines or  SpatialPolygons
#' \code{shud.att}
#' @param tri Triangles defination
#' @param r.soil raster of soil classification
#' @param r.geol raster of geology layer
#' @param r.lc raster of land cover, LAI and Roughness Length
#' @param r.forc raster of forcing data
#' @param r.mf raster of melt factor
#' @param r.BC raster of boundary condition
#' @param r.SS raster of Source/Sink
#' @param sp.lake SpatialPolygon of lake layers.
#' @return data.frame of SHUD .att
#' @noRd
shud.att_legacy <- function(tri, r.soil = NULL, r.geol = NULL, r.lc = NULL, r.forc = NULL,
                     r.mf = NULL, r.BC = NULL, r.SS = NULL, sp.lake = NULL){
  shud.att(
    tri = tri,
    r.soil = r.soil,
    r.geol = r.geol,
    r.lc = r.lc,
    r.forc = r.forc,
    r.mf = r.mf,
    r.BC = r.BC,
    r.SS = r.SS,
    sp.lake = sp.lake
  )
}

# 
# pa=shud.att(tri, 
#             r.soil = r.soil, r.geol = r.geol, 
#             r.lc = rlc.idx, 
#             r.forc = sp.forc, 
#             r.BC = 0, 
#             sp.lake = sp.lake)
