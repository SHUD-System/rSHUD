# Legacy/internal compatibility file.
#
# These helpers are kept out of the public API and only forward older internal
# names to the modern implementations in R/MeshDomain.R. They remain in the
# package build so archived scripts that source these symbols internally can be
# migrated without changing runtime behavior.

#' Legacy internal wrapper for SHUD mesh generation
#' \code{shud.mesh_legacy}
#' @param tri Triangles defination
#' @param dem Elevation as a \code{terra::SpatRaster}; CRS must match the
#'   watershed boundary.
#' @param AqDepth Aquifer depth, numeric.
#' @param r.aq Aquifer thickness as a \code{terra::SpatRaster}
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

#' Legacy internal wrapper for triangulation centroids
#' \code{Tri2Centroid_legacy}
#' @param tri Triangles defination
#' @return Centroids of the triangles, m x 2;
#' @noRd
Tri2Centroid_legacy <- function(tri){
  Tri2Centroid(tri)
}

#' Legacy internal wrapper for mesh attribute extraction
#' \code{shud.att_legacy}
#' @param tri Triangles defination
#' @param r.soil Soil classification as a \code{terra::SpatRaster}
#' @param r.geol Geology layer as a \code{terra::SpatRaster}
#' @param r.lc Land cover, LAI, and roughness length as a
#'   \code{terra::SpatRaster}
#' @param r.forc Forcing data zones as a \code{terra::SpatRaster} or
#'   \code{sf} object
#' @param r.mf Melt factor as a \code{terra::SpatRaster}
#' @param r.BC Boundary condition as a \code{terra::SpatRaster} or numeric
#' @param r.SS Source/sink as a \code{terra::SpatRaster} or numeric
#' @param sp.lake Lake layers as an \code{sf} object.
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
