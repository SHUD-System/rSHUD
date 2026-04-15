#' Convert spatial object to PSLG (Planar Straight Line Graph)
#'
#' Internal helper function to convert sf or sp objects to PSLG format
#' for RTriangle triangulation. Automatically handles sp to sf conversion.
#'
#' @param sp_obj sf or sp object (POLYGON or LINESTRING)
#' @return list with P (points matrix) and S (segments matrix)
#' @keywords internal
sf_to_pslg <- function(sp_obj) {
  # Auto-convert sp to sf if needed
  if (inherits(sp_obj, "Spatial")) {
    sp_obj <- sf::st_as_sf(sp_obj)
  }
  
  sf_obj <- sp_obj
  if (!inherits(sf_obj, "sf")) {
    stop("Input must be an sf object", call. = FALSE)
  }
  
  # Extract coordinates from sf geometry
  coords_list <- sf::st_coordinates(sf_obj)
  
  # Get unique points
  unique_coords <- unique(coords_list[, c("X", "Y")])
  P <- as.matrix(unique_coords)

  
  # Build segments
  S <- NULL
  H <- NULL  # Holes
  
  # Process each feature
  geom_type <- as.character(sf::st_geometry_type(sf_obj, by_geometry = FALSE))
  
  if (geom_type %in% c("POLYGON", "MULTIPOLYGON")) {
    # For polygons, create segments connecting consecutive points
    for (i in seq_len(nrow(sf_obj))) {
      feature_coords <- sf::st_coordinates(sf_obj[i, ])
      
      # Get unique polygon IDs (L2 for MULTIPOLYGON, L1 for POLYGON)
      if ("L2" %in% colnames(feature_coords)) {
        poly_ids <- unique(feature_coords[, "L2"])
      } else {
        poly_ids <- unique(feature_coords[, "L1"])
      }
      
      for (poly_id in poly_ids) {
        if ("L2" %in% colnames(feature_coords)) {
          poly_coords <- feature_coords[feature_coords[, "L2"] == poly_id, c("X", "Y")]
          ring_ids <- feature_coords[feature_coords[, "L2"] == poly_id, "L1"]
        } else {
          poly_coords <- feature_coords[feature_coords[, "L1"] == poly_id, c("X", "Y")]
          ring_ids <- rep(1, nrow(poly_coords))
        }
        
        # First ring is exterior, others are holes
        unique_rings <- unique(ring_ids)
        for (ring_idx in seq_along(unique_rings)) {
          ring_id <- unique_rings[ring_idx]
          ring_coords <- poly_coords[ring_ids == ring_id, , drop = FALSE]
          
          # Remove duplicate last point if it equals first point
          if (nrow(ring_coords) > 1 && 
              all(ring_coords[1, ] == ring_coords[nrow(ring_coords), ])) {
            ring_coords <- ring_coords[-nrow(ring_coords), , drop = FALSE]
          }
          
          # Find indices in P
          indices <- apply(ring_coords, 1, function(pt) {
            which(P[, 1] == pt[1] & P[, 2] == pt[2])[1]
          })
          
          # Create segments
          n_pts <- length(indices)
          if (n_pts > 1) {
            ring_segments <- cbind(indices, c(indices[-1], indices[1]))
            S <- rbind(S, ring_segments)
            
            # If this is a hole (not the first ring), add hole point
            if (ring_idx > 1) {
              hole_pt <- colMeans(ring_coords)
              H <- rbind(H, hole_pt)
            }
          }
        }
      }
    }
  } else if (geom_type %in% c("LINESTRING", "MULTILINESTRING")) {
    # For lines, create segments connecting consecutive points
    for (i in seq_len(nrow(sf_obj))) {
      feature_coords <- sf::st_coordinates(sf_obj[i, ])
      
      # Get unique line IDs
      if ("L1" %in% colnames(feature_coords)) {
        line_ids <- unique(feature_coords[, "L1"])
      } else {
        line_ids <- 1
      }
      
      for (line_id in line_ids) {
        if ("L1" %in% colnames(feature_coords)) {
          line_coords <- feature_coords[feature_coords[, "L1"] == line_id, c("X", "Y")]
        } else {
          line_coords <- feature_coords[, c("X", "Y")]
        }
        
        # Find indices in P
        # We use a vectorized approach for performance
        indices <- apply(line_coords, 1, function(pt) {
          which(P[, 1] == pt[1] & P[, 2] == pt[2])[1]
        })
        
        # Create segments
        n_pts <- length(indices)
        if (n_pts > 1) {
          line_segments <- cbind(indices[-n_pts], indices[-1])
          S <- rbind(S, line_segments)
        }
      }
    }
  }
  
  # Remove duplicate segments
  if (!is.null(S)) {
    S <- unique(S)
  }
  
  list(P = P, S = S, H = H)
}


#' Generate triangular mesh domain
#'
#' Creates an unstructured triangular mesh using constrained Delaunay
#' triangulation. Modernized implementation using sf and terra internally,
#' but maintains backward compatibility with sp/raster objects.
#'
#' @param wb Watershed boundary - sf object with POLYGON geometry, or legacy
#'   sp object (SpatialPolygons, auto-converted to sf)
#' @param riv River network - sf object with LINESTRING geometry, or legacy
#'   sp object (SpatialLines, auto-converted to sf) (optional)
#' @param dem Elevation data - terra SpatRaster or legacy raster object
#'   (auto-converted to SpatRaster) (optional, for future use)
#' @param hole Holes/lakes - sf object with POLYGON geometry, or legacy
#'   sp object (auto-converted to sf) (optional)
#' @param pts Additional constraint points - matrix (n x 2) (optional)
#' @param q Minimum angle constraint for triangles (degrees, default 30)
#' @param ... Additional arguments passed to RTriangle::triangulate()
#'   (e.g., a = max_area)
#'
#' @return A triangulation object from RTriangle with components:
#'   \item{P}{Matrix of point coordinates (n x 2)}
#'   \item{T}{Matrix of triangle vertex indices (m x 3)}
#'   \item{E}{Matrix of edge indices}
#'   \item{NB}{Matrix of neighbor triangle indices}
#'
#' @details
#' This function has been modernized to use sf and terra internally for
#' better performance and maintainability. However, it automatically converts
#' legacy sp and raster objects, so existing code will continue to work.
#'
#' The function uses RTriangle for constrained Delaunay triangulation with
#' quality constraints. The minimum angle parameter (q) controls triangle
#' quality - higher values produce better-shaped triangles but may increase
#' computation time. Values above 35 degrees may cause issues.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' 
#' # Create simple boundary
#' boundary <- sf::st_as_sf(
#'   sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
#'                                          ncol=2, byrow=TRUE))))
#' )
#' 
#' # Generate mesh
#' tri <- shud.triangle(wb = boundary, q = 30, a = 2)
#' plot(tri, asp = 1)
#' 
#' # With rivers
#' rivers <- sf::st_as_sf(
#'   sf::st_sfc(sf::st_linestring(matrix(c(2,2, 8,8), ncol=2, byrow=TRUE)))
#' )
#' tri <- shud.triangle(wb = boundary, riv = rivers, q = 30)
#' }
shud.triangle <- function(wb, dem = NULL, riv = NULL, hole = NULL,
                         pts = NULL, q = 30, ...) {
  
  # Validate inputs
  if (missing(wb)) {
    stop("Parameter 'wb' (watershed boundary) is required", call. = FALSE)
  }
  
  # Auto-convert sp to sf
  if (inherits(wb, "Spatial")) {
    wb <- sf::st_as_sf(wb)
  }
  if (!is.null(riv) && inherits(riv, "Spatial")) {
    riv <- sf::st_as_sf(riv)
  }
  if (!is.null(hole) && inherits(hole, "Spatial")) {
    hole <- sf::st_as_sf(hole)
  }
  
  # Auto-convert raster to terra
  if (!is.null(dem) && inherits(dem, "Raster")) {
    dem <- terra::rast(dem)
  }
  
  # Validate wb is now sf
  if (!inherits(wb, "sf")) {
    stop(
      "Parameter 'wb' must be an sf or sp object, but received ",
      class(wb)[1],
      call. = FALSE
    )
  }
  
  # Check geometry type
  geom_type <- as.character(sf::st_geometry_type(wb, by_geometry = FALSE))
  if (!geom_type %in% c("POLYGON", "MULTIPOLYGON")) {
    stop(
      "Parameter 'wb' must have POLYGON or MULTIPOLYGON geometry, ",
      "but has ", geom_type,
      call. = FALSE
    )
  }
  
  # Validate q (min_angle)
  if (q <= 0 || q > 35) {
    if (q > 35) {
      warning(
        "q > 35 may cause triangulation issues, setting to 35",
        call. = FALSE
      )
      q <- 35
    } else {
      stop("q (minimum angle) must be positive", call. = FALSE)
    }
  }
  
  # Handle holes if provided
  if (!is.null(hole)) {
    if (!inherits(hole, "sf")) {
      stop("Parameter 'hole' must be an sf object", call. = FALSE)
    }
    # Subtract holes from boundary using sf
    wb <- sf::st_difference(wb, sf::st_union(hole))
  }
  
  # Convert boundary to PSLG
  pslg_boundary <- sf_to_pslg(wb)
  P <- pslg_boundary$P
  S <- pslg_boundary$S
  H <- pslg_boundary$H
  
  # Add rivers if provided
  if (!is.null(riv)) {
    if (!inherits(riv, "sf")) {
      stop("Parameter 'riv' must be an sf object", call. = FALSE)
    }
    
    pslg_rivers <- sf_to_pslg(riv)
    n_boundary_pts <- nrow(P)
    
    # Combine points and segments
    P <- rbind(P, pslg_rivers$P)
    S <- rbind(S, pslg_rivers$S + n_boundary_pts)
  }
  
  # Add extra constraint points if provided
  if (!is.null(pts)) {
    if (!is.matrix(pts) && !is.data.frame(pts)) {
      stop("Parameter 'pts' must be a matrix or data.frame", call. = FALSE)
    }
    if (ncol(pts) != 2) {
      stop("Parameter 'pts' must have 2 columns (x, y)", call. = FALSE)
    }
    P <- rbind(P, as.matrix(pts))
  }
  
  # Create PSLG object for RTriangle
  if (is.null(H)) {
    pslg <- RTriangle::pslg(P = P, S = S)
  } else {
    pslg <- RTriangle::pslg(P = P, S = S, H = H)
  }
  
  # Perform triangulation
  tri_args <- list(p = pslg, q = q)
  tri_args <- c(tri_args, list(...))
  
  tri <- do.call(RTriangle::triangulate, tri_args)
  
  return(tri)
}


#' Convert mesh triangulation to sf shapefile
#'
#' Converts a triangular mesh to an sf object with POLYGON geometry.
#' Modernized implementation using sf, with automatic sp conversion.
#'
#' @param pm SHUD.MESH S4 object, or triangulation object from RTriangle/shud.triangle
#' @param dbf data.frame of attributes to attach to polygons (optional)
#' @param crs Coordinate reference system (CRS) as EPSG code, proj4string, or WKT
#'
#' @return sf object with POLYGON geometry and attributes
#'
#' @details
#' This function has been modernized to return sf objects instead of sp
#' SpatialPolygonsDataFrame. The sf format is more efficient and integrates
#' better with modern R spatial packages.
#'
#' Area is automatically calculated using sf::st_area() and added to the
#' attributes table.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Create simple mesh
#' boundary <- sf::st_as_sf(
#'   sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
#'                                          ncol=2, byrow=TRUE))))
#' )
#' tri <- shud.triangle(wb = boundary, q = 30, a = 2)
#'
#' # Convert to sf
#' mesh_sf <- sp.mesh2Shape(pm = tri, crs = 4326)
#' plot(mesh_sf)
#' }
#'
#' @name sp.mesh2Shape
#' @rdname mesh_to_sf
NULL

#' Convert SHUD mesh to sf object
#'
#' Converts a SHUD mesh structure (from \code{\link{readmesh}} or \code{\link{shud.mesh}})
#' to an sf object with POLYGON geometry.
#'
#' @param pm A \code{SHUD.MESH} object or list with \code{T} and \code{P} elements.
#' @param dbf A \code{data.frame} of attributes to attach to polygons. Default is NULL.
#' @param crs Coordinate reference system (e.g., EPSG code or proj4string). Default is NULL.
#'
#' @return An \code{sf} object with POLYGON geometry.
#' @export
#' @aliases mesh_to_sf
mesh_to_sf <- function(pm = readmesh(), dbf = NULL, crs = NULL) {
  
  mesh <- pm
  attributes <- dbf
  
  # Handle different input types
  if (inherits(mesh, "SHUD.MESH") || methods::is(mesh, "Untructure Domain")) {
    # Extract from SHUD.MESH S4 object
    triangles <- as.matrix(mesh@mesh[, c("Node1", "Node2", "Node3")])
    points <- as.matrix(mesh@point[, c("X", "Y")])
    
    # If no attributes provided, use mesh attributes
    if (is.null(attributes)) {
      # Calculate average aquifer depth and elevation for each triangle
      aqd <- (mesh@point$AqDepth[triangles[, 1]] + 
              mesh@point$AqDepth[triangles[, 2]] + 
              mesh@point$AqDepth[triangles[, 3]]) / 3
      zs <- (mesh@point$Elevation[triangles[, 1]] + 
             mesh@point$Elevation[triangles[, 2]] + 
             mesh@point$Elevation[triangles[, 3]]) / 3
      
      attributes <- data.frame(
        ID = mesh@mesh$ID,
        AqDepth = aqd,
        Zsurf = zs
      )
    }
  } else if (is.list(mesh) && all(c("T", "P") %in% names(mesh))) {
    # RTriangle triangulation object
    triangles <- mesh$T
    points <- mesh$P
    
    if (is.null(attributes)) {
      attributes <- data.frame(ID = seq_len(nrow(triangles)))
    }
  } else {
    stop(
      "Parameter 'mesh' must be either a SHUD.MESH object or ",
      "a triangulation object with T (triangles) and P (points) components",
      call. = FALSE
    )
  }
  
  # Validate dimensions
  if (ncol(points) != 2) {
    stop("Points matrix must have 2 columns (x, y)", call. = FALSE)
  }
  if (ncol(triangles) != 3) {
    stop("Triangles matrix must have 3 columns (vertex indices)", call. = FALSE)
  }
  
  # Create polygon geometries for each triangle
  n_triangles <- nrow(triangles)
  polygon_list <- vector("list", n_triangles)
  
  for (i in seq_len(n_triangles)) {
    # Get vertex indices
    v1 <- triangles[i, 1]
    v2 <- triangles[i, 2]
    v3 <- triangles[i, 3]
    
    # Get coordinates and close the polygon
    coords <- rbind(
      points[v1, ],
      points[v2, ],
      points[v3, ],
      points[v1, ]  # Close the polygon
    )
    
    # Create polygon
    polygon_list[[i]] <- sf::st_polygon(list(coords))
  }
  
  # Create sf object
  geom <- sf::st_sfc(polygon_list)
  
  # Set CRS if provided
  if (!is.null(crs)) {
    geom <- sf::st_set_crs(geom, crs)
  }
  
  # Create sf data frame
  mesh_sf <- sf::st_sf(attributes, geometry = geom)
  
  # Calculate area
  mesh_sf$Area <- as.numeric(sf::st_area(mesh_sf))
  
  return(mesh_sf)
}

#' @rdname mesh_to_sf
#' @export
sp.mesh2Shape <- function(pm = readmesh(), dbf = NULL, crs = NULL) {
  warning("sp.mesh2Shape is deprecated. Please use mesh_to_sf instead.")
  mesh_to_sf(pm, dbf, crs)
}


#' Convert triangulation to shapefile
#'
#' Converts a triangular mesh to an sf object with POLYGON geometry.
#' Modernized implementation using sf.
#'
#' @param tri Triangulation object from RTriangle/shud.triangle
#' @param dbf data.frame of attributes to attach to polygons (optional)
#' @param crs Coordinate reference system (optional)
#'
#' @return sf object with POLYGON geometry
#' @export
sp.Tri2Shape <- function(tri, dbf = NULL, crs = NA) {
  if (is.na(crs)) {
    crs <- NULL
  }
  
  # Use sp.mesh2Shape implementation
  sp.mesh2Shape(pm = tri, dbf = dbf, crs = crs)
}
