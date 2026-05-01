#' Low-level compatibility SHUD mesh domain builder
#'
#' Creates a SHUD.MESH S4 object from a triangulation, extracting elevation
#' and aquifer depth from \code{terra::SpatRaster} data. This is a low-level
#' compatibility/advanced API used by the automated model builder; most users
#' should call \code{\link{auto_build_model}}.
#'
#' @param tri Triangulation object from RTriangle or shud.triangle()
#' @param dem Elevation data as a \code{terra::SpatRaster}
#' @param AqDepth Numeric value for uniform aquifer depth (default: 10)
#' @param r.aq Aquifer thickness as a \code{terra::SpatRaster} (optional,
#'   overrides AqDepth if provided)
#'
#' @return SHUD.MESH S4 object with mesh and point data
#'
#' @details
#' This function has been modernized to use terra for raster operations,
#' providing better performance. It is retained for advanced manual workflows
#' and compatibility with older scripts. Deprecated low-level callers may still
#' pass legacy objects, which are converted for compatibility; new code should
#' pass \code{terra::SpatRaster} directly or use
#' \code{\link{auto_build_model}} for the full workflow.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' # Create boundary and generate mesh
#' boundary <- sf::st_as_sf(
#'   sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
#'                                          ncol=2, byrow=TRUE))))
#' )
#' tri <- shud.triangle(wb = boundary, q = 30, a = 2)
#'
#' # Create DEM
#' dem <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
#' terra::values(dem) <- runif(100, 100, 200)
#'
#' # Generate mesh domain
#' mesh <- shud.mesh(tri, dem = dem, AqDepth = 10)
#' }
shud.mesh <- function(tri, dem, AqDepth = 10, r.aq = dem * 0 + AqDepth) {
  
  triangulation <- tri
  
  # Compatibility conversion for older low-level callers.
  if (inherits(dem, "Raster")) {
    dem <- terra::rast(dem)
  }
  
  # Handle aquifer depth
  if (!missing(r.aq) && !is.null(r.aq)) {
    if (inherits(r.aq, "Raster")) {
      r.aq <- terra::rast(r.aq)
    }
    aquifer_depth <- r.aq
  } else {
    aquifer_depth <- AqDepth
  }
  
  # Validate inputs
  if (missing(triangulation)) {
    stop("Parameter 'tri' (triangulation) is required", call. = FALSE)
  }
  if (missing(dem)) {
    stop("Parameter 'dem' (elevation) is required", call. = FALSE)
  }
  
  if (!inherits(dem, "SpatRaster")) {
    stop(
      "Parameter 'dem' must be a terra SpatRaster object, but received ",
      class(dem)[1],
      call. = FALSE
    )
  }
  
  # Extract triangulation components
  if (!is.list(triangulation) || !all(c("T", "P", "NB") %in% names(triangulation))) {
    stop(
      "Parameter 'tri' must be a triangulation object with ",
      "T (triangles), P (points), and NB (neighbors) components",
      call. = FALSE
    )
  }
  
  topo <- triangulation$NB
  topo[topo < 0] <- 0
  
  pt <- triangulation$P
  pid <- triangulation$T
  
  # Calculate triangle centroids
  x <- (pt[pid[, 1], 1] + pt[pid[, 2], 1] + pt[pid[, 3], 1]) / 3
  y <- (pt[pid[, 1], 2] + pt[pid[, 2], 2] + pt[pid[, 3], 2]) / 3
  cxy <- cbind(x, y)

  # Extract elevation at centroids using terra
  zc <- .extract_spatraster_values(dem, cxy, method = "bilinear")
  
  # Create mesh data frame
  mesh_df <- data.frame(
    ID = seq_len(nrow(topo)),
    Node1 = pid[, 1],
    Node2 = pid[, 2],
    Node3 = pid[, 3],
    Nabr1 = topo[, 1],
    Nabr2 = topo[, 2],
    Nabr3 = topo[, 3],
    Zmax = zc
  )
  
  # Extract elevation at vertices
  zs <- .extract_spatraster_values(dem, pt, method = "bilinear")
  
  # Handle aquifer depth
  if (is.numeric(aquifer_depth) && length(aquifer_depth) == 1) {
    # Uniform aquifer depth
    aq <- rep(aquifer_depth, nrow(pt))
  } else if (inherits(aquifer_depth, "SpatRaster")) {
    # Spatially variable aquifer depth
    aq <- .extract_spatraster_values(aquifer_depth, pt, method = "bilinear")
  } else {
    stop(
      "Parameter 'r.aq' must be either a single numeric value or ",
      "a terra SpatRaster object",
      call. = FALSE
    )
  }
  
  # Create point data frame
  point_df <- data.frame(
    ID = seq_len(nrow(pt)),
    X = pt[, 1],
    Y = pt[, 2],
    AqDepth = aq,
    Elevation = zs
  )
  
  # Create SHUD.MESH S4 object
  mesh_obj <- SHUD.MESH(mesh = mesh_df, point = point_df)
  
  return(mesh_obj)
}


#' Low-level compatibility mesh attribute builder
#'
#' Extracts attribute values from \code{terra::SpatRaster} layers and
#' \code{sf} polygons to mesh triangle centroids. This is a low-level
#' compatibility/advanced API used by the automated model builder; most users
#' should call \code{\link{auto_build_model}}.
#'
#' @param tri Triangulation object from RTriangle or shud.triangle()
#' @param r.soil Soil classification as a \code{terra::SpatRaster} (optional)
#' @param r.geol Geology classification as a \code{terra::SpatRaster} (optional)
#' @param r.lc Land cover classification as a \code{terra::SpatRaster}
#'   (optional)
#' @param r.forc Forcing zones as a \code{terra::SpatRaster} or \code{sf}
#'   polygon/point object (optional)
#' @param r.mf Melt factor as a \code{terra::SpatRaster} (optional)
#' @param r.BC Boundary condition raster or numeric (optional)
#' @param r.SS Source/sink raster or numeric (optional)
#' @param sp.lake Lake polygons as an \code{sf} object (optional)
#'
#' @return data.frame with attribute indices for each mesh triangle
#'
#' @details
#' This function has been modernized to use terra for raster operations,
#' providing 2-5x better performance. It is retained for advanced manual
#' workflows and compatibility with older scripts. Deprecated low-level callers
#' may still pass legacy objects, which are converted for compatibility; new
#' code should pass \code{terra::SpatRaster} and \code{sf} objects directly or
#' use \code{\link{auto_build_model}} for the full workflow.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' # Generate mesh
#' boundary <- sf::st_as_sf(
#'   sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
#'                                          ncol=2, byrow=TRUE))))
#' )
#' tri <- shud.triangle(wb = boundary, q = 30, a = 2)
#'
#' # Create attribute rasters
#' soil <- terra::rast(ncol=10, nrow=10, xmin=0, xmax=10, ymin=0, ymax=10)
#' terra::values(soil) <- sample(1:5, 100, replace=TRUE)
#'
#' # Calculate attributes
#' att <- shud.att(tri, r.soil = soil)
#' }
shud.att <- function(tri, r.soil = NULL, r.geol = NULL, r.lc = NULL,
                    r.forc = NULL, r.mf = NULL, r.BC = NULL,
                    r.SS = NULL, sp.lake = NULL) {
  
  triangulation <- tri
  
  # Compatibility conversion for older low-level callers.
  convert_to_terra <- function(r) {
    if (is.null(r)) return(NULL)
    if (inherits(r, "Raster")) return(terra::rast(r))
    return(r)
  }
  
  soil <- convert_to_terra(r.soil)
  geology <- convert_to_terra(r.geol)
  landcover <- convert_to_terra(r.lc)
  forcing <- convert_to_terra(r.forc)
  melt_factor <- convert_to_terra(r.mf)
  boundary_condition <- convert_to_terra(r.BC)
  source_sink <- convert_to_terra(r.SS)
  
  # Compatibility conversion for older low-level callers.
  if (!is.null(sp.lake) && inherits(sp.lake, "Spatial")) {
    sp.lake <- sf::st_as_sf(sp.lake)
  }
  lake <- sp.lake
  
  # Validate triangulation
  if (!is.list(triangulation) || !all(c("T", "P") %in% names(triangulation))) {
    stop(
      "Parameter 'triangulation' must be a triangulation object with ",
      "T (triangles) and P (points) components",
      call. = FALSE
    )
  }
  
  # Calculate centroids
  centroids <- Tri2Centroid(triangulation)
  n_cells <- nrow(centroids)
  
  # Initialize attribute data frame
  att <- data.frame(
    INDEX = seq_len(n_cells),
    SOIL = 1L,
    GEOL = 1L,
    LC = 1L,
    FORC = 1L,
    MF = 1L,
    BC = 0L,
    SS = 0L,
    LAKE = 0L
  )
  
  # Helper function to extract raster values
  extract_raster_values <- function(raster_obj, coords, default_value) {
    if (is.null(raster_obj)) {
      return(rep(default_value, nrow(coords)))
    }
    
    if (is.numeric(raster_obj) && length(raster_obj) == 1) {
      return(rep(raster_obj, nrow(coords)))
    }
    
    if (!inherits(raster_obj, "SpatRaster")) {
      stop(
        "Raster parameter must be a terra SpatRaster object or numeric value",
        call. = FALSE
      )
    }
    
    values <- .extract_spatraster_values(raster_obj, coords, method = "simple")
    
    # Replace NA with default
    values[is.na(values)] <- default_value
    
    return(as.integer(values))
  }
  
  # Helper function to extract from sf polygons
  extract_sf_values <- function(sf_obj, coords, default_value) {
    if (is.null(sf_obj)) {
      return(rep(default_value, nrow(coords)))
    }
    
    if (!inherits(sf_obj, "sf")) {
      stop("Spatial polygon parameter must be an sf object", call. = FALSE)
    }
    
    # Convert coordinates to sf points
    points_sf <- sf::st_as_sf(
      data.frame(x = coords[, 1], y = coords[, 2]),
      coords = c("x", "y"),
      crs = sf::st_crs(sf_obj)
    )
    
    # Spatial join to find which polygon each point is in
    joined <- sf::st_join(points_sf, sf_obj, join = sf::st_within)
    
    # Extract ID or use default
    if ("ID" %in% colnames(joined)) {
      values <- joined$ID
    } else {
      # Use row number as ID
      values <- as.integer(rownames(joined))
    }
    
    values[is.na(values)] <- default_value
    
    return(as.integer(values))
  }
  
  # Extract attributes (vectorized operations)
  att$SOIL <- extract_raster_values(soil, centroids, 1L)
  att$GEOL <- extract_raster_values(geology, centroids, 1L)
  att$LC <- extract_raster_values(landcover, centroids, 1L)
  
  if (!is.null(r.forc) && (inherits(r.forc, "sf") || inherits(r.forc, "Spatial"))) {
    if (inherits(r.forc, "Spatial")) r.forc <- sf::st_as_sf(r.forc)
    att$FORC <- extract_sf_values(r.forc, centroids, 1L)
  } else {
    att$FORC <- extract_raster_values(forcing, centroids, 1L)
  }
  
  att$MF <- extract_raster_values(melt_factor, centroids, 1L)
  att$BC <- extract_raster_values(boundary_condition, centroids, 0L)
  att$SS <- extract_raster_values(source_sink, centroids, 0L)
  
  # Handle lake (sf polygon)
  if (!is.null(lake)) {
    att$LAKE <- extract_sf_values(lake, centroids, 0L)
  }
  
  return(att)
}


#' Calculate triangle centroids
#'
#' Computes the centroids of triangles in a triangulation.
#'
#' @param tri Triangulation object with T (triangles) and P (points)
#'
#' @return Matrix of centroid coordinates (n x 2)
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#'
#' boundary <- sf::st_as_sf(
#'   sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
#'                                          ncol=2, byrow=TRUE))))
#' )
#' tri <- shud.triangle(wb = boundary, q = 30, a = 2)
#' centroids <- Tri2Centroid(tri)
#' }
Tri2Centroid <- function(tri) {
  if (!is.list(tri) || !all(c("T", "P") %in% names(tri))) {
    stop(
      "Parameter 'tri' must have T (triangles) and P (points) components",
      call. = FALSE
    )
  }
  
  tt <- tri$T
  pt <- tri$P
  
  # Vectorized centroid calculation
  xc <- (pt[tt[, 1], 1] + pt[tt[, 2], 1] + pt[tt[, 3], 1]) / 3
  yc <- (pt[tt[, 1], 2] + pt[tt[, 2], 2] + pt[tt[, 3], 2]) / 3
  
  return(cbind(xc, yc))
}
