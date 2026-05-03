#' River Processing Functions
#'
#' Functions for processing river networks including order calculation,
#' path analysis, and topology determination using modern sf library.
#'
#' @name river_processing
NULL

#' Calculate River Order (Strahler Classification)
#'
#' Calculates the Strahler stream order for a river network. This function
#' uses sf for spatial operations and implements the Strahler algorithm
#' for stream classification.
#'
#' @param rivers sf object with LINESTRING geometry representing river segments
#' @return Integer vector of stream orders for each river segment
#'
#' @details
#' The Strahler stream order is calculated as follows:
#' - First-order streams are headwater streams with no tributaries
#' - When two streams of order i join, they form a stream of order i+1
#' - When streams of different orders join, the resulting stream has the
#'   order of the higher-order stream
#'
#' @section Migration Note:
#' This function replaces deprecated \code{sp.RiverOrder()} and uses sf.
#' The algorithm and results are equivalent, but the input must be an sf object.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' # Load river network as sf object
#' rivers <- st_read("rivers.shp")
#' # Calculate stream order
#' order <- calc_river_order(rivers)
#' # Plot by order
#' plot(rivers["geometry"], col = order)
#' }
calc_river_order <- function(rivers) {
  msg <- "calc_river_order::"
  
  # Validate input
  if (!inherits(rivers, "sf")) {
    stop("'rivers' must be an sf object with LINESTRING geometry.\n",
         "To migrate from sp: rivers_sf <- sf::st_as_sf(rivers_sp)",
         call. = FALSE)
  }
  
  # Check geometry type
  geom_types <- unique(as.character(sf::st_geometry_type(rivers)))
  if (!all(geom_types %in% c("LINESTRING", "MULTILINESTRING"))) {
    stop("'rivers' must contain LINESTRING or MULTILINESTRING geometry",
         call. = FALSE)
  }
  
  # Convert MULTILINESTRING to LINESTRING if needed
  if (any(geom_types == "MULTILINESTRING")) {
    rivers <- suppressWarnings(sf::st_cast(rivers, "LINESTRING"))
  }
  
  # Extract coordinates
  coords <- get_coords(rivers)
  
  # Get from-to nodes
  ft <- get_from_to_nodes(rivers, coords)
  
  # Helper function to identify first-order streams
  get_first_order <- function(ft_matrix) {
    fr <- ft_matrix[, 2]
    to <- ft_matrix[, 3]
    
    # Count occurrences of each point
    tb <- table(c(fr, to))
    pid <- as.numeric(names(tb))
    ncount <- as.numeric(tb)
    
    # Outlet points (to points not in from)
    p_out <- to[!to %in% fr]
    
    # Junction points (appear more than twice)
    p_jnt <- pid[ncount > 2]
    
    # Key points (outlets and junctions)
    p_key <- c(p_out, p_jnt)
    
    # Target from-id (from points not in to - headwaters)
    tid <- fr[!fr %in% to]
    
    nd <- length(tid)
    ret <- NULL
    
    # Trace from each headwater to outlet or junction
    for (i in seq_len(nd)) {
      cr <- tid[i]
      for (j in seq_len(100000)) {
        rid <- which(fr == cr)
        ret <- c(ret, rid)
        sid <- to[rid]
        
        # Stop if reached key point
        if (any(sid %in% p_key)) {
          break
        } else {
          cr <- sid
        }
      }
    }
    
    ret
  }
  
  # Build from-to matrix
  nsp <- nrow(rivers)
  x <- ft
  y <- x
  x_ord <- rep(0, nsp)
  
  # Iteratively assign orders
  for (i in seq_len(10000)) {
    ids <- y[, 1]
    message(msg, "Order = ", i)
    
    id <- get_first_order(y)
    x_ord[ids[id]] <- i
    
    y <- y[-id, , drop = FALSE]
    ny <- nrow(y)
    
    if (ny <= 0) {
      break
    }
  }
  
  x_ord
}

#' Get river coordinates
#'
#' Extract unique coordinates from river network
#'
#' @param x \code{sf} object with LINESTRING geometry. Deprecated river
#'   wrappers may pass legacy line objects, which are converted internally for
#'   compatibility.
#' @return Matrix of unique coordinates
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' # Create a simple river network
#' line1 <- st_linestring(matrix(c(0,0, 1,1), ncol=2, byrow=TRUE))
#' line2 <- st_linestring(matrix(c(1,1, 2,0), ncol=2, byrow=TRUE))
#' rivers <- st_sfc(line1, line2)
#' 
#' # Extract unique coordinates
#' coords <- get_coords(rivers)
#' print(coords)
#' }
get_coords <- function(x) {
  if (inherits(x, "sf") || inherits(x, "sfc")) {
    # Extract all coordinates
    coords_list <- lapply(sf::st_geometry(x), function(line) {
      coords <- sf::st_coordinates(line)
      if (nrow(coords) == 0) {
        return(matrix(numeric(0), ncol = 2))
      }
      coords[, 1:2, drop = FALSE]
    })
    
    # Combine and get unique
    all_coords <- do.call(rbind, coords_list)
  } else {
    x_sf <- sf::st_as_sf(x)
    coords_list <- lapply(sf::st_geometry(x_sf), function(line) {
      coords <- sf::st_coordinates(line)
      if (nrow(coords) == 0) {
        return(matrix(numeric(0), ncol = 2))
      }
      coords[, 1:2, drop = FALSE]
    })
    all_coords <- do.call(rbind, coords_list)
  }
  
  unique_coords <- unique(all_coords)
  unique_coords
}

#' Get from-to nodes for river segments
#'
#' Identifies the start and end node indices for each river segment
#'
#' @param sf_line \code{sf} object with LINESTRING geometry. Deprecated river
#'   wrappers may pass legacy line objects, which are converted internally for
#'   compatibility.
#' @param coords Matrix of unique coordinates (default: get_coords(sf_line))
#' @return Matrix with columns: ID, FrNode, ToNode. Invalid line geometries
#'   with fewer than two coordinates, non-finite coordinates, or identical start
#'   and end coordinates return `NA` for `FrNode` and `ToNode`.
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' # Create a simple river network with 2 segments
#' line1 <- st_linestring(matrix(c(0,0, 1,1), ncol=2, byrow=TRUE))
#' line2 <- st_linestring(matrix(c(1,1, 2,0), ncol=2, byrow=TRUE))
#' rivers <- st_sfc(line1, line2)
#' 
#' # Get from-to node indices
#' ft_nodes <- get_from_to_nodes(rivers)
#' print(ft_nodes)
#' }
get_from_to_nodes <- function(sf_line, coords = get_coords(sf_line)) {
  if (inherits(sf_line, "sf") || inherits(sf_line, "sfc")) {
    nsp <- length(sf::st_geometry(sf_line))
    # Get all coordinates matrices
    coords_list <- lapply(sf::st_geometry(sf_line), function(line) {
      coords <- sf::st_coordinates(line)
      if (nrow(coords) == 0) {
        return(matrix(numeric(0), ncol = 2))
      }
      coords[, 1:2, drop = FALSE]
    })
  } else {
    sf_line <- sf::st_as_sf(sf_line)
    nsp <- nrow(sf_line)
    coords_list <- lapply(sf::st_geometry(sf_line), function(line) {
      coords <- sf::st_coordinates(line)
      if (nrow(coords) == 0) {
        return(matrix(numeric(0), ncol = 2))
      }
      coords[, 1:2, drop = FALSE]
    })
  }

  invalid_line <- vapply(coords_list, function(x) {
    nrow(x) < 2 || any(!is.finite(x)) ||
      identical(unname(x[1, ]), unname(x[nrow(x), ]))
  }, logical(1), USE.NAMES = FALSE)
  
  # Get first and last points of each line segment
  fr_pts <- do.call(rbind, lapply(coords_list, function(x) {
    if (nrow(x) < 2 || any(!is.finite(x)) ||
        identical(unname(x[1, ]), unname(x[nrow(x), ]))) {
      return(matrix(c(NA_real_, NA_real_), nrow = 1))
    }
    x[1, , drop = FALSE]
  }))
  to_pts <- do.call(rbind, lapply(coords_list, function(x) {
    if (nrow(x) < 2 || any(!is.finite(x)) ||
        identical(unname(x[1, ]), unname(x[nrow(x), ]))) {
      return(matrix(c(NA_real_, NA_real_), nrow = 1))
    }
    x[nrow(x), , drop = FALSE]
  }))
  
  # Use string hashing for fast exact matching
  fr_str <- paste(fr_pts[, 1], fr_pts[, 2])
  to_str <- paste(to_pts[, 1], to_pts[, 2])
  coord_str <- paste(coords[, 1], coords[, 2])
  
  # Vectorized match
  fr_idx <- match(fr_str, coord_str)
  to_idx <- match(to_str, coord_str)
  invalid_line <- invalid_line | (!is.na(fr_idx) & !is.na(to_idx) & fr_idx == to_idx)
  fr_idx[invalid_line] <- NA_integer_
  to_idx[invalid_line] <- NA_integer_
  
  # Build return matrix
  result <- cbind(seq_len(nsp), fr_idx, to_idx)
  colnames(result) <- c("ID", "FrNode", "ToNode")
  
  result
}

#' Calculate Downstream Relationships
#'
#' Determines the downstream segment index for each river segment in a network.
#' A segment's downstream is the segment whose start node matches this segment's
#' end node.
#'
#' @param sf_line sf object with LINESTRING geometry representing river segments
#' @param coords Matrix of unique coordinates (optional, will be calculated if NULL)
#' @return Integer vector of downstream indices (-1 for outlets, -3 for no downstream)
#'
#' @details
#' The function identifies downstream relationships by matching end nodes to
#' start nodes. Outlets (segments with no downstream) are marked with -1.
#' If multiple downstream segments are found (unusual), the first is selected
#' and a message is displayed.
#'
#' @section Migration Note:
#' This function replaces deprecated \code{sp.RiverDown()} and uses sf.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' rivers <- st_read("rivers.shp")
#' downstream <- calc_river_downstream(rivers)
#' # Identify outlets
#' outlets <- which(downstream < 0)
#' }
calc_river_downstream <- function(sf_line, coords = NULL) {
  msg <- "calc_river_downstream::"
  
  # Validate input
  if (!inherits(sf_line, "sf")) {
    stop("'sf_line' must be an sf object with LINESTRING geometry.\n",
         "To migrate from sp: sf_line_sf <- sf::st_as_sf(rivers_sp)",
         call. = FALSE)
  }
  
  # Check geometry type
  geom_types <- unique(as.character(sf::st_geometry_type(sf_line)))
  if (!all(geom_types %in% c("LINESTRING", "MULTILINESTRING"))) {
    stop("'sf_line' must contain LINESTRING or MULTILINESTRING geometry",
         call. = FALSE)
  }
  
  # Convert MULTILINESTRING to LINESTRING if needed
  if (any(geom_types == "MULTILINESTRING")) {
    sf_line <- suppressWarnings(sf::st_cast(sf_line, "LINESTRING"))
  }
  
  # Get coordinates if not provided
  if (is.null(coords)) {
    coords <- get_coords(sf_line)
  }
  
  # Get from-to nodes
  ft <- get_from_to_nodes(sf_line, coords)
  
  nsp <- nrow(sf_line) 
  idown <- rep(-3, nsp)
  
  # Find downstream for each segment
  for (i in seq_len(nsp)) {
    pto <- ft[i, 3]  # To node of current segment
    id <- which(ft[, 2] == pto)  # Segments starting at this node
    nid <- length(id)
    
    if (nid == 1) {
      idown[i] <- id
    } else if (nid > 1) {
      message(msg, i, "/", nsp, "\t", nid, " downstream segments found")
      message("Downstream IDs: ", paste(id, collapse = ", "))
      idown[i] <- id[1]
    } else {
      # No downstream - this is an outlet
      idown[i] <- -1
    }
  }
  
  idown
}

#' Calculate River Paths
#'
#' Dissolves river network into continuous paths from headwaters to outlets
#' or junctions. Each path represents a continuous flow line.
#'
#' @param rivers sf object with LINESTRING geometry representing river segments
#' @param coords Matrix of unique coordinates (optional, will be calculated if NULL)
#' @param downstream Integer vector of downstream indices (optional, will be calculated if NULL)
#' @return List with components:
#'   \item{seg_ids}{List of segment IDs for each path}
#'   \item{point_ids}{List of point IDs for each path}
#'   \item{paths}{sf object with dissolved paths as LINESTRING}
#'
#' @details
#' The function traces upstream from outlets and junctions to create continuous
#' flow paths. This is useful for simplifying river networks and identifying
#' main stems vs. tributaries.
#'
#' @section Migration Note:
#' This function replaces deprecated \code{sp.RiverPath()} and uses sf.
#' The output structure is similar but uses sf objects.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' rivers <- st_read("rivers.shp")
#' paths <- calc_river_path(rivers)
#' # Plot dissolved paths
#' plot(paths$paths)
#' }
calc_river_path <- function(rivers, coords = NULL, downstream = NULL) {
  msg <- "calc_river_path::"
  
  # Validate input
  if (!inherits(rivers, "sf")) {
    stop("'rivers' must be an sf object with LINESTRING geometry.\n",
         "To migrate from sp: rivers_sf <- sf::st_as_sf(rivers_sp)",
         call. = FALSE)
  }
  
  # Check geometry type
  geom_types <- unique(as.character(sf::st_geometry_type(rivers)))
  if (!all(geom_types %in% c("LINESTRING", "MULTILINESTRING"))) {
    stop("'rivers' must contain LINESTRING or MULTILINESTRING geometry",
         call. = FALSE)
  }
  
  # Convert MULTILINESTRING to LINESTRING if needed
  if (any(geom_types == "MULTILINESTRING")) {
    rivers <- sf::st_cast(rivers, "LINESTRING")
  }
  
  # Get coordinates if not provided
  if (is.null(coords)) {
    coords <- get_coords(rivers)
  }
  
  # Get from-to nodes
  ft <- get_from_to_nodes(rivers, coords)
  
  # Calculate downstream if not provided
  if (is.null(downstream)) {
    message(msg, "Calculating downstream relationships...")
    downstream <- calc_river_downstream(rivers, coords = coords)
  }
  
  # Get node IDs for each segment
  node_ids_list <- lapply(sf::st_geometry(rivers), function(line) {
    line_coords <- sf::st_coordinates(line)[, 1:2, drop = FALSE]
    
    # Find indices in unique coords for all points
    ids <- apply(line_coords, 1, function(pt) {
      which(coords[, 1] == pt[1] & coords[, 2] == pt[2])[1]
    })
    
    ids
  })
  
  # Identify key segments (outlets and segments ending at junctions)
  # Count frequency of each point
  all_points <- unlist(node_ids_list)
  p_frq <- table(all_points)
  
  # Points appearing more than twice are junctions
  p_jnt <- as.numeric(names(p_frq[p_frq > 2]))
  
  # Key segments: outlets or segments ending at junctions
  riv_keyid <- c(which(downstream < 0), which(ft[, 3] %in% p_jnt))
  riv_keyid <- unique(riv_keyid)
  
  # Function to trace upstream
  go_up <- function(id0) {
    uid <- which(downstream == id0)
    if (length(uid) == 1) {
      c(go_up(uid), id0)
    } else {
      id0
    }
  }
  
  # Trace upstream from each key segment
  message(msg, "Tracing upstream paths...")
  stream_paths <- lapply(riv_keyid, function(x) go_up(x))
  
  # Get point IDs for each path
  message(msg, "Extracting point sequences...")
  nstr <- length(stream_paths)
  point_lists <- lapply(seq_len(nstr), function(i) {
    sid <- stream_paths[[i]]
    unique(unlist(node_ids_list[sid]))
  })
  
  # Build sf object with dissolved paths
  message(msg, "Building dissolved paths...")
  path_geoms <- lapply(seq_len(nstr), function(i) {
    pid <- point_lists[[i]]
    sf::st_linestring(coords[pid, , drop = FALSE])
  })
  
  paths_sf <- sf::st_sf(
    path_id = seq_len(nstr),
    n_segments = sapply(stream_paths, length),
    geometry = sf::st_sfc(path_geoms, crs = sf::st_crs(rivers))
  )
  
  list(
    seg_ids = stream_paths,
    point_ids = point_lists,
    paths = paths_sf
  )
}
