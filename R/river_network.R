#' River Network Functions
#'
#' Functions for building and analyzing river networks, including
#' integrated property calculations.
#'
#' @name river_network
NULL

#' Calculate River Properties
#'
#' Integrated function to calculate multiple river properties including
#' stream order, downstream relationships, length, and slope. This function
#' provides a unified interface to avoid redundant calculations.
#'
#' @param rivers sf object with LINESTRING geometry representing river segments
#' @param dem SpatRaster object with elevation data (required for slope calculation)
#' @param properties Character vector specifying which properties to calculate.
#'   Options: "order", "downstream", "length", "slope", "all" (default: "all")
#' @return sf object with calculated properties added as columns
#'
#' @details
#' This function integrates multiple river property calculations to avoid
#' redundant processing. Properties are calculated efficiently by reusing
#' intermediate results (e.g., coordinates, from-to nodes).
#'
#' Available properties:
#' \itemize{
#'   \item \strong{order}: Strahler stream order
#'   \item \strong{downstream}: Index of downstream segment (-1 for outlets)
#'   \item \strong{length}: Length of each segment in map units
#'   \item \strong{slope}: Bed slope calculated from DEM
#' }
#'
#' @section Migration Note:
#' This function integrates functionality from multiple legacy functions:
#' \code{sp.RiverOrder()}, \code{sp.RiverDown()}, \code{RiverLength()},
#' and \code{RiverSlope()}. It uses sf and terra instead of sp and raster.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' 
#' # Load data
#' rivers <- st_read("rivers.shp")
#' dem <- rast("dem.tif")
#' 
#' # Calculate all properties
#' rivers_props <- calc_river_properties(rivers, dem)
#' 
#' # Calculate only specific properties
#' rivers_order <- calc_river_properties(rivers, properties = "order")
#' }
calc_river_properties <- function(rivers, dem = NULL,
                                  properties = "all") {
  msg <- "calc_river_properties::"
  
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
  
  # Expand "all" to all available properties
  if ("all" %in% properties) {
    properties <- c("order", "downstream", "length", "slope")
  }
  
  # Validate properties
  valid_props <- c("order", "downstream", "length", "slope")
  invalid <- setdiff(properties, valid_props)
  if (length(invalid) > 0) {
    stop("Invalid properties: ", paste(invalid, collapse = ", "), "\n",
         "Valid options: ", paste(valid_props, collapse = ", "),
         call. = FALSE)
  }
  
  # Check if DEM is required
  if ("slope" %in% properties && is.null(dem)) {
    stop("'dem' is required for slope calculation", call. = FALSE)
  }
  
  # Validate DEM if provided
  if (!is.null(dem)) {
    if (!inherits(dem, "SpatRaster")) {
      stop("'dem' must be a SpatRaster object.\n",
           "To migrate from raster: dem_terra <- terra::rast(dem_raster)",
           call. = FALSE)
    }
  }
  
  # Convert MULTILINESTRING to LINESTRING if needed
  if (any(geom_types == "MULTILINESTRING")) {
    rivers <- suppressWarnings(sf::st_cast(rivers, "LINESTRING"))
  }
  
  # Make a copy to avoid modifying original
  result <- rivers
  
  # Calculate shared intermediate results
  message(msg, "Extracting coordinates...")
  coords <- get_coords(rivers)
  
  message(msg, "Identifying from-to nodes...")
  ft <- get_from_to_nodes(rivers, coords)
  
  # Calculate downstream if needed for order or explicitly requested
  downstream <- NULL
  if ("order" %in% properties || "downstream" %in% properties) {
    message(msg, "Calculating downstream relationships...")
    downstream <- calc_river_downstream(rivers, coords = coords)
    if ("downstream" %in% properties) {
      result$downstream <- downstream
    }
  }
  
  # Calculate order if requested
  if ("order" %in% properties) {
    message(msg, "Calculating stream order...")
    order <- calc_river_order(rivers)
    result$order <- order
  }
  
  # Calculate length if requested
  if ("length" %in% properties) {
    message(msg, "Calculating segment lengths...")
    # Use sf to calculate length
    lengths <- as.numeric(sf::st_length(rivers))
    result$length <- lengths
  }
  
  # Calculate slope if requested
  if ("slope" %in% properties) {
    message(msg, "Calculating slopes from DEM...")
    
    # Extract elevations at from and to nodes
    from_coords <- coords[ft[, 2], , drop = FALSE]
    to_coords <- coords[ft[, 3], , drop = FALSE]
    
    # Convert to terra vect for extraction
    from_pts <- terra::vect(from_coords, crs = terra::crs(dem))
    to_pts <- terra::vect(to_coords, crs = terra::crs(dem))
    
    # Extract elevations
    z_from <- .extract_spatraster_values(dem, from_pts)
    z_to <- .extract_spatraster_values(dem, to_pts)
    
    # Calculate lengths if not already done
    if (!"length" %in% names(result)) {
      lengths <- as.numeric(sf::st_length(rivers))
    } else {
      lengths <- result$length
    }
    
    # Calculate slope
    slopes <- (z_from - z_to) / lengths
    
    # Handle negative or zero slopes (set to small positive value)
    slopes[slopes <= 0] <- 1e-5
    
    result$slope <- slopes
  }
  
  result
}

#' Generate River Type Parameters
#'
#' Creates a data frame of hydraulic parameters for different river types/orders.
#' This is used to assign width, depth, and other hydraulic properties based
#' on stream order or other classification.
#'
#' @param n Integer, number of river types
#' @param width Numeric vector of widths for each type (default: 2 * (1:n))
#' @param depth Numeric vector of depths (default: 5.0 + 0.5 * (1:n))
#' @param manning Numeric, Manning's roughness coefficient (default: 0.04)
#' @return data.frame with hydraulic parameters for each river type
#'
#' @details
#' The function generates default hydraulic parameters that scale with river
#' type/order. Users can customize width and depth, while other parameters
#' use reasonable defaults.
#'
#' Parameters included:
#' \itemize{
#'   \item Index: River type index
#'   \item Depth: Channel depth (m)
#'   \item BankSlope: Bank slope (default: 0)
#'   \item Width: Channel width (m)
#'   \item Sinuosity: Channel sinuosity (default: 1.0)
#'   \item Manning: Manning's n (s/m^(1/3))
#'   \item Cwr: Coefficient for width-discharge relationship (default: 0.6)
#'   \item KsatH: Horizontal hydraulic conductivity (default: 0.1)
#'   \item BedThick: Bed sediment thickness (default: 0.1)
#' }
#'
#' @section Migration Note:
#' This function replaces \code{RiverType()} with the same interface.
#'
#' @export
#' @examples
#' # Generate parameters for 5 river types
#' river_types <- generate_river_types(5)
#' 
#' # Custom widths based on watershed area
#' widths <- c(2, 5, 10, 20, 40)
#' river_types <- generate_river_types(5, width = widths)
generate_river_types <- function(n, width = 2 * (1:n),
                                 depth = 5.0 + 0.5 * (1:n),
                                 manning = 0.04) {
  # Validate inputs
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("'n' must be a positive integer", call. = FALSE)
  }
  
  n <- as.integer(n)
  
  if (length(width) == 1) {
    width <- rep(width, n)
  } else if (length(width) != n) {
    stop("'width' must have length 1 or n", call. = FALSE)
  }
  
  if (length(depth) == 1) {
    depth <- rep(depth, n)
  } else if (length(depth) != n) {
    stop("'depth' must have length 1 or n", call. = FALSE)
  }
  
  if (length(manning) == 1) {
    manning <- rep(manning, n)
  } else if (length(manning) != n) {
    stop("'manning' must have length 1 or n", call. = FALSE)
  }
  
  # Build data frame
  river_types <- data.frame(
    Index = 1:n,
    Depth = depth,
    BankSlope = 0,
    Width = width,
    Sinuosity = 1.0,
    Manning = manning,
    Cwr = 0.6,
    KsatH = 0.1,
    BedThick = 0.1
  )
  
  river_types
}

#' Calculate River Width from Watershed Area
#'
#' Estimates river width for different stream orders based on watershed area
#' using an empirical relationship.
#'
#' @param area Numeric, watershed area (in appropriate units)
#' @param n_types Integer, number of river types/orders
#' @return Numeric vector of widths for each type
#'
#' @details
#' Uses an empirical power-law relationship to estimate channel width
#' based on drainage area. The relationship accounts for decreasing width
#' with increasing stream order (higher order = larger rivers).
#'
#' @export
#' @examples
#' # Calculate widths for 5 stream orders in a 1000 km² watershed
#' widths <- calc_river_width_from_area(1000, 5)
calc_river_width_from_area <- function(area, n_types = 10) {
  check_positive(area, "area")
  
  if (!is.numeric(n_types) || length(n_types) != 1 || n_types < 1) {
    stop("'n_types' must be a positive integer", call. = FALSE)
  }
  
  n_types <- as.integer(n_types)
  
  # Empirical relationship
  dd <- (1 / (1:n_types)) ^ 0.8
  widths <- rev(8 * log10(area + 1) * area ^ 0.25 * dd)
  
  round(widths, 2)
}

#' Build River Network
#'
#' Constructs a complete river network with all properties calculated,
#' including stream order, downstream relationships, lengths, slopes,
#' and hydraulic parameters. This is the main function for river network
#' preparation in SHUD models.
#'
#' @param rivers sf object with LINESTRING geometry representing river segments
#' @param dem SpatRaster object with elevation data
#' @param area Numeric, watershed area for width estimation (optional)
#' @param river_order Integer vector of stream orders (optional, will be calculated if NULL)
#' @param downstream Integer vector of downstream indices (optional, will be calculated if NULL)
#' @return List with components:
#'   \item{network}{sf object with all river properties}
#'   \item{river_types}{data.frame with hydraulic parameters for each type}
#'   \item{points}{data.frame with from/to node coordinates and elevations}
#'
#' @details
#' This function integrates all river processing steps:
#' \enumerate{
#'   \item Calculate stream order (Strahler classification)
#'   \item Identify downstream relationships
#'   \item Calculate segment lengths
#'   \item Extract elevations and calculate slopes
#'   \item Generate hydraulic parameters for each river type
#' }
#'
#' The output is ready for use in SHUD model input files.
#'
#' @section Migration Note:
#' This function replaces \code{shud.river()} and uses sf/terra instead of
#' legacy spatial packages.
#' The output structure is similar but uses modern spatial formats.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' 
#' # Load data
#' rivers <- st_read("rivers.shp")
#' dem <- rast("dem.tif")
#' 
#' # Build complete river network
#' river_network <- build_river_network(rivers, dem, area = 1000)
#' 
#' # Access components
#' network_sf <- river_network$network
#' river_params <- river_network$river_types
#' node_info <- river_network$points
#' }
build_river_network <- function(rivers, dem, area = NULL,
                                river_order = NULL, downstream = NULL) {
  msg <- "build_river_network::"
  
  # Validate inputs
  if (!inherits(rivers, "sf")) {
    stop("'rivers' must be an sf object with LINESTRING geometry.\n",
         "To migrate from sp: rivers_sf <- sf::st_as_sf(rivers_sp)",
         call. = FALSE)
  }
  
  if (!inherits(dem, "SpatRaster")) {
    stop("'dem' must be a SpatRaster object.\n",
         "To migrate from raster: dem_terra <- terra::rast(dem_raster)",
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
  
  nsp <- nrow(rivers)
  
  # Extract coordinates
  message(msg, "Extracting coordinates...")
  coords <- get_coords(rivers)
  
  # Get from-to nodes
  message(msg, "Identifying from/to nodes...")
  ft <- get_from_to_nodes(rivers, coords)
  
  # Calculate river order if not provided
  if (is.null(river_order)) {
    message(msg, "Calculating river order...")
    river_order <- calc_river_order(rivers)
  } else {
    if (length(river_order) != nsp) {
      stop("'river_order' must have length equal to number of river segments",
           call. = FALSE)
    }
  }
  
  # Calculate downstream if not provided
  if (is.null(downstream)) {
    message(msg, "Identifying downstream relationships...")
    downstream <- calc_river_downstream(rivers, coords = coords)
  } else {
    if (length(downstream) != nsp) {
      stop("'downstream' must have length equal to number of river segments",
           call. = FALSE)
    }
  }
  
  # Calculate lengths
  message(msg, "Calculating segment lengths...")
  lengths <- as.numeric(sf::st_length(rivers))
  
  # Extract elevations and calculate slopes
  message(msg, "Extracting elevations and calculating slopes...")
  from_coords <- coords[ft[, 2], , drop = FALSE]
  to_coords <- coords[ft[, 3], , drop = FALSE]
  
  # Convert to terra vect for extraction
  from_pts <- terra::vect(from_coords, crs = terra::crs(dem))
  to_pts <- terra::vect(to_coords, crs = terra::crs(dem))
  
  # Extract elevations
  z_from <- .extract_spatraster_values(dem, from_pts)
  z_to <- .extract_spatraster_values(dem, to_pts)
  
  # Calculate slopes
  slopes <- (z_from - z_to) / lengths
  
  # Handle negative or zero slopes
  slopes[slopes <= 0] <- 1e-5
  
  # Build river data frame
  river_df <- data.frame(
    Index = 1:nsp,
    Down = downstream,
    Type = river_order,
    Slope = slopes,
    Length = lengths,
    BC = 0  # Boundary condition (0 = no BC)
  )
  
  # Add to sf object
  result_sf <- rivers
  result_sf$Index <- river_df$Index
  result_sf$Down <- river_df$Down
  result_sf$Type <- river_df$Type
  result_sf$Slope <- river_df$Slope
  result_sf$Length <- river_df$Length
  result_sf$BC <- river_df$BC
  
  # Build point data frame
  point_df <- data.frame(
    From.x = from_coords[, 1],
    From.y = from_coords[, 2],
    From.z = z_from,
    To.x = to_coords[, 1],
    To.y = to_coords[, 2],
    To.z = z_to
  )
  rownames(point_df) <- 1:nsp
  
  # Generate river types
  message(msg, "Generating river type parameters...")
  ntype <- max(river_order)
  
  if (is.null(area)) {
    river_types <- generate_river_types(ntype)
  } else {
    widths <- calc_river_width_from_area(area, ntype)
    river_types <- generate_river_types(ntype, width = widths)
  }
  
  # Return list
  list(
    network = result_sf,
    river_types = river_types,
    points = point_df
  )
}

#' Get River Outlets
#'
#' Identifies outlet segments in a river network (segments with no downstream).
#'
#' @param rivers sf object with river network, or integer vector of downstream indices
#' @return Integer vector of outlet segment indices
#'
#' @details
#' Outlets are identified as segments where the downstream index is negative
#' (typically -1).
#'
#' @section Migration Note:
#' This function replaces \code{getOutlets()} with a simpler interface.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' rivers <- st_read("rivers.shp")
#' downstream <- calc_river_downstream(rivers)
#' outlets <- get_river_outlets(downstream)
#' }
get_river_outlets <- function(rivers) {
  if (inherits(rivers, "sf")) {
    # Extract downstream from sf object
    if (!"Down" %in% names(rivers) && !"downstream" %in% names(rivers)) {
      stop("River network must have 'Down' or 'downstream' column",
           call. = FALSE)
    }
    
    downstream <- if ("Down" %in% names(rivers)) {
      rivers$Down
    } else {
      rivers$downstream
    }
  } else if (is.numeric(rivers)) {
    downstream <- rivers
  } else {
    stop("'rivers' must be an sf object or numeric vector of downstream indices",
         call. = FALSE)
  }
  
  which(downstream < 0)
}

#' Convert River Network to SHUD.RIVER S4 Object
#'
#' Converts a river network (from build_river_network) to the SHUD.RIVER
#' S4 class format for compatibility with existing SHUD workflows.
#'
#' @param network_list List output from \code{build_river_network()}
#' @return SHUD.RIVER S4 object
#'
#' @details
#' This function creates a SHUD.RIVER object that includes both the SHUD tabular
#' river fields and the modern sf network representation stored in the
#' \code{network} slot.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' rivers <- st_read("rivers.shp")
#' dem <- rast("dem.tif")
#' network_list <- build_river_network(rivers, dem)
#' shud_river <- as_shud_river(network_list)
#' }
as_shud_river <- function(network_list) {
  # Validate input
  if (!is.list(network_list)) {
    stop("'network_list' must be a list", call. = FALSE)
  }
  
  required <- c("network", "river_types", "points")
  missing <- setdiff(required, names(network_list))
  if (length(missing) > 0) {
    stop("Missing required components: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  
  # Extract components
  network_sf <- network_list$network
  river_types <- network_list$river_types
  points <- network_list$points
  
  # Build SHUD tabular river fields.
  river_df <- data.frame(
    Index = network_sf$Index,
    Down = network_sf$Down,
    Type = network_sf$Type,
    Slope = network_sf$Slope,
    Length = network_sf$Length,
    BC = network_sf$BC
  )
  
  # Get CRS
  crs_str <- as.character(sf::st_crs(network_sf)$input)
  if (is.na(crs_str)) {
    crs_str <- ""
  }
  
  # Create SHUD.RIVER object
  SHUD.RIVER(
    river = river_df,
    rivertype = river_types,
    point = points,
    network = network_sf,
    crs = crs_str
  )
}

#' Convert SHUD.RIVER to sf Object
#'
#' Extracts the sf network from a SHUD.RIVER object, or reconstructs it from
#' tabular river fields if needed.
#'
#' @param shud_river SHUD.RIVER S4 object
#' @return sf object with river network
#'
#' @details
#' If the SHUD.RIVER object has a network slot, it is returned directly.
#' Otherwise, the function attempts to reconstruct an sf object from the point
#' coordinates.
#'
#' @export
#' @examples
#' \dontrun{
#' shud_river <- read_river("model.riv")
#' # Convert to sf
#' rivers_sf <- shud_river_to_sf(shud_river)
#' }
shud_river_to_sf <- function(shud_river) {
  # Validate input
  if (!inherits(shud_river, "SHUD River")) {
    stop("'shud_river' must be a SHUD.RIVER object", call. = FALSE)
  }
  
  # Check if modern format with network slot
  if (!is.null(shud_river@network) && inherits(shud_river@network, "sf")) {
    return(shud_river@network)
  }
  
  # Reconstruct from SHUD tabular river fields.
  message("Reconstructing sf object from SHUD tabular river fields...")
  
  river_df <- shud_river@river
  point_df <- shud_river@point
  
  nriv <- nrow(river_df)
  
  # Build linestrings from point coordinates
  geoms <- lapply(seq_len(nriv), function(i) {
    coords <- matrix(c(
      point_df[i, "From.x"], point_df[i, "From.y"],
      point_df[i, "To.x"], point_df[i, "To.y"]
    ), ncol = 2, byrow = TRUE)
    sf::st_linestring(coords)
  })
  
  # Get CRS if available
  crs_val <- if (length(shud_river@crs) > 0 && shud_river@crs != "") {
    shud_river@crs
  } else {
    NA
  }
  
  # Create sf object
  rivers_sf <- sf::st_sf(
    Index = river_df$Index,
    Down = river_df$Down,
    Type = river_df$Type,
    Slope = river_df$Slope,
    Length = river_df$Length,
    BC = river_df$BC,
    geometry = sf::st_sfc(geoms, crs = crs_val)
  )
  
  rivers_sf
}

#' Check if SHUD.RIVER Uses Modern Format
#'
#' Determines if a SHUD.RIVER object uses the modern sf-based format
#' or tabular fields only.
#'
#' @param shud_river SHUD.RIVER S4 object
#' @return Logical, TRUE if modern format
#'
#' @export
is_modern_river <- function(shud_river) {
  if (!inherits(shud_river, "SHUD River")) {
    stop("'shud_river' must be a SHUD.RIVER object", call. = FALSE)
  }
  
  !is.null(shud_river@network) && inherits(shud_river@network, "sf")
}
