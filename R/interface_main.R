#' Main Interface Functions
#'
#' High-level interface functions for automated SHUD model building.
#' These functions provide a streamlined workflow from raw geospatial data
#' to complete SHUD model input files.
#'
#' @name interface_main
NULL

#' Automatically Build SHUD Model
#'
#' Automated workflow to build a complete SHUD hydrological model from
#' geospatial data. This modernized function uses terra and sf exclusively
#' for all spatial operations.
#'
#' @param project_name Character, project name for output files
#' @param domain sf object with POLYGON geometry defining watershed boundary
#' @param rivers sf object with LINESTRING geometry for river network (optional)
#' @param dem SpatRaster object with elevation data
#' @param soil SpatRaster object with soil data (optional)
#' @param geology SpatRaster object with geology data (optional)
#' @param landcover SpatRaster object with land cover data (optional)
#' @param forcing_sites sf object with POINT geometry for forcing locations (optional)
#' @param forcing_files Character vector of forcing data filenames
#' @param output_dir Character, output directory path (default: current directory)
#' @param mesh_options List of mesh generation options:
#'   \itemize{
#'     \item max_area: Maximum triangle area (default: 2e5 m²)
#'     \item min_angle: Minimum triangle angle (default: 33 degrees)
#'     \item simplify_boundary: Tolerance for boundary simplification (default: 200 m)
#'     \item simplify_rivers: Tolerance for river simplification (default: 0 m)
#'   }
#' @param river_options List of river processing options:
#'   \itemize{
#'     \item min_length: Minimum river segment length (default: 0 m)
#'   }
#' @param aquifer_depth Numeric, depth of aquifer bottom below surface (default: 10 m)
#' @param years Integer vector, years for LAI/melt factor time series
#' @param clean Logical, whether to clean existing files in output directory
#' @param parameters List or data.frame, model parameters (default: shud.para())
#' @param calibration List or data.frame, calibration parameters (default: shud.calib())
#' @param melt_factor xts object, melt factor time series (optional)
#' @param remove_outliers Logical, remove outliers in soil/geology data
#' @param verbose Logical, display progress messages (default: TRUE)
#'
#' @details
#' Inputs used for meter-based buffering, simplification, cropping, and mesh
#' generation must have a defined projected CRS in metres/meters. Longitude/
#' latitude CRS inputs and projected CRSs in non-metre units are rejected;
#' transform inputs to an appropriate metre-based projected CRS before calling
#' this workflow. All supplied spatial inputs must use the same CRS as
#' \code{domain}; this workflow does not auto-reproject.
#'
#' @return List with components:
#'   \item{mesh}{SHUD.MESH object}
#'   \item{river}{SHUD.RIVER object}
#'   \item{attributes}{Attribute data frame}
#'   \item{files}{Named vector of output file paths}
#'
#' @section Migration Note:
#' This function replaces \code{autoBuildModel()} with modern spatial libraries.
#' Key differences:
#' \itemize{
#'   \item Only accepts terra::SpatRaster and sf objects
#'   \item Rejects legacy spatial objects with clear error messages
#'   \item Legacy inputs belong in deprecated wrapper workflows such as
#'     \code{\link{autoBuildModel}}, which may perform limited conversion before
#'     forwarding to this modern API
#'   \item Uses snake_case parameter names (e.g., project_name vs prjname)
#'   \item Returns structured list instead of just mesh object
#'   \item Provides verbose progress reporting
#' }
#'
#' To migrate from legacy code:
#' \preformatted{
#' # Modern code (sf/terra):
#' wbd <- sf::st_read("boundary.shp")
#' dem <- terra::rast("dem.tif")
#'
#' # Deprecated autoBuildModel() remains the migration compatibility entry point
#' # for limited legacy conversion/forwarding.
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' # Load data using modern libraries
#' domain <- st_read("watershed.shp")
#' rivers <- st_read("rivers.shp")
#' dem <- rast("dem.tif")
#' soil <- rast("soil.tif")
#' lc <- rast("landcover.tif")
#'
#' # Build model
#' model <- auto_build_model(
#'   project_name = "my_watershed",
#'   domain = domain,
#'   rivers = rivers,
#'   dem = dem,
#'   soil = soil,
#'   landcover = lc,
#'   output_dir = "model_output",
#'   years = 2000:2010
#' )
#'
#' # Access results
#' mesh <- model$mesh
#' river_network <- model$river
#' }
auto_build_model <- function(
  project_name,
  domain,
  rivers = NULL,
  dem,
  soil = NULL,
  geology = NULL,
  landcover = NULL,
  forcing_sites = NULL,
  forcing_files = NULL,
  output_dir = getwd(),
  mesh_options = list(),
  river_options = list(),
  aquifer_depth = 10,
  years = 2000:2010,
  clean = TRUE,
  parameters = NULL,
  calibration = NULL,
  melt_factor = NULL,
  remove_outliers = TRUE,
  verbose = TRUE) {

  msg <- paste0("auto_build_model(", project_name, ")::")

  if (inherits(dem, "PackedSpatRaster")) dem <- terra::unwrap(dem)
  if (inherits(soil, "PackedSpatRaster")) soil <- terra::unwrap(soil)
  if (inherits(geology, "PackedSpatRaster")) geology <- terra::unwrap(geology)
  if (inherits(landcover, "PackedSpatRaster")) landcover <- terra::unwrap(landcover)
  if (inherits(domain, "PackedSpatVector")) domain <- terra::unwrap(domain)
  if (inherits(rivers, "PackedSpatVector")) rivers <- terra::unwrap(rivers)
  if (inherits(forcing_sites, "PackedSpatVector")) forcing_sites <- terra::unwrap(forcing_sites)

  if (verbose) {
    message(msg, " Starting automated SHUD model building")
  }

  # ============================================================================
  # 1. VALIDATE INPUTS
  # ============================================================================

  if (verbose) message(msg, " Validating inputs...")

  # Check required parameters
  if (missing(project_name) || is.null(project_name)) {
    stop("Parameter 'project_name' is required", call. = FALSE)
  }
  if (missing(domain) || is.null(domain)) {
    stop("Parameter 'domain' (watershed boundary) is required", call. = FALSE)
  }
  if (missing(dem) || is.null(dem)) {
    stop("Parameter 'dem' (elevation data) is required", call. = FALSE)
  }

  # Reject legacy raster/sp objects with migration guidance
  if (inherits(domain, c("Spatial", "SpatialPolygons", "SpatialPolygonsDataFrame"))) {
    stop(
      "Legacy sp objects are not supported.\n\n",
      "To migrate your code:\n",
      "  NEW: domain <- sf::st_read('boundary.shp')\n\n",
      "See inst/MIGRATION_GUIDE.md for complete migration instructions.",
      call. = FALSE
    )
  }

  if (inherits(dem, c("Raster", "RasterLayer", "RasterStack", "RasterBrick"))) {
    stop(
      "Legacy raster objects are not supported.\n\n",
      "To migrate your code:\n",
      "  NEW: dem <- terra::rast('dem.tif')\n\n",
      "See inst/MIGRATION_GUIDE.md for complete migration instructions.",
      call. = FALSE
    )
  }

  # Validate domain is sf
  if (!inherits(domain, "sf")) {
    stop(
      "Parameter 'domain' must be an sf object with POLYGON geometry.\n",
      "Convert using: domain <- sf::st_as_sf(domain)",
      call. = FALSE
    )
  }

  # Validate dem is SpatRaster
  if (!inherits(dem, "SpatRaster")) {
    stop(
      "Parameter 'dem' must be a terra SpatRaster object.\n",
      "Convert using: dem <- terra::rast(dem)",
      call. = FALSE
    )
  }

  # Validate optional spatial inputs
  if (!is.null(rivers)) {
    if (inherits(rivers, c("Spatial", "SpatialLines", "SpatialLinesDataFrame"))) {
      stop(
        "Legacy sp river objects are not supported.\n",
        "Convert using: rivers <- sf::st_read('rivers.shp')",
        call. = FALSE
      )
    }
    if (!inherits(rivers, "sf")) {
      stop("Parameter 'rivers' must be an sf object", call. = FALSE)
    }
  }

  if (!is.null(soil) && !inherits(soil, "SpatRaster")) {
    stop("Parameter 'soil' must be a SpatRaster object", call. = FALSE)
  }
  if (!is.null(geology) && !inherits(geology, "SpatRaster")) {
    stop("Parameter 'geology' must be a SpatRaster object", call. = FALSE)
  }
  if (!is.null(landcover) && !inherits(landcover, "SpatRaster")) {
    stop("Parameter 'landcover' must be a SpatRaster object", call. = FALSE)
  }
  if (!is.null(forcing_sites) && !inherits(forcing_sites, "sf")) {
    stop("Parameter 'forcing_sites' must be an sf object", call. = FALSE)
  }

  # Meter-based distances in this workflow are unsafe in missing or geographic CRS.
  check_sf_projected_crs(domain, "domain")
  check_raster_projected_crs(dem, "dem")
  if (!is.null(rivers)) {
    check_sf_projected_crs(rivers, "rivers")
  }
  if (!is.null(forcing_sites)) {
    check_sf_projected_crs(forcing_sites, "forcing_sites")
  }
  if (!is.null(soil)) {
    check_raster_projected_crs(soil, "soil")
  }
  if (!is.null(geology)) {
    check_raster_projected_crs(geology, "geology")
  }
  if (!is.null(landcover)) {
    check_raster_projected_crs(landcover, "landcover")
  }

  check_model_builder_crs_match(dem, domain, "dem", "domain")
  if (!is.null(rivers)) {
    check_model_builder_crs_match(rivers, domain, "rivers", "domain")
  }
  if (!is.null(forcing_sites)) {
    check_model_builder_crs_match(forcing_sites, domain,
                                  "forcing_sites", "domain")
  }
  if (!is.null(soil)) {
    check_model_builder_crs_match(soil, domain, "soil", "domain")
  }
  if (!is.null(geology)) {
    check_model_builder_crs_match(geology, domain, "geology", "domain")
  }
  if (!is.null(landcover)) {
    check_model_builder_crs_match(landcover, domain, "landcover", "domain")
  }

  # Validate numeric parameters
  check_positive(aquifer_depth, "aquifer_depth")

  # ============================================================================
  # 2. SETUP OUTPUT DIRECTORIES
  # ============================================================================

  if (verbose) message(msg, " Setting up output directories...")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- file.path(output_dir, "fig")
  gis_dir <- file.path(output_dir, "gis")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(gis_dir, showWarnings = FALSE, recursive = TRUE)

  # Setup file paths
  file_paths <- shud.filein(project_name, inpath = output_dir,
                            outpath = output_dir)

  # ============================================================================
  # 3. PROCESS MESH OPTIONS
  # ============================================================================

  # Set default mesh options
  default_mesh_opts <- list(
    max_area = 2e5,
    min_angle = 33,
    simplify_boundary = 200,
    simplify_rivers = 0
  )
  mesh_options <- utils::modifyList(default_mesh_opts, mesh_options)

  # Set default river options
  default_river_opts <- list(
    min_length = 0
  )
  river_options <- utils::modifyList(default_river_opts, river_options)

  # ============================================================================
  # 4. CROP AND BUFFER DATA
  # ============================================================================

  if (verbose) message(msg, " Buffering and cropping spatial data...")

  # Buffer domain for data extraction
  buffer_dist <- max(c(2000, mesh_options$simplify_boundary))
  domain_buffered <- sf::st_buffer(domain, dist = buffer_dist)

  # Crop DEM to buffered domain
  dem_cropped <- terra::crop(dem, terra::vect(domain_buffered))

  # Save initial plot
  if (verbose) {
    png(filename = file.path(fig_dir, "data_0_input.png"),
        height = 11, width = 11, res = 100, units = "in")
    terra::plot(dem_cropped, main = "Input Data")
    plot(sf::st_geometry(domain), add = TRUE, border = "red", lwd = 2)
    if (!is.null(rivers)) {
      plot(sf::st_geometry(rivers), add = TRUE, col = "blue", lwd = 2)
    }
    dev.off()
  }

  # ============================================================================
  # 5. SIMPLIFY BOUNDARY AND RIVERS
  # ============================================================================

  if (verbose) message(msg, " Simplifying boundary and river network...")

  # Simplify boundary
  domain_union <- sf::st_union(domain)
  if (mesh_options$simplify_boundary > 0) {
    domain_simplified <- sf::st_simplify(domain_union,
                                         dTolerance = mesh_options$simplify_boundary)
  } else {
    domain_simplified <- domain_union
  }
  domain_simplified <- sf::st_sf(geometry = domain_simplified)

  # Simplify rivers if provided
  rivers_simplified <- NULL
  if (!is.null(rivers)) {
    if (mesh_options$simplify_rivers > 0) {
      rivers_simplified <- sf::st_simplify(rivers,
                                          dTolerance = mesh_options$simplify_rivers)
    } else {
      rivers_simplified <- rivers
    }

    # Further simplify by length if specified
    if (river_options$min_length > 0) {
      lengths <- as.numeric(sf::st_length(rivers_simplified))
      rivers_simplified <- rivers_simplified[lengths >= river_options$min_length, ]
    }
  }

  # Save simplified plot
  if (verbose) {
    png(filename = file.path(fig_dir, "data_1_simplified.png"),
        height = 11, width = 11, res = 100, units = "in")
    terra::plot(dem_cropped, main = "Simplified Boundary and Rivers")
    plot(sf::st_geometry(domain_simplified), add = TRUE, border = "red", lwd = 2)
    if (!is.null(rivers_simplified)) {
      plot(sf::st_geometry(rivers_simplified), add = TRUE, col = "blue", lwd = 2)
    }
    dev.off()
  }

  # ============================================================================
  # 6. GENERATE TRIANGULAR MESH
  # ============================================================================

  if (verbose) message(msg, " Generating triangular mesh...")

  tri <- shud.triangle(
    wb = domain_simplified,
    riv = rivers_simplified,
    q = mesh_options$min_angle,
    a = mesh_options$max_area
  )

  n_triangles <- nrow(tri$T)
  if (verbose) {
    message(msg, " Number of triangles: ", n_triangles)
  }

  # ============================================================================
  # 7. CREATE SHUD MESH OBJECT
  # ============================================================================

  if (verbose) message(msg, " Creating SHUD mesh object...")

  mesh_obj <- shud.mesh(tri, dem = dem_cropped, AqDepth = aquifer_depth)

  # Convert to shapefile and save
  mesh_sf <- mesh_to_sf(mesh_obj, crs = sf::st_crs(domain))
  sf::st_write(mesh_sf, file.path(gis_dir, "domain.shp"),
               delete_dsn = TRUE, quiet = TRUE)

  # Save mesh plot
  if (verbose) {
    png(filename = file.path(fig_dir, "data_2_mesh.png"),
        height = 11, width = 11, res = 100, units = "in")
    plot(sf::st_geometry(mesh_sf), main = "Generated Mesh")
    dev.off()
  }

  # ============================================================================
  # 8. CALCULATE MESH ATTRIBUTES
  # ============================================================================

  if (verbose) message(msg, " Calculating mesh attributes...")

  # Prepare attribute rasters
  attr_soil <- if (!is.null(soil)) soil else 1:n_triangles
  attr_geol <- if (!is.null(geology)) geology else 1:n_triangles

  # Calculate attributes
  mesh_att <- shud.att(
    tri,
    r.soil = attr_soil,
    r.geol = attr_geol,
    r.lc = landcover,
    r.forc = forcing_sites
  )

  # ============================================================================
  # 9. PROCESS RIVER NETWORK
  # ============================================================================

  river_obj <- NULL
  river_seg <- NULL

  if (!is.null(rivers_simplified)) {
    if (verbose) message(msg, " Processing river network...")

    # Calculate watershed area for width estimation
    watershed_area <- as.numeric(sf::st_area(domain)) * 1e-6  # Convert to km²

    # Build river network
    river_network <- build_river_network(
      rivers = rivers_simplified,
      dem = dem_cropped,
      area = watershed_area
    )

    # Convert to SHUD.RIVER object
    river_obj <- as_shud_river(river_network)

    # Save river shapefile
    sf::st_write(river_network$network, file.path(gis_dir, "river.shp"),
                 delete_dsn = TRUE, quiet = TRUE)

    # Check for multiple outlets
    outlets <- get_river_outlets(river_network$network)
    if (length(outlets) > 1) {
      message(msg, " WARNING: ", length(outlets), " outlets detected in river network")
    }

    # Cut rivers with mesh triangles
    if (verbose) message(msg, " Cutting rivers with mesh triangles...")
    river_seg <- shud.rivseg(mesh_sf, river_network$network)
    sf::st_write(river_seg, file.path(gis_dir, "seg.shp"),
                 delete_dsn = TRUE, quiet = TRUE)
  }

  # ============================================================================
  # 10. GENERATE INITIAL CONDITIONS
  # ============================================================================

  if (verbose) message(msg, " Generating initial conditions...")

  n_rivers <- if (!is.null(river_obj)) nrow(river_obj@river) else 0
  initial_conditions <- shud.ic(nrow(mesh_obj@mesh), n_rivers,
                                AqD = aquifer_depth)

  # ============================================================================
  # 11. GENERATE TIME SERIES DATA
  # ============================================================================

  if (verbose) message(msg, " Generating time series data...")

  # LAI and roughness length
  if (!is.null(landcover)) {
    lc_values <- unique(terra::values(landcover, na.rm = TRUE))
    lai_rl <- fun.lairl(lc_values, years = years)

    # Save LAI/RL plot
    if (verbose) {
      png(filename = file.path(fig_dir, "data_3_lai_rl.png"),
          height = 11, width = 11, res = 100, units = "in")
      par(mfrow = c(2, 1))
      zoo::plot.zoo(lai_rl$LAI, main = "Leaf Area Index")
      zoo::plot.zoo(lai_rl$RL, main = "Roughness Length")
      dev.off()
    }

    write_tsd(lai_rl$LAI, file = file_paths["md.lai"])
    write_tsd(lai_rl$RL, file = file_paths["md.rl"])
  }

  # Melt factor
  if (is.null(melt_factor)) {
    melt_factor <- MeltFactor(years = years)
  }
  write_tsd(melt_factor, file = file_paths["md.mf"])

  # ============================================================================
  # 12. WRITE OUTPUT FILES
  # ============================================================================

  if (verbose) message(msg, " Writing SHUD model input files...")

  # Write mesh and river
  write_mesh(mesh_obj, file = file_paths["md.mesh"])
  if (!is.null(river_obj)) {
    write_river(river_obj, file = file_paths["md.riv"])
  }

  # Write initial conditions
  write_ic(initial_conditions, file = file_paths["md.ic"])

  # Write attributes
  write_df(mesh_att, file = file_paths["md.att"])
  if (!is.null(river_seg)) {
    write_df(river_seg, file = file_paths["md.rivseg"])
  }

  # Write parameters and calibration
  if (is.null(parameters)) {
    n_days <- 365 * length(years) + round(length(years) / 4)
    parameters <- shud.para(nday = n_days)
  }
  if (is.null(calibration)) {
    calibration <- shud.calib()
  }
  write_config(parameters, file_paths["md.para"])
  write_config(calibration, file_paths["md.calib"])

  # Write forcing file list if provided
  if (!is.null(forcing_files)) {
    start_date <- paste0(min(years), "0101")
    write_forc(forcing_files, file = file_paths["md.forc"],
               startdate = start_date)
  }

  # ============================================================================
  # 13. SUMMARY
  # ============================================================================

  if (verbose) {
    message(msg, " Model building complete!")
    message(msg, " Summary:")
    message(msg, "   Mesh cells: ", nrow(mesh_obj@mesh))
    if (!is.null(river_obj)) {
      message(msg, "   River segments: ", nrow(river_obj@river))
    }
    if (!is.null(river_seg)) {
      message(msg, "   River-mesh segments: ", nrow(river_seg))
    }
    message(msg, "   Output directory: ", output_dir)
  }

  # Return results
  invisible(list(
    mesh = mesh_obj,
    river = river_obj,
    attributes = mesh_att,
    river_segments = river_seg,
    files = file_paths
  ))
}


#' @title Deprecated: Build SHUD Model (Legacy Interface)
#' @description
#' This function is deprecated. Use \code{\link{auto_build_model}} instead.
#' It is retained only as a migration compatibility entry point. It performs
#' limited conversion of legacy \code{sp}/\code{raster} inputs to
#' \code{sf}/\code{terra} before forwarding to \code{auto_build_model()}; new
#' code should pass \code{sf} and \code{terra::SpatRaster} objects directly.
#'
#' @param indata Named list with components: wbd (domain), riv (rivers), dem, rsoil, rgeol, rlc, forc
#' @param forcfiles Forcing file paths or data.frame
#' @param prjname Project name string
#' @param ... Additional arguments passed to auto_build_model()
#' @return Same as auto_build_model()
#' @export
autoBuildModel <- function(indata, forcfiles = NULL, prjname = NULL, ...) {
  .Deprecated("auto_build_model",
    msg = paste0("autoBuildModel() is deprecated, use auto_build_model().\n",
                 "Will be removed in v3.1.0."))

  safe_convert_sf <- function(x) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "Spatial")) return(sf::st_as_sf(x))
    x
  }
  safe_convert_terra <- function(x) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "RasterLayer") || inherits(x, "RasterStack") || inherits(x, "RasterBrick"))
      return(terra::rast(x))
    x
  }

  dots <- list(...)
  out_dir     <- dots$outdir %||% dots$output_dir %||% getwd()
  aq_depth    <- dots$AqDepth %||% dots$aquifer_depth %||% 10
  yrs         <- dots$years %||% 2000:2010
  do_clean    <- dots$clean %||% TRUE
  rm_out      <- dots$rm.outlier %||% dots$remove_outliers %||% TRUE
  be_quiet    <- dots$quiet
  verbose_val <- if (!is.null(be_quiet)) !be_quiet else TRUE
  mesh_opts   <- dots$mesh_options %||% list()
  if (!is.null(dots$a.max))    mesh_opts$max_area   <- dots$a.max
  if (!is.null(dots$q.min))    mesh_opts$min_angle  <- dots$q.min
  if (!is.null(dots$tol.wb))   mesh_opts$tol_domain <- dots$tol.wb
  if (!is.null(dots$tol.riv))  mesh_opts$tol_river  <- dots$tol.riv
  if (!is.null(dots$tol.len))  mesh_opts$tol_length <- dots$tol.len

  auto_build_model(
    project_name    = if (is.null(prjname)) "unnamed" else prjname,
    domain          = safe_convert_sf(indata$wbd),
    rivers          = safe_convert_sf(indata$riv),
    dem             = safe_convert_terra(indata$dem),
    soil            = safe_convert_terra(indata$rsoil),
    geology         = safe_convert_terra(indata$rgeol),
    landcover       = safe_convert_terra(indata$rlc),
    forcing_sites   = safe_convert_sf(indata$forc),
    forcing_files   = forcfiles,
    output_dir      = out_dir,
    mesh_options    = mesh_opts,
    aquifer_depth   = aq_depth,
    years           = yrs,
    clean           = do_clean,
    parameters      = dots$cfg.para,
    calibration     = dots$cfg.calib,
    melt_factor     = dots$mf,
    remove_outliers = rm_out,
    verbose         = verbose_val
  )
}


#' Quick Model Building with Default Parameters
#'
#' Simplified interface for building SHUD models with sensible defaults.
#' This function is designed for quick prototyping and common use cases.
#'
#' @param project_name Character, project name for output files
#' @param domain sf object with POLYGON geometry defining watershed boundary
#' @param dem SpatRaster object with elevation data
#' @param rivers sf object with LINESTRING geometry for river network (optional)
#' @param landcover SpatRaster object with land cover data (optional)
#' @param output_dir Character, output directory path (default: current directory)
#' @param years Integer vector, years for time series (default: 2000:2010)
#' @param verbose Logical, display progress messages (default: TRUE)
#'
#' @return List with components:
#'   \item{mesh}{SHUD.MESH object}
#'   \item{river}{SHUD.RIVER object (if rivers provided)}
#'   \item{files}{Named vector of output file paths}
#'
#' @details
#' This function uses default parameters optimized for typical watersheds:
#' \itemize{
#'   \item Maximum triangle area: 200,000 m² (20 ha)
#'   \item Minimum triangle angle: 30 degrees
#'   \item Boundary simplification: 200 m
#'   \item Aquifer depth: 10 m
#' }
#'
#' Input CRS requirements are inherited from \code{\link{auto_build_model}}:
#' spatial inputs must have a defined projected CRS in metres/meters because
#' the workflow uses meter-based distances.
#'
#' For more control over parameters, use \code{\link{auto_build_model}}.
#'
#' @section Migration Note:
#' This is a new function with no legacy equivalent. It provides a streamlined
#' interface for users who want to quickly build models without specifying
#' detailed options.
#'
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#'
#' # Load minimal required data
#' domain <- st_read("watershed.shp")
#' dem <- rast("dem.tif")
#'
#' # Quick model build
#' model <- quick_model(
#'   project_name = "test_model",
#'   domain = domain,
#'   dem = dem
#' )
#'
#' # With rivers and land cover
#' rivers <- st_read("rivers.shp")
#' lc <- rast("landcover.tif")
#'
#' model <- quick_model(
#'   project_name = "full_model",
#'   domain = domain,
#'   dem = dem,
#'   rivers = rivers,
#'   landcover = lc,
#'   output_dir = "quick_model"
#' )
#' }
quick_model <- function(
  project_name,
  domain,
  dem,
  rivers = NULL,
  landcover = NULL,
  output_dir = getwd(),
  years = 2000:2010,
  verbose = TRUE) {

  if (verbose) {
    message("quick_model: Building SHUD model with default parameters...")
  }

  # Validate required inputs
  if (missing(project_name) || is.null(project_name)) {
    stop("Parameter 'project_name' is required", call. = FALSE)
  }
  if (missing(domain) || is.null(domain)) {
    stop("Parameter 'domain' (watershed boundary) is required", call. = FALSE)
  }
  if (missing(dem) || is.null(dem)) {
    stop("Parameter 'dem' (elevation data) is required", call. = FALSE)
  }

  # Validate spatial object types
  if (!inherits(domain, "sf")) {
    stop(
      "Parameter 'domain' must be an sf object.\n",
      "Convert using: domain <- sf::st_read('boundary.shp')",
      call. = FALSE
    )
  }
  if (!inherits(dem, "SpatRaster")) {
    stop(
      "Parameter 'dem' must be a terra SpatRaster object.\n",
      "Convert using: dem <- terra::rast('dem.tif')",
      call. = FALSE
    )
  }
  if (!is.null(rivers) && !inherits(rivers, "sf")) {
    stop(
      "Parameter 'rivers' must be an sf object.\n",
      "Convert using: rivers <- sf::st_read('rivers.shp')",
      call. = FALSE
    )
  }
  if (!is.null(landcover) && !inherits(landcover, "SpatRaster")) {
    stop(
      "Parameter 'landcover' must be a terra SpatRaster object.\n",
      "Convert using: landcover <- terra::rast('landcover.tif')",
      call. = FALSE
    )
  }

  # Set default mesh options optimized for typical watersheds
  mesh_opts <- list(
    max_area = 2e5,        # 20 hectares
    min_angle = 30,        # Good quality triangles
    simplify_boundary = 200,  # Moderate simplification
    simplify_rivers = 0    # No river simplification
  )

  # Set default river options
  river_opts <- list(
    min_length = 0
  )

  # Call auto_build_model with defaults
  result <- auto_build_model(
    project_name = project_name,
    domain = domain,
    rivers = rivers,
    dem = dem,
    soil = NULL,
    geology = NULL,
    landcover = landcover,
    forcing_sites = NULL,
    forcing_files = NULL,
    output_dir = output_dir,
    mesh_options = mesh_opts,
    river_options = river_opts,
    aquifer_depth = 10,
    years = years,
    clean = TRUE,
    parameters = NULL,
    calibration = NULL,
    melt_factor = NULL,
    remove_outliers = TRUE,
    verbose = verbose
  )

  if (verbose) {
    message("quick_model: Model building complete!")
    message("quick_model: Output directory: ", output_dir)
  }

  # Return simplified result (exclude some internal details)
  list(
    mesh = result$mesh,
    river = result$river,
    files = result$files
  )
}
