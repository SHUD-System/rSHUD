#' SHUD File I/O Functions
#' 
#' This module provides modern functions for reading and writing SHUD model files.
#' All functions use snake_case naming convention and provide comprehensive
#' error handling and validation.
#' 
#' @name shud-io
NULL

# =============================================================================
# Modern SHUD File Reading Functions (snake_case naming)
# =============================================================================

#' Read SHUD mesh file
#' 
#' Reads a SHUD mesh file (.mesh) and returns a SHUD.MESH object.
#' 
#' @param file Character. Full path to the .mesh file
#' @return SHUD.MESH object containing mesh and point data
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' mesh <- read_mesh("model.mesh")
#' }
read_mesh <- function(file = shud.filein()['md.mesh']) {
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("Mesh file does not exist: ", file, call. = FALSE)
  }
  
  d <- read_df(file = file)
  ret <- SHUD.MESH(mesh = data.frame(d[[1]]),
                   point = data.frame(d[[2]]))
  return(ret)
}

#' Read SHUD river file
#' 
#' Reads a SHUD river file (.riv) and returns a SHUD.RIVER object.
#' 
#' @param file Character. Full path to the .riv file
#' @return SHUD.RIVER object containing river network data
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' river <- read_river("model.riv")
#' }
read_river <- function(file = shud.filein()['md.riv']) {
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("River file does not exist: ", file, call. = FALSE)
  }
  
  d <- read_df(file = file)
  ret <- SHUD.RIVER(river = data.frame(d[[1]]), 
                    rivertype = data.frame(d[[2]]),
                    point = data.frame())
  return(ret)
}

#' Read SHUD attributes file
#' 
#' Reads a SHUD attributes file (.att) and returns attribute data.
#' 
#' @param file Character. Full path to the .att file
#' @return Matrix containing attribute data
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' att <- read_att("model.att")
#' }
read_att <- function(file = shud.filein()['md.att']) {
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("Attributes file does not exist: ", file, call. = FALSE)
  }
  
  x <- read_df(file)
  ret <- x[[1]]
  return(ret)
}

#' Read SHUD river segments file
#' 
#' Reads a SHUD river segments file (.rivseg) and returns segment data.
#' 
#' @param file Character. Full path to the .rivseg file
#' @return Matrix containing river segment data
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' rivseg <- read_rivseg("model.rivseg")
#' }
read_rivseg <- function(file = shud.filein()['md.rivseg']) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("River segments file does not exist: ", file, call. = FALSE)
  }
  
  d <- read_df(file = file)
  ret <- d[[1]]
  return(ret)
}

#' Read SHUD parameters file
#' 
#' Reads a SHUD parameters file (.para) and returns parameter data.
#' 
#' @param file Character. Full path to the .para file
#' @return Data frame containing parameter data
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' para <- read_para("model.para")
#' }
read_para <- function(file = shud.filein()['md.para']) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  return(read_config(file))
}

#' Read SHUD calibration file
#' 
#' Reads a SHUD calibration file (.calib) and returns calibration data.
#' 
#' @param file Character. Full path to the .calib file
#' @return Data frame containing calibration data
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' calib <- read_calib("model.calib")
#' }
read_calib <- function(file = shud.filein()['md.calib']) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  return(read_config(file))
}

#' Read SHUD configuration file
#' 
#' Reads a SHUD configuration file (.para or .calib) and returns configuration data.
#' 
#' @param file Character. Full path to the configuration file
#' @return Data frame containing configuration data
#' @importFrom utils type.convert read.table
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' config <- read_config("model.para")
#' }
read_config <- function(file = shud.filein()['md.para']) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("Configuration file does not exist: ", file, call. = FALSE)
  }
  
  tline <- readLines(file, skipNul = TRUE)
  tline <- tline[!grepl('^#', tline)]
  x <- utils::read.table(text = tline, header = FALSE, stringsAsFactors = FALSE)
  xdf <- data.frame(rbind(t(x[, -1]), NULL), stringsAsFactors = FALSE)
  colnames(xdf) <- toupper(as.character(as.matrix(x[, 1])))
  ret <- xdf
  return(ret)
}

#' Read SHUD initial conditions file
#' 
#' Reads a SHUD initial conditions file (.ic) and returns initial condition data.
#' 
#' @param file Character. Full path to the .ic file
#' @return List containing initial condition data
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' ic <- read_ic("model.ic")
#' }
read_ic <- function(file = shud.filein()['md.ic']) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("Initial conditions file does not exist: ", file, call. = FALSE)
  }
  
  x <- read_df(file)
  cn <- c('minit', 'rinit', 'linit')
  names(x) <- cn[seq_along(x)]
  return(x)
}

#' Read SHUD soil file
#' 
#' Reads a SHUD soil file (.soil) and returns soil data.
#' 
#' @param file Character. Full path to the .soil file
#' @return Matrix containing soil data
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' soil <- read_soil("model.soil")
#' }
read_soil <- function(file = shud.filein()['md.soil']) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("Soil file does not exist: ", file, call. = FALSE)
  }
  
  x <- read_df(file)
  ret <- x[[1]]
  return(ret)
}

#' Read SHUD geology file
#' 
#' Reads a SHUD geology file (.geol) and returns geology data.
#' 
#' @param file Character. Full path to the .geol file
#' @return Matrix containing geology data
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' geol <- read_geol("model.geol")
#' }
read_geol <- function(file = shud.filein()['md.geol']) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("Geology file does not exist: ", file, call. = FALSE)
  }
  
  x <- read_df(file)
  ret <- x[[1]]
  return(ret)
}

#' Read SHUD land cover file
#' 
#' Reads a SHUD land cover file (.lc) and returns land cover data.
#' 
#' @param file Character. Full path to the .lc file
#' @return Matrix containing land cover data
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' lc <- read_lc("model.lc")
#' }
read_lc <- function(file = shud.filein()['md.lc']) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("Land cover file does not exist: ", file, call. = FALSE)
  }
  
  x <- read_df(file)
  ret <- x[[1]]
  return(ret)
}

#' Read SHUD forcing file list
#' 
#' Reads a SHUD forcing file (.forc) and returns forcing site information.
#' 
#' @param file Character. Full path to the .forc file
#' @return List containing start time and forcing sites data frame
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' forc <- read_forc_fn("model.forc")
#' }
read_forc_fn <- function(file = shud.filein()['md.forc']) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop("Forcing file does not exist: ", file, call. = FALSE)
  }
  
  txt <- readLines(file)
  hd <- read.table(text = txt[1])
  path <- txt[2]
  df <- read.table(text = txt[-1 * 1:2], header = TRUE)
  nc <- ncol(df)
  df[, nc] <- file.path(path, df[, nc])
  ret <- list('StartTime' = as.character(hd[2]),
              Sites = df)
  return(ret)
}

#' Read SHUD forcing CSV files
#' 
#' Reads the CSV files referenced in a SHUD forcing file (.forc).
#' 
#' @param file Character. Full path to the .forc file
#' @param id Numeric vector. Index of forcing sites to read. If NULL, returns average of all sites
#' @return Forcing data as xts object or list of xts objects
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' forc_data <- read_forc_csv("model.forc", id = 1:3)
#' }
read_forc_csv <- function(file = shud.filein()['md.forc'], id = NULL) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  msg <- 'read_forc_csv::'
  xf <- read_forc_fn(file = file)
  fns <- xf$Sites[, ncol(xf$Sites)]
  
  if (is.null(id)) {
    fn <- fns
    RID <- seq_along(fns)
  } else {
    fn <- fns[id]
    RID <- id
  }
  
  nf <- length(fn)
  xl <- list()
  for (i in 1:nf) {
    message(msg, i, '/', nf, '\t', fn[i])
    tsd <- read_tsd(fns[RID[i]])[[1]]
    xl[[i]] <- tsd
  }
  names(xl) <- basename(fns[RID])
  
  if (is.null(id)) {
    mx <- xl[[1]] * 0
    for (i in 1:nf) {
      mx <- mx + xl[[i]]
    }
    return(mx / nf)
  } else {
    return(xl)
  }
}

#' Read matrix or data frame file with header
#' 
#' Reads a file containing matrix or data frame data with (nrow, ncol) in header.
#' This is a core utility function used by other SHUD file readers.
#' 
#' @param file Character. Full path to file
#' @param text Character vector. Text content from readLines(file). If provided, file is ignored
#' @param sep Character. Separator for the table (default: tab)
#' @return List of matrices/data frames
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' data <- read_df("data.txt")
#' }
read_df <- function(file, text = readLines(file), sep = '\t') {
  if (missing(file) && missing(text)) {
    stop("Either 'file' or 'text' parameter is required", call. = FALSE)
  }
  
  if (!missing(file) && !file.exists(file)) {
    stop("File does not exist: ", file, call. = FALSE)
  }
  
  r0 <- 1
  nrow <- length(text)
  xl <- list()
  
  for (i in 1:100) {
    ndim <- as.numeric(utils::read.table(text = text[r0]))
    nr <- ndim[1]
    
    if (nr <= 0) {
      nr <- nrow - 2
    }
    
    xl[[i]] <- rbind(utils::read.table(text = text[0:nr + 1 + r0], header = TRUE, sep = sep))
    r0 <- r0 + nr + 2
    
    if (r0 + 1 > nrow) {
      break
    }
  }
  
  return(xl)
}

# =============================================================================
# Spatial Data Reading Functions
# =============================================================================

#' Read river shapefile
#' 
#' Reads a river shapefile and returns a SpatialLines object.
#' Note: This function uses legacy spatial libraries and will be updated
#' in future versions to use sf.
#' 
#' @param file Character. Full path to river shapefile
#' @return SpatialLines object
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' river_sp <- read_river_sp("rivers.shp")
#' }
read_river_sp <- function(file = file.path(shud.filein()['inpath'], 'gis', 'river.shp')) {
  if (missing(file) || is.null(file) || is.na(file)) {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!file.exists(file)) {
    stop('River shapefile does not exist: ', file, call. = FALSE)
  }
  
  spr <- as(sf::st_read(file, quiet = TRUE), "Spatial")
  return(spr)
}

# =============================================================================
# Modern SHUD File Writing Functions (snake_case naming)
# =============================================================================

#' Write SHUD mesh file
#' 
#' Writes a SHUD.MESH object to a mesh file (.mesh).
#' 
#' @param pm SHUD.MESH object containing mesh and point data
#' @param file Character. Full path to the output .mesh file
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' write_mesh(mesh_object, "model.mesh")
#' }
write_mesh <- function(pm, file) {
  msg <- "write_mesh::"
  
  if (missing(pm)) {
    stop("Parameter 'pm' (SHUD.MESH object) is required", call. = FALSE)
  }
  
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!inherits(pm, "Untructure Domain")) {
    stop("Parameter 'pm' must be a SHUD.MESH object", call. = FALSE)
  }
  
  # Ensure directory exists
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  
  message(msg, "Writing to file ", file)
  write_df(pm@mesh, file = file, append = FALSE, quiet = TRUE)
  write_df(pm@point, file = file, append = TRUE, quiet = TRUE)
}

#' Write SHUD river file
#' 
#' Writes a SHUD.RIVER object to a river file (.riv).
#' 
#' @param pr SHUD.RIVER object containing river network data
#' @param file Character. Full path to the output .riv file
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' write_river(river_object, "model.riv")
#' }
write_river <- function(pr, file) {
  msg <- "write_river::"
  
  if (missing(pr)) {
    stop("Parameter 'pr' (SHUD.RIVER object) is required", call. = FALSE)
  }
  
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!inherits(pr, "SHUD River")) {
    stop("Parameter 'pr' must be a SHUD.RIVER object", call. = FALSE)
  }
  
  # Ensure directory exists
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  
  message(msg, "Writing to file ", file)
  write_df(pr@river, file = file, append = FALSE, quiet = TRUE)
  write_df(pr@rivertype, file = file, append = TRUE, quiet = TRUE)
  if (length(pr@point) > 0) {
    write_df(pr@point, file = file, append = TRUE, quiet = TRUE)
  }
}

#' Write SHUD initial conditions file
#' 
#' Writes initial condition data to a SHUD initial conditions file (.ic).
#' 
#' @param x List containing initial condition data
#' @param file Character. Full path to the output .ic file
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' write_ic(ic_data, "model.ic")
#' }
write_ic <- function(x, file) {
  msg <- "write_ic::"
  
  if (missing(x)) {
    stop("Parameter 'x' (initial condition data) is required", call. = FALSE)
  }
  
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!is.list(x) || length(x) < 2) {
    stop("Parameter 'x' must be a list with at least 2 elements", call. = FALSE)
  }
  
  # Ensure directory exists
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  
  message(msg, "Writing to file ", file)
  write_df(x[[1]], file = file, append = FALSE, quiet = TRUE)
  write_df(x[[2]], file = file, append = TRUE, quiet = TRUE)
}

#' Write SHUD configuration file
#' 
#' Writes SHUD model configuration parameters to a file (.para, .calib, etc.).
#' 
#' @param x SHUD model configuration parameters or calibration data
#' @param file Character. Full path to the output configuration file
#' @importFrom utils write.table
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' write_config(config_data, "model.para")
#' }
write_config <- function(x, file) {
  msg <- "write_config::"
  
  if (missing(x)) {
    stop("Parameter 'x' (configuration data) is required", call. = FALSE)
  }
  
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  # Ensure directory exists
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  
  message(msg, "Writing to file ", file)
  out <- cbind(names(x), t(x))
  utils::write.table(out, file, append = FALSE, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
}

#' Write SHUD forcing file
#' 
#' Writes SHUD forcing site information to a forcing file (.forc).
#' 
#' @param x Data frame of forcing sites information with columns:
#'   (id, lon, lat, x, y, z, filename)
#' @param file Character. Full path to the output .forc file
#' @param path Character. Common path of the forcing files (default: "")
#' @param startdate Character. Start date in format YYYYMMDD (default: "20000101")
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' write_forc(forc_sites, "model.forc", path = "forcing/", startdate = "20100101")
#' }
write_forc <- function(x, file, path = "", startdate = "20000101") {
  msg <- "write_forc::"
  
  if (missing(x)) {
    stop("Parameter 'x' (forcing sites data) is required", call. = FALSE)
  }
  
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (!is.data.frame(x)) {
    stop("Parameter 'x' must be a data frame", call. = FALSE)
  }
  
  # Ensure directory exists
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  
  nf <- nrow(x)
  nc <- ncol(x)
  message(msg, "Writing to file ", file)
  write(paste(nf, startdate), file = file, append = FALSE)
  write(path, file = file, append = TRUE, ncolumns = 1)
  write(colnames(x), file = file, ncolumns = nc, append = TRUE, sep = "\t")
  write(t(x), file = file, ncolumns = nc, append = TRUE, sep = "\t")
}

#' Write data frame to file with header
#' 
#' Writes a data frame to file with (nrow, ncol) header.
#' This is a core utility function used by other SHUD file writers.
#' 
#' @param x Data frame or matrix to write
#' @param file Character. Full path to output file
#' @param append Logical. Whether to append to existing file (default: FALSE)
#' @param quiet Logical. Whether to suppress messages (default: FALSE)
#' @param header Numeric vector. Custom header (default: c(nrow, ncol))
#' @family shud-io
#' @export
#' @examples
#' \dontrun{
#' write_df(data_matrix, "output.txt")
#' }
write_df <- function(x, file, append = FALSE, quiet = FALSE, header = NULL) {
  msg <- "write_df::"
  
  if (missing(x)) {
    stop("Parameter 'x' (data to write) is required", call. = FALSE)
  }
  
  if (missing(file) || is.null(file) || is.na(file) || file == "") {
    stop("Parameter 'file' is required", call. = FALSE)
  }
  
  if (inherits(x, "sf")) {
    x <- as.data.frame(sf::st_drop_geometry(x), stringsAsFactors = FALSE)
  }
  
  x <- as.matrix(rbind(x))
  nr <- nrow(x)
  nc <- ncol(x)
  
  if (is.null(header)) {
    header <- c(nr, nc)
  }
  
  if (!quiet) {
    message(msg, "Writing to file ", file)
  }
  
  write(header, file = file, append = append, sep = "\t", 
        ncolumns = length(header))
  write(colnames(x), file = file, ncolumns = nc, append = TRUE, sep = "\t")
  write(t(x), file = file, ncolumns = nc, append = TRUE, sep = "\t")
}


# =============================================================================
# Deprecated Function Aliases (for backward compatibility)
# =============================================================================

#' @rdname read_mesh
#' @export
#' @section Deprecated:
#' \code{readmesh()} is deprecated. Use \code{read_mesh()} instead.
readmesh <- function(file = shud.filein()["md.mesh"]) {
  .Deprecated("read_mesh", package = "rSHUD",
              msg = "readmesh() is deprecated. Use read_mesh() instead.")
  read_mesh(file)
}

#' @rdname read_river
#' @export
#' @section Deprecated:
#' \code{readriv()} is deprecated. Use \code{read_river()} instead.
readriv <- function(file = shud.filein()["md.riv"]) {
  .Deprecated("read_river", package = "rSHUD",
              msg = "readriv() is deprecated. Use read_river() instead.")
  read_river(file)
}

#' @rdname read_att
#' @export
#' @section Deprecated:
#' \code{readatt()} is deprecated. Use \code{read_att()} instead.
readatt <- function(file = shud.filein()["md.att"]) {
  .Deprecated("read_att", package = "rSHUD",
              msg = "readatt() is deprecated. Use read_att() instead.")
  read_att(file)
}

#' @rdname read_rivseg
#' @export
#' @section Deprecated:
#' \code{readrivseg()} is deprecated. Use \code{read_rivseg()} instead.
readrivseg <- function(file = shud.filein()["md.rivseg"]) {
  .Deprecated("read_rivseg", package = "rSHUD",
              msg = "readrivseg() is deprecated. Use read_rivseg() instead.")
  read_rivseg(file)
}

#' @rdname read_para
#' @export
#' @section Deprecated:
#' \code{readpara()} is deprecated. Use \code{read_para()} instead.
readpara <- function(file = shud.filein()["md.para"]) {
  .Deprecated("read_para", package = "rSHUD",
              msg = "readpara() is deprecated. Use read_para() instead.")
  read_para(file)
}

#' @rdname read_calib
#' @export
#' @section Deprecated:
#' \code{readcalib()} is deprecated. Use \code{read_calib()} instead.
readcalib <- function(file = shud.filein()["md.calib"]) {
  .Deprecated("read_calib", package = "rSHUD",
              msg = "readcalib() is deprecated. Use read_calib() instead.")
  read_calib(file)
}

#' @rdname read_config
#' @export
#' @section Deprecated:
#' \code{readconfig()} is deprecated. Use \code{read_config()} instead.
readconfig <- function(file = shud.filein()["md.para"]) {
  .Deprecated("read_config", package = "rSHUD",
              msg = "readconfig() is deprecated. Use read_config() instead.")
  read_config(file)
}

#' @rdname read_ic
#' @export
#' @section Deprecated:
#' \code{readic()} is deprecated. Use \code{read_ic()} instead.
readic <- function(file = shud.filein()["md.ic"]) {
  .Deprecated("read_ic", package = "rSHUD",
              msg = "readic() is deprecated. Use read_ic() instead.")
  read_ic(file)
}

#' @rdname read_soil
#' @export
#' @section Deprecated:
#' \code{readsoil()} is deprecated. Use \code{read_soil()} instead.
readsoil <- function(file = shud.filein()["md.soil"]) {
  .Deprecated("read_soil", package = "rSHUD",
              msg = "readsoil() is deprecated. Use read_soil() instead.")
  read_soil(file)
}

#' @rdname read_geol
#' @export
#' @section Deprecated:
#' \code{readgeol()} is deprecated. Use \code{read_geol()} instead.
readgeol <- function(file = shud.filein()["md.geol"]) {
  .Deprecated("read_geol", package = "rSHUD",
              msg = "readgeol() is deprecated. Use read_geol() instead.")
  read_geol(file)
}

#' @rdname read_lc
#' @export
#' @section Deprecated:
#' \code{readlc()} is deprecated. Use \code{read_lc()} instead.
readlc <- function(file = shud.filein()["md.lc"]) {
  .Deprecated("read_lc", package = "rSHUD",
              msg = "readlc() is deprecated. Use read_lc() instead.")
  read_lc(file)
}

#' @rdname read_forc_fn
#' @export
#' @section Deprecated:
#' \code{readforc.fn()} is deprecated. Use \code{read_forc_fn()} instead.
readforc.fn <- function(file = shud.filein()["md.forc"]) {
  .Deprecated("read_forc_fn", package = "rSHUD",
              msg = "readforc.fn() is deprecated. Use read_forc_fn() instead.")
  read_forc_fn(file)
}

#' @rdname read_forc_csv
#' @export
#' @section Deprecated:
#' \code{readforc.csv()} is deprecated. Use \code{read_forc_csv()} instead.
readforc.csv <- function(file = shud.filein()["md.forc"], id = NULL) {
  .Deprecated("read_forc_csv", package = "rSHUD",
              msg = "readforc.csv() is deprecated. Use read_forc_csv() instead.")
  read_forc_csv(file, id)
}

#' @rdname read_df
#' @export
#' @section Deprecated:
#' \code{read.df()} is deprecated. Use \code{read_df()} instead.
read.df <- function(file, text = readLines(file), sep = "\t") {
  .Deprecated("read_df", package = "rSHUD",
              msg = "read.df() is deprecated. Use read_df() instead.")
  read_df(file, text, sep)
}

#' @rdname read_river_sp
#' @export
#' @section Deprecated:
#' \code{readriv.sp()} is deprecated. Use \code{read_river_sp()} instead.
readriv.sp <- function(file = file.path(shud.filein()["inpath"], "gis", "river.shp")) {
  .Deprecated("read_river_sp", package = "rSHUD",
              msg = "readriv.sp() is deprecated. Use read_river_sp() instead.")
  read_river_sp(file)
}
