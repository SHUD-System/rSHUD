
#' Class of SHUD.MESH
#' 
#' S4 class representing an unstructured triangular mesh domain for SHUD models.
#' The class stores mesh topology, point coordinates, and associated attributes.
#' 
#' @slot mesh data.frame with mesh topology (triangles and neighbors)
#' @slot point data.frame with point coordinates and attributes
#' 
#' @section Mesh Data Frame:
#' The mesh slot contains:
#' \itemize{
#'   \item ID: Triangle ID
#'   \item Node1, Node2, Node3: Vertex indices
#'   \item Nabr1, Nabr2, Nabr3: Neighbor triangle indices
#'   \item Zmax: Maximum elevation
#' }
#' 
#' @section Point Data Frame:
#' The point slot contains:
#' \itemize{
#'   \item ID: Point ID
#'   \item X, Y: Coordinates
#'   \item AqDepth: Aquifer depth
#'   \item Elevation: Surface elevation
#' }
#' 
#' @return Class of SHUD.MESH
#' @importFrom methods new setMethod callNextMethod
#' @export
SHUD.MESH <- methods::setClass("Untructure Domain",
                               slots = c(mesh = "data.frame", 
                                         point = "data.frame"))

#' Class of SHUD.RIVER
#'
#' S4 class representing a river network for SHUD models.
#' The class stores river segment properties, hydraulic parameters,
#' and node information. Now supports both legacy data.frame format
#' and modern sf spatial format.
#'
#' @slot river data.frame with river segment properties
#' @slot rivertype data.frame with hydraulic parameters for each river type
#' @slot point data.frame with from/to node coordinates and elevations
#' @slot network sf object with river network geometry (optional, for modern format)
#' @slot crs character string with coordinate reference system (optional)
#'
#' @section River Data Frame:
#' The river slot contains:
#' \itemize{
#'   \item Index: Segment ID
#'   \item Down: Downstream segment index (-1 for outlets)
#'   \item Type: River type/order
#'   \item Slope: Bed slope
#'   \item Length: Segment length
#'   \item BC: Boundary condition flag
#' }
#'
#' @section River Type Data Frame:
#' The rivertype slot contains hydraulic parameters:
#' \itemize{
#'   \item Index: Type index
#'   \item Depth: Channel depth (m)
#'   \item BankSlope: Bank slope
#'   \item Width: Channel width (m)
#'   \item Sinuosity: Channel sinuosity
#'   \item Manning: Manning's roughness coefficient
#'   \item Cwr: Width-discharge coefficient
#'   \item KsatH: Horizontal hydraulic conductivity
#'   \item BedThick: Bed sediment thickness
#' }
#'
#' @section Point Data Frame:
#' The point slot contains node information:
#' \itemize{
#'   \item From.x, From.y, From.z: Start node coordinates and elevation
#'   \item To.x, To.y, To.z: End node coordinates and elevation
#' }
#'
#' @section Modern Format:
#' When created with modern functions (e.g., \code{build_river_network()}),
#' the network slot contains an sf object with LINESTRING geometry and
#' all river properties as attributes. The crs slot stores the coordinate
#' reference system.
#'
#' @return Class of SHUD.RIVER
#' @importFrom methods new
#' @export
SHUD.RIVER <- methods::setClass("SHUD River",
                                slots = c(river = "data.frame",
                                          rivertype = "data.frame",
                                          point = "data.frame",
                                          network = "ANY",
                                          crs = "character"),
                                prototype = list(network = NULL,
                                                crs = ""))

#' @importFrom methods setMethod callNextMethod
setMethod("initialize", "SHUD River", function(.Object, ...) {
  .Object <- callNextMethod()
  if (is.null(.Object@network)) .Object@network <- NULL
  if (length(.Object@crs) == 0) .Object@crs <- ""
  .Object
})

#' Upgrade legacy SHUD.RIVER objects to v3 format
#'
#' Adds missing slots (network, crs) to v2 serialized SHUD.RIVER objects.
#'
#' @param old_river A SHUD.RIVER object loaded from v2 serialization
#' @return Updated SHUD.RIVER object with all required slots
#' @export
upgrade_shud_river <- function(old_river) {
  if (!methods::is(old_river, "SHUD River")) {
    stop("Input must be a SHUD.RIVER object", call. = FALSE)
  }
  if (!methods::.hasSlot(old_river, "network")) {
    methods::slot(old_river, "network") <- NULL
  }
  if (!methods::.hasSlot(old_river, "crs")) {
    methods::slot(old_river, "crs") <- ""
  }
  methods::validObject(old_river)
  old_river
}

#' This is data to be included in my package
#' @name sh
#' @docType data
#' @keywords data
NULL

#' This is data of one of sub-catchmetn in Sacramento Watershed  to be included in my package
#' @name sac
#' @docType data
#' @keywords data
NULL

#' This is water-balance example data bundled with the package
#' @name wb
#' @docType data
#' @keywords data
NULL

#' This is Waerma Watershed example data bundled with the package
#' @name waerma
#' @docType data
#' @keywords data
NULL

#' This is data of Lake Sunappe to be included in my package
#' @name sunapee
#' @docType data
#' @keywords sunapee
NULL


#' This is .shud environment.
#' @name .shud
.shud <- new.env()

