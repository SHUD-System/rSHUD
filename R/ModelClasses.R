
#' Class of SHUD.MESH
#' \code{SHUD.MESH} 
#' @slot mesh data.frame or matrix
#' @slot point data.frame or matrix
#' @return Class of SHUD.MESH
#' @importFrom methods new
#' @export
SHUD.MESH <- methods::setClass("Untructure Domain",
                      slots = c(mesh="data.frame", 
                                point="data.frame"))
#' Class of SHUD.RIVER
#' \code{SHUD.RIVER} 
#' @slot river data.frame or matrix
#' @slot rivertype data.frame or matrix
#' @slot point data.frame or matrix
#' @return Class of SHUD.RIVER
#' @importFrom methods new
#' @export
SHUD.RIVER <- methods::setClass("SHUD River", 
                      slots = c(river="data.frame", 
                                rivertype="data.frame",
                                point = 'data.frame'))

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

#' This is data of Lake Sunappe to be included in my package
#' @name sunapee
#' @docType data
#' @keywords sunapee
NULL


#' This is .shud environment.
#' @name .shud
.shud<-new.env()

