
#' Return NLCD colors.
#' \code{NLCD.colors} 
#' @param  x, which is the integer number defined in NLCD 2001,2006, and 2011.  Any value out of (0,100) will be assigned as black.
#' @param def.col Default color, if the landuse code is missing.
#' @keywords NLCD, Landuse, colormap
#' @return HEX color strings
#' @source See the link of NLCD: \url{http://www.mrlc.gov/nlcd11_leg.php} 
NLCD.colors <- function(x,def.col = '#EEEEEE'){
  # 2001, 2006, 2011 version 
  #reference http://www.mrlc.gov/nlcd11_leg.php
  lccol = matrix(def.col, nrow=100)
  lccol[1]=	"#00fa00"
  lccol[11]=	"#476ba1"
  lccol[12]=	"#d1defa"
  lccol[21]=	"#decaca"
  lccol[22]=	"#d99482"
  lccol[23]=	"#ee0000"
  lccol[24]=	"#ab0000"
  lccol[31]=	"#b3aea3"
  lccol[32]=	"#fafafa"
  lccol[41]=	"#68ab63"
  lccol[42]=	"#1c6330"
  lccol[43]=	"#b5ca8f"
  lccol[51]=	"#a68c30"
  lccol[52]=	"#ccba7d"
  lccol[71]=	"#e3e3c2"
  lccol[72]=	"#caca78"
  lccol[73]=	"#99c247"
  lccol[74]=	"#78ae94"
  lccol[81]=	"#dcd93d"
  lccol[82]=	"#ab7028"
  lccol[90]=	"#bad9eb"
  lccol[91]=	"#b5d4e6"
  lccol[92]=	"#b5d4e6"
  lccol[93]=	"#b5d4e6"
  lccol[94]=	"#b5d4e6"
  lccol[95]=	"#70a3ba"
  
  x=round(x)
  #x=sort(round(x), decreasing=FALSE);
  ret = rep(def.col, length(x));
  id = which( x %in% 1:100)
  ret[id] = lccol[x[id]]
  return(ret)
}

#' Return NLCD class names.
#' \code{NLCD.names} 
#' @param  x, which is the integer number defined in NLCD 2001,2006, and 2011.  Any value out of (0,100) will be assigned as black.
#' @param def.name Default name.
#' @keywords NLCD, Landuse
#' @return Characters
NLCD.names <- function(x,def.name = 'n/a'){
  # 2001, 2006, 2011 version 
  #reference http://www.mrlc.gov/nlcd11_leg.php
  
  nlcd.names = matrix(def.name, nrow=100)
  nlcd.names[11]=	"Open Water"
  nlcd.names[12]=	"Perennial Ice/Snow"
  nlcd.names[21]=	"Developed, Open Space"
  nlcd.names[22]=	"Developed, Low Intensity"
  nlcd.names[23]=	"Developed, Medium Intensity"
  nlcd.names[24]=	"Developed, High Intensity"
  nlcd.names[31]=	"Barren Land (Rock/Sand/Clay)"
  nlcd.names[41]=	"Deciduous Forest"
  nlcd.names[42]=	"Evergreen Forest"
  nlcd.names[43]=	"Mixed Forest"
  nlcd.names[51]=	"Dwarf Scrub"
  nlcd.names[52]=	"Shrub/Scrub"
  nlcd.names[71]=	"Grassland/Herbaceous"
  nlcd.names[72]=	"Sedge/Herbaceous"
  nlcd.names[73]=	"Lichens"
  nlcd.names[74]=	"Moss"
  nlcd.names[81]=	"Pasture/Hay"
  nlcd.names[82]=	"Cultivated Crops"
  nlcd.names[90]=	"Woody Wetlands"
  nlcd.names[95]=	"Emergent Herbaceous Wetlands"
  x=round(x);
  ret = rep(def.name, length(x));
  id = which( x %in% 1:100)
  ret[id] = nlcd.names[x[id]]
  return(ret)
}
