
#' Deprecated: split legacy polylines by length
#' \code{SimplifybyLen} 
#' @param sp \code{sf} line object; legacy \code{SpatialLines} input is accepted
#'   for compatibility.
#' @param split_length the length of the segments to split the lines into, in
#'   map units.
#' @param plot.results TRUE/FALSE 
#' @importFrom grDevices dev.off graphics.off png rgb topo.colors 
#' @importFrom graphics grid hist lines par plot points 
#' @importFrom methods as 
#' @importFrom stats dist rnorm time 
#' @importFrom utils read.table 
#' @return Legacy \code{SpatialLinesDataFrame} for compatibility.
#' @source  https://stackoverflow.com/questions/38700246/how-do-i-split-divide-polyline-shapefiles-into-equally-length-smaller-segments
#' @export
SimplifybyLen <- function(sp, split_length = 20, plot.results = FALSE) {
  sp_sf = if (inherits(sp, "sf")) sp else sf::st_as_sf(sp)
  sp_legacy = methods::as(sp_sf, 'SpatialLines')
  #### Define support functions ####
  # Compatibility helper for the legacy SpatialLines representation used below.
  SpatialLines2df = function(sp) {
    df = data.frame(
      id = character(),
      mline_id = character(),
      segment_id = character(),
      fx = numeric(),
      fy = numeric(),
      tx = numeric(),
      ty = numeric(),
      stringsAsFactors = FALSE
    )
    for (i in 1:length(sp)) {
      coords = sp@lines[[i]]@Lines[[1]]@coords # For each line takes the coords of the vertex
      row_nums = 1:(nrow(coords) - 1)
      mline_id = formatC(i, width = 9, flag = '0') # Creates id for the line
      segment_id = formatC(row_nums, width = 3, flag = '0') # Creates id for each single segment belonging to the line
      id = paste0(mline_id, '_', segment_id) # Creates a composite id
      for (j in row_nums) {
        # For each segment stores ids and coordinates
        df[nrow(df) + 1, ] = c(id[j],
                                mline_id,
                                segment_id[j],
                                coords[j, 1],
                                coords[j, 2],
                                coords[j + 1, 1],
                                coords[j + 1, 2])
      }
    }
    row.names(df) = NULL
    df$fx = as.numeric(df$fx)
    df$fy = as.numeric(df$fy)
    df$tx = as.numeric(df$tx)
    df$ty = as.numeric(df$ty)
    return(df)
  }
  
  # Convert a data frame of IDs and coordinates into sf line features.
  # linedf: a data.frame with columns as:
  #         id = generic ids of the lines,
  #         fx = coordinates x of the first point of the line
  #         fy = coordinates y of the first point of the line
  #         tx = coordinates x of the second point of the line
  #         tx = coordinates y of the second point of the line
  

  linedf2sf = function(linedf, crs) {
    geoms = vector("list", nrow(linedf))
    for (i in 1:nrow(linedf)) {
      c1 = cbind(rbind(linedf$fx[i], linedf$tx[i]),
                 rbind(linedf$fy[i], linedf$ty[i]))
      geoms[[i]] = sf::st_linestring(c1)
    }
    sf::st_sf(
      id = linedf$id,
      geometry = sf::st_sfc(geoms, crs = crs)
    )
  }
  
  
  #### Split the lines ####
  # Convert the compatibility line object into a data frame.
  linedf = SpatialLines2df(sp_legacy)
  df = data.frame(
    id = character(),
    fx = numeric(),
    fy = numeric(),
    tx = numeric(),
    ty = numeric(),
    stringsAsFactors = FALSE
  )
  
  
  for (i in 1:nrow(linedf)) {
    # For each line of the dataframe, corresponding to a single line of the spatial object
    # skips if length is less then split_length
    v_seg = linedf[i, ]
    seg_length = sqrt((v_seg$fx - v_seg$tx) ^ 2 + (v_seg$fy - v_seg$ty) ^ 2) # segment length
    if (seg_length <= split_length) {
      df[nrow(df) + 1, ] = c(paste0(v_seg$id, '_', '0000'),
                              v_seg$fx,
                              v_seg$fy,
                              v_seg$tx,
                              v_seg$ty)
      next()
    }
    # Create a vector of direction as the line and unit length
    # vector v corresponding to the line
    v = c(v_seg$tx  -  v_seg$fx,
          v_seg$ty  -  v_seg$fy)
    # vector of direction v and unit length
    u = c(v[1]  /  sqrt(v[1]  ^  2 + v[2]  ^  2), v[2]  /  sqrt(v[1]  ^  2 + v[2]  ^ 2))
    # Calculates how many segment the line is split into and the leftover
    num_seg = floor(seg_length  /  split_length)
    seg_left = seg_length - (num_seg  *  split_length)
    
    # Add to the output dataframe each segment plus the leftover
    for (i in 0:(num_seg  -  1)) {
      # Add num_seg segments
      df[nrow(df)  +  1,] = c(
        paste0(v_seg$id, '_', formatC(i, width = 4, flag = '0')),
        v_seg$fx + u[1]  *  split_length  *  i,
        v_seg$fy + u[2]  *  split_length  *  i,
        v_seg$fx + u[1]  *  split_length  *  (i  +  1),
        v_seg$fy + u[2]  *  split_length  *  (i  +  1)
      )
    }
    df[nrow(df) + 1,] = c(
      paste0(v_seg$id, '_', formatC(
        num_seg, width = 4, flag = '0'
      )),
      # Add leftover segment
      v_seg$fx + u[1] * split_length * num_seg,
      v_seg$fy + u[2] * split_length * num_seg,
      v_seg$tx,
      v_seg$ty
    )
    
  }
  
  #### Visualise the results to check ####
  if (plot.results) {
    plot(sf::st_geometry(sp_sf))
    coords = cbind(as.numeric(df$fx), as.numeric(df$fy))
    coords = rbind(coords, as.numeric(df$tx[nrow(df)]), as.numeric(df$ty)[nrow(df)])
    pt_sf = sf::st_as_sf(data.frame(x = coords[, 1], y = coords[, 2]), coords = c("x", "y"),
                         crs = sf::st_crs(sp_sf))
    plot(sf::st_geometry(pt_sf), col = 'red', add = TRUE)
  }
  
  #### Output ####
  df$fx = as.numeric(df$fx)
  df$fy = as.numeric(df$fy)
  df$tx = as.numeric(df$tx)
  df$ty = as.numeric(df$ty)
  sl = linedf2sf(df, sf::st_crs(sp_sf))
  
  att=data.frame('Index'=1:nrow(df), df[,-1], 'Length' = as.numeric(sf::st_length(sl)))
  rownames(att) = df[,1]
  sld=sf::st_sf(att, geometry = sf::st_geometry(sl))
  return(methods::as(sld, "Spatial"))
}
