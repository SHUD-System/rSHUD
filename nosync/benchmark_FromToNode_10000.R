# benchmark_FromToNode_10000.R
suppressPackageStartupMessages({
  library(sp)
  devtools::load_all(".")
})

cat("Generating SpatialLines object with 10000 points...\n")

# Create a large SpatialLines object
# Let's create 1000 lines, each with 10 points -> 10000 points total.
n_lines <- 1000
pts_per_line <- 10
lines_list <- vector("list", n_lines)

for (i in 1:n_lines) {
  x <- seq(i, i + 1, length.out = pts_per_line)
  y <- seq(0, 1, length.out = pts_per_line) + rnorm(pts_per_line, 0, 0.01)
  lines_list[[i]] <- sp::Lines(list(sp::Line(cbind(x, y))), ID = as.character(i))
}

sl <- sp::SpatialLines(lines_list, proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
cat("Total lines:", length(sl), "\n")
cat("Total points:", sum(sapply(sl@lines, function(x) nrow(x@Lines[[1]]@coords))), "\n")

cat("Extracting unique coordinates...\n")
coord <- get_coords(sl)

cat("Benchmarking get_from_to_nodes...\n")
t0 <- Sys.time()
ft <- get_from_to_nodes(sl, coords = coord)
t1 <- Sys.time()

time_diff <- as.numeric(difftime(t1, t0, units = "secs"))
cat(sprintf("Calculation Speed: %.4f seconds\n", time_diff))

cat("Verifying correctness...\n")
# Check if FrNode and ToNode point to the correct coordinates
correct <- TRUE
for (i in 1:length(sl)) {
  # Get actual first and last points of the line
  pts <- sl@lines[[i]]@Lines[[1]]@coords
  fr_pt_actual <- pts[1, , drop = FALSE]
  to_pt_actual <- pts[nrow(pts), , drop = FALSE]
  
  # Get points from FromToNode result
  fr_idx <- as.numeric(ft[i, "FrNode"])
  to_idx <- as.numeric(ft[i, "ToNode"])
  
  fr_pt_calc <- unname(coord[fr_idx, , drop = FALSE])
  to_pt_calc <- unname(coord[to_idx, , drop = FALSE])
  fr_pt_actual <- unname(fr_pt_actual)
  to_pt_actual <- unname(to_pt_actual)
  
  # Check if they match
  if (is.na(fr_idx) || is.na(to_idx) || any(fr_pt_actual != fr_pt_calc) || any(to_pt_actual != to_pt_calc)) {
    correct <- FALSE
    cat(sprintf("Mismatch found at line ID %s\n", ft[i, "ID"]))
    cat("Actual From:\n")
    print(fr_pt_actual)
    cat("Calculated From Index:", fr_idx, "\n")
    cat("Calculated From Point:\n")
    print(fr_pt_calc)
    break
  }
}

if (correct) {
  cat("Correctness Verification: PASSED. All From/To nodes perfectly match the start/end coordinates of the lines.\n")
} else {
  cat("Correctness Verification: FAILED.\n")
}
