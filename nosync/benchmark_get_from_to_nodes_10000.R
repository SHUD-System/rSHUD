suppressPackageStartupMessages({
  library(sf)
  devtools::load_all(".")
})

cat("Generating sf object with 10000 points...\n")

n_lines <- 1000
pts_per_line <- 10
lines_list <- vector("list", n_lines)

for (i in 1:n_lines) {
  x <- seq(i, i + 1, length.out = pts_per_line)
  y <- seq(0, 1, length.out = pts_per_line) + rnorm(pts_per_line, 0, 0.01)
  lines_list[[i]] <- sf::st_linestring(cbind(x, y))
}

sfc <- sf::st_sfc(lines_list, crs = 4326)
rivers <- sf::st_sf(ID = 1:n_lines, geometry = sfc)

cat("Total lines:", nrow(rivers), "\n")

cat("Extracting unique coordinates...\n")
coords <- get_coords(rivers)

cat("Benchmarking get_from_to_nodes...\n")
t0 <- Sys.time()
ft <- get_from_to_nodes(rivers, coords)
t1 <- Sys.time()

time_diff <- as.numeric(difftime(t1, t0, units = "secs"))
cat(sprintf("Calculation Speed: %.4f seconds\n", time_diff))

cat("Verifying correctness...\n")
correct <- TRUE
for (i in 1:nrow(rivers)) {
  pts <- sf::st_coordinates(sf::st_geometry(rivers)[[i]])[, 1:2, drop = FALSE]
  fr_pt_actual <- unname(pts[1, , drop = FALSE])
  to_pt_actual <- unname(pts[nrow(pts), , drop = FALSE])
  
  fr_idx <- as.numeric(ft[i, "FrNode"])
  to_idx <- as.numeric(ft[i, "ToNode"])
  
  fr_pt_calc <- unname(coords[fr_idx, , drop = FALSE])
  to_pt_calc <- unname(coords[to_idx, , drop = FALSE])
  
  if (is.na(fr_idx) || is.na(to_idx) || any(fr_pt_actual != fr_pt_calc) || any(to_pt_actual != to_pt_calc)) {
    correct <- FALSE
    cat(sprintf("Mismatch found at line ID %s\n", ft[i, "ID"]))
    break
  }
}

if (correct) {
  cat("Correctness Verification: PASSED.\n")
} else {
  cat("Correctness Verification: FAILED.\n")
}
