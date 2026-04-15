suppressPackageStartupMessages({
  library(sf)
  devtools::load_all(".")
})

n_lines <- 100
pts_per_line <- 10
lines_list <- vector("list", n_lines)

for (i in 1:n_lines) {
  x <- seq(i, i + 1, length.out = pts_per_line)
  y <- seq(0, 1, length.out = pts_per_line) + rnorm(pts_per_line, 0, 0.01)
  lines_list[[i]] <- sf::st_linestring(cbind(x, y))
}

sfc <- sf::st_sfc(lines_list, crs = 4326)
rivers <- sf::st_sf(ID = 1:n_lines, geometry = sfc)

cat("Testing get_from_to_nodes with sf object...\n")
ft <- get_from_to_nodes(rivers)

cat("Success! First few rows:\n")
print(head(ft))
