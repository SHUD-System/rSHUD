#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rSHUD)
  library(sp)
  library(raster)
  library(rgeos)
})

repo_dir <- "/Volumes/CloudDisk/CloudDrive/Development/rSHUD"
data_file <- file.path(repo_dir, "data", "sh.rda")
out_dir <- file.path(repo_dir, "nosync", "river_fromto_points_check_sh_only")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_lookup <- function(coords, ids) {
  out <- matrix(NA_real_, nrow = length(ids), ncol = 2)
  colnames(out) <- colnames(coords)
  ok <- !is.na(ids) & ids >= 1 & ids <= nrow(coords)
  if (any(ok)) out[ok, ] <- coords[ids[ok], , drop = FALSE]
  out
}

coord_key <- function(xy) {
  if (nrow(xy) == 0) return(character())
  apply(signif(xy, 12), 1, function(z) paste(z[1], z[2], sep = "|"))
}

topology_pairs <- function(fr_xy, to_xy) {
  fr_key <- coord_key(fr_xy)
  to_key <- coord_key(to_xy)
  n <- length(fr_key)
  if (n == 0) return(character())
  pairs <- character()
  for (i in seq_len(n)) {
    j <- which(fr_key == to_key[i])
    j <- j[j != i]
    if (length(j) > 0) pairs <- c(pairs, paste(i, j, sep = "->"))
  }
  unique(pairs)
}


riv <- sh_data$riv
coords_orig <- get_coords(riv)
ft_truth <- get_from_to_nodes(riv, coords = coords_orig)
ft_current <- get_from_to_nodes(riv)
ft_simpl_with_coord <- get_from_to_nodes(riv, coords = coords_orig)

fr_truth <- safe_lookup(coords_orig, ft_truth[, "FrNode"])
to_truth <- safe_lookup(coords_orig, ft_truth[, "ToNode"])
fr_current <- safe_lookup(coords_orig, ft_current[, "FrNode"])
to_current <- safe_lookup(coords_orig, ft_current[, "ToNode"])
fr_simpl_with_coord <- safe_lookup(coords_orig, ft_simpl_with_coord[, "FrNode"])
to_simpl_with_coord <- safe_lookup(coords_orig, ft_simpl_with_coord[, "ToNode"])

fr_mis <- rowSums(abs(fr_truth - fr_current), na.rm = TRUE) > 0 |
  is.na(rowSums(fr_truth + fr_current))
to_mis <- rowSums(abs(to_truth - to_current), na.rm = TRUE) > 0 |
  is.na(rowSums(to_truth + to_current))
any_mis <- fr_mis | to_mis

truth_pairs <- topology_pairs(fr_truth, to_truth)
current_pairs <- topology_pairs(fr_current, to_current)
miss_pairs <- setdiff(truth_pairs, current_pairs)
extra_pairs <- setdiff(current_pairs, truth_pairs)

png(file.path(out_dir, "sh_fromto_compare.png"), width = 2200, height = 1000, res = 160)
par(mfrow = c(1, 2), mar = c(4, 4, 4, 1))
plot(riv, col = "grey55", lwd = 1.2, main = "sh$riv | truth simplify=FALSE")
points(fr_truth, pch = 16, cex = 0.45, col = "#1f77b4")
points(to_truth, pch = 1, cex = 0.55, col = "#d62728")
legend("topleft", legend = c("From", "To"), pch = c(16, 1), col = c("#1f77b4", "#d62728"), bty = "n")
plot(riv, col = "grey55", lwd = 1.2, main = "sh$riv | current simplify=TRUE")
points(fr_current, pch = 16, cex = 0.45, col = "#1f77b4")
points(to_current, pch = 1, cex = 0.55, col = "#d62728")
bad <- which(any_mis & !is.na(fr_current[, 1]) & !is.na(to_current[, 1]))
if (length(bad) > 0) {
  segments(fr_truth[bad, 1], fr_truth[bad, 2], fr_current[bad, 1], fr_current[bad, 2], col = "#1f77b488", lwd = 1)
  segments(to_truth[bad, 1], to_truth[bad, 2], to_current[bad, 1], to_current[bad, 2], col = "#d6272888", lwd = 1)
}
dev.off()

detail <- data.frame(
  dataset = "sh__sh",
  reach_id = ft_truth[, "ID"],
  fr_truth_x = fr_truth[, 1],
  fr_truth_y = fr_truth[, 2],
  to_truth_x = to_truth[, 1],
  to_truth_y = to_truth[, 2],
  fr_current_x = fr_current[, 1],
  fr_current_y = fr_current[, 2],
  to_current_x = to_current[, 1],
  to_current_y = to_current[, 2],
  fr_mismatch = fr_mis,
  to_mismatch = to_mis,
  any_mismatch = any_mis,
  stringsAsFactors = FALSE
)
write.csv(detail, file.path(out_dir, "sh_detail.csv"), row.names = FALSE)

summary_df <- data.frame(
  dataset = "sh__sh",
  n_reach = nrow(ft_truth),
  n_fr_mismatch = sum(fr_mis, na.rm = TRUE),
  n_to_mismatch = sum(to_mis, na.rm = TRUE),
  n_any_mismatch = sum(any_mis, na.rm = TRUE),
  truth_topology_pairs = length(truth_pairs),
  current_topology_pairs = length(current_pairs),
  missing_topology_pairs = length(miss_pairs),
  extra_topology_pairs = length(extra_pairs),
  n_any_mismatch_simplify_with_coord = sum(
    rowSums(abs(fr_truth - fr_simpl_with_coord), na.rm = TRUE) > 0 |
      rowSums(abs(to_truth - to_simpl_with_coord), na.rm = TRUE) > 0,
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(out_dir, "sh_summary.csv"), row.names = FALSE)
print(summary_df)
