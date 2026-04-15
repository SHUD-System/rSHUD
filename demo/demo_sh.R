# =============================================================================
# Demo: SHUD model construction for the 'sh' (Shale Hills) watershed
# =============================================================================
# Manual workflow: simplify boundary/rivers, triangulation, mesh, river
# prep, forcing metadata, and example I/O.
#
# Extra packages (not all are hard Imports of rSHUD): sp, raster, rgeos,
# fields, colorspace. Install with install.packages(...).
# =============================================================================

set.seed(123)

pkg_demo <- c("sp", "raster", "rgeos", "fields", "colorspace", "xts")
miss <- pkg_demo[!vapply(pkg_demo, requireNamespace, NA, quietly = TRUE)]
if (length(miss)) {
  stop(
    "demo_sh needs: ", paste(miss, collapse = ", "),
    "\nInstall, e.g. install.packages(c(\"",
    paste(miss, collapse = "\", \""), "\"))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(sp)
  library(raster)
  library(rgeos)
  library(fields)
  library(xts)
})
library(rSHUD)

# Resolve writers from namespace (works when installed package predates export())
.get_rshud <- function(fname) {
  if (!fname %in% names(asNamespace("rSHUD"))) {
    stop("rSHUD is missing ", fname, "(); install current rSHUD from source.", call. = FALSE)
  }
  get(fname, envir = asNamespace("rSHUD"))
}
write_forc <- .get_rshud("write_forc")
write_tsd <- .get_rshud("write_tsd")

prjname <- "sh"
model.in <- file.path(tempdir(), "demo_sh", "input", prjname)
model.out <- file.path(tempdir(), "demo_sh", "output", paste0(prjname, ".out"))
fin <- shud.filein(prjname, inpath = model.in, outpath = model.out)

dir.create(model.in, showWarnings = FALSE, recursive = TRUE)
dir.create(model.out, showWarnings = FALSE, recursive = TRUE)

message("Input directory: ", model.in)
message("Output directory: ", model.out)

utils::data("sh", package = "rSHUD")

wbd <- sh[["wbd"]]
riv <- sh[["riv"]]
dem <- sh[["dem"]]

# Forcing series: bundle may omit 'forc'; then use a short synthetic daily xts
if (!is.null(sh[["forc"]])) {
  tsd.forc <- sh[["forc"]]
} else {
  message("Dataset 'sh' has no 'forc' component; using synthetic daily forcing for demo I/O.")
  nd <- 365L
  dts <- seq(as.Date("2010-01-01"), by = "day", length.out = nd)
  tsd.forc <- xts(
    cbind(
      Precip_mm.d = runif(nd, 0, 0.01) * 86400,
      Temp_C = runif(nd, 5, 25),
      RH_1 = runif(nd, 0.4, 0.9),
      Wind_m.s = runif(nd, 1, 5),
      RN_w.m2 = runif(nd, 100, 400),
      Pres_pa = runif(nd, 95000, 102000)
    ),
    order.by = dts
  )
}

a.max <- 200
q.min <- 33
tol.riv <- 5
tol.wb <- 5
aqd <- 3
NX <- 800

years <- seq(
  as.numeric(format(min(time(tsd.forc)), "%Y")),
  as.numeric(format(max(time(tsd.forc)), "%Y"))
)
ndays <- days_in_year(years)

message("Model period: ", min(years), " to ", max(years))
message("Total days: ", sum(ndays))

# Simplify geometries (sf); convert back to sp for legacy helpers
wbd_sf <- sf::st_as_sf(wbd)
wb_u <- sf::st_union(wbd_sf)
wb_simp_sf <- sf::st_simplify(sf::st_as_sf(wb_u), dTolerance = tol.wb)
wb.simp <- sf::as_Spatial(wb_simp_sf)

riv.simp <- sf::as_Spatial(
  sf::st_simplify(sf::st_as_sf(riv), dTolerance = tol.riv)
)
riv.simp <- sp.CutSptialLines(sl = riv.simp, tol = 20)

tri <- shud.triangle(wb = wb.simp, q = q.min, a = a.max, S = NX)
message("Generated ", nrow(tri$T), " triangles")

pm <- shud.mesh(tri, dem = dem, AqDepth = aqd)
mesh_crs <- sf::st_crs(wbd_sf)
if (inherits(mesh_crs, "crs") && is.na(mesh_crs)) {
  mesh_crs <- raster::crs(wbd)
}
mesh_sf <- rSHUD::mesh_to_sf(pm = pm, crs = mesh_crs)
ncell <- nrow(pm@mesh)
message("Mesh cells: ", ncell)

pa <- shud.att(tri)
message("Mesh attribute table rows: ", nrow(pa))

pr <- shud.river(riv.simp, dem = dem)

spr <- riv.simp
sp.seg <- shud.rivseg(mesh_sf, sf::st_as_sf(spr))
prs <- sp.seg
message("River segments after mesh intersection: ", nrow(prs))

pic <- shud.ic(nrow(pm@mesh), nrow(pr@river), AqD = aqd, p1 = 0.2, p2 = 0.2)
message("IC object prepared (cells x rivers); see ?write_ic to export.")

go.plot <- function() {
  z <- getElevation(pm = pm)
  idx.ord <- order(z)
  col <- colorspace::diverge_hcl(length(z))
  plot(mesh_sf[idx.ord, ], axes = TRUE, col = col, lwd = 0.5,
       main = "SHUD Mesh with Elevation")
  plot(spr, col = "blue", add = TRUE, lwd = 3)
  image.plot(
    legend.only = TRUE, zlim = range(z), col = col, horizontal = FALSE,
    legend.args = list(text = "Elevation (m)", side = 3, line = 0.05, font = 2, adj = 0.2),
    smallplot = c(0.79, 0.86, 0.20, 0.4)
  )
}

ia <- getArea(pm = pm)

png.filename <- file.path(model.out, "sh_mesh.png")
grDevices::png(filename = png.filename, height = 9, width = 6, res = 400, units = "in")
graphics::par(mfrow = c(2, 1), mar = c(3, 3.5, 1.5, 1))
go.plot()
graphics::par(mfrow = c(1, 1))
graphics::mtext(side = 3, text = "(a) Mesh with Elevation")
graphics::hist(ia, xlab = "", nclass = 20, main = "", ylab = "")
graphics::mtext(side = 3, text = "(b) Area Distribution")
graphics::mtext(side = 2, text = "Frequency", line = 2)
graphics::mtext(side = 1, text = expression(paste("Area (", m^2, ")")), line = 2)
graphics::box()
grDevices::dev.off()

message("Visualization saved to: ", png.filename)

wb_cent <- sf::st_centroid(sf::st_geometry(wb_simp_sf))
sp.c <- sf::as_Spatial(
  sf::st_sf(data.frame(ID = "forcing", stringsAsFactors = FALSE), geometry = wb_cent)
)
sp.forc <- ForcingCoverage(
  sp.meteoSite = sp.c,
  pcs = raster::crs(wb.simp),
  dem = dem,
  wbd = wbd
)

write_forc(
  sp.forc@data,
  file = fin[["md.forc"]],
  path = file.path(model.in, prjname),
  startdate = format(min(time(tsd.forc)), "%Y%m%d")
)

write_tsd(tsd.forc, file = file.path(fin[["inpath"]], "forcing.csv"))

message("Demo completed successfully!")
message("All files written under: ", model.in)
