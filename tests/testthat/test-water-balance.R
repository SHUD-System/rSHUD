test_that("wb.DS respects supplied initial conditions", {
  dates <- as.Date("2020-01-01") + 0:1
  xl <- list(
    eleysurf = xts::xts(rbind(c(1, 2), c(4, 8)), order.by = dates),
    eleyunsat = xts::xts(rbind(c(2, 3), c(8, 15)), order.by = dates),
    eleygw = xts::xts(rbind(c(3, 4), c(15, 24)), order.by = dates),
    rivystage = xts::xts(rbind(1, 7), order.by = dates)
  )
  supplied_ic <- list(
    minit = data.frame(Surface = c(1, 2), Unsat = c(2, 3), GW = c(3, 4)),
    rinit = data.frame(Stage = 1)
  )
  fallback_ic <- list(
    minit = data.frame(Surface = c(100, 200), Unsat = c(100, 200), GW = c(100, 200)),
    rinit = data.frame(Stage = 100)
  )
  river <- SHUD.RIVER(
    river = data.frame(Type = 1, Length = 10),
    rivertype = data.frame(Width = 2),
    point = data.frame()
  )

  testthat::local_mocked_bindings(
    readic = function(...) fallback_ic,
    read_geol = function(...) data.frame(ThetaS.m3_m3. = 0.5),
    read_att = function(...) data.frame(GEOL = c(1, 1)),
    read_calib = function(...) data.frame(GEOL_THETAS = 1),
    read_river = function(...) river,
    readgeol = function(...) data.frame(ThetaS.m3_m3. = 0.5),
    readatt = function(...) data.frame(GEOL = c(1, 1)),
    readcalib = function(...) data.frame(GEOL_THETAS = 1),
    readriv = function(...) river,
    getArea = function(...) c(10, 10),
    .package = "rSHUD"
  )

  ds <- wb.DS(xl = xl, ic = supplied_ic)

  expect_equal(as.numeric(ds$Ele["ds.sf", ]), c(3, 6))
  expect_equal(as.numeric(ds$Ele["ds.us", ]), c(3, 6))
  expect_equal(as.numeric(ds$Ele["ds.gw", ]), c(6, 10))
  expect_equal(as.numeric(ds$Riv["ds.riv", ]), 6)
})

test_that("wb.riv uses modern read_river without deprecated self-call", {
  dates <- as.Date("2020-01-01") + 0:1
  xl <- list(
    rivqdown = xts::xts(matrix(c(10, 20), ncol = 1), order.by = dates),
    rivqsurf = xts::xts(matrix(c(1, 2), ncol = 1), order.by = dates),
    rivqsub = xts::xts(matrix(c(3, 4), ncol = 1), order.by = dates)
  )
  river <- SHUD.RIVER(
    river = data.frame(Down = -1, Type = 1, Length = 10),
    rivertype = data.frame(Width = 2),
    point = data.frame()
  )

  testthat::local_mocked_bindings(
    read_river = function(...) river,
    readriv = function(...) stop("deprecated readriv() called"),
    getArea = function(...) 10,
    .package = "rSHUD"
  )

  wb <- wb.riv(xl = xl, plot = FALSE)

  expect_equal(colnames(wb), c("DH", "Qout", "Qin_sf", "Qin_gw"))
  expect_equal(as.numeric(wb[, "Qout"]), c(1, 2))
})

test_that("wb.all uses modern read_river without deprecated self-call", {
  dates <- as.Date("2020-01-01") + 0:1
  xl <- list(
    elevprcp = xts::xts(matrix(c(1, 2, 3, 4), nrow = 2), order.by = dates),
    eleveta = xts::xts(matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2), order.by = dates),
    elevetp = xts::xts(matrix(c(0.5, 0.6, 0.7, 0.8), nrow = 2), order.by = dates),
    rivqdown = xts::xts(matrix(c(10, 20), ncol = 1), order.by = dates)
  )
  river <- SHUD.RIVER(
    river = data.frame(Down = -1, Type = 1, Length = 10),
    rivertype = data.frame(Width = 2),
    point = data.frame()
  )

  summarize_all <- function(x, FUN) {
    xts::xts(matrix(FUN(x), ncol = 1), order.by = as.Date("2020-01-01"))
  }

  testthat::local_mocked_bindings(
    read_river = function(...) river,
    readriv = function(...) stop("deprecated readriv() called"),
    getArea = function(...) c(10, 10),
    wb.DS = function(...) list(
      Ele = matrix(0, nrow = 3, ncol = 2),
      Riv = matrix(0, nrow = 1, ncol = 1)
    ),
    .package = "rSHUD"
  )

  wb <- wb.all(xl = xl, ic = list(), fun = summarize_all, plot = FALSE)

  expect_equal(colnames(wb), c("D_pqe", "P", "Q", "AET", "PET"))
  expect_equal(as.numeric(wb[, "Q"]), 1.5)
})

test_that("datafilter.riv uses modern readers without deprecated self-calls", {
  river <- SHUD.RIVER(
    river = data.frame(Type = c(1, 1)),
    rivertype = data.frame(Depth = 1),
    point = data.frame()
  )

  testthat::local_mocked_bindings(
    read_river = function(...) river,
    read_calib = function(...) c(RIV_DPTH = 0),
    readriv = function(...) stop("deprecated readriv() called"),
    readcalib = function(...) stop("deprecated readcalib() called"),
    .package = "rSHUD"
  )

  filtered <- datafilter.riv(
    matrix(c(0.1, 0.2, 1.1, 1.2), nrow = 2),
    plot = FALSE
  )

  expect_equal(filtered$ID, 2)
})
