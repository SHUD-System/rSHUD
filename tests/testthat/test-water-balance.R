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
