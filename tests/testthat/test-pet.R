# Tests for PET calculations

test_that("PET_PM returns positive finite scalar PET", {
  expect_silent(pet <- PET_PM(
    Wind = 1,
    Temp = 25,
    RH = 0.5,
    RadNet = 400,
    Press = NULL,
    Elevation = 100
  ))

  expect_true(is.numeric(pet))
  expect_length(pet, 1)
  expect_true(is.finite(pet))
  expect_gt(pet, 0)
})

test_that("PET_PM matches FAO-56 daily reference ET with W to MJ conversion", {
  wind <- 1.8
  temp <- 24
  rh <- 0.55
  radnet <- 360
  elevation <- 250
  wind_height <- 10

  press <- 101.325 * ((293 - 0.0065 * elevation) / 293)^5.26
  lambda <- 2.501 - 0.002361 * temp
  gamma <- 0.0016286 * press / lambda
  es <- 0.6108 * exp(17.27 * temp / (temp + 237.3))
  ed <- es * (1 - rh)
  delta <- 4098 * es / (temp + 237.3)^2
  u2 <- wind * 4.87 / log(67.8 * wind_height - 5.42)
  gamma_wind <- gamma * (1 + 0.33 * u2)
  rg_mj_day <- radnet * 0.9 * 86400 * 1e-6
  expected <- (0.408 * delta * rg_mj_day + gamma * 900 / (temp + 273) * u2 * ed) /
    (delta + gamma_wind)

  expect_silent(pet <- PET_PM(
    Wind = wind,
    Temp = temp,
    RH = rh,
    RadNet = radnet,
    Press = NULL,
    Elevation = elevation
  ))

  expect_equal(pet, expected, tolerance = 1e-12)
  expect_gt(abs(pet - expected / (86400 * 1e-6)), 1)
})

test_that("PET_PM returns positive finite vector PET", {
  wind <- c(1, 1.5, 2)
  expect_silent(pet <- PET_PM(
    Wind = wind,
    Temp = c(20, 25, 30),
    RH = c(0.6, 0.5, 0.4),
    RadNet = c(300, 400, 500),
    Press = NULL,
    Elevation = 100
  ))

  expect_true(is.numeric(pet))
  expect_length(pet, length(wind))
  expect_true(all(is.finite(pet)))
  expect_true(all(pet > 0))
})

test_that("PET_PM broadcasts scalar inputs across vector meteorology", {
  temp <- c(20, 25, 30)
  expect_silent(pet <- PET_PM(
    Wind = 1,
    Temp = temp,
    RH = 0.5,
    RadNet = 400,
    Press = 101.3
  ))

  expect_true(is.numeric(pet))
  expect_length(pet, length(temp))
  expect_true(all(is.finite(pet)))
})

test_that("PET_PM accepts matching vector elevation when pressure is computed", {
  elevation <- c(0, 100, 200)
  expect_silent(pet <- PET_PM(
    Wind = c(1, 1.5, 2),
    Temp = c(20, 25, 30),
    RH = 0.5,
    RadNet = 400,
    Press = NULL,
    Elevation = elevation
  ))

  expect_true(is.numeric(pet))
  expect_length(pet, length(elevation))
  expect_true(all(is.finite(pet)))
})

test_that("PET_PM rejects incompatible vector lengths", {
  expect_error(
    PET_PM(
      Wind = c(1, 1.5, 2),
      Temp = c(20, 25),
      RH = 0.5,
      RadNet = 400,
      Press = 101.3
    ),
    "length 1 or share one common length"
  )

  expect_error(
    PET_PM(
      Wind = c(1, 1.5, 2),
      Temp = c(20, 25, 30),
      RH = 0.5,
      RadNet = 400,
      Press = NULL,
      Elevation = c(100, 200)
    ),
    "length 1 or share one common length"
  )
})
