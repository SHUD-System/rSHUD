# Tests for PET calculations

test_that("PET_PM returns positive finite scalar PET", {
  pet <- suppressMessages(PET_PM(
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

test_that("PET_PM returns positive finite vector PET", {
  wind <- c(1, 1.5, 2)
  pet <- suppressMessages(PET_PM(
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
  pet <- suppressMessages(PET_PM(
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
  pet <- suppressMessages(PET_PM(
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
    suppressMessages(PET_PM(
      Wind = c(1, 1.5, 2),
      Temp = c(20, 25),
      RH = 0.5,
      RadNet = 400,
      Press = 101.3
    )),
    "length 1 or share one common length"
  )

  expect_error(
    suppressMessages(PET_PM(
      Wind = c(1, 1.5, 2),
      Temp = c(20, 25, 30),
      RH = 0.5,
      RadNet = 400,
      Press = NULL,
      Elevation = c(100, 200)
    )),
    "length 1 or share one common length"
  )
})
