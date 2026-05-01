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
