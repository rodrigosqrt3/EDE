test_that("solow_roberts2003 implements the published p-value", {
  d <- sighting_data(data.frame(
    year = c(1900, 1905, 1910),
    sightings = 1
  ))
  curve <- solow_roberts2003(d, test_year = 1990, data_out = TRUE)
  expect_equal(curve$chance[curve$time == 1990], 5 / 85)
})

test_that("solow_roberts2003 returns the first alpha crossing", {
  d <- sighting_data(data.frame(
    year = c(1900, 1905, 1910),
    sightings = 1
  ))
  res <- solow_roberts2003(d, alpha = 0.05, test_year = 2010)
  expect_equal(res$estimate, 2005)
})

test_that("solow_roberts2003 validates its inputs", {
  d <- sighting_data(data.frame(year = c(1900, 1910), sightings = 1))
  expect_error(solow_roberts2003(d, alpha = 0, test_year = 2000), "must be in")
  expect_error(solow_roberts2003(d), "test_year")
  expect_error(solow_roberts2003(d, test_year = 1910), "later than")

  one <- sighting_data(data.frame(year = 1900, sightings = 1))
  expect_error(solow_roberts2003(one, test_year = 2000), "at least 2")
})

test_that("solow_roberts2003 warns when alpha is not reached", {
  d <- sighting_data(data.frame(year = c(1900, 1905, 1910), sightings = 1))
  expect_warning(
    res <- solow_roberts2003(d, alpha = 0.01, test_year = 1911),
    "never falls"
  )
  expect_true(is.na(res$estimate))
})
