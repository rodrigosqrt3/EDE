test_that("solow_roberts2003 implements the published p-value", {
  d <- sighting_data(data.frame(
    year = c(1900, 1905, 1910),
    sightings = 1
  ))
  curve <- solow_roberts2003(d, test_year = 1990, data_out = TRUE)
  expect_equal(curve$chance[curve$time == 1990], 5 / 85)
})

test_that("solow_roberts2003 reproduces the Mauritian orchid example", {
  # A. aff. rutenbergianum: last two collections 1976 and 1978, p = 0.10 in 1996
  orchid <- sighting_data(data.frame(
    year = c(1933, 1941, 1952, 1962, 1965, 1966, 1967, 1971, 1974, 1976, 1978),
    sightings = 1
  ))
  curve <- solow_roberts2003(orchid, test_year = 1996, data_out = TRUE)
  expect_equal(curve$chance[curve$time == 1996], 0.10)
  expect_equal(solow_roberts2003(orchid, test_year = 2100)$upper, 2016)
})

test_that("solow_roberts2003 returns the Robson-Whitlock upper bound", {
  d <- sighting_data(data.frame(
    year = c(1900, 1905, 1910),
    sightings = 1
  ))
  res <- solow_roberts2003(d, alpha = 0.05, test_year = 2010)
  expect_true(is.na(res$estimate))
  expect_equal(res$interval_type, "one-sided")
  expect_equal(res$upper, 2005)
  expect_equal(res$upper, robson1964(d)$upper)
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
  expect_true(is.na(res$upper))
})

test_that("solow_roberts2003 never tests beyond test_year", {
  d <- sighting_data(data.frame(year = c(1900, 1905, 1910), sightings = 1))
  curve <- solow_roberts2003(d, test_year = 1910.5, data_out = TRUE)
  expect_equal(nrow(curve), 0L)
})
