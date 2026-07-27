test_that("solow1993 reproduces the reference point estimate", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  res <- solow1993(sighting_data(d), alpha = 0.05, test_year = 2000)
  expect_equal(res$estimate, 1936)
})

test_that("solow2005 reproduces the reference point estimate", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  res <- solow2005(sighting_data(d), alpha = 0.05, test_year = 2000)
  expect_equal(res$estimate, 1941)
})

test_that("solow1993/solow2005 require test_year past the last sighting", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  expect_error(solow1993(sighting_data(d), 0.05, test_year = 1905), "later than")
  expect_error(solow2005(sighting_data(d), 0.05, test_year = 1905), "later than")
})

test_that("solow1993/solow2005 validate alpha", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  sd <- sighting_data(d)
  expect_error(solow1993(sd, alpha = 0, test_year = 2000), "must be in \\(0, 1\\)")
  expect_error(solow1993(sd, alpha = 1, test_year = 2000), "must be in \\(0, 1\\)")
  expect_error(solow2005(sd, alpha = 0, test_year = 2000), "must be in \\(0, 1\\)")
  expect_error(solow2005(sd, alpha = 1, test_year = 2000), "must be in \\(0, 1\\)")
})

test_that("solow1993/solow2005 require test_year to be supplied and numeric", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  sd <- sighting_data(d)
  expect_error(solow1993(sd, 0.05), "test_year")
  expect_error(solow1993(sd, 0.05, test_year = "2000"), "test_year")
  expect_error(solow2005(sd, 0.05), "test_year")
  expect_error(solow2005(sd, 0.05, test_year = "2000"), "test_year")
})

test_that("solow1993/solow2005 reject test_year equal to the last sighting", {
  # expand_record() allows end_time == last sighting, but there is then no
  # candidate year strictly after it, which must surface as this error
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  sd <- sighting_data(d)
  expect_error(solow1993(sd, 0.05, test_year = 1910), "later than")
  expect_error(solow2005(sd, 0.05, test_year = 1910), "later than")
})

test_that("solow1993/solow2005 data_out = TRUE returns the full chance curve", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  sd <- sighting_data(d)

  curve1 <- solow1993(sd, 0.05, test_year = 2000, data_out = TRUE)
  expect_named(curve1, c("time", "chance"))
  expect_equal(nrow(curve1), 69L)

  curve2 <- solow2005(sd, 0.05, test_year = 2000, data_out = TRUE)
  expect_named(curve2, c("time", "chance"))
  expect_equal(nrow(curve2), 69L)
})

test_that("solow1993/solow2005 warn and return NA when persistence is never rejected", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  sd <- sighting_data(d)

  expect_warning(res1 <- solow1993(sd, alpha = 0.5, test_year = 1911), "never falls")
  expect_true(is.na(res1$estimate))

  expect_warning(res2 <- solow2005(sd, alpha = 0.5, test_year = 1911), "never falls")
  expect_true(is.na(res2$estimate))
})

test_that("solow2005_survival_series handles y >= 1 (empty inclusion-exclusion sum)", {
  expect_equal(solow2005_survival_series(1.5, 3), 1)
  expect_equal(solow2005_survival_series(0.3, 3), 0, tolerance = 1e-10)
})
