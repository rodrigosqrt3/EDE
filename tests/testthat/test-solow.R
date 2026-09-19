test_that("solow1993 reproduces the reference point estimate", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  res <- solow1993(sighting_data(d), alpha = 0.05, test_year = 2000)
  expect_equal(res$estimate, 1944)
})

test_that("solow1993 reproduces the Caribbean monk seal example", {
  d <- sighting_data(data.frame(
    year = c(1915, 1922, 1932, 1948, 1952),
    sightings = 1
  ))
  curve <- solow1993(d, test_year = 1992, data_out = TRUE)
  # Solow (1993), equation 2: the first sighting defines the origin and is
  # omitted, leaving n = 4, t[n] = 37 and T = 77.
  expect_equal(curve$chance[curve$time == 1992], (37 / 77)^4)
})

test_that("solow2005 implements the Weibull p-value in equation 16", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  res <- solow2005(sighting_data(d), alpha = 0.05, test_year = 2000)
  expect_equal(res$estimate, 1939)

  curve <- solow2005(sighting_data(d), test_year = 1940, data_out = TRUE)
  times <- sort(rep(d$years, d$sightings), decreasing = TRUE)
  shape <- EDE:::ole_shape(times)
  expected <- exp(-length(times) * ((1940 - times[1]) /
                                      (1940 - times[length(times)]))^(1 / shape))
  expect_equal(curve$chance[curve$time == 1940], expected)
})

test_that("solow1993/solow2005 require test_year past the last sighting", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  expect_error(solow1993(sighting_data(d), 0.05, test_year = 1905), "later than")
  d3 <- data.frame(years = c(1900, 1905, 1910), sightings = c(1, 1, 1))
  expect_error(solow2005(sighting_data(d3), 0.05, test_year = 1905), "later than")
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
  sd3 <- sighting_data(data.frame(years = c(1900, 1905, 1910), sightings = 1))
  expect_error(solow2005(sd3, 0.05, test_year = 1910), "later than")
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

  sd3 <- sighting_data(data.frame(years = c(1900, 1905, 1910), sightings = 1))
  expect_warning(res2 <- solow2005(sd3, alpha = 0.001, test_year = 1911), "never falls")
  expect_true(is.na(res2$estimate))
})

test_that("solow2005 validates k and the number of events", {
  d <- sighting_data(data.frame(year = c(1900, 1905, 1910, 1915), sightings = 1))
  expect_error(solow2005(d, test_year = 2000, k = 2), "between 3")
  expect_error(solow2005(d, test_year = 2000, k = 5), "number of sighting events")

  d2 <- sighting_data(data.frame(year = c(1900, 1910), sightings = 1))
  expect_error(solow2005(d2, test_year = 2000), "at least 3")
})

test_that("solow2005 rejects tied most recent sighting events", {
  d <- sighting_data(data.frame(
    year = c(1900, 1910, 1920),
    sightings = c(1, 1, 2)
  ))
  expect_error(
    solow2005(d, test_year = 2000),
    "strictly greater than the second most recent"
  )
})

test_that("solow1993 internal helper functions and edge cases handle single sightings", {
  d_single <- data.frame(time = 1900, count = 1)
  expect_error(solow1993(sighting_data(d_single), alpha = 0.05, test_year = 2000), "at least 2")
  expect_equal(EDE:::solow1993_chance(d_single), 1.0)
})
