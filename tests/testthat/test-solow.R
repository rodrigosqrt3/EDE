test_that("solow1993 returns the inverted-test upper bound", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  res <- solow1993(sighting_data(d), alpha = 0.05, test_year = 2000)
  # one sighting at the 1907 origin is omitted: n = 15 events, t[n] = 24
  expect_true(is.na(res$estimate))
  expect_equal(res$interval_type, "one-sided")
  expect_equal(res$upper, 1907 + 24 / 0.05^(1 / 15))
})

test_that("solow1993 reproduces the Caribbean monk seal example", {
  d <- sighting_data(data.frame(
    year = c(1915, 1922, 1932, 1948, 1952),
    sightings = 1
  ))
  curve <- solow1993(d, test_year = 1992, data_out = TRUE)
  # Solow (1993a), equation 2: the first sighting defines the origin and is
  # omitted, leaving n = 4, t[n] = 37 and T = 77.
  expect_equal(curve$chance[curve$time == 1992], (37 / 77)^4)
  expect_equal(round(curve$chance[curve$time == 1992], 3), 0.053)
  expect_equal(solow1993(d, test_year = 2100)$upper, 1915 + 37 / 0.05^(1 / 4))
})

test_that("solow1993 counts every sighting event", {
  single <- sighting_data(data.frame(year = c(1900, 1910, 1920), sightings = 1))
  double <- sighting_data(data.frame(year = c(1900, 1910, 1920), sightings = c(1, 2, 2)))
  p1 <- solow1993(single, test_year = 1930, data_out = TRUE)$chance
  p2 <- solow1993(double, test_year = 1930, data_out = TRUE)$chance
  expect_equal(p1[10], (20 / 30)^2)
  expect_equal(p2[10], (20 / 30)^4)

  # extra sightings in the origin time unit are events, the first is not
  origin_tie <- sighting_data(data.frame(year = c(1900, 1920), sightings = c(2, 1)))
  p3 <- solow1993(origin_tie, test_year = 1930, data_out = TRUE)$chance
  expect_equal(p3[10], (20 / 30)^2)
})

test_that("solow2005 implements the Weibull p-value in equation 16", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  curve <- solow2005(sighting_data(d), test_year = 1940, data_out = TRUE)
  times <- sort(rep(d$years, d$sightings), decreasing = TRUE)
  shape <- EDE:::ole_shape(times)
  expected <- exp(-length(times) * ((1940 - times[1]) /
                                      (1940 - times[length(times)]))^(1 / shape))
  expect_equal(curve$chance[curve$time == 1940], expected)
})

test_that("solow2005 reproduces the monk seal example (Section 5)", {
  seal <- sighting_data(data.frame(year = c(1915, 1922, 1932, 1948, 1952), sightings = 1))
  times <- c(1952, 1948, 1932, 1922, 1915)
  expect_equal(round(EDE:::ole_shape(times), 2), 0.76)
  curve <- solow2005(seal, test_year = 2004, data_out = TRUE)
  expect_lt(abs(curve$chance[curve$time == 2004] - 0.085), 1e-3)
  expect_lt(abs(ole(seal)$estimate - 1964), 1)
  # the paper reports an upper bound "around 2027" with the rounded m = 0.76
  expect_lt(abs(solow2005(seal, test_year = 2100)$upper - 2027), 3)
})

test_that("solow2005 upper bound is equation 17 and matches ole() at 2 * alpha", {
  d <- sighting_data(data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  ))
  res <- solow2005(d, alpha = 0.05, test_year = 2000)
  times <- sort(rep(d$time, d$count), decreasing = TRUE)
  k <- length(times)
  c_alpha <- (-log(0.05) / k)^EDE:::ole_shape(times)
  expect_true(is.na(res$estimate))
  expect_equal(res$upper, (times[1] - c_alpha * times[k]) / (1 - c_alpha))
  expect_equal(res$upper, ole(d, alpha = 0.10)$upper)
})

test_that("solow2005 reports a missing bound when -log(alpha) >= k", {
  d <- sighting_data(data.frame(year = c(1900, 1905, 1910, 1915), sightings = 1))
  expect_warning(res <- solow2005(d, alpha = 0.01, test_year = 1e6, k = 3), "never falls")
  expect_true(is.na(res$upper))
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
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  sd <- sighting_data(d)
  expect_error(solow1993(sd, 0.05, test_year = 1910), "later than")
  sd3 <- sighting_data(data.frame(years = c(1900, 1905, 1910), sightings = 1))
  expect_error(solow2005(sd3, 0.05, test_year = 1910), "later than")
})

test_that("solow1993/solow2005 data_out = TRUE returns the full p-value curve", {
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
  expect_true(is.na(res1$upper))

  sd3 <- sighting_data(data.frame(years = c(1900, 1905, 1910), sightings = 1))
  expect_warning(res2 <- solow2005(sd3, alpha = 0.001, test_year = 1911), "never falls")
  expect_true(is.na(res2$upper))
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

test_that("solow1993 needs a sighting after the observation origin", {
  d_single <- data.frame(time = 1900, count = 1)
  expect_error(solow1993(sighting_data(d_single), alpha = 0.05, test_year = 2000), "at least 1")
  expect_equal(EDE:::solow1993_chance(d_single), 1.0)

  d_zero <- sighting_data(data.frame(a = c(1900, 1910), b = c(0, 0)))
  expect_error(solow1993(d_zero, 0.05, test_year = 1920), "no sightings")
})
