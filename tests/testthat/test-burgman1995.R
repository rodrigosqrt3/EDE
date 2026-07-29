test_that("burgman1995 reproduces the reference point estimate", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  res <- burgman1995(sighting_data(d), alpha = 0.05, test_year = 1945)
  expect_equal(res$estimate, 1942)
})

test_that("burgman1995 validates alpha and test_year", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  sd <- sighting_data(d)
  expect_error(burgman1995(sd, alpha = 0, test_year = 2000), "must be in \\(0, 1\\)")
  expect_error(burgman1995(sd, alpha = 1, test_year = 2000), "must be in \\(0, 1\\)")
  expect_error(burgman1995(sd, alpha = 0.05), "test_year")
  expect_error(burgman1995(sd, alpha = 0.05, test_year = "2000"), "test_year")
})

test_that("burgman1995 rejects test_year equal to the last sighting", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  expect_error(burgman1995(sighting_data(d), 0.05, test_year = 1931), "later than")
})

test_that("burgman1995 data_out = TRUE returns the full chance curve", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  curve <- burgman1995(sighting_data(d), alpha = 0.05, test_year = 1935, data_out = TRUE)
  expect_named(curve, c("time", "chance"))
  expect_equal(nrow(curve), 4L)
})

test_that("burgman1995 warns and returns NA with a single candidate year", {
  # with only one candidate year, diff(chance) has length 0 and the leading
  # `decreasing` flag is hardcoded FALSE, so rejection is never possible
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  expect_warning(
    res <- burgman1995(sighting_data(d), alpha = 0.05, test_year = 1932),
    "never falls"
  )
  expect_true(is.na(res$estimate))
})

test_that("burgman internal helper functions handle edge cases", {
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 0, r = 2), 1.0)
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 5, r = 0), 1.0)
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 5, r = 10), 1.0)
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 5, r = 7), 0.00876, tolerance = 1e-5)
  d_zero <- data.frame(time = 1900:1905, count = rep(0, 6))
  expect_equal(EDE:::burgman_chance(d_zero), 1.0)
})
