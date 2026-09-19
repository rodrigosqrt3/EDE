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

test_that("burgman1995 reproduces equation 4 values in Table 2", {
  seal <- sighting_data(data.frame(
    year = c(1915, 1922, 1932, 1948, 1952),
    sightings = 1
  ))
  seal_curve <- burgman1995(seal, test_year = 1992, data_out = TRUE)
  expect_equal(seal_curve$chance[seal_curve$time == 1992], 0.258, tolerance = 5e-4)
})

test_that("burgman1995 can reject persistence at the first candidate year", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  curve <- burgman1995(sighting_data(d), alpha = 0.05, test_year = 1932,
                       data_out = TRUE)
  if (curve$chance[1] <= 0.05) {
    expect_equal(burgman1995(sighting_data(d), alpha = 0.05,
                             test_year = 1932)$estimate, 1932)
  } else {
    expect_warning(burgman1995(sighting_data(d), alpha = 0.05,
                               test_year = 1932), "never falls")
  }
})

test_that("burgman internal helper functions handle edge cases", {
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 0, r = 2), 1.0)
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 5, r = 0), 1.0)
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 5, r = 10), 1.0)
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 5, r = 7), 0.00876, tolerance = 1e-5)
  d_zero <- data.frame(time = 1900:1905, count = rep(0, 6))
  expect_equal(EDE:::burgman_chance(d_zero), 1.0)
})
