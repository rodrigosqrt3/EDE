test_that("ole reproduces the reference point estimate and confidence interval", {
  # reference values from the closed-form OLE formula (Roberts & Solow 2003),
  # worked example, alpha = 0.05
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  sd <- sighting_data(d)
  res <- ole(sd, alpha = 0.05)

  expect_equal(res$estimate, 1933.019, tolerance = 1e-3)
  expect_equal(res$lower, 1931.079, tolerance = 1e-3)
  expect_equal(res$upper, 1939.964, tolerance = 1e-3)
})

test_that("ole requires at least 3 nonzero sightings", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  expect_error(ole(sighting_data(d), 0.05), "at least 3")
})

test_that("ole validates alpha", {
  d <- data.frame(years = c(1900, 1905, 1910), sightings = c(1, 1, 1))
  sd <- sighting_data(d)
  expect_error(ole(sd, alpha = 0), "must be in \\(0, 1\\)")
  expect_error(ole(sd, alpha = 1), "must be in \\(0, 1\\)")
})

test_that("ole rejects datasets where T1 equals T2", {
  d_dup <- data.frame(time = c(1900, 1910, 1920), count = c(1, 1, 2))
  expect_error(ole(sighting_data(d_dup)), "strictly greater")
})

test_that("ole reproduces the dodo point estimate and unrounded CI formula", {
  d <- data.frame(
    year = c(1598, 1601, 1602, 1607, 1611, 1628, 1631, 1638, 1662),
    sightings = c(1, 1, 1, 1, 1, 2, 1, 1, 1)
  )
  res <- ole(sighting_data(d))
  expect_equal(res$estimate, 1690.42, tolerance = 0.02)
  expect_equal(res$lower, 1669.13, tolerance = 0.02)
  expect_equal(res$upper, 1798.88, tolerance = 0.02)
})

test_that("ole validates k", {
  sd <- sighting_data(data.frame(year = c(1900, 1905, 1910, 1915), sightings = 1))
  expect_error(ole(sd, k = 2), "between 3")
  expect_error(ole(sd, k = 5), "number of sighting events")
  expect_s3_class(ole(sd, k = 3), "ede_estimate")
})

test_that("ole_shape handles an insufficient number of times", {
  expect_true(is.na(EDE:::ole_shape(c(1900, 1910))))
})
