test_that("ole reproduces the reference point estimate and confidence interval", {
  # reference values from the closed-form OLE formula (Roberts & Solow 2003),
  # worked example, alpha = 0.05
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  sd <- sighting_data(d)
  res <- ole(sd, alpha = 0.05)

  expect_equal(res$estimate, 1934.661, tolerance = 1e-3)
  expect_equal(res$lower, 1931.13, tolerance = 1e-2)
  expect_equal(res$upper, 1954.551, tolerance = 1e-2)
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
