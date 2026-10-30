test_that("ole keeps its point estimate and confidence interval", {
  # regression values from the package on a synthetic record; the published
  # examples are tested below
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

test_that("ole reproduces the dodo example (Roberts & Solow 2003)", {
  d <- data.frame(
    year = c(1598, 1601, 1602, 1607, 1611, 1628, 1631, 1638, 1662),
    sightings = c(1, 1, 1, 1, 1, 2, 1, 1, 1)
  )
  res <- ole(sighting_data(d))
  # published: 1690 (95% CI 1669-1797)
  expect_equal(round(res$estimate), 1690)
  expect_equal(round(res$lower), 1669)
  # the paper's 1797 follows from the shape estimate rounded to 0.39; the
  # printed formula gives 0.385 and an upper limit of 1798.9
  expect_lt(abs(res$upper - 1798.88), 0.01)
})

test_that("ole validates k", {
  sd <- sighting_data(data.frame(year = c(1900, 1905, 1910, 1915), sightings = 1))
  expect_error(ole(sd, k = 2), "between 3")
  expect_error(ole(sd, k = 5), "number of sighting events")
})

test_that("ole returns Inf when the upper confidence limit does not exist", {
  sd <- sighting_data(data.frame(year = c(1900, 1905, 1910, 1915), sightings = 1))
  # -log(0.025) > 3, so the upper limit does not exist for k = 3
  expect_warning(res <- ole(sd, k = 3), "does not exist")
  expect_equal(res$upper, Inf)
  expect_gt(res$lower, 1915)
  expect_gt(res$estimate, 1915)
})

test_that("ole_shape handles an insufficient number of times", {
  expect_true(is.na(EDE:::ole_shape(c(1900, 1910))))
})
