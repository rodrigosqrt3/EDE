seal <- function() {
  sighting_data(data.frame(year = c(1915, 1922, 1932, 1948, 1952), sightings = 1))
}

test_that("solow1993_bayes reproduces the monk seal example (Solow 1993a)", {
  res <- solow1993_bayes(seal(), test_year = 1992)
  # n = 4, t[n] = 37, T = 77: B = 0.37 and p(H0 | t) = 0.27 for prior 0.5
  expect_equal(res$bayes_factor, 3 / ((77 / 37)^3 - 1))
  expect_equal(round(res$bayes_factor, 2), 0.37)
  expect_equal(round(res$posterior, 2), 0.27)
  expect_s3_class(res, "ede_bayes")
})

test_that("solow1993_bayes uses the n = 1 limit", {
  d <- sighting_data(data.frame(year = c(1900, 1910), sightings = 1))
  res <- solow1993_bayes(d, test_year = 1930)
  expect_equal(res$bayes_factor, 1 / log(3))
  near_one <- (1 + 1e-6 - 1) / (3^(1e-6) - 1)
  expect_equal(res$bayes_factor, near_one, tolerance = 1e-5)
})

test_that("solow1993_bayes posterior follows the prior", {
  res <- solow1993_bayes(seal(), test_year = 1992, prior = 0.8)
  b <- res$bayes_factor
  expect_equal(res$posterior, 1 / (1 + 0.2 / (0.8 * b)))
})

test_that("solow1993_bayes data_out returns a decreasing curve", {
  curve <- solow1993_bayes(seal(), test_year = 2000, data_out = TRUE)
  expect_named(curve, c("time", "bayes_factor", "posterior"))
  expect_equal(curve$time, 1953:2000)
  expect_true(all(diff(curve$bayes_factor) < 0))
  expect_equal(curve$bayes_factor[curve$time == 1992],
               solow1993_bayes(seal(), test_year = 1992)$bayes_factor)
})

test_that("solow1993_bayes validates its inputs", {
  expect_error(solow1993_bayes(seal(), test_year = 1992, prior = 1), "prior")
  expect_error(solow1993_bayes(seal()), "test_year")
  expect_error(solow1993_bayes(seal(), test_year = 1950), "later than")
  one <- sighting_data(data.frame(year = 1900, sightings = 1))
  expect_error(solow1993_bayes(one, test_year = 2000), "at least 1")
})

test_that("print.ede_bayes shows the Bayes factor and posterior", {
  out <- capture.output(print(solow1993_bayes(seal(), test_year = 1992)))
  expect_match(out[1], "Bayes factor")
  expect_match(out[2], "0\\.374")
  expect_match(out[3], "posterior")
})
