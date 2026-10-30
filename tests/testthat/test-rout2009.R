# Helenium amarum (Rout et al. 2009, Case study): S = 178 surveys, last
# sighting at survey s[n] = 169, seen in n = 142 surveys. Only these summaries
# are published, and the Bayes factor depends on the record only through
# them, so the record is rebuilt with survey 0 as the origin and the 142
# sightings in surveys 28-169.
amarum <- function() {
  sighting_data(data.frame(survey = c(0, 28:169), seen = c(0, rep(1, 142))))
}

test_that("rout2009 reproduces the H. amarum stopping times (Rout et al. 2009)", {
  # Results: rule of thumb with R = 0.00282 and prior 0.5 gives d* = 13 for
  # the constant rate model (Fig. 2) and d* = 14 with a = 0.1
  # ("Declining sighting rate model"); Fig. 6 / text: d* = 23 when a = 0.5.
  expect_equal(rout2009_absences(amarum(), cost_ratio = 0.00282)$optimal_absent, 13)
  expect_equal(rout2009_absences(amarum(), cost_ratio = 0.00282, a = 0.1)$optimal_absent, 14)
  expect_equal(rout2009_absences(amarum(), cost_ratio = 0.00282, a = 0.5)$optimal_absent, 23)
  # Approximation (equation 5): d* = 13.15 (Results; Fig. S3)
  res <- rout2009_absences(amarum(), cost_ratio = 0.00282)
  expect_equal(round(res$approximation, 2), 13.15)
  expect_equal(res$stop_time, 169 + 13)
  expect_equal(res$n, 142)
  expect_equal(res$sn, 169)
  expect_true(is.na(rout2009_absences(amarum(), cost_ratio = 0.00282, a = 0.1)$approximation))
})

test_that("rout2009 has the Bayes factor of Rout et al. (2009)", {
  res <- rout2009(amarum(), test_year = 178, a = 0.1)
  k <- 142 * 0.9
  expect_equal(res$bayes_factor, (k - 1) / ((178 / 169)^(k - 1) - 1))
  expect_equal(res$posterior, 1 / (1 + 1 / res$bayes_factor))
  expect_s3_class(res, "ede_bayes")
  expect_named(res, c("bayes_factor", "posterior", "prior", "test_year", "method", "a"))
})

test_that("rout2009 with a = 0 equals solow1993_bayes (Rout et al. equation 2)", {
  seal <- sighting_data(data.frame(year = c(1915, 1922, 1932, 1948, 1952), count = 1))
  r <- rout2009(seal, test_year = 1992)
  s <- solow1993_bayes(seal, test_year = 1992)
  expect_equal(r$bayes_factor, s$bayes_factor)
  expect_equal(r$posterior, s$posterior)
  expect_equal(rout2009(seal, test_year = 2000, data_out = TRUE)$bayes_factor,
               solow1993_bayes(seal, test_year = 2000, data_out = TRUE)$bayes_factor)
})

test_that("a declining rate gives a higher probability of persistence (Fig. 5)", {
  const <- rout2009(amarum(), test_year = 189, data_out = TRUE)
  decl <- rout2009(amarum(), test_year = 189, a = 0.1, data_out = TRUE)
  expect_equal(const$time, 170:189)
  expect_true(all(decl$posterior > const$posterior))
  expect_true(all(diff(decl$posterior) < 0))
})

test_that("rout2009 uses the limit when n(1 - a) = 1", {
  d <- sighting_data(data.frame(t = c(0, 5, 10), count = c(0, 1, 1)))
  res <- rout2009(d, test_year = 30, a = 0.5)
  expect_equal(res$bayes_factor, 1 / log(3))
  near <- rout2009(d, test_year = 30, a = 0.5 - 1e-7)$bayes_factor
  expect_lt(abs(near - 1 / log(3)), 1e-5)
})

test_that("rout2009 and rout2009_absences validate their inputs", {
  expect_error(rout2009(amarum(), test_year = 178, a = 1), "`a`")
  expect_error(rout2009(amarum(), test_year = 178, a = -0.1), "`a`")
  expect_error(rout2009(amarum(), test_year = 178, prior = 0), "prior")
  expect_error(rout2009(amarum()), "test_year")
  expect_error(rout2009(amarum(), test_year = 160), "later than")
  expect_error(rout2009_absences(amarum()), "cost_ratio")
  expect_error(rout2009_absences(amarum(), cost_ratio = -1), "cost_ratio")
  expect_error(rout2009_absences(amarum(), cost_ratio = 0.01, max_absent = 0), "max_absent")
  one <- sighting_data(data.frame(t = 1, count = 1))
  expect_error(rout2009(one, test_year = 5), "at least 1")
})

test_that("print.ede_bayes works for rout2009", {
  out <- capture.output(print(rout2009(amarum(), test_year = 178, a = 0.1)))
  expect_match(out[1], "Rout et al. \\(2009\\)")
  expect_match(out[3], "posterior")
})
