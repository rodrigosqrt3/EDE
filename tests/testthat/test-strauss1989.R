test_that("strauss1989 returns the unbiased point estimate and confidence bound", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  res <- strauss1989(sighting_data(d), alpha = 0.05)
  expect_equal(res$estimate, 1934.429, tolerance = 1e-3)
  expect_equal(res$lower, 1931)
  expect_equal(res$upper, 1943.819, tolerance = 1e-3)
})

test_that("strauss1989 reproduces Table 1 (one-parameter case)", {
  # c = (z_upper - z) / range for 95% confidence, from the paper's Table 1
  table1 <- c(`2` = 19.000, `5` = 1.115, `8` = 0.534, `16` = 0.221, `50` = 0.063)
  for (n in as.integer(names(table1))) {
    sd <- sighting_data(data.frame(time = seq(0, 1, length.out = n), count = 1))
    expect_lt(abs(strauss1989(sd)$upper - 1 - table1[[as.character(n)]]), 1e-3)
  }
})

test_that("strauss1989 reproduces the Section 5 example", {
  # taxon 5: 16 finds from 800 to 1147; taxon 10: 8 finds from 979 to 1091
  taxon5 <- sighting_data(data.frame(time = seq(800, 1147, length.out = 16), count = 1))
  taxon10 <- sighting_data(data.frame(time = seq(979, 1091, length.out = 8), count = 1))
  expect_equal(round(strauss1989(taxon5)$estimate), 1170)
  expect_equal(round(strauss1989(taxon10)$estimate), 1107)
})

test_that("strauss1989 validates alpha and minimum n", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  sd <- sighting_data(d)
  expect_error(strauss1989(sd, alpha = 0), "must be in \\(0, 1\\)")
  expect_error(strauss1989(sd, alpha = 1), "must be in \\(0, 1\\)")

  d1 <- data.frame(years = 1900, sightings = 1)
  expect_error(strauss1989(sighting_data(d1), 0.05), "at least 2")
})

test_that("strauss1989_curve requires at least 2 sightings", {
  d1 <- data.frame(years = 1900, sightings = 1)
  expect_error(strauss1989_curve(sighting_data(d1)), "at least 2")
})

test_that("strauss1989_curve reproduces the upper bound at matching alpha", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  sd <- sighting_data(d)
  curve <- strauss1989_curve(sd)

  expect_equal(nrow(curve), 100L)
  expect_named(curve, c("time", "chance"))

  row95 <- curve[abs(curve$chance - 0.95) < 1e-9, ]
  pt <- strauss1989(sd, alpha = 0.05)
  expect_equal(row95$time, pt$upper, tolerance = 1e-3)
})

ammonite_taxon <- function(n, from, to) {
  # Strauss & Sadler (1989), Section 5: only n and the range are given, so
  # the finds are spaced evenly (the Bayesian results depend only on these)
  sighting_data(data.frame(time = seq(from, to, length.out = n), count = 1))
}

test_that("strauss1989_bayes reproduces the posterior mean difference of Fig. 4", {
  # taxa 5 (16 finds, 800-1147 m) and 10 (8 finds, 979-1091 m), prior range
  # (0, 1500) m: posterior mean of theta2 - theta2' is 60 m (p. 425)
  t5 <- strauss1989_bayes(ammonite_taxon(16, 800, 1147), prior_range = c(0, 1500))
  t10 <- strauss1989_bayes(ammonite_taxon(8, 979, 1091), prior_range = c(0, 1500))
  expect_equal(round(t5$estimate - t10$estimate), 60)
  expect_lt(abs(t5$estimate - 1173.672), 1e-3)
  expect_lt(abs(t10$estimate - 1113.351), 1e-3)

  # 95% prediction interval for the difference by numerical convolution.
  # The paper reports -24 to 150 m; the equal-tailed interval from the
  # posteriors of eq. 26 is -20.3 to 144.8 m (the paper's integration is
  # presumably coarser), so only the package value is checked here.
  cdf_diff <- function(d) {
    stats::integrate(function(t) t5$density(t) * (1 - t10$cdf(t - d)),
                     1147, 1500, subdivisions = 2000L, rel.tol = 1e-10)$value
  }
  q <- function(p) stats::uniroot(function(d) cdf_diff(d) - p, c(-353, 409), tol = 1e-8)$root
  expect_lt(abs(q(0.025) - (-20.29)), 0.05)
  expect_lt(abs(q(0.975) - 144.76), 0.05)
})

test_that("strauss1989_bayes matches equations 25-27 and numerical integration", {
  sd <- ammonite_taxon(16, 800, 1147)
  res <- strauss1989_bayes(sd, alpha = 0.05, prior_range = c(0, 1500))
  # eq. 27 on the paper's (0, 1) scale
  n <- 16; y <- 800 / 1500; z <- 1147 / 1500
  u <- function(k) (z - y)^(-k + 2) - (1 - y)^(-k + 2) - z^(-k + 2) + 1
  eq27 <- (n - 2) * u(n - 1) / ((n - 3) * u(n)) +
    y * ((z - y)^(-n + 2) - (1 - y)^(-n + 2)) / u(n)
  expect_lt(abs(res$estimate - 1500 * eq27), 1e-8)
  # eq. 26 on the (0, 1) scale, transformed back
  t <- c(1150, 1200, 1300)
  h26 <- (n - 2) * ((t / 1500 - y)^(-n + 1) - (t / 1500)^(-n + 1)) / u(n) / 1500
  expect_lt(max(abs(res$density(t) - h26)), 1e-12)

  expect_lt(abs(stats::integrate(res$density, 1147, 1500)$value - 1), 1e-8)
  expect_lt(abs(stats::integrate(function(t) t * res$density(t), 1147, 1500)$value -
                  res$estimate), 1e-6)
  expect_equal(res$lower, 1147)
  expect_lt(abs(res$cdf(res$upper) - 0.95), 1e-8)
  expect_equal(res$interval_type, "one-sided")
  expect_equal(res$density(c(1000, 1600)), c(0, 0))
})

test_that("strauss1989_bayes handles n = 3 and an unbounded prior range", {
  sd3 <- sighting_data(data.frame(time = c(0, 4, 10), count = 1))
  r3 <- strauss1989_bayes(sd3, prior_range = c(-5, 40))
  m3 <- stats::integrate(function(t) t * r3$density(t), 10, 40)$value
  expect_lt(abs(r3$estimate - m3), 1e-6)
  expect_equal(strauss1989_bayes(sd3, prior_range = c(-5, Inf))$estimate, Inf)

  sd <- sighting_data(data.frame(time = c(1, 2, 3, 7, 10), count = 1))
  r_inf <- strauss1989_bayes(sd, prior_range = c(0, Inf))
  m <- stats::integrate(function(t) t * r_inf$density(t), 10, Inf)$value
  expect_lt(abs(r_inf$estimate - m), 1e-6)
  expect_lt(abs(r_inf$cdf(r_inf$upper) - 0.95), 1e-8)
  # a wider prior range moves the posterior to later times
  r_narrow <- strauss1989_bayes(sd, prior_range = c(0, 20))
  expect_lt(r_narrow$upper, r_inf$upper)
})

test_that("strauss1989_bayes ignores counts and returns the posterior curve", {
  sd <- sighting_data(data.frame(time = c(1, 3, 6, 9), count = c(1, 4, 2, 1)))
  sd1 <- sighting_data(data.frame(time = c(1, 3, 6, 9), count = 1))
  expect_equal(strauss1989_bayes(sd, prior_range = c(0, 30))$estimate,
               strauss1989_bayes(sd1, prior_range = c(0, 30))$estimate)
  curve <- strauss1989_bayes(sd, prior_range = c(0, 30), data_out = TRUE)
  expect_named(curve, c("time", "density", "chance"))
  expect_equal(nrow(curve), 200L)
  expect_equal(curve$chance[1], 0)
  expect_lt(abs(curve$chance[200] - 0.999), 1e-8)
  expect_true(all(diff(curve$density) <= 0))
})

test_that("strauss1989_bayes validates its inputs", {
  sd <- sighting_data(data.frame(time = c(1, 3, 6, 9), count = 1))
  expect_error(strauss1989_bayes(sd), "prior_range")
  expect_error(strauss1989_bayes(sd, prior_range = 0), "prior_range")
  expect_error(strauss1989_bayes(sd, prior_range = c(2, 20)), "earlier than")
  expect_error(strauss1989_bayes(sd, prior_range = c(-Inf, 20)), "finite")
  expect_error(strauss1989_bayes(sd, prior_range = c(0, 9)), "later than")
  expect_error(strauss1989_bayes(sd, alpha = 0, prior_range = c(0, 20)), "must be in \\(0, 1\\)")
  expect_error(strauss1989_bayes(sd, alpha = 1, prior_range = c(0, 20)), "must be in \\(0, 1\\)")
  sd2 <- sighting_data(data.frame(time = c(1, 3, 5), count = c(1, 2, 0)))
  expect_error(strauss1989_bayes(sd2, prior_range = c(0, 20)), "at least 3")
})
