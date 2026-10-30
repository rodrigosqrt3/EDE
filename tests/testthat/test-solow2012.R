# Ivory-billed Woodpecker record (Roberts et al. 2010) as listed in Solow &
# Beet (2014, p. 3): one sighting per year. "Physical" sightings are certain
# in both papers (1897 marks the start of the observation period). Solow et al.
# (2012) also treated the 17 expert sightings up to 1944 as certain, leaving 29
# uncertain sightings, all after 1944; Solow & Beet (2014) treat all 46
# non-physical sightings as uncertain.
ibw_physical <- function() {
  c(1897, 1898:1902, 1904:1910, 1913, 1914, 1917, 1924, 1925, 1932, 1935,
    1938, 1939)
}
ibw_other <- function() {
  c(1911, 1916, 1920, 1921, 1923, 1926, 1929:1931, 1933, 1934, 1936, 1937,
    1941:1944, 1946, 1948:1952, 1955, 1958, 1959, 1962, 1966:1969,
    1971:1974, 1976, 1981, 1982, 1985:1988, 1999, 2004:2006)
}
years <- function(y) sighting_data(data.frame(year = y, sightings = 1))
ibw2012 <- function() {
  other <- ibw_other()
  list(certain = years(c(ibw_physical(), other[other <= 1944])),
       uncertain = years(other[other > 1944]))
}

test_that("woodpecker record has the counts given in the papers", {
  expect_length(ibw_physical(), 22)
  expect_length(ibw_other(), 17 + 29)
  r <- ibw2012()
  expect_equal(sum(r$certain$count), 39)
  expect_equal(sum(r$uncertain$count), 29)
})

test_that("solow2012_bayes reproduces the omega posterior (Solow et al. 2012, Fig. 2)", {
  r <- ibw2012()
  om <- solow2012_bayes(r$certain, r$uncertain, 2010, tau_prior = "exponential",
                        exp_quantile = c(1940, 0.5), data_out = "omega")
  # density is normalised by p(t | E) itself, so it integrates to ~1
  expect_lt(abs(sum(om$density) * 0.005 - 1), 0.005)
  cdf <- cumsum(om$density) * 0.005
  q <- stats::approx(cdf, om$omega, c(0.25, 0.75))$y
  # paper p. 183: mode 0.67, quartiles 0.64 and 0.77
  expect_lt(abs(om$omega[which.max(om$density)] - 0.67), 0.0051)  # 0.665 grid
  expect_equal(round(q[1], 2), 0.64)
  expect_lt(abs(q[2] - 0.77), 0.011)
})

test_that("solow2012_bayes extinction-time posterior has the paper's modes (Fig. 1)", {
  r <- ibw2012()
  ex <- solow2012_bayes(r$certain, r$uncertain, 2010, tau_prior = "exponential",
                        exp_quantile = c(1940, 0.5), data_out = TRUE)
  expect_named(ex, c("time", "likelihood_ratio", "density"))
  expect_equal(range(ex$time), c(1944, 2010))
  # main modes in 1944 (boundary), 1952 and 1988 (paper p. 182)
  top <- ex$time[order(-ex$density)][1]
  expect_equal(floor(top), 1988)
  around <- function(a, b) max(ex$density[ex$time >= a & ex$time < b])
  expect_gt(around(1952, 1953), around(1953.5, 1954.5))
  expect_gt(ex$density[1], around(1944.5, 1945.5))
})

test_that("solow2012_bayes supports extinction for the woodpecker (Solow et al. 2012)", {
  r <- ibw2012()
  e <- solow2012_bayes(r$certain, r$uncertain, 2010, tau_prior = "exponential",
                       exp_quantile = c(1940, 0.5))
  u <- solow2012_bayes(r$certain, r$uncertain, 2010)
  # Paper: B = 4.03 (exponential) and 4.97 (uniform). Not reproduced exactly
  # (see ?solow2012_bayes); the package values are pinned as regression values.
  expect_equal(round(e$bayes_factor_extinction, 2), 3.28)
  expect_equal(round(u$bayes_factor_extinction, 2), 3.87)
  expect_gt(e$bayes_factor_extinction, 3)  # "substantial" evidence, as in paper
  expect_gt(u$bayes_factor_extinction, e$bayes_factor_extinction)
  expect_equal(e$bayes_factor, 1 / e$bayes_factor_extinction)
  expect_equal(e$posterior, 1 / (1 + 1 / e$bayes_factor))
  expect_s3_class(e, "ede_bayes")
})

test_that("solow_beet2014_bayes gives the paper's conclusions (Solow & Beet 2014)", {
  cert <- years(ibw_physical())
  unc <- years(ibw_other())
  bf <- function(model, q) {
    tp <- if (is.null(q)) "uniform" else "exponential"
    solow_beet2014_bayes(cert, unc, 2010, model = model, tau_prior = tp,
                         exp_quantile = q)$bayes_factor_extinction
  }
  m1 <- c(bf(1, c(1940, 0.5)), bf(1, NULL), bf(1, c(1940, 0.9)))
  m2 <- c(bf(2, c(1940, 0.5)), bf(2, NULL), bf(2, c(1940, 0.9)))
  # Paper p. 3-4: model 1 0.13, 0.19, 0.007; model 2 ~4e6, 4e6, 3e6.
  # Not reproduced exactly; package values pinned as regression values.
  expect_equal(signif(m1, 2), c(0.24, 0.36, 0.013))
  expect_equal(signif(m2, 2), c(3.3e6, 3.1e6, 1.8e6))
  expect_true(all(m1 < 1))       # model 1: evidence against extinction
  expect_true(all(m2 > 1e6))     # model 2: overwhelming evidence for it
  expect_true(m1[3] < m1[1] && m1[1] < m1[2])
})

test_that("model 2 without uncertain sightings reduces to solow1993_bayes", {
  cert <- years(ibw_physical())
  sb <- solow_beet2014_bayes(cert, NULL, 2010, model = 2)
  s93 <- solow1993_bayes(cert, 2010)
  expect_equal(sb$bayes_factor, s93$bayes_factor, tolerance = 1e-6)
})

test_that("the omega/lambda integration matches brute force on equation 5", {
  # toy record: origin 0, certain at 2, 4, 5; uncertain at 6, 7, 9; T = 12
  cert <- sighting_data(data.frame(t = c(0, 2, 4, 5), n = 1))
  unc <- sighting_data(data.frame(t = c(6, 7, 9), n = 1))
  rec <- uncertain_record(cert, unc, 12, "test")
  brute <- function(tau, form) {
    k <- sum(rec$u_time > tau)
    nc <- rec$nc; nu <- rec$nu; n <- nc + nu; Tt <- rec$T; tL <- rec$tL
    inner <- function(w) vapply(w, function(wi) {
      f <- function(l) {
        if (form == "s2012") {
          # Solow et al. (2012) eq 5, times the 1 / lambda prior
          l^(n - 1) * wi^(nc - n) * (1 - wi)^k *
            exp(-l * (tau + (1 - wi) / wi * (Tt - tL)))
        } else {
          # Solow & Beet (2014) eq 5 (binomial sum form), times 1 / lambda
          nb <- nu - k
          s <- sum(choose(nb, 0:nb) * ((1 - wi) / wi)^(n - (nc + 0:nb)))
          l^(n - 1) * exp(-l * (tau + (1 - wi) / wi * Tt)) * s
        }
      }
      stats::integrate(f, 0, Inf, rel.tol = 1e-12)$value
    }, numeric(1))
    stats::integrate(inner, 0, 1, rel.tol = 1e-10)$value
  }
  s12 <- list(a = rec$nc, n = rec$nc + rec$nu, L = rec$T - rec$tL, certain_power = 0)
  m1 <- list(a = rec$nc, n = rec$nc + rec$nu, L = rec$T, certain_power = 0)
  for (tau in c(5.5, 6.5, 8, 10)) {
    expect_equal(
      uncertain_loglik(tau, rec, s12) - uncertain_loglik(12, rec, s12),
      log(brute(tau, "s2012") / brute(12, "s2012")), tolerance = 1e-6)
    expect_equal(
      uncertain_loglik(tau, rec, m1) - uncertain_loglik(12, rec, m1),
      log(brute(tau, "sb1") / brute(12, "sb1")), tolerance = 1e-6)
  }
})

test_that("exponential prior matches its quantile and allows increasing densities", {
  r <- ibw2012()
  rec <- uncertain_record(r$certain, r$uncertain, 2010, "test")
  for (p in c(0.5, 0.9, 0.2)) {
    pri <- uncertain_tau_prior("exponential", c(1940, p), rec)
    expect_equal(stats::integrate(pri$density, 0, 43)$value, p, tolerance = 1e-8)
    expect_equal(stats::integrate(pri$density, 0, 113)$value, 1, tolerance = 1e-8)
  }
  # a prior probability equal to the uniform one gives the uniform prior
  a <- solow2012_bayes(r$certain, r$uncertain, 2010)
  b <- solow2012_bayes(r$certain, r$uncertain, 2010, tau_prior = "exponential",
                       exp_quantile = c(1940, 43 / 113))
  expect_equal(a$bayes_factor, b$bayes_factor, tolerance = 1e-6)
  # a user-supplied prior density, normalised internally
  f <- solow2012_bayes(r$certain, r$uncertain, 2010,
                       tau_prior = function(x) rep(5, length(x)))
  expect_equal(f$bayes_factor, a$bayes_factor, tolerance = 1e-6)
})

test_that("uncertain-sighting methods validate their inputs", {
  r <- ibw2012()
  expect_error(solow2012_bayes(r$certain, r$uncertain, 2010, prior = 1), "prior")
  expect_error(solow2012_bayes(r$certain, r$uncertain), "test_year")
  expect_error(solow2012_bayes(r$certain, r$uncertain, 2006), "later than")
  expect_error(solow2012_bayes(data.frame(1), r$uncertain, 2010), "sighting_data")
  expect_error(solow2012_bayes(r$certain, r$uncertain, 2010,
                               tau_prior = "exponential"), "exp_quantile")
  expect_error(solow2012_bayes(r$certain, r$uncertain, 2010,
                               tau_prior = "exponential",
                               exp_quantile = c(1800, 0.5)), "inside")
  expect_error(solow2012_bayes(r$certain, r$uncertain, 2010, data_out = "x"),
               "data_out")
  # Solow et al. (2012) needs uncertain sightings after the last certain one
  expect_error(solow2012_bayes(years(ibw_physical()), years(ibw_other()), 2010),
               "solow_beet2014_bayes")
  expect_error(solow_beet2014_bayes(r$certain, r$uncertain, 2010), "model")
  expect_error(solow_beet2014_bayes(r$certain, r$uncertain, 2010, model = 3),
               "model")
  early <- years(c(1890, 1950))
  expect_error(solow_beet2014_bayes(r$certain, early, 2010, model = 1), "origin")
  one <- years(1900)
  expect_error(solow2012_bayes(one, NULL, 2010), "at least 1")
})

test_that("print.ede_bayes works for uncertain-sighting results", {
  r <- ibw2012()
  out <- capture.output(print(solow2012_bayes(r$certain, r$uncertain, 2010)))
  expect_match(out[1], "Solow et al. \\(2012\\)")
  expect_match(out[2], "Bayes factor")
})
