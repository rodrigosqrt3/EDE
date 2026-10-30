sd_years <- function(year, count = 1) {
  sighting_data(data.frame(year = year, count = count))
}
near <- function(x, y, tol) expect_lt(max(abs(x - y)), tol)

# Greater stick-nest rat Leporillus conditor, Lee et al. (2014) Table 2
# (from Burgman et al. 1995); 1835 is the start of observation
rat_years <- c(1844, 1857, 1879, 1896, 1897, 1907, 1920, 1921, 1930, 1933,
               1938, 1986)

test_that("L. conditor reproduces the 0.697 of the text and Fig. 3", {
  # Fig. 3 convention: origin sighting not counted, 2011 - t0 periods
  rat <- sd_years(c(1835, rat_years), c(0, rep(1, 12)))
  res <- lee2014(rat, test_year = 2010)
  expect_s3_class(res, "ede_bayes")
  # p. 256: P(extant) = 0.697 (MCMC); exact value 0.6955
  near(res$posterior, 0.697, 0.002)
  expect_equal(res$periods, 176)
  expect_equal(sum(res$extinction_time$probability) + res$posterior, 1)
  expect_equal(range(res$extinction_time$time), c(1987, 2010))
  expect_equal(res$posterior,
               1 / (1 + (1 - res$prior) / (res$prior * res$bayes_factor)))
})

test_that("L. conditor with the 1986 sighting uncertain is close to 0.025", {
  rat_c <- sd_years(c(1835, rat_years[-12]), c(0, rep(1, 11)))
  rat_u <- sd_years(1986)
  res <- lee2014(rat_c, rat_u, test_year = 2010)
  # p. 256: 0.697 drops to 0.025 (OpenBUGS MCMC); exact value 0.0316
  # (0.0295 with observation through 2011). Not reproduced to MCMC precision.
  near(res$posterior, 0.025, 0.007)
  expect_match(res$method, "uncertain")
})

test_that("Fig. 3 black bars are reproduced (bar heights read from figure)", {
  # Heights measured from the rendered figure (+-0.003); Fig. 3 convention
  trop <- sd_years(c(1915, 1922, 1932, 1948, 1952), c(0, 1, 1, 1, 1))
  near(lee2014(trop, test_year = 2010)$posterior, 0.078, 0.003)
  nig <- sd_years(setdiff(1972:1984, c(1978, 1980)), c(0, rep(1, 10)))
  expect_lt(lee2014(nig, test_year = 2010)$posterior, 0.005)
  api <- sd_years(c(1835, 1853, 1857, 1872, 1873, 1903, 1929, 1933),
                  c(0, rep(1, 7)))
  # bar 0.106, exact 0.100
  near(lee2014(api, test_year = 2010)$posterior, 0.106, 0.008)
})

test_that("Fig. 4 bars (ii) and (iii) are reproduced (read from figure)", {
  # Birdlife records of Table 3, observation from t0 to 2011; bar heights
  # measured from the rendered figure (+-0.005). (ii) certain only,
  # (iii) uncertain sightings grouped together.
  ruf_c <- sd_years(c(1929, 1960, 1963, 1969))
  ruf_u <- sd_years(c(1947, 1970, 1971, 1972, 1982, 1985, 1986, 1988))
  near(lee2014(ruf_c, test_year = 2011, prior = 0.9)$posterior, 0.665, 0.01)
  near(lee2014(ruf_c, test_year = 2011)$posterior, 0.177, 0.01)
  near(lee2014(ruf_c, ruf_u, test_year = 2011, prior = 0.9)$posterior,
       0.63, 0.01)
  near(lee2014(ruf_c, ruf_u, test_year = 2011)$posterior, 0.16, 0.01)

  car_c <- sd_years(c(1789, 1829, 1866, 1879), c(0, 1, 1, 1))
  car_u <- sd_years(c(1789, 1847, 1891))
  near(lee2014(car_c, test_year = 2011, prior = 0.825)$posterior, 0.504, 0.01)
  near(lee2014(car_c, test_year = 2011)$posterior, 0.181, 0.01)
  near(lee2014(car_c, car_u, test_year = 2011, prior = 0.825)$posterior,
       0.187, 0.01)
  near(lee2014(car_c, car_u, test_year = 2011)$posterior, 0.048, 0.01)

  pel_c <- sd_years(c(1930, 1995))
  pel_u <- sd_years(2008)
  near(lee2014(pel_c, test_year = 2011, prior = 0.36)$posterior, 0.668, 0.01)
  near(lee2014(pel_c, test_year = 2011)$posterior, 0.776, 0.01)
  near(lee2014(pel_c, pel_u, test_year = 2011, prior = 0.36)$posterior,
       0.758, 0.01)
})

test_that("exact computation matches brute-force integration (Appendix S1 data)", {
  # Appendix S1 example: y = 1,0,0,0,0,1,0,0,0; z = 0,1,0,1,0,1,0,1,0; T = 9
  y <- c(1, 0, 0, 0, 0, 1, 0, 0, 0)
  z <- c(0, 1, 0, 1, 0, 1, 0, 1, 0)
  M <- 10
  res <- lee2014(sd_years(which(y == 1)), sd_years(c(1, which(z == 1)),
                                                   c(0, rep(1, 4))),
                 test_year = 9, rate_max = M)
  lik <- function(alive) {
    pc <- function(m) {
      vapply(m, function(mm) {
        p <- ifelse(seq_along(y) <= alive, 1 - exp(-mm), 0)
        prod(ifelse(y == 1, p, 1 - p))
      }, numeric(1))
    }
    pu <- function(m, f) {
      p <- 1 - exp(-(ifelse(seq_along(z) <= alive, m, 0) + f))
      prod(ifelse(z == 1, p, 1 - p))
    }
    lc <- stats::integrate(pc, 0, M, rel.tol = 1e-8)$value / M
    inner <- function(f) {
      vapply(f, function(ff) {
        stats::integrate(function(m) vapply(m, pu, numeric(1), f = ff), 0, M,
                         rel.tol = 1e-8)$value
      }, numeric(1))
    }
    lu <- stats::integrate(inner, 0, M, rel.tol = 1e-8)$value / M^2
    lc * lu
  }
  l0 <- lik(9)
  l1 <- mean(vapply(0:8, lik, numeric(1)))
  near(res$bayes_factor, l0 / l1, 1e-5 * l0 / l1)
})

test_that("certain-only likelihood matches the incomplete Beta form", {
  # (1/M) int_0^M (1 - e^-m)^S e^(-m (k - S)) dm, eqns 9-11
  f <- function(S, k, M) {
    stats::integrate(function(m) (1 - exp(-m))^S * exp(-m * (k - S)), 0, M,
                     rel.tol = 1e-12)$value / M
  }
  near(exp(lee2014_lcert(3, 10, 5)), f(3, 10, 5), 1e-12)
  near(exp(lee2014_lcert(3, 3, 5)), f(3, 3, 5), 1e-12)
  near(exp(lee2014_lcert(0, 0, 5)), 1, 1e-12)
  expect_equal(lee2014_lcert(4, 3, 5), -Inf)
  # cancellation fallback (small rate bound, many sightings)
  near(exp(lee2014_lcert(40, 40, 0.5)), f(40, 40, 0.5),
       1e-6 * f(40, 40, 0.5))
})

test_that("records with no certain sightings and the persistence limit work", {
  cert <- sd_years(2000, 0)
  unc <- sd_years(c(2001, 2003, 2004))
  res <- lee2014(cert, unc, test_year = 2010)
  expect_true(res$posterior > 0 && res$posterior < 1)
  expect_equal(res$extinction_time$time, 2000:2010)
  # after a sighting in every period, one empty period is strong evidence of
  # extinction under the wide uniform rate prior (high rates are plausible)
  every <- lee2014(sd_years(2000:2009), test_year = 2010)
  expect_lt(every$posterior, 0.05)
  expect_equal(every$extinction_time$time, 2010)
})

test_that("input validation", {
  s <- sd_years(c(2000, 2005))
  expect_error(lee2014(data.frame(), test_year = 2010), "sighting_data")
  expect_error(lee2014(s, data.frame(), test_year = 2010), "sighting_data")
  expect_error(lee2014(s, test_year = 2010, prior = 1), "prior")
  expect_error(lee2014(s, test_year = 2010, rate_max = 0), "rate_max")
  expect_error(lee2014(sd_years(c(2000, 2005.5)), test_year = 2010),
               "whole-number")
  expect_error(lee2014(sd_years(2000, 0), test_year = 2010), "no sightings")
  expect_error(lee2014(s, test_year = 2005), "later than")
  expect_error(lee2014(s, sd_years(2008), test_year = 2008), "later than")
  expect_error(lee2014(s, test_year = 2010.5), "whole number")
})
